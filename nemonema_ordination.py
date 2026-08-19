"""
nemonema_ordination.py — constrained ordination of nematode communities on
soil variables.

WHY THIS EXISTS
Correlating every index against every soil variable produces hundreds of tests
and a near-certain false positive. Redundancy analysis asks the question once:
how much of the variation in community composition is explained by the soil
variables together?

  rda()                  redundancy analysis with a permutation test
  variance_partition()   how much each block of variables explains, and how much
                         they share
  forward_selection()    which soil variables earn their place, with correction
  soil_pca()             reduce many correlated soil variables to a few axes

METHOD
Community data are Hellinger-transformed before analysis, which makes ordinary
linear methods appropriate for species abundances (Legendre & Gallagher 2001,
Oecologia 129:271-280). Explanatory variables are standardised. Significance
comes from permuting rows, so no distributional assumption is made about the
response.

Legendre, P. & Legendre, L. (2012) Numerical Ecology, 3rd English edition.
Elsevier. Chapter 11.
Peres-Neto, P.R. et al. (2006) Ecology 87:2614-2625 (adjusted R2 for RDA).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _prep(counts: pd.DataFrame, env: pd.DataFrame):
    """Align samples, Hellinger-transform the community, standardise the soil
    variables, and drop anything constant."""
    common = counts.columns.intersection(env.index)
    Y = counts[common].T.astype(float)
    X = env.loc[common].select_dtypes(include=[np.number])
    X = X.loc[:, X.notna().sum() >= max(3, 0.6 * len(X))]
    X = X.fillna(X.mean())
    X = X.loc[:, X.std(ddof=1) > 0]
    tot = Y.sum(axis=1).replace(0, np.nan)
    Y = np.sqrt(Y.div(tot, axis=0)).fillna(0)       # Hellinger
    Y = Y - Y.mean(axis=0)
    Xs = (X - X.mean()) / X.std(ddof=1)
    return Y, Xs, list(common)


def _rda_core(Y: pd.DataFrame, X: pd.DataFrame):
    """Fitted and residual variation from a least-squares fit of Y on X."""
    Xm = np.column_stack([np.ones(len(X)), X.values])
    B, *_ = np.linalg.lstsq(Xm, Y.values, rcond=None)
    fitted = Xm @ B
    ss_total = float((Y.values ** 2).sum())
    ss_fit = float((fitted ** 2).sum())
    return ss_fit, ss_total, fitted


def rda(counts: pd.DataFrame, env: pd.DataFrame, permutations: int = 999,
        seed: int = 0) -> dict:
    """Redundancy analysis of community composition on soil variables.

    R2_adj is the honest figure: raw R2 rises mechanically with the number of
    explanatory variables, and with many variables and few samples it can reach
    1.0 while explaining nothing (Peres-Neto et al. 2006).
    """
    Y, X, common = _prep(counts, env)
    n, m = len(Y), X.shape[1]
    if n < 4:
        return {"error": f"only {n} samples match between the two tables"}
    if m == 0:
        return {"error": "no usable soil variables (all constant or missing)"}
    if m >= n - 1:
        return {"error": f"{m} soil variables and only {n} samples. RDA needs "
                         f"fewer variables than samples minus one, otherwise the "
                         f"fit is saturated and R2 is meaningless. Reduce the "
                         f"variables first (see soil_pca) or add samples."}

    ss_fit, ss_tot, fitted = _rda_core(Y, X)
    R2 = ss_fit / ss_tot if ss_tot else np.nan
    R2adj = 1 - (1 - R2) * (n - 1) / (n - m - 1)
    F = (ss_fit / m) / ((ss_tot - ss_fit) / (n - m - 1))

    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(permutations):
        Yp = Y.iloc[rng.permutation(n)].reset_index(drop=True)
        Yp.index = Y.index
        sf, stt, _ = _rda_core(Yp, X)
        Fp = (sf / m) / ((stt - sf) / (n - m - 1)) if stt > sf else np.inf
        if Fp >= F:
            cnt += 1

    # constrained axes, for plotting
    U, S, Vt = np.linalg.svd(fitted - fitted.mean(axis=0), full_matrices=False)
    k = min(2, (S > 1e-9).sum())
    scores = pd.DataFrame(U[:, :k] * S[:k], index=Y.index,
                          columns=[f"RDA{i+1}" for i in range(k)])
    axis_pct = (S ** 2 / (S ** 2).sum() * 100)[:k] if S.sum() else np.array([np.nan] * k)
    # variable loadings: correlation of each soil variable with each axis
    load = pd.DataFrame(
        {c: [np.corrcoef(X[c].values, scores.iloc[:, a])[0, 1] for a in range(k)]
         for c in X.columns}, index=scores.columns).T

    return {"n": n, "n_variables": m, "R2": float(R2), "R2_adjusted": float(R2adj),
            "F": float(F), "p": float((cnt + 1) / (permutations + 1)),
            "permutations": permutations, "scores": scores,
            "axis_pct": axis_pct, "loadings": load, "variables": list(X.columns)}


def forward_selection(counts: pd.DataFrame, env: pd.DataFrame,
                      permutations: int = 199, alpha: float = 0.05,
                      seed: int = 0) -> pd.DataFrame:
    """Add soil variables one at a time, keeping only those that improve the fit.

    Each candidate is tested by permutation, and the threshold is tightened by
    the number of candidates still available (a Bonferroni step) so that
    selecting from many variables does not manufacture significance.
    """
    Y, X, _ = _prep(counts, env)
    n = len(Y)
    remaining = list(X.columns)
    chosen, rows = [], []
    while remaining and len(chosen) < n - 2:
        best = None
        for cand in remaining:
            cols = chosen + [cand]
            sf, stt, _ = _rda_core(Y, X[cols])
            R2 = sf / stt
            if best is None or R2 > best[1]:
                best = (cand, R2)
        cand, R2 = best
        cols = chosen + [cand]
        m = len(cols)
        F = (R2 / m) / ((1 - R2) / (n - m - 1))
        rng = np.random.default_rng(seed + len(chosen))
        cnt = 0
        for _ in range(permutations):
            Yp = Y.iloc[rng.permutation(n)]
            Yp.index = Y.index
            sf, stt, _ = _rda_core(Yp, X[cols])
            R2p = sf / stt
            Fp = (R2p / m) / ((1 - R2p) / (n - m - 1)) if R2p < 1 else np.inf
            if Fp >= F:
                cnt += 1
        p = (cnt + 1) / (permutations + 1)
        thresh = alpha / max(1, len(remaining))
        rows.append({"step": len(chosen) + 1, "variable": cand,
                     "cumulative_R2": round(R2, 4),
                     "R2_adjusted": round(1 - (1 - R2) * (n - 1) / (n - m - 1), 4),
                     "F": round(F, 3), "p": round(p, 4),
                     "threshold": round(thresh, 5),
                     "retained": p <= thresh})
        if p > thresh:
            break
        chosen.append(cand)
        remaining.remove(cand)
    return pd.DataFrame(rows)


def variance_partition(counts: pd.DataFrame, blocks: dict,
                       env: pd.DataFrame) -> pd.DataFrame:
    """Partition explained variation between two sets of variables.

    blocks: {"label A": [col, col], "label B": [col, ...]}
    Reports each block alone, both together, the unique part of each, and the
    shared part they cannot be credited to either.
    """
    keys = list(blocks)
    if len(keys) != 2:
        return pd.DataFrame([{"error": "supply exactly two blocks"}])
    Y, X, _ = _prep(counts, env)
    a = [c for c in blocks[keys[0]] if c in X.columns]
    b = [c for c in blocks[keys[1]] if c in X.columns]
    if not a or not b:
        return pd.DataFrame([{"error": "a block has no usable variables"}])
    n = len(Y)

    def adj(cols):
        sf, stt, _ = _rda_core(Y, X[cols])
        R2 = sf / stt
        m = len(cols)
        return 1 - (1 - R2) * (n - 1) / (n - m - 1)

    Ra, Rb, Rab = adj(a), adj(b), adj(a + b)
    shared = Ra + Rb - Rab
    return pd.DataFrame([
        {"component": f"{keys[0]} alone (unique)", "adjusted_R2": round(Rab - Rb, 4)},
        {"component": f"{keys[1]} alone (unique)", "adjusted_R2": round(Rab - Ra, 4)},
        {"component": "shared between both", "adjusted_R2": round(shared, 4)},
        {"component": "unexplained", "adjusted_R2": round(1 - Rab, 4)},
        {"component": "both together", "adjusted_R2": round(Rab, 4)},
    ])


def soil_pca(env: pd.DataFrame, n_axes: int = 3):
    """Reduce many correlated soil variables to a few uncorrelated axes.

    Use when you have more soil variables than RDA can accept. Correlate against
    the axes instead of the raw variables: three tests rather than thirty.
    """
    X = env.select_dtypes(include=[np.number])
    X = X.loc[:, X.std(ddof=1) > 0].fillna(X.mean())
    Xs = (X - X.mean()) / X.std(ddof=1)
    U, S, Vt = np.linalg.svd(Xs.values - Xs.values.mean(axis=0),
                             full_matrices=False)
    k = min(n_axes, len(S))
    scores = pd.DataFrame(U[:, :k] * S[:k], index=X.index,
                          columns=[f"SoilPC{i+1}" for i in range(k)])
    var = (S ** 2) / (S ** 2).sum() * 100
    load = pd.DataFrame(Vt[:k].T, index=X.columns, columns=scores.columns)
    return scores, var[:k], load


def triplot(res: dict, groups=None, pal=None, ax=None):
    """RDA sample scores with soil-variable arrows."""
    import matplotlib.pyplot as plt
    if "error" in res:
        return None
    if ax is None:
        _, ax = plt.subplots(figsize=(6.6, 5.6))
    sc, load, pct = res["scores"], res["loadings"], res["axis_pct"]
    if groups is None:
        ax.scatter(sc.iloc[:, 0], sc.iloc[:, 1], s=70, color="#2E7D5B",
                   edgecolor="w", zorder=3)
    else:
        levs = list(pd.unique(groups))
        cols = pal(len(levs)) if callable(pal) else \
            ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7"][:len(levs)]
        for i, lv in enumerate(levs):
            m = groups.reindex(sc.index) == lv
            ax.scatter(sc.loc[m].iloc[:, 0], sc.loc[m].iloc[:, 1], s=70,
                       color=cols[i % len(cols)], label=str(lv), edgecolor="w",
                       zorder=3)
        ax.legend(frameon=False, fontsize=8, loc="upper center",
                  ncol=min(len(levs), 5), bbox_to_anchor=(0.5, -0.13))
    scale = 0.8 * max(abs(sc.values).max(), 1e-9)
    for v in load.index:
        x, y = load.loc[v, load.columns[0]] * scale, load.loc[v, load.columns[1]] * scale
        ax.arrow(0, 0, x, y, color="#A03030", width=0.002,
                 head_width=scale * 0.045, length_includes_head=True, zorder=4)
        ax.text(x * 1.09, y * 1.09, v, color="#A03030", fontsize=8, ha="center",
                va="center")
    ax.axhline(0, color="grey", lw=.6, ls=":")
    ax.axvline(0, color="grey", lw=.6, ls=":")
    ax.set_xlabel(f"RDA1 ({pct[0]:.1f}% of constrained variation)")
    ax.set_ylabel(f"RDA2 ({pct[1]:.1f}%)" if len(pct) > 1 else "RDA2")
    ax.set_title("Redundancy analysis: community constrained by soil variables",
                 loc="left", fontweight="bold", fontsize=10)
    return ax
