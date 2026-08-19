"""
nemonema_qc.py — data-quality checks and design-aware statistics.

Nothing here changes a computed value. It reports what the numbers rest on.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# Degeneracy and dominance
# --------------------------------------------------------------------------
def degeneracy(counts, taxa) -> pd.DataFrame:
    """Indices that compute but carry no information for this dataset.

    Worked example: with no fungivores present, CI = 100 x 0.8 Fu2 /
    (3.2 Ba1 + 0.8 Fu2) evaluates to 0. That is arithmetically correct, but it
    reflects the taxon list rather than the soil, and must not be reported as
    'decomposition is entirely bacterial'.
    """
    tr = taxa["trophic"].reindex(counts.index)
    cp = pd.to_numeric(taxa["cp"], errors="coerce").reindex(counts.index)
    tot = counts.fillna(0).astype(float).sum(axis=1)
    grp = {g: float(tot[tr == g].sum()) for g in ["PP", "BF", "FF", "OM", "PR"]}
    ba1 = float(tot[(tr == "BF") & (cp == 1)].sum())
    fu2 = float(tot[(tr == "FF") & (cp == 2)].sum())
    struct = float(tot[((tr.isin(["BF", "FF", "OM"])) & (cp >= 3))
                       | ((tr == "PR") & (cp >= 2))].sum())
    ncp = int(cp[tot > 0].nunique())
    rows = []

    def add(idx, ok, why, act):
        rows.append({"index": idx, "informative": ok, "reason": why, "action": act})

    add("CI", grp["FF"] > 0,
        "no fungivores in the taxon list" if grp["FF"] == 0 else "fungivores present",
        "computes as 0; do not interpret as a bacterial channel"
        if grp["FF"] == 0 else "usable")
    add("NCR", grp["FF"] > 0,
        "no fungivores" if grp["FF"] == 0 else "fungivores present",
        "computes as 1.0; carries no information" if grp["FF"] == 0 else "usable")
    add("EI", (ba1 + fu2) > 0,
        "no Ba1 or Fu2 guilds" if (ba1 + fu2) == 0 else "enrichment guilds present",
        "undefined — omit" if (ba1 + fu2) == 0 else "usable")
    add("SI", struct > 0,
        "no structure-indicator guilds" if struct == 0 else "structure guilds present",
        "undefined — omit" if struct == 0 else "usable")
    add("PPI", grp["PP"] > 0, "no plant parasites" if grp["PP"] == 0 else "present",
        "undefined — omit" if grp["PP"] == 0 else "usable")
    add("MI", grp["PP"] < sum(grp.values()),
        "every taxon is a plant parasite; MI uses free-living taxa only"
        if grp["PP"] == sum(grp.values()) else "free-living taxa present",
        "undefined — omit" if grp["PP"] == sum(grp.values()) else "usable")
    add("faunal indices overall", ncp >= 3, f"only {ncp} c-p class(es) represented",
        "interpret with strong caution — faunal analysis assumes a spread of "
        "c-p classes" if ncp < 3 else "usable")
    out = pd.DataFrame(rows)
    out.attrs["trophic_totals"] = grp
    return out


def dominance(counts) -> pd.DataFrame:
    c = counts.fillna(0).astype(float)
    rel = c.div(c.sum(axis=0).replace(0, np.nan), axis=1) * 100
    o = pd.DataFrame({"dominant_taxon": rel.idxmax(), "dominant_pct": rel.max(),
                      "taxa_present": (c > 0).sum(axis=0)})
    o["flag"] = np.where(o["dominant_pct"] > 70, "extreme dominance",
                  np.where(o["dominant_pct"] > 50, "high dominance", ""))
    return o


# --------------------------------------------------------------------------
# Design
# --------------------------------------------------------------------------
def replication(samples, group_col="group") -> dict:
    cols = {c.lower(): c for c in samples.columns}
    gc = cols.get(group_col)
    if gc is None:
        return {"ok": False, "message": "No `group` column in the samples sheet."}
    counts = samples[gc].value_counts().sort_index()
    smallest = int(counts.min())
    out = {"per_group": counts.to_dict(), "n_groups": int(len(counts)),
           "smallest": smallest, "balanced": bool(counts.nunique() == 1)}
    if smallest < 2:
        out["ok"] = False
        out["message"] = (f"A group has n={smallest}. No dispersion can be "
                          "estimated and no test can run.")
    elif smallest < 3:
        out["ok"] = True
        out["message"] = (f"Smallest group n={smallest}. Tests run but are very "
                          "weakly powered; three is the usual minimum.")
    else:
        out["ok"] = True
        out["message"] = (f"n={smallest} or more per group"
                          + ("; balanced design."
                             if out["balanced"] else "; unbalanced — state this."))
    return out


def detect_blocking(samples, group_col="group") -> dict:
    """Find a column whose levels each appear exactly once per group.

    That structure means the design is paired or blocked. Analysing it as
    independent samples violates the independence assumption — the assumption
    that matters most.
    """
    cols = {c.lower(): c for c in samples.columns}
    gc = cols.get(group_col)
    if gc is None:
        return {"blocked": False}
    g = samples[gc]
    hits = []
    for c in samples.columns:
        if c == gc or samples[c].nunique() < 2:
            continue
        tab = pd.crosstab(samples[c], g)
        if tab.shape[1] > 1 and (tab.values == 1).all():
            hits.append(c)
    if not hits:
        return {"blocked": False,
                "message": "No column appears exactly once per group; samples "
                           "look independent."}
    return {"blocked": True, "column": hits[0], "candidates": hits,
            "n_blocks": int(samples[hits[0]].nunique()),
            "message": (f"'{hits[0]}' appears exactly once in every group. This is "
                        "a PAIRED or BLOCKED design. One-way ANOVA treats these as "
                        "independent samples. Use a paired test, or fit "
                        f"{hits[0]} as a random effect.")}


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
def compare(table, groups) -> pd.DataFrame:
    from scipy import stats as st
    g = groups.reindex(table.index)
    rows = []
    for col in table.select_dtypes(include=[np.number]).columns:
        arr = [table.loc[g == l, col].dropna().values for l in pd.unique(g.dropna())]
        arr = [a for a in arr if len(a) > 1]
        if len(arr) < 2:
            rows.append({"index": col, "n": np.nan, "df_between": np.nan,
                         "df_within": np.nan, "anova_F": np.nan, "anova_p": np.nan,
                         "eta_squared": np.nan, "kruskal_H": np.nan,
                         "kruskal_p": np.nan})
            continue
        N, k = sum(len(a) for a in arr), len(arr)
        try:
            F, pa = st.f_oneway(*arr)
        except Exception:
            F, pa = np.nan, np.nan
        try:
            H, pk = st.kruskal(*arr)
        except Exception:
            H, pk = np.nan, np.nan
        allv = np.concatenate(arr)
        gm = allv.mean()
        ssb = sum(len(a) * (a.mean() - gm) ** 2 for a in arr)
        sst = ((allv - gm) ** 2).sum()
        rows.append({"index": col, "n": N, "df_between": k - 1, "df_within": N - k,
                     "anova_F": F, "anova_p": pa,
                     "eta_squared": ssb / sst if sst else np.nan,
                     "kruskal_H": H, "kruskal_p": pk})
    out = pd.DataFrame(rows).set_index("index")
    out["anova_p_BH"] = benjamini_hochberg(out["anova_p"].values)
    out["kruskal_p_BH"] = benjamini_hochberg(out["kruskal_p"].values)
    return out


def benjamini_hochberg(pvalues) -> np.ndarray:
    """FDR adjustment. Testing 21 indices at alpha 0.05 without correction gives
    roughly a 66% chance of at least one false positive."""
    p = np.asarray(pvalues, dtype=float)
    ok = ~np.isnan(p)
    out = np.full(p.shape, np.nan)
    q = p[ok]
    m = len(q)
    if m == 0:
        return out
    order = np.argsort(q)
    adj = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, q[i] * m / (rank + 1))
        adj[i] = prev
    out[ok] = np.minimum(adj, 1.0)
    return out


def diagnostics(values, groups) -> dict:
    """Everything a referee expects beside an F statistic."""
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    levs = list(pd.unique(df["g"]))
    arr = [df.loc[df["g"] == l, "v"].values for l in levs]
    arr = [a for a in arr if len(a) > 1]
    if len(arr) < 2:
        return {"error": "fewer than two groups with n > 1"}
    # Shapiro on RESIDUALS, not raw values — testing raw values tests the wrong thing
    resid = np.concatenate([a - a.mean() for a in arr])
    k, N = len(arr), sum(len(a) for a in arr)
    F, p = st.f_oneway(*arr)
    allv = np.concatenate(arr)
    gm = allv.mean()
    ssb = sum(len(a) * (a.mean() - gm) ** 2 for a in arr)
    sst = ((allv - gm) ** 2).sum()
    sw = st.shapiro(resid).pvalue if 3 <= len(resid) <= 5000 else np.nan
    lv = st.levene(*arr).pvalue
    normal = not (sw == sw and sw < 0.05)
    equal = not (lv == lv and lv < 0.05)
    rec = ("ANOVA assumptions met — report the F test." if normal and equal else
           "Variances unequal — use Welch's ANOVA." if normal else
           "Residuals non-normal — use Kruskal-Wallis, a permutation test, or a "
           "GLM with an appropriate error family. For raw counts a negative "
           "binomial GLM is the principled choice.")
    return {"n": N, "n_groups": k, "df_between": k - 1, "df_within": N - k,
            "group_n": {str(l): int(len(a)) for l, a in zip(levs, arr)},
            "F": float(F), "p": float(p),
            "eta_squared": float(ssb / sst) if sst else np.nan,
            "shapiro_p_residuals": float(sw) if sw == sw else None,
            "levene_p": float(lv) if lv == lv else None,
            "residuals_normal": bool(normal), "variances_equal": bool(equal),
            "recommendation": rec}


def paired(values, groups, blocks) -> pd.DataFrame:
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index),
                       "b": blocks.reindex(values.index)}).dropna()
    levs = list(pd.unique(df["g"]))
    if len(levs) != 2:
        return pd.DataFrame([{"error": f"paired test needs 2 groups, found {len(levs)}"}])
    a = df[df.g == levs[0]].set_index("b")["v"]
    b = df[df.g == levs[1]].set_index("b")["v"]
    common = a.index.intersection(b.index)
    if len(common) < 3:
        return pd.DataFrame([{"error": "fewer than three complete pairs"}])
    x, y = a[common].values, b[common].values
    d = x - y
    t, pt = st.ttest_rel(x, y)
    try:
        W, pw = st.wilcoxon(x, y)
    except Exception:
        W, pw = np.nan, np.nan
    F, pf = st.f_oneway(x, y)
    sd = d.std(ddof=1)
    return pd.DataFrame([{
        "group_1": levs[0], "group_2": levs[1], "n_pairs": len(common),
        "mean_difference": float(d.mean()), "sd_of_differences": float(sd),
        "cohens_dz": float(d.mean() / sd) if sd else np.nan,
        "paired_t": float(t), "p_paired": float(pt),
        "wilcoxon_W": float(W), "p_wilcoxon": float(pw),
        "p_unpaired_for_comparison": float(pf)}])


# --------------------------------------------------------------------------
# Multivariate
# --------------------------------------------------------------------------
def permanova(dist, groups, permutations=999, seed=0) -> dict:
    """Anderson (2001) Austral Ecology 26:32-46."""
    rng = np.random.default_rng(seed)
    g = groups.reindex(dist.index).dropna()
    d2 = dist.loc[g.index, g.index].values.astype(float) ** 2
    n = len(g)
    lab = g.values
    levs = pd.unique(lab)
    a = len(levs)
    if a < 2 or n < 3:
        return {"error": "need at least two groups and three samples"}
    sst = d2.sum() / (2 * n)

    def within(L):
        return sum(d2[np.ix_(np.where(L == lv)[0], np.where(L == lv)[0])].sum()
                   / (2 * (L == lv).sum())
                   for lv in levs if (L == lv).sum() > 1)

    ssw = within(lab)
    F = ((sst - ssw) / (a - 1)) / (ssw / (n - a)) if ssw > 0 else np.nan
    cnt = 0
    for _ in range(permutations):
        L = rng.permutation(lab)
        w = within(L)
        Fp = ((sst - w) / (a - 1)) / (w / (n - a)) if w > 0 else np.nan
        if Fp >= F:
            cnt += 1
    return {"pseudo_F": float(F), "R2": float((sst - ssw) / sst) if sst else np.nan,
            "p": float((cnt + 1) / (permutations + 1)), "n": n, "n_groups": a,
            "permutations": permutations}


def permdisp(dist, groups, permutations=999, seed=0) -> dict:
    """Anderson (2006) Biometrics 62:245-253.

    Run WITH permanova, never instead of it: a significant PERMANOVA can mean
    the centroids differ, or merely that one group is more variable."""
    from scipy import stats as st
    rng = np.random.default_rng(seed)
    g = groups.reindex(dist.index).dropna()
    d = dist.loc[g.index, g.index].values.astype(float)
    n = len(g)
    J = np.eye(n) - np.ones((n, n)) / n
    G = J @ (-0.5 * d ** 2) @ J
    w, V = np.linalg.eigh(G)
    keep = w > 1e-9
    coords = V[:, keep] * np.sqrt(w[keep])
    lab = g.values
    resid = np.empty(n)
    for lv in pd.unique(lab):
        idx = np.where(lab == lv)[0]
        resid[idx] = np.sqrt(((coords[idx] - coords[idx].mean(axis=0)) ** 2).sum(axis=1))
    arr = [resid[lab == lv] for lv in pd.unique(lab)]
    arr = [x for x in arr if len(x) > 1]
    if len(arr) < 2:
        return {"error": "need two groups with n > 1"}
    F = st.f_oneway(*arr).statistic
    cnt = 0
    for _ in range(permutations):
        L = rng.permutation(lab)
        a2 = [resid[L == lv] for lv in pd.unique(lab)]
        a2 = [x for x in a2 if len(x) > 1]
        try:
            if st.f_oneway(*a2).statistic >= F:
                cnt += 1
        except Exception:
            pass
    return {"F": float(F), "p": float((cnt + 1) / (permutations + 1)),
            "mean_dispersion": {str(lv): float(resid[lab == lv].mean())
                                for lv in pd.unique(lab)}}


# --------------------------------------------------------------------------
def nsh_diagnostics(nsh_tbl) -> pd.DataFrame:
    cols = [c for c in nsh_tbl.columns if c.endswith("_score")]
    o = pd.DataFrame({"min": nsh_tbl[cols].min(), "max": nsh_tbl[cols].max(),
                      "mean": nsh_tbl[cols].mean().round(2),
                      "sd": nsh_tbl[cols].std().round(3),
                      "distinct": nsh_tbl[cols].nunique()})
    o["status"] = np.where(o["distinct"] == 1, "CONSTANT — no discrimination",
                    np.where(o["max"] - o["min"] == 1, "near-constant (1 step)",
                             "varies"))
    if "NSH" in nsh_tbl:
        v = nsh_tbl["NSH"].dropna()
        o.attrs["scale_pct"] = float((v.max() - v.min()) / 24 * 100) if len(v) else np.nan
        o.attrs["n_constant"] = int((o["distinct"] == 1).sum())
    return o.sort_values("sd", ascending=False)


def redundancy(table, threshold=0.98) -> pd.DataFrame:
    num = table.select_dtypes(include=[np.number])
    keep = [c for c in num.columns if num[c].nunique() > 1]
    if len(keep) < 2:
        return pd.DataFrame()
    corr = num[keep].corr(method="spearman")
    rows = [{"index_1": keep[i], "index_2": keep[j],
             "spearman_r": round(float(corr.iloc[i, j]), 4),
             "note": "near-duplicate — report one"}
            for i in range(len(keep)) for j in range(i + 1, len(keep))
            if abs(corr.iloc[i, j]) >= threshold]
    return pd.DataFrame(rows) if rows else pd.DataFrame()
