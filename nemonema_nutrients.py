"""
nemonema_nutrients.py — relate nematode indices to soil chemical and physical data.

WHY THIS MATTERS MORE THAN IT LOOKS
Nematode-based indices are validated largely against other nematode-based
indices, which is circular. Correlating them with variables measured
independently of the nematode community — organic carbon, mineral nitrogen,
available phosphorus, pH, yield — is the non-circular evidence. Ghaderi et al.
(2025) rest their case for the NSH index on exactly this: its correlations with
organic C, NO3-, NH4+, PO4 and potentially mineralisable N.

WHAT THIS MODULE REFUSES TO DO
It does not compute correlations between NSH and its own component indices
without flagging them. NSH is constructed from MI, PPI, EI, SI, BI, CI and
abundance, so a correlation between NSH and MI is arithmetic, not evidence.
Those pairs are labelled CIRCULAR and excluded from the headline results.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Indices that enter the NSH sum. A correlation between NSH and any of these is
# a property of the formula, not a finding.
NSH_COMPONENTS = {"MI", "PPI", "EI", "SI", "BI", "CI", "total_N", "NSH"}

# Common soil variables, for recognising a nutrient sheet automatically.
KNOWN_NUTRIENTS = {
    "ph", "ec", "organic_c", "organic_carbon", "oc", "soc", "organic_m",
    "organic_matter", "om", "labile_c", "no3", "no3_n", "nitrate", "nh4",
    "nh4_n", "ammonium", "total_n", "available_n", "pmn",
    "potentially_mineralisable_n", "p", "olsen_p", "colwell_p", "available_p",
    "po4", "k", "available_k", "ca", "mg", "s", "na", "cec", "pbi",
    "bulk_density", "moisture", "clay", "silt", "sand", "yield",
    "microbial_biomass_c", "mbc", "respiration",
}


def find_nutrient_columns(samples: pd.DataFrame) -> list:
    """Numeric columns in the samples sheet that look like soil measurements."""
    out = []
    for c in samples.columns:
        cl = str(c).strip().lower().replace(" ", "_").replace("-", "_")
        if not pd.api.types.is_numeric_dtype(samples[c]):
            continue
        if cl in {"replicate", "rep", "plot", "block", "n", "sample"}:
            continue
        if cl in KNOWN_NUTRIENTS or samples[c].nunique() > 2:
            out.append(c)
    return out


def load_nutrients(xl_or_samples, samples: pd.DataFrame | None = None):
    """Nutrients come either from a dedicated 'nutrients' sheet or from extra
    numeric columns already in the samples sheet."""
    if isinstance(xl_or_samples, pd.DataFrame):
        s = xl_or_samples
        cols = find_nutrient_columns(s)
        return (s[cols].apply(pd.to_numeric, errors="coerce") if cols
                else pd.DataFrame(index=s.index))
    try:
        nut = pd.read_excel(xl_or_samples, "nutrients", index_col=0)
        nut.index = nut.index.astype(str).str.strip()
        nut = nut.apply(pd.to_numeric, errors="coerce")
        return nut.dropna(axis=1, how="all")
    except Exception:
        return (load_nutrients(samples) if samples is not None
                else pd.DataFrame())


def correlate(indices: pd.DataFrame, nutrients: pd.DataFrame,
              method: str = "spearman", min_n: int = 5) -> pd.DataFrame:
    """Pairwise correlation of every index against every nutrient variable.

    Spearman by default: it does not assume linearity or normality, and
    nematode indices are frequently neither.
    """
    from scipy import stats as st
    common = indices.index.intersection(nutrients.index)
    if len(common) < min_n:
        return pd.DataFrame([{"error": f"only {len(common)} samples match between "
                                       f"the index and nutrient tables; at least "
                                       f"{min_n} are needed"}])
    I = indices.loc[common].select_dtypes(include=[np.number])
    N = nutrients.loc[common].select_dtypes(include=[np.number])
    rows = []
    for ic in I.columns:
        for nn in N.columns:
            x, y = I[ic], N[nn]
            ok = x.notna() & y.notna()
            n = int(ok.sum())
            if n < min_n or x[ok].nunique() < 3 or y[ok].nunique() < 3:
                continue
            f = st.spearmanr if method == "spearman" else st.pearsonr
            r, p = f(x[ok], y[ok])
            rows.append({"index": ic, "soil_variable": nn, "n": n,
                         "r": float(r), "p": float(p),
                         "r_squared": float(r) ** 2,
                         "circular": ic in NSH_COMPONENTS and nn in NSH_COMPONENTS})
    if not rows:
        return pd.DataFrame([{"error": "no pair had enough varying, matching data"}])
    out = pd.DataFrame(rows)
    out["p_BH"] = _bh(out["p"].values)
    out["strength"] = pd.cut(out["r"].abs(), [-0.01, .2, .4, .6, .8, 1.0],
                             labels=["negligible", "weak", "moderate", "strong",
                                     "very strong"])
    return out.sort_values("p_BH")


def _bh(p):
    p = np.asarray(p, float)
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


def matrix(corr: pd.DataFrame, value: str = "r") -> pd.DataFrame:
    if "error" in corr.columns:
        return pd.DataFrame()
    return corr.pivot(index="index", columns="soil_variable", values=value)


def heatmap(corr: pd.DataFrame, ax=None, annotate=True, alpha=0.05):
    """Correlation heatmap. Cells significant after BH correction are marked;
    unmarked cells should not be interpreted."""
    import matplotlib.pyplot as plt
    m = matrix(corr, "r")
    if m.empty:
        return None
    pm = matrix(corr, "p_BH")
    if ax is None:
        _, ax = plt.subplots(figsize=(1.1 * m.shape[1] + 3, 0.36 * m.shape[0] + 2))
    im = ax.imshow(m.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels(m.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels(m.index, fontsize=8)
    if annotate:
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = m.values[i, j]
                if pd.isna(v):
                    continue
                sig = pm.values[i, j] < alpha if not pd.isna(pm.values[i, j]) else False
                ax.text(j, i, f"{v:.2f}" + ("*" if sig else ""), ha="center",
                        va="center", fontsize=7,
                        color="white" if abs(v) > 0.55 else "black",
                        fontweight="bold" if sig else "normal")
    plt.colorbar(im, ax=ax, label="Spearman r", shrink=.8)
    ax.set_title("Nematode indices versus soil variables", loc="left",
                 fontweight="bold")
    return ax


def interpret(corr: pd.DataFrame, alpha: float = 0.05) -> dict:
    """Summarise, separating non-circular evidence from arithmetic."""
    if "error" in corr.columns:
        return {"error": corr["error"].iloc[0]}
    real = corr[~corr["circular"]]
    sig = real[real["p_BH"] < alpha]
    n = int(corr["n"].max()) if len(corr) else 0
    notes = []
    if n < 10:
        notes.append(f"Only {n} samples. Correlation estimates at this n are very "
                     "unstable — a single point can change r substantially. Treat "
                     "these as exploratory.")
    if int(corr["circular"].sum()):
        notes.append(f"{int(corr['circular'].sum())} pairs were flagged CIRCULAR "
                     "and excluded from the headline count: NSH is constructed "
                     "from MI, PPI, EI, SI, BI, CI and abundance, so correlations "
                     "among them are arithmetic, not evidence.")
    notes.append("Correlation is not causation, and a soil variable and a "
                 "nematode index can both respond to a third factor such as "
                 "management history.")
    return {"n_pairs": len(real), "n_significant": len(sig), "max_n": n,
            "top": (sig.head(8)[["index", "soil_variable", "n", "r", "p_BH",
                                 "strength"]] if len(sig) else pd.DataFrame()),
            "notes": notes}


TEMPLATE_NOTE = """\
Add a fourth sheet named `nutrients` to your workbook, or simply put extra
numeric columns in the `samples` sheet.

    Sample | pH  | organic_C | NO3_N | Olsen_P | ...
    S01    | 6.5 | 1.88      | 63    | 49      | ...
    S02    | 5.6 | 1.72      | 77    | 40      | ...

Sample names must match the counts sheet exactly. Any numeric column is
treated as a soil variable. Units are yours to record in the methods.
"""
