"""
nemacom_plots.py — figures and extended statistics for NemaCom.

Kept separate from the calculation engine so plotting choices can never change
a computed value. Import alongside nemacom_all.

Additions here are deliberately few. Du Preez et al. (2022) Soil Biol Biochem
169:108640 warn that index proliferation is a real problem in this field: every
extra index is another chance to report whichever one happened to come out
significant. Each item below either (a) is a standard visualisation that already
exists in the primary literature, or (b) directly serves the calibration gap
identified for the NSH index.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# --------------------------------------------------------------------------
# Palettes
# --------------------------------------------------------------------------
# The first three are colour-vision-deficiency safe. Roughly 8% of men of
# northern European descent, and a smaller proportion of other populations, have
# some form of red-green colour deficiency, so a CVD-safe default is not a
# nicety. Journals increasingly ask for it.
PALETTES = {
    "Okabe-Ito (CVD-safe)": ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                             "#0072B2", "#D55E00", "#CC79A7", "#000000"],
    "Set2 (CVD-safe)":      ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                             "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"],
    "Viridis (CVD-safe)":   ["#440154", "#414487", "#2A788E", "#22A884",
                             "#7AD151", "#FDE725", "#3B528B", "#5EC962"],
    "Dark2":                ["#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                             "#66A61E", "#E6AB02", "#A6761D", "#666666"],
    "Vivid":                ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                             "#FF7F00", "#FFFF33", "#A65628", "#F781BF"],
    "Earth":                ["#8C6D31", "#BD9E39", "#637939", "#8CA252",
                             "#843C39", "#AD494A", "#7B4173", "#A55194"],
    "Grayscale (print)":    ["#1A1A1A", "#4D4D4D", "#808080", "#A6A6A6",
                             "#CCCCCC", "#E6E6E6", "#333333", "#999999"],
}

TROPHIC_ORDER = ["PP", "BF", "FF", "OM", "PR"]


def palette(name, n=8):
    cols = PALETTES.get(name, PALETTES["Okabe-Ito (CVD-safe)"])
    return [cols[i % len(cols)] for i in range(n)]


def trophic_colors(name):
    p = palette(name, 5)
    return dict(zip(TROPHIC_ORDER, p))


def apply_style(font_size=10, spines=False):
    plt.rcParams.update({
        "font.size": font_size,
        "axes.spines.top": spines,
        "axes.spines.right": spines,
        "axes.grid": False,
        "figure.autolayout": True,
    })


# --------------------------------------------------------------------------
# Error bars
# --------------------------------------------------------------------------
def group_summary(values: pd.Series, groups: pd.Series, err="SE"):
    """Mean and error per group. err = 'SD', 'SE' or 'CI95'."""
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    out = df.groupby("g")["v"].agg(["mean", "std", "count"])
    out["sem"] = out["std"] / np.sqrt(out["count"])
    if err == "SD":
        out["err"] = out["std"]
    elif err == "CI95":
        # t-based, not 1.96 — matters at the small n typical of field trials
        tcrit = out["count"].apply(
            lambda n: st.t.ppf(0.975, n - 1) if n > 1 else np.nan)
        out["err"] = tcrit * out["sem"]
    else:
        out["err"] = out["sem"]
    out.attrs["err_type"] = err
    return out


def bar_with_error(values, groups, err="SE", pal="Okabe-Ito (CVD-safe)",
                   ylabel="", title="", letters=None, ax=None):
    s = group_summary(values, groups, err)
    cols = palette(pal, len(s))
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(s))
    ax.bar(x, s["mean"], yerr=s["err"], capsize=5, color=cols,
           edgecolor="black", linewidth=0.6,
           error_kw=dict(ecolor="black", lw=1.1))
    # overlay the raw points: a bar chart alone hides n and spread
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    for i, lev in enumerate(s.index):
        v = df.loc[df["g"] == lev, "v"].values
        ax.scatter(np.full(len(v), i) + np.random.uniform(-.09, .09, len(v)),
                   v, color="black", s=16, zorder=3, alpha=.75)
    if letters is not None:
        for i, lev in enumerate(s.index):
            if lev in letters:
                ax.text(i, s["mean"].iloc[i] + s["err"].iloc[i],
                        f"  {letters[lev]}", ha="center", va="bottom",
                        fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(s.index, rotation=0)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}  (error bars = {err})", loc="left", fontweight="bold")
    return ax


# --------------------------------------------------------------------------
# Post-hoc tests
# --------------------------------------------------------------------------
def tukey_posthoc(values, groups):
    """Tukey HSD pairwise comparisons. Returns a tidy p-value table."""
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    levs = list(df["g"].unique())
    arrays = [df.loc[df["g"] == l, "v"].values for l in levs]
    if len(arrays) < 2 or any(len(a) < 2 for a in arrays):
        return pd.DataFrame()
    res = st.tukey_hsd(*arrays)
    rows = []
    for i in range(len(levs)):
        for j in range(i + 1, len(levs)):
            rows.append({"group_1": levs[i], "group_2": levs[j],
                         "difference": res.statistic[i, j],
                         "p_value": res.pvalue[i, j]})
    return pd.DataFrame(rows)


def dunn_posthoc(values, groups, adjust="bonferroni"):
    """Dunn's test after Kruskal-Wallis, with Bonferroni or Holm adjustment."""
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    levs = list(df["g"].unique())
    if len(levs) < 2:
        return pd.DataFrame()
    df["r"] = st.rankdata(df["v"])
    N = len(df)
    g = df.groupby("g")["r"].agg(["mean", "count"])
    # tie correction
    _, counts = np.unique(df["v"], return_counts=True)
    ties = (counts ** 3 - counts).sum()
    sigma2 = (N * (N + 1) / 12) - ties / (12 * (N - 1))
    rows = []
    for i in range(len(levs)):
        for j in range(i + 1, len(levs)):
            a, b = levs[i], levs[j]
            se = np.sqrt(sigma2 * (1 / g.loc[a, "count"] + 1 / g.loc[b, "count"]))
            z = (g.loc[a, "mean"] - g.loc[b, "mean"]) / se if se else np.nan
            rows.append({"group_1": a, "group_2": b, "z": z,
                         "p_unadjusted": 2 * (1 - st.norm.cdf(abs(z)))})
    out = pd.DataFrame(rows)
    m = len(out)
    if adjust == "bonferroni":
        out["p_adjusted"] = np.minimum(out["p_unadjusted"] * m, 1.0)
    else:  # holm
        order = out["p_unadjusted"].rank(method="first").astype(int)
        adj = out["p_unadjusted"] * (m - order + 1)
        out["p_adjusted"] = np.minimum.accumulate(
            adj.sort_values(ascending=False))[out.index].clip(upper=1.0)
    return out


def compact_letters(posthoc, alpha=0.05, pcol="p_value"):
    """Compact letter display from a pairwise p-value table.

    Greedy assignment: adequate for the small designs typical of field trials.
    For complex designs verify against a dedicated package before publishing."""
    if posthoc.empty:
        return {}
    levs = sorted(set(posthoc["group_1"]) | set(posthoc["group_2"]))
    diff = {(r["group_1"], r["group_2"]): r[pcol] < alpha
            for _, r in posthoc.iterrows()}

    def differs(a, b):
        return diff.get((a, b), diff.get((b, a), False))

    letters = {l: set() for l in levs}
    groups_of_letter = []
    for lev in levs:
        placed = False
        for k, members in enumerate(groups_of_letter):
            if not any(differs(lev, m) for m in members):
                members.append(lev)
                letters[lev].add(k)
                placed = True
        if not placed:
            groups_of_letter.append([lev])
            letters[lev].add(len(groups_of_letter) - 1)
    return {l: "".join(sorted("abcdefgh"[k] for k in v)) for l, v in letters.items()}


# --------------------------------------------------------------------------
# Ferris (2010) metabolic footprint faunal profile
# --------------------------------------------------------------------------
def footprint_profile(faunal_tbl, fp_tbl, groups=None, k=None,
                      pal="Okabe-Ito (CVD-safe)", ax=None):
    """The signature Ferris (2010) figure: rhomboid footprints centred on the
    EI/SI intersection.

    Per Ferris (2010) sec. 2.1.3, x-coordinates are SI +/- 0.5*Fs/k and
    y-coordinates EI +/- 0.5*Fe/k, where k is a scalar chosen for legibility and
    held constant across every footprint on one figure. Footprints are comparable
    within a figure, not between figures drawn with different k."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 6))
    Fe = fp_tbl["enrichment_footprint"]
    Fs = fp_tbl["structure_footprint"]
    if k is None:
        k = max(float(np.nanmax([Fe.max(), Fs.max()])) / 40.0, 1e-9)

    levs = list(pd.unique(groups)) if groups is not None else ["all"]
    cmap = dict(zip(levs, palette(pal, len(levs))))
    seen = set()
    for s in faunal_tbl.index:
        if s not in fp_tbl.index:
            continue
        EI, SI = faunal_tbl.loc[s, "EI"], faunal_tbl.loc[s, "SI"]
        fe, fs = Fe.get(s, 0) / k, Fs.get(s, 0) / k
        lev = groups.get(s, "all") if groups is not None else "all"
        pts = [(SI - .5 * fs, EI), (SI, EI + .5 * fe),
               (SI + .5 * fs, EI), (SI, EI - .5 * fe)]
        ax.add_patch(Polygon(pts, closed=True, facecolor=cmap[lev],
                             alpha=.35, edgecolor=cmap[lev], lw=1.3,
                             label=lev if lev not in seen else None))
        seen.add(lev)
        ax.plot(SI, EI, "o", color=cmap[lev], ms=4)

    ax.axhline(50, ls=":", c="grey", lw=1)
    ax.axvline(50, ls=":", c="grey", lw=1)
    for x, y, t in [(6, 94, "A Maturing"), (74, 94, "B Structured"),
                    (6, 3, "D Degraded"), (74, 3, "C Disturbed")]:
        ax.text(x, y, t, color="grey", fontsize=8)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Structure index (SI)")
    ax.set_ylabel("Enrichment index (EI)")
    ax.set_title("Metabolic footprint faunal profile (Ferris 2010)",
                 loc="left", fontweight="bold")
    if groups is not None:
        ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=len(levs))
    ax.text(.99, .01, f"scalar k = {k:.3g}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7, color="grey")
    return ax


# --------------------------------------------------------------------------
# NSH calibration explorer
# --------------------------------------------------------------------------
# This is the piece aimed squarely at the identified research gap. Ghaderi et al.
# (2025) state the NSH index "requires further calibration, testing and
# standardisation on more nematode community datasets", and their abundance bands
# derive from global biome ranges (van den Hoogen et al. 2019) in which tropical
# cropland is thinly represented. Recalibrating for a region means changing band
# edges and seeing what happens to discrimination between treatments.
def nsh_discrimination(nsh_tbl, groups):
    """How well does NSH separate the groups? Returns effect size and spread.

    eta_squared is the proportion of NSH variance explained by group. A low
    value with a significant p is the signature of an index that is
    'working' statistically but not discriminating usefully."""
    from scipy import stats as st
    df = pd.DataFrame({"n": nsh_tbl["NSH"],
                       "g": groups.reindex(nsh_tbl.index)}).dropna()
    levs = list(df["g"].unique())
    arrays = [df.loc[df["g"] == l, "n"].values for l in levs]
    arrays = [a for a in arrays if len(a) > 1]
    if len(arrays) < 2:
        return {}
    F, p = st.f_oneway(*arrays)
    grand = df["n"].mean()
    ss_b = sum(len(a) * (a.mean() - grand) ** 2 for a in arrays)
    ss_t = ((df["n"] - grand) ** 2).sum()
    return {
        "n_groups": len(arrays),
        "NSH_min": float(df["n"].min()),
        "NSH_max": float(df["n"].max()),
        "NSH_range_used": float(df["n"].max() - df["n"].min()),
        "range_available": 24.0,           # the index runs 8-32
        "pct_of_scale_used": float((df["n"].max() - df["n"].min()) / 24 * 100),
        "anova_F": float(F),
        "anova_p": float(p),
        "eta_squared": float(ss_b / ss_t) if ss_t else np.nan,
        "n_categories_occupied": int(nsh_tbl["interpretation"].nunique()),
    }


def subscore_contributions(nsh_tbl):
    """Which subscores actually vary? A subscore that is constant across every
    sample contributes nothing to discrimination and is a candidate for
    recalibration in a new region."""
    cols = [c for c in nsh_tbl.columns if c.endswith("_score")]
    out = pd.DataFrame({
        "mean": nsh_tbl[cols].mean(),
        "sd": nsh_tbl[cols].std(),
        "min": nsh_tbl[cols].min(),
        "max": nsh_tbl[cols].max(),
    })
    out["range"] = out["max"] - out["min"]
    out["varies"] = out["range"] > 0
    return out.sort_values("sd", ascending=False)


def nsh_band_sweep(faunal_tbl, abundance, groups, nsh_func,
                   index_name="EI", shifts=(-15, -10, -5, 0, 5, 10, 15)):
    """Shift one index's band edges and watch discrimination change.

    A crude but honest first pass at recalibration: if shifting a band edge
    markedly improves separation between treatments, that band is mis-set for
    the system being studied. Report this as exploratory - shifting edges to
    maximise a p-value is circular if you then test on the same data."""
    import copy
    import nemacom_all as nc
    rows = []
    for sh in shifts:
        bands = copy.deepcopy(nc.NSH_BANDS)
        if index_name in bands:
            bands[index_name] = [(lo + sh if np.isfinite(lo) else lo,
                                  hi + sh if np.isfinite(hi) else hi, sc)
                                 for lo, hi, sc in bands[index_name]]
        saved = nc.NSH_BANDS
        try:
            nc.NSH_BANDS = bands
            tbl = nsh_func(faunal_tbl, abundance)
            d = nsh_discrimination(tbl, groups)
        finally:
            nc.NSH_BANDS = saved
        rows.append({"shift": sh, **{k: d.get(k) for k in
                                     ("eta_squared", "anova_p",
                                      "NSH_range_used")}})
    return pd.DataFrame(rows).set_index("shift")
