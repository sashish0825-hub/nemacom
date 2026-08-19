"""nemonema_plots.py — figures and post-hoc tests for NEMO-NEMA."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# The first three are safe for colour-vision deficiency, which affects roughly
# 8% of men of northern European descent. Several journals now require it.
PALETTES = {
    "Okabe-Ito (CVD-safe)": ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                             "#0072B2", "#D55E00", "#CC79A7", "#000000"],
    "Set2 (CVD-safe)":      ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                             "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"],
    "Viridis (CVD-safe)":   ["#440154", "#414487", "#2A788E", "#22A884",
                             "#7AD151", "#FDE725", "#3B528B", "#5EC962"],
    "Dark2":  ["#1B9E77", "#D95F02", "#7570B3", "#E7298A",
               "#66A61E", "#E6AB02", "#A6761D", "#666666"],
    "Vivid":  ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
               "#FF7F00", "#FFC300", "#A65628", "#F781BF"],
    "Earth":  ["#8C6D31", "#BD9E39", "#637939", "#8CA252",
               "#843C39", "#AD494A", "#7B4173", "#A55194"],
    "Grayscale (print)": ["#1A1A1A", "#4D4D4D", "#808080", "#A6A6A6",
                          "#CCCCCC", "#E6E6E6", "#333333", "#999999"],
}
TROPHIC_ORDER = ["PP", "BF", "FF", "OM", "PR"]


def palette(name, n=8):
    cols = PALETTES.get(name, PALETTES["Okabe-Ito (CVD-safe)"])
    return [cols[i % len(cols)] for i in range(n)]


def trophic_colors(name):
    return dict(zip(TROPHIC_ORDER, palette(name, 5)))


def style(font_size=10):
    plt.rcParams.update({"font.size": font_size, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.autolayout": True})


# --------------------------------------------------------------------------
def group_stats(values, groups, err="SE"):
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    o = df.groupby("g")["v"].agg(["mean", "std", "count"])
    o["sem"] = o["std"] / np.sqrt(o["count"])
    if err == "SD":
        o["err"] = o["std"]
    elif err == "CI95":
        # t-based, not 1.96 — at n = 3 the normal approximation is far too narrow
        o["err"] = o["count"].apply(
            lambda n: st.t.ppf(0.975, n - 1) if n > 1 else np.nan) * o["sem"]
    else:
        o["err"] = o["sem"]
    return o


def bar_with_error(values, groups, err="SE", pal="Okabe-Ito (CVD-safe)",
                   ylabel="", title="", letters=None, ax=None):
    """Bars are group means; raw points are overlaid so n and spread stay visible.

    Letter placement is computed from the highest visible element in each column
    (bar top, error bar, or the highest data point), so significance letters
    never sit on top of the dots.
    """
    s = group_stats(values, groups, err)
    cols = palette(pal, len(s))
    if ax is None:
        _, ax = plt.subplots(figsize=(max(6.4, 1.15 * len(s) + 2), 4.6))
    x = np.arange(len(s))
    ax.bar(x, s["mean"], yerr=s["err"], capsize=5, color=cols,
           edgecolor="black", linewidth=0.6,
           error_kw=dict(ecolor="black", lw=1.1), width=0.72)

    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    rng = np.random.default_rng(0)
    tops = []
    for i, lev in enumerate(s.index):
        v = df.loc[df["g"] == lev, "v"].values
        ax.scatter(np.full(len(v), i) + rng.uniform(-.075, .075, len(v)), v,
                   color="black", s=15, zorder=3, alpha=.8)
        e = s["err"].iloc[i] if s["err"].iloc[i] == s["err"].iloc[i] else 0
        tops.append(max(s["mean"].iloc[i] + e, v.max() if len(v) else 0))

    span = (max(tops) - min(0, df["v"].min())) or 1.0
    if letters:
        for i, lev in enumerate(s.index):
            if lev in letters:
                ax.text(i, tops[i] + span * 0.045, letters[lev], ha="center",
                        va="bottom", fontweight="bold",
                        fontsize=plt.rcParams["font.size"] + 1)
    ax.set_ylim(top=max(tops) + span * 0.16)

    # long category names collide when set horizontally; rotate only when needed
    labs = [f"{l}\n(n={int(s.loc[l, 'count'])})" for l in s.index]
    longest = max(len(str(l)) for l in s.index)
    rot = 0 if longest <= 8 and len(s) <= 6 else (30 if longest <= 14 else 45)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=rot,
                       ha="right" if rot else "center",
                       rotation_mode="anchor" if rot else None)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}   error bars = {err}", loc="left", fontweight="bold")
    return ax


def faunal_profile(faunal_tbl, groups=None, pal="Okabe-Ito (CVD-safe)", ax=None):
    """Ferris et al. (2001) faunal profile.

    Quadrat labels verified against Fig. 1 and Table 5 of that paper, which
    describes conventional annual cropping as 'D-A-B: Degraded-disturbed-
    structured' and low-input perennial cropping as 'C: Stable'.
        A upper-left  Disturbed     B upper-right Structured
        D lower-left  Degraded      C lower-right Stable
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5.8))
    levs = list(pd.unique(groups)) if groups is not None else ["all"]
    cols = palette(pal, len(levs))
    if groups is None:
        ax.scatter(faunal_tbl["SI"], faunal_tbl["EI"], s=80, color=cols[0],
                   edgecolor="w", zorder=3)
    else:
        for i, lev in enumerate(levs):
            m = groups.reindex(faunal_tbl.index) == lev
            ax.scatter(faunal_tbl.loc[m, "SI"], faunal_tbl.loc[m, "EI"], s=80,
                       color=cols[i], label=str(lev), edgecolor="w", zorder=3)
        ax.legend(frameon=False, fontsize=8, loc="upper center",
                  ncol=min(len(levs), 5), bbox_to_anchor=(0.5, -0.13))
    ax.axhline(50, ls=":", c="grey", lw=1)
    ax.axvline(50, ls=":", c="grey", lw=1)
    for x, y, t in [(2, 95, "A  Disturbed"), (66, 95, "B  Structured"),
                    (2, 2, "D  Degraded"), (66, 2, "C  Stable")]:
        ax.text(x, y, t, color="#6B6B6B", fontsize=8.5, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Structure index (SI)")
    ax.set_ylabel("Enrichment index (EI)")
    ax.set_title("Faunal profile (Ferris et al. 2001)", loc="left",
                 fontweight="bold")
    return ax


def footprint_profile(faunal_tbl, fp_tbl, groups=None, k=None,
                      pal="Okabe-Ito (CVD-safe)", ax=None):
    """Ferris (2010) rhomboid footprints centred on the EI/SI intersection.

    Per sec. 2.1.3 the x-coordinates are SI ± 0.5 Fs/k and the y-coordinates
    EI ± 0.5 Fe/k, where k is a legibility scalar held constant across the
    figure. Footprints are comparable within a figure, not between figures."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5.8))
    Fe, Fs = fp_tbl["enrichment_footprint"], fp_tbl["structure_footprint"]
    if k is None:
        k = max(float(np.nanmax([Fe.max(), Fs.max()])) / 40.0, 1e-9)
    levs = list(pd.unique(groups)) if groups is not None else ["all"]
    cmap = dict(zip(levs, palette(pal, len(levs))))
    seen = set()
    for s in faunal_tbl.index:
        if s not in fp_tbl.index:
            continue
        EI, SI = faunal_tbl.loc[s, "EI"], faunal_tbl.loc[s, "SI"]
        if pd.isna(EI) or pd.isna(SI):
            continue
        fe, fs = Fe.get(s, 0) / k, Fs.get(s, 0) / k
        lev = groups.get(s, "all") if groups is not None else "all"
        ax.add_patch(Polygon([(SI - .5 * fs, EI), (SI, EI + .5 * fe),
                              (SI + .5 * fs, EI), (SI, EI - .5 * fe)],
                             closed=True, facecolor=cmap[lev], alpha=.35,
                             edgecolor=cmap[lev], lw=1.3,
                             label=str(lev) if lev not in seen else None))
        seen.add(lev)
        ax.plot(SI, EI, "o", color=cmap[lev], ms=4)
    ax.axhline(50, ls=":", c="grey")
    ax.axvline(50, ls=":", c="grey")
    for x, y, t in [(2, 95, "A  Disturbed"), (66, 95, "B  Structured"),
                    (2, 2, "D  Degraded"), (66, 2, "C  Stable")]:
        ax.text(x, y, t, color="#6B6B6B", fontsize=8.5, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Structure index (SI)")
    ax.set_ylabel("Enrichment index (EI)")
    ax.set_title("Metabolic footprint profile (Ferris 2010)", loc="left",
                 fontweight="bold")
    if groups is not None:
        ax.legend(frameon=False, fontsize=8, loc="upper center",
                  ncol=min(len(levs), 5), bbox_to_anchor=(0.5, -0.13))
    ax.text(.99, .01, f"scalar k = {k:.3g}", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=7, color="grey")
    return ax


LONG = {"PP": "plant parasite", "BF": "bacterivore", "FF": "fungivore",
        "OM": "omnivore", "PR": "predator"}


def trophic_bars(faunal_tbl, pal="Okabe-Ito (CVD-safe)", ax=None):
    """Stacked trophic composition.

    The legend is placed BELOW the axes. Placed inside, it sits on top of the
    bars and hides data — which defeats the purpose of the figure."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    cols = trophic_colors(pal)
    comp = faunal_tbl[[c for c in faunal_tbl.columns if c.startswith("pct_")]]
    bot = np.zeros(len(comp))
    for c in comp.columns:
        k = c.split("_")[1]
        ax.bar(range(len(comp)), comp[c], bottom=bot, color=cols[k],
               label=f"{k} — {LONG.get(k, k)}", width=0.82)
        bot += comp[c].fillna(0).values
    ax.set_xticks(range(len(comp)))
    ax.set_xticklabels(comp.index, rotation=90, fontsize=7)
    ax.set_ylabel("% of nematode community")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.6, len(comp) - 0.4)
    ax.set_title("Trophic composition", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8, ncol=5, loc="upper center",
              bbox_to_anchor=(0.5, -0.32), handlelength=1.2, columnspacing=1.4)
    return ax


# --------------------------------------------------------------------------
def tukey(values, groups):
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    levs = list(pd.unique(df["g"]))
    arr = [df.loc[df["g"] == l, "v"].values for l in levs]
    if len(arr) < 2 or any(len(a) < 2 for a in arr):
        return pd.DataFrame()
    r = st.tukey_hsd(*arr)
    return pd.DataFrame([{"group_1": levs[i], "group_2": levs[j],
                          "difference": r.statistic[i, j],
                          "p_value": r.pvalue[i, j]}
                         for i in range(len(levs)) for j in range(i + 1, len(levs))])


def dunn(values, groups, adjust="bonferroni"):
    """Dunn's test after Kruskal-Wallis, with tie correction."""
    from scipy import stats as st
    df = pd.DataFrame({"v": values, "g": groups.reindex(values.index)}).dropna()
    levs = list(pd.unique(df["g"]))
    if len(levs) < 2:
        return pd.DataFrame()
    df["r"] = st.rankdata(df["v"])
    N = len(df)
    g = df.groupby("g")["r"].agg(["mean", "count"])
    _, cnt = np.unique(df["v"], return_counts=True)
    sigma2 = (N * (N + 1) / 12) - (cnt ** 3 - cnt).sum() / (12 * (N - 1))
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
    out["p_adjusted"] = np.minimum(out["p_unadjusted"] * m, 1.0) \
        if adjust == "bonferroni" else out["p_unadjusted"]
    return out


def letters(posthoc, alpha=0.05, pcol="p_value"):
    """Compact letter display. Greedy assignment — adequate for the small
    designs typical of field trials; verify against a dedicated package for
    complex ones."""
    if posthoc.empty:
        return {}
    levs = sorted(set(posthoc["group_1"]) | set(posthoc["group_2"]))
    diff = {(r["group_1"], r["group_2"]): r[pcol] < alpha
            for _, r in posthoc.iterrows()}

    def differs(a, b):
        return diff.get((a, b), diff.get((b, a), False))

    assigned = {l: set() for l in levs}
    groups_of = []
    for lev in levs:
        placed = False
        for k, members in enumerate(groups_of):
            if not any(differs(lev, m) for m in members):
                members.append(lev)
                assigned[lev].add(k)
                placed = True
        if not placed:
            groups_of.append([lev])
            assigned[lev].add(len(groups_of) - 1)
    return {l: "".join(sorted("abcdefgh"[k] for k in v))
            for l, v in assigned.items()}
