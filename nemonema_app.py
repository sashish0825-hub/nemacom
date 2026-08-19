"""
NEMO-NEMA — Nematode Ecological Metrics and Ordination for
             Nematode Ecosystem Monitoring and Assessment

Run:  streamlit run nemonema_app.py

Developed by
  Ashish Kumar Singh, Kavita Jain, Vishal Singh Somvanshi, Rashid Pervez,
  Anil Sirohi, Pankaj
  Division of Nematology, ICAR-Indian Agricultural Research Institute,
  New Delhi 110012
"""

import io
import os

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

import nemonema_core as nc
import nemonema_plots as npl
import nemonema_qc as qc
import nemonema_summary as nsum
import nemonema_autofill as af
import nemonema_report as rpt
import nemonema_nutrients as nut
import nemonema_ordination as ordn

st.set_page_config(page_title="NEMO-NEMA", layout="wide",
                   page_icon="logo_mark.png" if os.path.exists("logo_mark.png")
                   else "\U0001F52C")

st.markdown("""<style>
  .block-container {padding-top: 2rem; max-width: 1500px;}
  .stTabs [data-baseweb="tab-list"] {gap:2px; border-bottom:1px solid #DED6C6;}
  .stTabs [data-baseweb="tab"] {height:42px; padding:0 14px; background:transparent;
     border-radius:6px 6px 0 0; font-size:0.9rem;}
  .stTabs [aria-selected="true"] {background:#F1ECE1; border-bottom:2px solid #2E7D5B;
     color:#17453A;}
  .nn-rule {height:3px; width:74px; background:linear-gradient(90deg,#17453A,#8CBB72);
     border-radius:2px; margin:2px 0 12px 0;}
  .nn-sub {color:#5B6B62; font-size:0.9rem; line-height:1.5;}
  .nn-cp {display:inline-block; width:9px; height:9px; border-radius:50%;
     margin-right:5px; vertical-align:middle;}
</style>""", unsafe_allow_html=True)

_h = st.columns([1, 9])
with _h[0]:
    if os.path.exists("icar_logo.png"):
        st.image("icar_logo.png", width=88)
with _h[1]:
    st.markdown(
        "<h1 style='margin-bottom:0'>NEMO-NEMA</h1><div class='nn-rule'></div>"
        "<div class='nn-sub'><b>N</b>ematode <b>E</b>cological <b>M</b>etrics and "
        "<b>O</b>rdination for <b>N</b>ematode <b>E</b>cosystem <b>M</b>onitoring "
        "and <b>A</b>ssessment<br><span style='color:#8A9690'>"
        + "".join(f"<span class='nn-cp' style='background:{c}'></span>"
                  for c in ["#D55E00", "#E69F00", "#F0E442", "#56B4E9", "#17453A"])
        + f"c-p 1 to 5 &nbsp;&middot;&nbsp; {nc.VERSION}</span></div>",
        unsafe_allow_html=True)
st.write("")


def footer():
    st.divider()
    st.markdown(
        "<div style='font-size:0.86rem; line-height:1.6; color:#3E4A44'>"
        "<b>Developed by</b><br>Ashish Kumar Singh &middot; Kavita Jain &middot; "
        "Vishal Singh Somvanshi &middot; Rashid Pervez &middot; Anil Sirohi "
        "&middot; Pankaj<br><span style='color:#5B6B62'>Division of Nematology, "
        "ICAR-Indian Agricultural Research Institute, New Delhi 110012"
        "</span></div>", unsafe_allow_html=True)


def fig_bytes(fig, dpi=300):
    """bbox_inches='tight' is essential: legends now sit below the axes and
    would otherwise be cropped out of the exported file."""
    b = io.BytesIO()
    fig.savefig(b, format="png", dpi=dpi, bbox_inches="tight")
    return b.getvalue()


def show(fig_fn, table, key, caption="", dpi=300, default="Figure"):
    """Every analysis can be shown as a figure or a table. The toggle is per
    analysis, so a user can take the table for a supplement and the figure for
    the paper without re-running anything."""
    mode = st.radio("View as", ["Figure", "Table"], horizontal=True,
                    key=key, index=0 if default == "Figure" else 1,
                    label_visibility="collapsed")
    if mode == "Table":
        st.dataframe(table.round(3) if hasattr(table, "round") else table,
                     use_container_width=True)
        st.download_button("Download table (CSV)", table.to_csv().encode(),
                           f"{key}.csv", "text/csv", key=f"csv_{key}")
    else:
        fig = fig_fn()
        if fig is None:
            st.info("No figure available for this analysis.")
            return
        st.pyplot(fig, use_container_width=False)
        if caption:
            st.caption(caption)
        st.download_button("Download figure (PNG)", fig_bytes(fig, dpi),
                           f"{key}.png", "image/png", key=f"png_{key}")
        plt.close(fig)


# ==========================================================================
with st.sidebar:
    st.header("1. Your data")
    up = st.file_uploader("Excel workbook (.xlsx)", type=["xlsx"])
    st.caption("Three sheets: counts, taxa, samples. "
               "Only `counts` is essential — the Auto-assign tab can propose "
               "trophic group and c-p from taxon names.")

    st.header("2. Figure style")
    pal = st.selectbox("Colour palette", list(npl.PALETTES))
    st.markdown(" ".join(
        f'<span style="display:inline-block;width:22px;height:14px;background:{c};'
        f'border:1px solid #999;margin-right:2px"></span>'
        for c in npl.palette(pal, 8)), unsafe_allow_html=True)
    errtype = st.radio("Error bars", ["SE", "SD", "CI95"], horizontal=True)
    fsize = st.slider("Font size", 8, 20, 11)
    dpi = st.select_slider("Figure resolution (dpi)", [100, 150, 300, 600], 300)

    st.header("3. Options")
    pvchoice = st.radio("Prominence value", ["D x sqrt(F)", "D x sqrt(F) / 10"],
                        help="Both variants are attributed to Norton (1978). "
                             "State which you used in your methods.")
    basis = st.selectbox("Counts are expressed per",
                         ["100 g dry soil", "200 cc soil", "250 g soil", "other"],
                         help="The NSH abundance band assumes per 100 g dry soil.")

if up is None:
    c1, c2 = st.columns([3, 2])
    with c2:
        if os.path.exists("micrograph.jpg"):
            st.image("micrograph.jpg",
                     caption="Soil nematode, 4x. Division of Nematology, ICAR-IARI.")
    with c1:
        st.info("Upload a workbook to begin.")
        st.markdown("""
**Minimum input — one sheet**

| Sheet | Contents |
|---|---|
| `counts` | taxa in rows, samples in columns, individuals recovered |

**Recommended — three sheets**

| Sheet | Contents |
|---|---|
| `taxa` | `trophic` (PP/BF/FF/OM/PR), `cp` (1–5), `source` |
| `samples` | `group` = the factor being compared; `field`/`plot` if available |

If the `taxa` sheet is missing or incomplete, the **Auto-assign** tab proposes
trophic group and c-p from the taxon names, with the source and a confidence
flag on every row. Nothing is applied without your confirmation.
""")
    footer()
    st.stop()

npl.style(fsize)
try:
    counts, taxa, samples = nc.load(up)
except Exception as exc:
    st.error(f"Could not read that workbook: {exc}")
    footer()
    st.stop()

if "taxa_edit" in st.session_state and st.session_state.get("taxa_src") == up.name:
    taxa = st.session_state["taxa_edit"]

problems = nc.validate(counts, taxa, samples)
divisor = 10.0 if pvchoice.endswith("/ 10") else 1.0
grp = samples["group"] if samples is not None and "group" in samples else None

st.success(f"{counts.shape[0]} taxa · {counts.shape[1]} samples · "
           f"{int(counts.fillna(0).values.sum()):,} individuals")

TABS = ["Auto-assign", "Data check", "Community", "Diversity", "Faunal analysis",
        "Soil health", "Footprints", "Summary", "Statistics", "Multivariate",
        "Soil nutrients", "Reference", "Validation", "Report", "Cite",
        "Manual"]
tabs = st.tabs(TABS)

# ---------------------------------------------------------------- auto-assign
with tabs[0]:
    st.subheader("Propose trophic group and c-p from taxon names")
    st.caption("Assignments come from Yeates et al. (1993) and Bongers & Bongers "
               "(1998) as tabulated in Ferris (2010) Table 1. They are at FAMILY "
               "level, so every proposal needs your confirmation: a fuzzy match "
               "on a misspelling, applied silently, becomes a wrong maturity "
               "index that looks entirely normal.")
    prop = af.propose(list(counts.index))
    s = af.summarise(prop)
    a, b, c_, d_ = st.columns(4)
    a.metric("Taxa", s["n_taxa"])
    b.metric("Exact match", s["exact"])
    c_.metric("Approximate", s["fuzzy"])
    d_.metric("No match", s["unmatched"])
    (st.success if s["ready"] and not s["needs_attention"] else st.warning)(s["message"])
    st.dataframe(prop[["taxon", "matched_genus", "match_type", "family", "trophic",
                       "trophic_name", "cp", "confidence", "action"]],
                 use_container_width=True)
    if st.button("Apply proposals to blank cells only"):
        st.session_state["taxa_edit"] = af.apply_proposal(
            taxa.reindex(counts.index), prop, overwrite=False)
        st.session_state["taxa_src"] = up.name
        st.success("Applied. Values you had already entered were not overwritten. "
                   "Re-check the other tabs.")
        st.rerun()
    st.download_button("Download proposals (CSV)", prop.to_csv(index=False).encode(),
                       "trait_proposals.csv", "text/csv")

if problems:
    with tabs[1]:
        st.error("Input problems — nothing has been computed:")
        for p in problems:
            st.write("•", p)
        st.info("If `trophic` or `cp` are missing, use the Auto-assign tab.")
    footer()
    st.stop()

# ---- computations
nrt = nc.norton(counts, divisor).join(taxa[["trophic", "cp"]])
dv = nc.diversity(counts)
fn = nc.faunal(counts, taxa)
ex = nc.extra_indices(counts, taxa)
fn_q = fn.assign(quadrant=[nc.quadrant(fn.EI[i], fn.SI[i]) for i in fn.index])
nsh = nc.nsh_index(fn[["MI", "PPI", "EI", "SI", "BI", "CI"]], counts.sum(axis=0))
fps = None
if "length_um" in taxa.columns and taxa["length_um"].notna().any():
    try:
        fps = nc.footprints(counts, taxa)
    except Exception as e:
        st.sidebar.warning(f"Footprints unavailable: {e}")

# ---------------------------------------------------------------- data check
with tabs[1]:
    st.subheader("Read this before interpreting anything")
    if samples is not None:
        rep = qc.replication(samples)
        (st.success if rep.get("ok") else st.error)(rep["message"])
        blk = qc.detect_blocking(samples)
        (st.error if blk.get("blocked") else st.success)(
            blk.get("message", "No blocking detected."))
    deg = qc.degeneracy(counts, taxa)
    nbad = int((~deg["informative"]).sum())
    if nbad:
        st.warning(f"{nbad} index/indices carry no information for this dataset. "
                   "They still appear in the tables — do not report them as findings.")
    st.dataframe(deg, use_container_width=True)
    dom = qc.dominance(counts)
    nx = int((dom["flag"] == "extreme dominance").sum())
    if nx:
        st.warning(f"{nx} of {len(dom)} samples exceed 70% dominance by one taxon. "
                   "Diversity indices then largely report that taxon's share.")
    st.dataframe(dom.round(1), use_container_width=True)
    if "source" not in taxa.columns or taxa["source"].isna().all():
        st.error("No `source` recorded for trophic and c-p. Every maturity and "
                 "food-web index is built from them; without sources these "
                 "results are not reproducible.")

# ---------------------------------------------------------------- community
with tabs[2]:
    st.subheader("Community structure — Norton (1978)")
    st.caption(f"Prominence value computed as **{pvchoice}**.")
    if nrt.attrs.get("pv_uninformative"):
        st.warning("Every taxon occurs in every sample, so absolute frequency is "
                   "100 throughout and prominence value = density × 10 exactly. "
                   "It adds nothing beyond mean density here.")

    def _f():
        d = nrt.sort_values("prominence_value").head(30)
        fig, ax = plt.subplots(figsize=(7.5, max(3, .3 * len(d) + 1)))
        cols = npl.trophic_colors(pal)
        ax.barh(range(len(d)), d["prominence_value"],
                color=[cols.get(t, "#999") for t in d["trophic"]])
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels(d.index, fontsize=8)
        ax.set_xlabel("Prominence value")
        ax.set_title("Ranked prominence value", loc="left", fontweight="bold")
        return fig
    show(_f, nrt, "community", "Bars coloured by trophic group.", dpi)

# ---------------------------------------------------------------- diversity
with tabs[3]:
    st.subheader("Diversity")
    m = st.selectbox("Metric", ["shannon_H", "simpson_1minusD", "pielou_J",
                                "margalef_d", "richness_S", "total_N"])

    def _f():
        if grp is None:
            fig, ax = plt.subplots(figsize=(7.5, 4))
            ax.bar(range(len(dv)), dv[m], color=npl.palette(pal, 1)[0])
            ax.set_xticks(range(len(dv)))
            ax.set_xticklabels(dv.index, rotation=90, fontsize=7)
            ax.set_ylabel(m)
            return fig
        fig, ax = plt.subplots(figsize=(7, 4.4))
        tk = npl.tukey(dv[m], grp)
        npl.bar_with_error(dv[m], grp, err=errtype, pal=pal, ylabel=m, title=m,
                           letters=npl.letters(tk) if not tk.empty else None, ax=ax)
        return fig
    show(_f, dv, "diversity",
         f"Bars are group means; error bars {errtype}; black points are individual "
         "samples. Bars sharing a letter do not differ (Tukey HSD, p > 0.05).", dpi)

# ---------------------------------------------------------------- faunal
with tabs[4]:
    st.subheader("Maturity and food-web indices")
    st.caption("MI, MI(2–5), σMI and PPI after Bongers (1990) and Bongers & "
               "Bongers (1998); EI, SI, CI and BI after Ferris et al. (2001); "
               "trophic groups after Yeates et al. (1993).")
    show(lambda: npl.faunal_profile(fn, grp, pal, plt.subplots(figsize=(6, 5.8))[1]).figure,
         fn_q, "faunal",
         "Each point is one sample, positioned by structure index (x) and "
         "enrichment index (y). Dashed lines mark the 50/50 split into the four "
         "quadrats of Ferris et al. (2001): A disturbed and enriched, "
         "B maturing to structured, C structured and stable, D degraded.", dpi)
    st.markdown("**Trophic composition**")
    show(lambda: npl.trophic_bars(fn, pal, plt.subplots(figsize=(9, 5))[1]).figure,
         fn[[c for c in fn.columns if c.startswith("pct_")]], "trophic",
         "Percentage of the nematode community in each trophic group.", dpi)
    st.markdown("**Additional indices**")
    st.dataframe(ex.round(3), use_container_width=True)

# ---------------------------------------------------------------- NSH
with tabs[5]:
    st.subheader("Nematode Soil Health index")
    st.caption("Ghaderi et al. (2025) Eur J Soil Sci 76:e70149. Seven subscores "
               "summed; range 8–32. Below 15 degraded, 15–24 moderate, "
               "above 25 well-functioning.")
    if basis != "100 g dry soil":
        st.warning(f"Your counts are per {basis}. The abundance subscore assumes "
                   "per 100 g dry soil, so that subscore — and every total — is "
                   "on the wrong scale until converted.")

    def _f():
        if grp is None:
            return None
        fig, ax = plt.subplots(figsize=(7, 4.4))
        tk = npl.tukey(nsh["NSH"], grp)
        npl.bar_with_error(nsh["NSH"], grp, err=errtype, pal=pal,
                           ylabel="NSH index", title="Soil health",
                           letters=npl.letters(tk) if not tk.empty else None, ax=ax)
        ax.axhline(15, ls="--", c="#C44E52", lw=1)
        ax.axhline(25, ls="--", c="#55A868", lw=1)
        return fig
    show(_f, nsh, "nsh", "Dashed lines mark the degraded (15) and "
         "well-functioning (25) thresholds.", dpi)

    st.markdown("**Calibration check**")
    dg = qc.nsh_diagnostics(nsh)
    nconst = int((dg["distinct"] == 1).sum())
    a, b = st.columns(2)
    a.metric("Scale used", f"{dg.attrs.get('scale_pct', float('nan')):.0f}%")
    b.metric("Constant subscores", f"{nconst} of 7")
    if nconst:
        st.warning(f"{nconst} subscores are identical in every sample. They add a "
                   "fixed amount to every NSH value and discriminate nothing. "
                   "This is worth reporting.")
    st.dataframe(dg, use_container_width=True)

# ---------------------------------------------------------------- footprints
with tabs[6]:
    st.subheader("Biomass and metabolic footprints")
    st.caption("Andrássy (1956) biomass; Ferris (2010) footprints.")
    if fps is None:
        st.info("Add `length_um` and `diameter_um` to the taxa sheet to compute "
                "these. Measure your own specimens where you can — published "
                "weights derive from adults, so juvenile-rich communities are "
                "overestimated.")
    else:
        miss = fps.attrs.get("missing") or []
        if miss:
            st.warning(f"Excluded for lack of measurements: {', '.join(miss)}")
        show(lambda: npl.footprint_profile(
                fn, fps, grp, None, pal, plt.subplots(figsize=(6, 5.8))[1]).figure,
             fps, "footprints",
             "Rhomboid width is the structure footprint, height the enrichment "
             "footprint, centred on each sample's SI/EI position (Ferris 2010). "
             "Comparable within this figure only.", dpi)

# ---------------------------------------------------------------- summary
with tabs[7]:
    st.subheader("Group means and dispersion")
    if grp is None:
        st.info("Add a `group` column to the samples sheet.")
    else:
        combined = dv.join(fn).join(nsh[["NSH"]])
        st.markdown("**Publication layout** — mean ± " + errtype + " (n)")
        wide = nsum.summary_wide(combined, grp, errtype)
        st.dataframe(wide, use_container_width=True)
        st.caption(wide.attrs.get("caption", ""))
        st.download_button("Download summary (CSV)", wide.to_csv().encode(),
                           "group_summary.csv", "text/csv")
        st.markdown("**Full detail** — mean, SD, SE, 95% CI, min, max, CV%")
        st.dataframe(nsum.group_summary(combined, grp, errtype),
                     use_container_width=True)

# ---------------------------------------------------------------- statistics
with tabs[8]:
    st.subheader("Group comparisons")
    if grp is None:
        st.info("Add a `group` column to the samples sheet.")
    else:
        combined = dv.join(fn).join(nsh[["NSH"]])
        corr_method = st.radio("Multiple-testing correction",
                               ["Benjamini-Hochberg (FDR)", "Holm (family-wise)"],
                               horizontal=True,
                               help="BH controls the false discovery rate and "
                                    "suits exploratory work. Holm controls the "
                                    "family-wise error rate and is stricter; use "
                                    "it when one false positive would be costly.")
        meth = "holm" if corr_method.startswith("Holm") else "BH"
        res = qc.compare(combined, grp)
        res["p_adjusted"] = qc.adjust(res["anova_p"].values, meth)
        disp = res.copy()
        for col in ("anova_p", "kruskal_p", "anova_p_BH", "kruskal_p_BH",
                    "p_adjusted"):
            if col in disp.columns:
                disp[col] = qc.format_p_column(res[col])
        for col in ("anova_F", "kruskal_H", "eta_squared"):
            if col in disp.columns:
                disp[col] = res[col].round(3)
        st.dataframe(disp, use_container_width=True)
        st.caption("p-values below 0.001 are shown in scientific notation. A "
                   "p-value is never exactly zero; full precision is retained in "
                   "the downloaded table.")
        nraw = int((res["anova_p"] < .05).sum())
        nadj = int((res["p_adjusted"] < .05).sum())
        st.info(f"{len(res)} indices tested. Significant at raw p < 0.05: "
                f"{nraw}. After {corr_method}: {nadj}. Report the corrected "
                "column, and state which correction you used.")
        st.warning("ANOVA and Kruskal–Wallis are both shown so that disagreement "
                   "between them is visible. They are diagnostic, not a menu — "
                   "choose one test before seeing the p-values, check its "
                   "assumptions, and report only that one.")
        st.markdown("**Assumption diagnostics**")
        which = st.selectbox("Index", list(combined.columns))
        d0 = qc.diagnostics(combined[which], grp)
        if "error" in d0:
            st.warning(d0["error"])
        else:
            a, b, c_, d_ = st.columns(4)
            a.metric("n", d0["n"])
            b.metric("df", f'{d0["df_between"]}, {d0["df_within"]}')
            c_.metric("F", f'{d0["F"]:.3f}')
            d_.metric("η²", f'{d0["eta_squared"]:.3f}')
            e, f_, g_ = st.columns(3)
            e.metric("Shapiro p (residuals)",
                     f'{d0["shapiro_p_residuals"]:.3f}'
                     if d0["shapiro_p_residuals"] else "n/a")
            f_.metric("Levene p", f'{d0["levene_p"]:.3f}' if d0["levene_p"] else "n/a")
            g_.metric("n per group", str(d0["group_n"]))
            (st.success if d0["residuals_normal"] and d0["variances_equal"]
             else st.warning)(d0["recommendation"])
        st.markdown("**Post-hoc**")
        p1, p2 = st.columns(2)
        with p1:
            st.write("Tukey HSD")
            st.dataframe(npl.tukey(combined[which], grp).round(4),
                         use_container_width=True)
        with p2:
            st.write("Dunn (Bonferroni)")
            st.dataframe(npl.dunn(combined[which], grp).round(4),
                         use_container_width=True)

# ---------------------------------------------------------------- multivariate
with tabs[9]:
    st.subheader("Community-level analysis")
    bc = nc.bray_curtis(counts)
    sc, var = nc.pca(counts)
    if grp is None:
        st.info("Add a `group` column for PERMANOVA.")
    else:
        st.caption("PERMANOVA (Anderson 2001) tests whether whole communities "
                   "differ between groups, in place of many univariate tests. "
                   "PERMDISP (Anderson 2006) checks whether any difference is a "
                   "location shift or merely unequal dispersion.")
        pm = qc.permanova(bc, grp, 999)
        pdz = qc.permdisp(bc, grp, 999)
        a, b, c_ = st.columns(3)
        a.metric("pseudo-F", f"{pm['pseudo_F']:.2f}")
        b.metric("R²", f"{pm['R2']:.3f}")
        c_.metric("p", f"{pm['p']:.3f}")
        (st.warning if pdz.get("p", 1) < .05 else st.success)(
            f"PERMDISP p = {pdz.get('p', float('nan')):.3f} — "
            + ("dispersions also differ, so the PERMANOVA result cannot be "
               "attributed to a location shift alone."
               if pdz.get("p", 1) < .05 else
               "dispersions homogeneous, so this reflects a genuine difference "
               "in composition."))

    def _f():
        fig, ax = plt.subplots(figsize=(6.4, 5.4))
        cols = npl.palette(pal, grp.nunique() if grp is not None else 1)
        if grp is None:
            ax.scatter(sc["PC1"], sc["PC2"], s=70, color=cols[0])
        else:
            for i, lv in enumerate(pd.unique(grp)):
                m = grp.reindex(sc.index) == lv
                ax.scatter(sc.loc[m, "PC1"], sc.loc[m, "PC2"], s=70,
                           color=cols[i], label=str(lv), edgecolor="w")
            ax.legend(frameon=False, fontsize=8)
        ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
        ax.set_title("PCA, Hellinger-transformed", loc="left", fontweight="bold")
        return fig
    show(_f, sc, "pca", "Principal coordinates of Hellinger-transformed "
         "relative abundances.", dpi)

    st.markdown("**Redundant indices**")
    red = qc.redundancy(dv.join(fn))
    if len(red):
        st.dataframe(red, use_container_width=True)
        st.caption("Pairs correlated above r = 0.98. Reporting both is reporting "
                   "one number twice and inflates the multiple-testing burden.")
    else:
        st.success("No near-duplicate pairs above r = 0.98.")

# ---------------------------------------------------------------- nutrients
with tabs[10]:
    st.subheader("Nematode indices versus soil variables")
    st.caption("Nematode indices are usually validated against other nematode "
               "indices, which is circular. Correlating them with variables "
               "measured independently of the nematode community — organic "
               "carbon, mineral nitrogen, available phosphorus, pH, yield — is "
               "the non-circular evidence. Ghaderi et al. (2025) rest their case "
               "for the NSH index on exactly this.")
    ndf = nut.load_nutrients(up, samples)
    if ndf is None or ndf.empty:
        st.info("No soil data found.")
        st.code(nut.TEMPLATE_NOTE, language="text")
    else:
        st.write(f"**{ndf.shape[1]} soil variables** found: "
                 + ", ".join(map(str, ndf.columns)))
        idx_tbl = dv.join(fn).join(nsh[["NSH"]])
        pick = st.multiselect("Indices to correlate", list(idx_tbl.columns),
                              default=[c for c in ["shannon_H", "MI", "PPI", "EI",
                                                   "SI", "BI", "NSH"]
                                       if c in idx_tbl.columns])
        meth = st.radio("Method", ["spearman", "pearson"], horizontal=True,
                        help="Spearman assumes neither linearity nor normality; "
                             "nematode indices are frequently neither.")
        if pick:
            cor = nut.correlate(idx_tbl[pick], ndf, method=meth)
            if "error" in cor.columns:
                st.warning(cor["error"].iloc[0])
            else:
                info = nut.interpret(cor)
                a, b, c_ = st.columns(3)
                a.metric("Pairs tested", info["n_pairs"])
                b.metric("Significant after BH", info["n_significant"])
                c_.metric("Samples", info["max_n"])
                for note in info["notes"]:
                    st.warning(note)
                show(lambda: nut.heatmap(
                        cor, plt.subplots(figsize=(1.1*ndf.shape[1]+3,
                                                   .36*len(pick)+2))[1]).figure,
                     cor, "nutrients",
                     "Spearman correlation coefficients. An asterisk marks pairs "
                     "significant after Benjamini-Hochberg correction; unmarked "
                     "cells should not be interpreted.", dpi)
                if len(info["top"]):
                    st.markdown("**Strongest non-circular relationships**")
                    st.dataframe(info["top"].round(4), use_container_width=True)

        st.divider()
        st.subheader("Redundancy analysis — one test instead of many")
        st.caption("RDA asks how much variation in whole-community composition "
                   "the soil variables explain together, rather than testing "
                   "every index against every variable. R² adjusted is the "
                   "honest figure: raw R² rises mechanically with each variable "
                   "added (Peres-Neto et al. 2006).")
        r = ordn.rda(counts, ndf, 499)
        if "error" in r:
            st.warning(r["error"])
            if ndf.shape[1] >= 4:
                st.markdown("**Reduce the soil variables first**")
                sc_, var_, load_ = ordn.soil_pca(ndf, 3)
                st.write("Variance explained: "
                         + ", ".join(f"{v:.1f}%" for v in var_))
                st.dataframe(load_.round(3), use_container_width=True)
                st.caption("Correlate against these axes instead of the raw "
                           "variables: three tests rather than many.")
        else:
            a, b, c_, d_ = st.columns(4)
            a.metric("R² adjusted", f"{r['R2_adjusted']:.3f}")
            b.metric("R² raw", f"{r['R2']:.3f}")
            c_.metric("F", f"{r['F']:.2f}")
            d_.metric("p", f"{r['p']:.3f}")
            if r["R2_adjusted"] < 0.05:
                st.warning("Adjusted R² is near zero: the soil variables explain "
                           "little community variation once the number of "
                           "variables is accounted for.")
            show(lambda: ordn.triplot(
                    r, grp, ax=plt.subplots(figsize=(7, 5.8))[1]).figure,
                 r["loadings"].round(3), "rda",
                 "Points are samples; arrows are soil variables. Arrow direction "
                 "shows the gradient, length its strength on these two axes.", dpi)

            st.markdown("**Which variables earn their place?**")
            st.caption("Forward selection adds variables one at a time. The "
                       "threshold tightens with the number of candidates, so "
                       "choosing from many does not manufacture significance.")
            st.dataframe(ordn.forward_selection(counts, ndf, 199),
                         use_container_width=True)

# ---------------------------------------------------------------- reference
with tabs[11]:
    st.subheader("Reference table — trophic group, c-p and morphometrics")
    st.caption("Family-level assignments from Yeates et al. (1993) for feeding "
               "habit and Bongers & Bongers (1998) for c-p, as tabulated in "
               "Ferris (2010) Table 1. Look a taxon up here to check what the "
               "software will assign it.")
    _ref = af.reference_table()
    q = st.text_input("Search genus or family")
    c1, c2, c3 = st.columns(3)
    tsel = c1.multiselect("Trophic group", sorted(_ref["trophic"].unique()))
    csel = c2.multiselect("c-p value", sorted(_ref["cp"].unique()))
    only_v = c3.checkbox("Show only rows needing verification")
    view = _ref.copy()
    if q:
        m = view["genus"].str.contains(q, case=False, na=False) | \
            view["family"].str.contains(q, case=False, na=False)
        view = view[m]
    if tsel:
        view = view[view["trophic"].isin(tsel)]
    if csel:
        view = view[view["cp"].isin(csel)]
    if only_v:
        view = view[view["trait_confidence"] == "VERIFY"]
    a, b, c_ = st.columns(3)
    a.metric("Genera shown", len(view))
    b.metric("Families", view["family"].nunique())
    c_.metric("Needing verification",
              int((view["trait_confidence"] == "VERIFY").sum()))
    st.dataframe(view, use_container_width=True, height=430)
    st.download_button("Download reference table (CSV)",
                       _ref.to_csv(index=False).encode(),
                       "NEMO-NEMA_reference_table.csv", "text/csv")
    st.divider()
    st.markdown("""
**Length and width are deliberately blank.**

No verified genus-level compilation of nematode body dimensions was used to
build this table. Filling those columns with plausible numbers would be
fabrication, and the consequence is specific: biomass follows
W = (L x D²)/1,600,000, so an error in L or D propagates directly into every
metabolic footprint, and the output looks entirely normal while being wrong.

**Three ways to obtain them, best first**

1. **Measure your own specimens.** Published dimensions are almost always adult
   females; field samples contain juveniles, so published values overestimate a
   mixed-stage community. Measure 20–30 individuals per taxon across the stages
   actually present.
2. **Taxonomic descriptions.** Species descriptions print de Man's indices:
   *L* (body length) and *a* (length ÷ greatest diameter). Enter `length_um` and
   `a_ratio`; diameter is derived. Sources include Andrássy, *Free-living
   Nematodes of Hungary* I–III; Bongers, *De Nematoden van Nederland*; Goodey,
   *Soil and Freshwater Nematodes*; Jairajpuri & Khan, *Predatory Nematodes*;
   and the *Indian Journal of Nematology* for Indian populations.
3. **Family average weights.** Ferris (2010) Table 1 gives these for about a
   hundred families. Note his own caveat: Dorylaimidae has a mean of 7.46 µg
   with SD 14.72 — a coefficient of variation near 200% — and he states that
   family-level resolution is problematic for biomass.

Every index except biomass and the metabolic footprints works without these
columns. An empty Footprints tab is a correct result; an inflated one is not.
""")

# ---------------------------------------------------------------- validation
with tabs[12]:
    st.subheader("Analysis validation")
    st.caption("Every check recomputes from the raw counts rather than reading a "
               "displayed value, so a disagreement means the pipeline is "
               "inconsistent — not that a number was typed wrongly.")
    try:
        _stats = (qc.compare(dv.join(fn).join(nsh[["NSH"]]), grp)
                  if grp is not None else None)
        vres = qc.validate_analysis(counts, taxa, samples, nrt, dv, fn, nsh,
                                    _stats, grp)
        verdict = vres.attrs["verdict"]
        (st.success if verdict == "PASS" else
         st.warning if verdict == "PASS WITH WARNINGS" else st.error)(
            f"**{verdict}** — {vres.attrs['n_fail']} failed, "
            f"{vres.attrs['n_warn']} warnings, "
            f"{int((vres['status']=='PASS').sum())} passed")
        for sec in vres["section"].unique():
            sub = vres[vres.section == sec]
            worst = ("FAIL" if (sub.status == "FAIL").any() else
                     "WARN" if (sub.status == "WARN").any() else "PASS")
            icon = {"PASS": "\u2705", "WARN": "\u26a0\ufe0f", "FAIL": "\u274c"}[worst]
            with st.expander(f"{icon}  {sec}", expanded=(worst != "PASS")):
                st.dataframe(sub[["check", "status", "detail"]],
                             use_container_width=True, hide_index=True)
        st.download_button("Download validation report (CSV)",
                           vres.to_csv(index=False).encode(),
                           "NEMO-NEMA_validation.csv", "text/csv")
    except Exception as exc:
        st.error(f"Validation could not run: {exc}")

# ---------------------------------------------------------------- report
with tabs[13]:
    st.subheader("Full report")
    st.write("One PDF containing every table, the figures with legends, the "
             "data-quality warnings, an interpretation of each index, a glossary "
             "of every parameter, the methods log and the references.")
    try:
        qcb = {"degeneracy": qc.degeneracy(counts, taxa),
               "dominance": qc.dominance(counts)}
        if samples is not None:
            qcb["replication"] = qc.replication(samples)
            qcb["blocking"] = qc.detect_blocking(samples)
        stats_out = (qc.compare(dv.join(fn).join(nsh[["NSH"]]), grp)
                     if grp is not None else None)
        pdf = rpt.build_report(
            counts, taxa, samples, nrt, dv, fn_q, ex, nsh, fps, stats_out,
            qc_blocks=qcb, plots=npl, palette=pal, err=errtype,
            settings={"prominence_value": pvchoice, "soil_basis": basis,
                      "error_bars": errtype},
            title="Nematode community analysis report")
        st.download_button("Download full report (PDF)", pdf,
                           "NEMO-NEMA_report.pdf", "application/pdf")
    except Exception as exc:
        st.warning(f"Report could not be built: {exc}")

# ---------------------------------------------------------------- cite
with tabs[14]:
    st.subheader("How to cite")
    st.markdown("""
**Cite the software**

> Singh, A.K., Jain, K., Somvanshi, V.S., Pervez, R., Sirohi, A. and Pankaj
> (2026) *NEMO-NEMA: Nematode Ecological Metrics and Ordination for Nematode
> Ecosystem Monitoring and Assessment.* Version 1.1. Division of Nematology,
> ICAR-Indian Agricultural Research Institute, New Delhi.

Add a DOI once the code is archived on Zenodo, and a URL if it is deployed
publicly.

**Also cite the primary sources of whatever you report.** The software
implements published methods; it does not replace them. Citing only the tool
leaves the actual authors of the indices uncredited.

| If you report | Cite |
|---|---|
| frequency, density, prominence value | Norton (1978) |
| Shannon | Shannon & Weaver (1949) |
| Simpson | Simpson (1949) |
| Pielou evenness | Pielou (1966) |
| Margalef richness | Margalef (1958) |
| trophic groups | Yeates et al. (1993) |
| MI, PPI | Bongers (1990) |
| c-p assignments, functional guilds | Bongers & Bongers (1998) |
| EI, SI, CI, BI, faunal profile | Ferris, Bongers & de Goede (2001) |
| biomass | Andrassy (1956) |
| metabolic footprints | Ferris (2010) |
| NSH index | Ghaderi et al. (2025) |
| PERMANOVA | Anderson (2001) |
| PERMDISP | Anderson (2006) |
| FDR correction | Benjamini & Hochberg (1995) |

**State these in your methods**

- which prominence value formula you used (two variants circulate)
- the soil basis of your counts
- your taxonomic resolution
- where each trophic and c-p assignment came from
- which statistical test you chose, and why, before seeing the p-values
""")
    st.download_button(
        "Download CITATION.cff",
        open("CITATION.cff").read().encode() if os.path.exists("CITATION.cff")
        else b"CITATION.cff not found in the app folder",
        "CITATION.cff", "text/plain")

# ---------------------------------------------------------------- manual
with tabs[15]:
    st.subheader("User manual")
    try:
        from nemonema_manual import MANUAL
        st.markdown(MANUAL)
    except Exception:
        st.info("nemonema_manual.py not found in this folder.")

# ---------------------------------------------------------------- download
st.divider()
buf = io.BytesIO()
with pd.ExcelWriter(buf, engine="openpyxl") as x:
    nrt.to_excel(x, sheet_name="Norton_community")
    dv.to_excel(x, sheet_name="Diversity")
    fn_q.to_excel(x, sheet_name="Faunal_indices")
    ex.to_excel(x, sheet_name="Extra_indices")
    nsh.to_excel(x, sheet_name="NSH_index")
    if fps is not None:
        fps.to_excel(x, sheet_name="Footprints")
    qc.degeneracy(counts, taxa).to_excel(x, sheet_name="QC_degeneracy", index=False)
    qc.dominance(counts).to_excel(x, sheet_name="QC_dominance")
    if grp is not None:
        qc.compare(dv.join(fn).join(nsh[["NSH"]]), grp).to_excel(
            x, sheet_name="Group_tests")
        nsum.summary_wide(dv.join(fn).join(nsh[["NSH"]]), grp, errtype).to_excel(
            x, sheet_name="Group_summary")
    taxa.to_excel(x, sheet_name="Trait_assignments")
    pd.DataFrame({"setting": ["software", "prominence_value", "soil_basis",
                              "error_bars", "n_samples", "n_taxa"],
                  "value": [nc.VERSION, pvchoice, basis, errtype,
                            counts.shape[1], counts.shape[0]]}).to_excel(
        x, sheet_name="Methods_log", index=False)
st.download_button("Download all results (.xlsx)", buf.getvalue(),
                   "NEMO-NEMA_results.xlsx",
                   "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
st.caption("Verify every trophic and c-p assignment against a primary source "
           "before publishing: all maturity and food-web indices depend on them.")
footer()
