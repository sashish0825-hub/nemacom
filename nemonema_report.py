"""
nemonema_report.py — one-file PDF report.

Every table, every figure, the data-quality warnings, an interpretation of each
index, the methods log and the references.

ON THE INTERPRETATION TEXT
Each sentence restates what an index means according to its source paper,
applied to the value computed. It is descriptive, never prescriptive: no
management recommendation is produced, because none can be justified from a
nematode community alone. Where an index is uninformative for the dataset, the
report says so rather than describing a meaningless number.
"""

from __future__ import annotations

import io
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm, inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, PageBreak)

VERSION = "NEMO-NEMA 1.1"
PAGE_W = 210 * mm - 40 * mm

GLOSSARY = [
 ("Absolute frequency %", "Proportion of samples containing the taxon",
  "Separates widespread taxa from patchy ones"),
 ("Mean density", "Individuals per sample unit",
  "The raw pressure measure; drives prominence value"),
 ("Prominence value", "Density x sqrt(frequency)",
  "Combines abundance and spread; heavily weighted towards abundant taxa"),
 ("Importance value", "Relative frequency + relative density",
  "Less dominated by sheer abundance than prominence value; report both"),
 ("Richness S", "Number of taxa present",
  "Simplest diversity measure, but blind to how evenly individuals are spread"),
 ("Shannon H'", "Diversity weighting richness and evenness",
  "The standard diversity index; falls when one taxon dominates"),
 ("Simpson 1-D", "Probability two random individuals differ",
  "More sensitive to dominant taxa than Shannon"),
 ("Pielou J'", "Evenness: H' divided by its maximum",
  "Isolates evenness from richness; undefined when only one taxon is present"),
 ("Margalef d", "Richness corrected for sample size",
  "Allows richness comparison across unequal sample totals"),
 ("MI", "Mean c-p value of free-living nematodes",
  "The primary disturbance indicator; low means colonisers dominate"),
 ("MI 2-5", "MI excluding c-p 1 taxa",
  "Detects perturbation unrelated to nutrient enrichment"),
 ("sigma MI", "MI including plant parasites",
  "Whole-assemblage maturity"),
 ("PPI", "Mean c-p value of plant parasites",
  "Characterises which parasite guilds dominate; rises with longidorids"),
 ("PPI/MI", "Ratio of the two maturity measures",
  "Above 1 suggests enrichment; the two respond oppositely to fertilisation"),
 ("EI", "Enrichment index, 0-100",
  "How much labile resource has recently entered the food web"),
 ("SI", "Structure index, 0-100",
  "Food-web complexity; high means higher trophic levels are intact"),
 ("BI", "Basal index, 0-100",
  "Share of the web in stress-tolerant basal taxa; high means depleted"),
 ("CI", "Channel index, 0-100",
  "Bacterial (low) versus fungal (high) decomposition; uninformative without "
  "fungivores"),
 ("NCR", "Bacterivores / (bacterivores + fungivores)",
  "Simpler decomposition-channel ratio; uninformative without fungivores"),
 ("Faunal quadrat", "Position on the EI x SI plane",
  "The most informative single food-web diagnostic: A disturbed, B structured, "
  "C stable, D degraded"),
 ("Metabolic footprint", "Carbon used in production and respiration",
  "Adds magnitude: two webs with identical EI and SI can move very different "
  "amounts of carbon"),
 ("NSH", "Seven nematode indices summed, 8-32",
  "A single soil health score; comparative, not absolute, until regionally "
  "calibrated"),
]


def _mi_text(v):
    if pd.isna(v):
        return "not computable"
    if v < 2.0:
        return "below 2.0 — a disturbed, coloniser-dominated community (Bongers 1990)"
    if v < 2.5:
        return "2.0–2.5 — low maturity, indicating recent enrichment or disturbance"
    if v < 3.0:
        return "2.5–3.0 — intermediate maturity"
    return "above 3.0 — a mature, relatively undisturbed community"


def _quad_text(ei, si):
    """Ferris et al. (2001) Fig. 1 and Table 5."""
    if pd.isna(ei) or pd.isna(si):
        return "quadrat not determinable"
    if ei >= 50 and si >= 50:
        return ("Quadrat B — maturing to structured: enriched and structured, "
                "low to moderate disturbance")
    if ei >= 50:
        return ("Quadrat A — disturbed: enriched but structurally simple, "
                "bacterial decomposition channel")
    if si >= 50:
        return ("Quadrat C — stable: undisturbed and mature, fungal "
                "decomposition channel")
    return "Quadrat D — degraded: depleted and structurally simple"


def _nsh_text(v):
    if pd.isna(v):
        return "not computable"
    if v < 15:
        return "below 15 — degraded, low biodiversity and disturbed food web"
    if v <= 24:
        return "15–24 — moderate: some structure but under stress or in transition"
    return "above 25 — well-functioning, diverse and structured"


def _t(df, rows=200, cols=9, fs=6.2):
    """Table fitted to the printable width; headers wrap rather than overflow."""
    n0 = len(df)
    d = df.copy().head(rows)
    if d.shape[1] > cols:
        d = d.iloc[:, :cols]
    if d.select_dtypes(include=[np.number]).shape[1]:
        d = d.round(2)
    hs = ParagraphStyle("hs", fontName="Helvetica-Bold", fontSize=fs,
                        leading=fs + 1.4, textColor=colors.white, alignment=1)
    cs = ParagraphStyle("cs", fontName="Helvetica", fontSize=fs, leading=fs + 1.4)
    header = [Paragraph(str(d.index.name or ""), hs)] + \
             [Paragraph(str(c).replace("_", " "), hs) for c in d.columns]
    body = [[Paragraph(str(i), cs)] +
            [Paragraph("" if pd.isna(v) else str(v), cs) for v in r]
            for i, r in zip(d.index, d.values)]
    n = len(header)
    first = min(PAGE_W * 0.26, PAGE_W / n * 2.2)
    rest = (PAGE_W - first) / (n - 1) if n > 1 else PAGE_W
    t = Table([header] + body, colWidths=[first] + [rest] * (n - 1),
              repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#17453A")),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#BCB09A")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#F8F6F1")]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5)]))
    t._n_hidden = max(0, n0 - rows)
    return t


def _im(fig):
    """Vector where svglib is available, so axis labels and legends remain
    editable text in the PDF. Falls back to raster otherwise."""
    w, h = fig.get_size_inches()
    try:
        from svglib.svglib import svg2rlg
        b = io.BytesIO()
        fig.savefig(b, format="svg", bbox_inches="tight")
        plt.close(fig)
        b.seek(0)
        d = svg2rlg(b)
        sc = min(PAGE_W / d.width, 1.0)
        d.width *= sc
        d.height *= sc
        d.scale(sc, sc)
        return d
    except Exception:
        pass
    b = io.BytesIO()
    fig.savefig(b, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    b.seek(0)
    sc = min(PAGE_W / (w * inch), 1.0)
    return Image(b, width=w * inch * sc, height=h * inch * sc)


def build_report(counts, taxa, samples, nrt, dv, fn, ex=None, nsh=None, fps=None,
                 stats=None, qc_blocks=None, plots=None, settings=None,
                 palette="Okabe-Ito (CVD-safe)", err="SE",
                 title="Nematode community analysis report"):
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, title=title, leftMargin=20 * mm,
                            rightMargin=20 * mm, topMargin=18 * mm,
                            bottomMargin=18 * mm)
    ss = getSampleStyleSheet()
    H1 = ParagraphStyle("a", parent=ss["Heading1"], fontSize=15,
                        textColor=colors.HexColor("#17453A"))
    H2 = ParagraphStyle("b", parent=ss["Heading2"], fontSize=11.5,
                        textColor=colors.HexColor("#2E7D5B"), spaceBefore=10)
    B = ParagraphStyle("c", parent=ss["Normal"], fontSize=9, leading=12.5)
    S = ParagraphStyle("d", parent=ss["Normal"], fontSize=7.6, leading=10,
                       textColor=colors.HexColor("#555555"))
    W = ParagraphStyle("e", parent=B, textColor=colors.HexColor("#A03030"),
                       leftIndent=8)
    g = samples["group"] if samples is not None and "group" in samples else None
    st = [Paragraph(title, H1),
          Paragraph(f"Generated {datetime.now():%d %B %Y, %H:%M} · {VERSION}", S),
          Spacer(1, 8),
          Paragraph(f"<b>{counts.shape[0]}</b> taxa · <b>{counts.shape[1]}</b> "
                    f"samples · <b>{int(counts.fillna(0).values.sum()):,}</b> "
                    "individuals"
                    + (f" · <b>{g.nunique()}</b> groups" if g is not None else ""),
                    B), Spacer(1, 10)]

    # ---- 1 data quality, deliberately first
    st += [Paragraph("1. Data quality and design", H2)]
    if qc_blocks:
        rep = qc_blocks.get("replication", {})
        if rep.get("message"):
            st += [Paragraph(rep["message"],
                             W if not rep.get("ok", True) else B)]
        blk = qc_blocks.get("blocking", {})
        if blk.get("blocked"):
            st += [Paragraph("<b>" + blk["message"] + "</b>", W)]
        deg = qc_blocks.get("degeneracy")
        if deg is not None:
            bad = deg[~deg["informative"]]
            if len(bad):
                st += [Spacer(1, 4),
                       Paragraph("<b>Indices carrying no information for this "
                                 "dataset:</b>", W)]
                for _, r in bad.iterrows():
                    st += [Paragraph(f"• <b>{r['index']}</b> — {r['reason']}. "
                                     f"{r['action']}.", W)]
        dom = qc_blocks.get("dominance")
        if dom is not None:
            n = int((dom["flag"] == "extreme dominance").sum())
            if n:
                st += [Spacer(1, 4),
                       Paragraph(f"{n} of {len(dom)} samples exceed 70% dominance "
                                 "by a single taxon. Diversity indices then "
                                 "largely report that taxon's share.", W)]
    if "source" not in taxa.columns or taxa["source"].isna().all():
        st += [Paragraph("No source recorded for trophic and c-p assignments. "
                         "Every maturity and food-web index is built from them; "
                         "without sources these results are not reproducible.", W)]

    # ---- 2 Norton
    st += [PageBreak(), Paragraph("2. Community structure — Norton (1978)", H2),
           Paragraph(f"Prominence value computed as "
                     f"<b>{(settings or {}).get('prominence_value', '')}</b>. Two "
                     "variants circulate, both attributed to Norton (1978).", S),
           Spacer(1, 4), _t(nrt)]
    if nrt.attrs.get("pv_uninformative"):
        st += [Spacer(1, 4),
               Paragraph("Every taxon occurs in every sample, so absolute "
                         "frequency is 100 throughout and prominence value equals "
                         "density × 10 exactly. It adds nothing beyond mean "
                         "density here.", W)]

    # ---- 3 diversity
    st += [PageBreak(), Paragraph("3. Diversity", H2), _t(dv)]
    if g is not None and plots is not None:
        try:
            fig, ax = plt.subplots(figsize=(6.2, 3.8))
            plots.bar_with_error(dv["shannon_H"], g, err=err, pal=palette,
                                 ylabel="Shannon H'", title="Shannon diversity",
                                 ax=ax)
            st += [Spacer(1, 6), _im(fig), Spacer(1, 3),
                   Paragraph(f"<b>Figure 1.</b> Shannon diversity by group. Bars "
                             f"are group means, error bars {err}, black points "
                             "individual samples so that n and spread remain "
                             "visible.", S)]
        except Exception:
            pass

    # ---- 4 faunal
    st += [PageBreak(), Paragraph("4. Maturity and food-web indices", H2), _t(fn)]
    mi, ei, si = fn["MI"].mean(), fn["EI"].mean(), fn["SI"].mean()
    st += [Spacer(1, 6), Paragraph("<b>Interpretation of dataset means</b>", B),
           Paragraph(f"• Maturity index {mi:.2f} — {_mi_text(mi)}.", B),
           Paragraph(f"• EI {ei:.1f}, SI {si:.1f} — {_quad_text(ei, si)} "
                     "(Ferris et al. 2001).", B)]
    if qc_blocks and qc_blocks.get("degeneracy") is not None:
        for _, r in qc_blocks["degeneracy"][
                ~qc_blocks["degeneracy"]["informative"]].iterrows():
            st += [Paragraph(f"• <b>{r['index']}</b> is not interpretable here "
                             f"({r['reason']}).", W)]
    if plots is not None:
        try:
            fig, ax = plt.subplots(figsize=(5.4, 5.2))
            plots.faunal_profile(fn, g, palette, ax)
            st += [Spacer(1, 6), _im(fig), Spacer(1, 3),
                   Paragraph("<b>Figure 2.</b> Faunal profile after Ferris et al. "
                             "(2001). Each point is one sample positioned by its "
                             "structure index (x) and enrichment index (y). Dashed "
                             "lines mark the 50/50 split into four quadrats: "
                             "A disturbed and enriched, B maturing to structured, "
                             "C structured and stable, D degraded.", S)]
        except Exception:
            pass
    if plots is not None:
        try:
            fig, ax = plt.subplots(figsize=(7.4, 4.4))
            plots.trophic_bars(fn, palette, ax)
            st += [PageBreak(), _im(fig), Spacer(1, 3),
                   Paragraph("<b>Figure 3.</b> Trophic composition of each sample, "
                             "as a percentage of total nematodes.", S)]
        except Exception:
            pass

    # ---- 5 NSH
    if nsh is not None:
        st += [PageBreak(), Paragraph("5. Nematode Soil Health index", H2),
               Paragraph("Ghaderi et al. (2025) Eur J Soil Sci 76:e70149. Seven "
                         "subscores summed; range 8–32.", S), Spacer(1, 4), _t(nsh)]
        v = nsh["NSH"].mean()
        vv = nsh["NSH"].dropna()
        used = (vv.max() - vv.min()) / 24 * 100 if len(vv) else np.nan
        st += [Spacer(1, 6), Paragraph(f"Mean NSH {v:.1f} — {_nsh_text(v)}.", B),
               Paragraph(f"Samples span {used:.0f}% of the 8–32 scale.", B)]
        sub = [c for c in nsh.columns if c.endswith("_score")]
        const = [c for c in sub if nsh[c].nunique() == 1]
        if const:
            st += [Paragraph(f"{len(const)} of {len(sub)} subscores are identical "
                             f"in every sample ({', '.join(const)}). They add a "
                             "fixed amount to every NSH value and contribute no "
                             "discrimination. This is worth reporting.", W)]
        if used < 35:
            st += [Paragraph("A narrow span suggests the index is not "
                             "discriminating strongly in this system. The authors "
                             "state it requires further calibration, and their "
                             "abundance bands derive from global biome ranges in "
                             "which tropical cropland is thinly represented.", W)]
        basis = (settings or {}).get("soil_basis", "")
        if basis and "100 g" not in str(basis):
            st += [Paragraph(f"Counts are expressed per {basis}. The abundance "
                             "subscore assumes per 100 g dry soil, so that "
                             "subscore — and every total above — is on the wrong "
                             "scale until converted.", W)]

    # ---- 6 footprints
    if fps is not None:
        st += [PageBreak(), Paragraph("6. Biomass and metabolic footprints", H2),
               Paragraph("Andrássy (1956) biomass; Ferris (2010) footprints.", S),
               Spacer(1, 4), _t(fps), Spacer(1, 6),
               Paragraph("Footprints assume all individuals of a taxon share one "
                         "body size and that life-course duration scales linearly "
                         "with c-p class. Published weights derive from adults, so "
                         "juvenile-rich communities are overestimated.", S)]

    # ---- 7 statistics
    if stats is not None and g is not None:
        st += [PageBreak(), Paragraph("7. Group comparisons", H2), _t(stats, rows=40),
               Spacer(1, 6),
               Paragraph("ANOVA and Kruskal–Wallis are both shown so that "
                         "disagreement between them is visible. They are "
                         "diagnostic, not a menu: choose one test before seeing "
                         "the p-values, check its assumptions, and report only "
                         "that one. Report the Benjamini–Hochberg corrected "
                         "columns, not the raw p-values.", S)]

    # ---- 8 glossary
    st += [PageBreak(), Paragraph("8. What each parameter means", H2),
           Paragraph("Every index carries assumptions; none is meaningful in "
                     "isolation.", S), Spacer(1, 5),
           _t(pd.DataFrame(GLOSSARY, columns=["Parameter", "What it measures",
                                              "Why it matters"]
                           ).set_index("Parameter"), rows=40, fs=6.6)]

    # ---- 9 methods
    st += [PageBreak(), Paragraph("9. Methods log", H2)]
    rows = [("Software", VERSION), ("Generated", f"{datetime.now():%Y-%m-%d %H:%M}"),
            ("Taxa", counts.shape[0]), ("Samples", counts.shape[1])]
    for k, v in (settings or {}).items():
        rows.append((str(k).replace("_", " "), v))
    if "source" in taxa.columns:
        rows.append(("Trait sources recorded",
                     f"{int(taxa['source'].notna().sum())} of {len(taxa)} taxa"))
    st += [_t(pd.DataFrame(rows, columns=["setting", "value"]).set_index("setting")),
           Spacer(1, 10), Paragraph("10. References", H2)]
    for r in ["Andrássy, I. (1956) Acta Zoologica 2: 1–15.",
              "Anderson, M.J. (2001) Austral Ecology 26: 32–46.",
              "Anderson, M.J. (2006) Biometrics 62: 245–253.",
              "Benjamini, Y. & Hochberg, Y. (1995) J R Stat Soc B 57: 289–300.",
              "Bongers, T. (1990) Oecologia 83: 14–19.",
              "Bongers, T. & Bongers, M. (1998) Applied Soil Ecology 10: 239–251.",
              "Ferris, H., Bongers, T. & de Goede, R.G.M. (2001) Applied Soil "
              "Ecology 18: 13–29.",
              "Ferris, H. (2010) European Journal of Soil Biology 46: 97–104.",
              "Ghaderi, R. et al. (2025) European Journal of Soil Science 76: e70149.",
              "Norton, D.C. (1978) Ecology of Plant-Parasitic Nematodes. Wiley.",
              "Yeates, G.W. et al. (1993) Journal of Nematology 25: 315–331."]:
        st += [Paragraph(r, S), Spacer(1, 2)]
    st += [Spacer(1, 12),
           Paragraph("Developed by Ashish Kumar Singh, Kavita Jain, Vishal Singh "
                     "Somvanshi, Rashid Pervez, Anil Sirohi and Pankaj. Division "
                     "of Nematology, ICAR-Indian Agricultural Research Institute, "
                     "New Delhi 110012.", S),
           Spacer(1, 6),
           Paragraph("This report describes nematode community structure. It "
                     "contains no management recommendation, because none can be "
                     "justified from a nematode community alone.", S)]
    doc.build(st)
    buf.seek(0)
    return buf.getvalue()
