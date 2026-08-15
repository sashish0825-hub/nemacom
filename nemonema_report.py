import io
from datetime import datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm, inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, PageBreak)
V = "NEMO-NEMA 1.0"
PAGE_W = 210*mm - 40*mm

def _mi(v):
    if np.isnan(v): return "not computable"
    if v < 2.0: return "below 2.0 - disturbed, coloniser-dominated"
    if v < 2.5: return "2.0-2.5 - low maturity, recent enrichment or disturbance"
    if v < 3.0: return "2.5-3.0 - intermediate maturity"
    return "above 3.0 - mature, relatively undisturbed"

def _q(ei, si):
    if np.isnan(ei) or np.isnan(si): return "quadrat not determinable"
    if ei >= 50 and si >= 50: return "Quadrat B - maturing to structured"
    if ei >= 50: return "Quadrat A - disturbed: enriched, structurally simple"
    if si >= 50: return "Quadrat C - structured and stable"
    return "Quadrat D - degraded"

def _n(v):
    if np.isnan(v): return "not computable"
    if v < 15: return "below 15 - degraded food web"
    if v <= 24: return "15-24 - moderate, under stress or in transition"
    return "above 25 - well-functioning"

def _t(df, rows=24, cols=9, fs=6.2):
    d = df.copy().head(rows)
    if d.shape[1] > cols: d = d.iloc[:, :cols]
    if d.select_dtypes(include=[np.number]).shape[1]: d = d.round(2)
    hs = ParagraphStyle("hs", fontName="Helvetica-Bold", fontSize=fs,
                        leading=fs+1.4, textColor=colors.white, alignment=1)
    cs = ParagraphStyle("cs", fontName="Helvetica", fontSize=fs, leading=fs+1.4)
    header = [Paragraph(str(d.index.name or ""), hs)] + \
             [Paragraph(str(c).replace("_", " "), hs) for c in d.columns]
    body = [[Paragraph(str(i), cs)] +
            [Paragraph("" if pd.isna(v) else str(v), cs) for v in r]
            for i, r in zip(d.index, d.values)]
    n = len(header)
    first = min(PAGE_W*0.26, PAGE_W/n*2.2)
    rest = (PAGE_W-first)/(n-1) if n > 1 else PAGE_W
    t = Table([header]+body, colWidths=[first]+[rest]*(n-1),
              repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#17453A")),
        ("GRID",(0,0),(-1,-1),0.25,colors.HexColor("#BCB09A")),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#F8F6F1")]),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("LEFTPADDING",(0,0),(-1,-1),3),("RIGHTPADDING",(0,0),(-1,-1),3),
        ("TOPPADDING",(0,0),(-1,-1),2.5),("BOTTOMPADDING",(0,0),(-1,-1),2.5)]))
    return t


def _im(fig):
    w, h = fig.get_size_inches()
    try:
        from svglib.svglib import svg2rlg
        b = io.BytesIO(); fig.savefig(b, format="svg", bbox_inches="tight")
        plt.close(fig); b.seek(0)
        d = svg2rlg(b)
        sc = min(PAGE_W/d.width, 1.0)
        d.width *= sc; d.height *= sc; d.scale(sc, sc)
        return d
    except Exception:
        pass
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig); b.seek(0)
    sc = min(PAGE_W/(w*inch), 1.0)
    return Image(b, width=w*inch*sc, height=h*inch*sc)

def build_report(counts, taxa, samples, nrt, dv, fn, ex=None, nsh=None, fps=None,
                 stats=None, plots=None, settings=None,
                 palette="Okabe-Ito (CVD-safe)", err="SE",
                 title="Nematode community analysis report"):
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, title=title, leftMargin=20*mm,
                            rightMargin=20*mm, topMargin=18*mm, bottomMargin=18*mm)
    ss = getSampleStyleSheet()
    H1 = ParagraphStyle("a", parent=ss["Heading1"], fontSize=15, textColor=colors.HexColor("#17453A"))
    H2 = ParagraphStyle("b", parent=ss["Heading2"], fontSize=11.5, textColor=colors.HexColor("#2E7D5B"), spaceBefore=10)
    B = ParagraphStyle("c", parent=ss["Normal"], fontSize=9, leading=12.5)
    S = ParagraphStyle("d", parent=ss["Normal"], fontSize=7.6, leading=10, textColor=colors.HexColor("#555555"))
    g = samples["group"] if samples is not None and "group" in samples else None
    st = [Paragraph(title, H1),
          Paragraph(f"Generated {datetime.now():%d %B %Y, %H:%M} - {V}", S), Spacer(1,8),
          Paragraph(f"<b>{counts.shape[0]}</b> taxa - <b>{counts.shape[1]}</b> samples - "
                    f"<b>{int(counts.fillna(0).values.sum()):,}</b> individuals", B), Spacer(1,10),
          Paragraph("1. Community structure - Norton (1978)", H2),
          Paragraph(f"Prominence value: {(settings or {}).get('prominence_value','')}. "
                    "Two variants circulate, both attributed to Norton (1978).", S),
          Spacer(1,4), _t(nrt), PageBreak(),
          Paragraph("2. Diversity", H2), _t(dv)]
    if g is not None and plots is not None:
        try:
            fig, ax = plt.subplots(figsize=(6,3.6))
            plots.bar_with_error(dv["shannon_H"], g, err=err, pal=palette,
                                 ylabel="Shannon H", title="Shannon diversity", ax=ax)
            st += [Spacer(1,6), _im(fig),
                   Paragraph(f"<b>Figure 1.</b> Shannon diversity by treatment group. Bars are "
                             f"group means, error bars {err}, black points individual "
                             "samples so n and spread stay visible.", S)]
        except Exception: pass
    st += [PageBreak(), Paragraph("3. Maturity and food-web indices", H2), _t(fn)]
    mi, ei, si = fn["MI"].mean(), fn["EI"].mean(), fn["SI"].mean()
    st += [Spacer(1,6), Paragraph("<b>Interpretation of dataset means</b>", B),
           Paragraph(f"MI {mi:.2f} - {_mi(mi)}.", B),
           Paragraph(f"EI {ei:.1f}, SI {si:.1f} - {_q(ei,si)} (Ferris et al. 2001).", B)]
    if plots is not None:
        try:
            fig, ax = plt.subplots(figsize=(5.2,5))
            cl = plots.palette(palette, g.nunique() if g is not None else 1)
            if g is None: ax.scatter(fn["SI"], fn["EI"], s=70, color=cl[0])
            else:
                for i, lv in enumerate(pd.unique(g)):
                    m = g.reindex(fn.index) == lv
                    ax.scatter(fn.loc[m,"SI"], fn.loc[m,"EI"], s=70, color=cl[i],
                               label=str(lv), edgecolor="w")
                ax.legend(frameon=False, fontsize=7)
            ax.axhline(50, ls=":", c="grey"); ax.axvline(50, ls=":", c="grey")
            ax.set_xlim(0,100); ax.set_ylim(0,100)
            ax.set_xlabel("Structure index (SI)"); ax.set_ylabel("Enrichment index (EI)")
            for x,y,t in [(4,94,"A Maturing"),(70,94,"B Structured"),(4,3,"D Degraded"),(70,3,"C Disturbed")]:
                ax.text(x,y,t,color="grey",fontsize=7)
            ax.set_title("Faunal profile (Ferris et al. 2001)", loc="left", fontweight="bold", fontsize=9)
            st += [Spacer(1,6), _im(fig), Spacer(1,3),
                   Paragraph("<b>Figure 2.</b> Faunal profile after Ferris et al. "
                             "(2001). Each point is one sample placed by structure "
                             "index (x) and enrichment index (y). Dashed lines mark "
                             "the 50/50 split: A disturbed, B structured, C stable, "
                             "D degraded.", S)]
        except Exception: pass
    if nsh is not None:
        st += [PageBreak(), Paragraph("4. Nematode Soil Health index", H2),
               Paragraph("Ghaderi et al. (2025) Eur J Soil Sci 76:e70149. Range 8-32.", S),
               Spacer(1,4), _t(nsh)]
        v = nsh["NSH"].mean(); used = (nsh["NSH"].max()-nsh["NSH"].min())/24*100
        st += [Spacer(1,6), Paragraph(f"Mean NSH {v:.1f} - {_n(v)}.", B),
               Paragraph(f"Samples span {used:.0f}% of the 8-32 scale. A narrow span "
                         "suggests limited discrimination; the authors state the index "
                         "requires further calibration.", B)]
    if fps is not None:
        st += [PageBreak(), Paragraph("5. Biomass and metabolic footprints", H2),
               Paragraph("Andrassy (1956) biomass; Ferris (2010) footprints.", S),
               Spacer(1,4), _t(fps), Spacer(1,6),
               Paragraph("Published weights derive from adults, so juvenile-rich "
                         "communities are overestimated.", S)]
    if stats is not None and g is not None:
        st += [PageBreak(), Paragraph("6. Group comparisons", H2), _t(stats, rows=30),
               Spacer(1,6),
               Paragraph("ANOVA and Kruskal-Wallis are both shown so disagreement is "
                         "visible. Choose one a priori and correct for multiple testing.", S)]
    st += [PageBreak(), Paragraph("7. Methods log", H2)]
    rows = [("Software", V), ("Generated", f"{datetime.now():%Y-%m-%d %H:%M}"),
            ("Taxa", counts.shape[0]), ("Samples", counts.shape[1])]
    for k, v in (settings or {}).items(): rows.append((k.replace("_"," "), v))
    st += [_t(pd.DataFrame(rows, columns=["setting","value"]).set_index("setting")),
           Spacer(1,10), Paragraph("8. References", H2)]
    for r in ["Andrassy, I. (1956) Acta Zoologica 2: 1-15.",
              "Bongers, T. (1990) Oecologia 83: 14-19.",
              "Bongers, T. & Bongers, M. (1998) Applied Soil Ecology 10: 239-251.",
              "Ferris, H., Bongers, T. & de Goede, R.G.M. (2001) Applied Soil Ecology 18: 13-29.",
              "Ferris, H. (2010) European Journal of Soil Biology 46: 97-104.",
              "Ghaderi, R. et al. (2025) European Journal of Soil Science 76: e70149.",
              "Norton, D.C. (1978) Ecology of Plant-Parasitic Nematodes. Wiley.",
              "Yeates, G.W. et al. (1993) Journal of Nematology 25: 315-331."]:
        st += [Paragraph(r, S), Spacer(1,2)]
    st += [Spacer(1,10),
           Paragraph("This report describes nematode community structure. It contains no "
                     "management recommendation, because none can be justified from a "
                     "nematode community alone.", S)]
    doc.build(st); buf.seek(0); return buf.getvalue()
