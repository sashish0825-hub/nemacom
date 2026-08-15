HELP = """
nemacom_all.py — NEMO-NEMA in one file.

  python3 nemacom_all.py template     -> writes NEMO-NEMA_input_template.xlsx
  python3 nemacom_all.py demo         -> writes NEMO-NEMA_demo_data.xlsx
  python3 nemacom_all.py run FILE.xlsx-> analyses, writes NEMO-NEMA_results.xlsx + figures
  streamlit run nemacom_all.py        -> interactive app

SOURCES (verify each before publishing):
  Norton (1978) Ecology of Plant-Parasitic Nematodes, Wiley  -> AF, RF, D, RD, PV
  Bongers (1990) Oecologia 83:14-19                          -> MI, PPI
  Bongers & Bongers (1998) Appl Soil Ecol 10:239-251         -> MI(2-5), sigma-MI
  Ferris, Bongers & de Goede (2001) Appl Soil Ecol 18:13-29  -> EI, SI, CI, BI
  Yeates et al. (1993) J Nematol 25:315-331                  -> trophic groups
  Andrassy (1956) Acta Zoologica 2:1-15                      -> biomass
  Ferris (2010) Eur J Soil Biol 46:97-104                    -> metabolic footprints
"""
import io
import sys
import numpy as np
import pandas as pd

TROPHIC = {"PP": "Plant parasite", "BF": "Bacterivore", "FF": "Fungivore",
           "OM": "Omnivore", "PR": "Predator"}

# Ferris et al. (2001) guild weights — TRANSCRIBED, NOT VERIFIED.
# Check every value against Table 1 of the paper before publishing.
W = {"enrichment": {("BF", 1): 3.2, ("FF", 2): 0.8},
     "basal": {("BF", 2): 0.8, ("FF", 2): 0.8},
     "structure": {("BF", 3): 1.8, ("BF", 4): 3.2, ("BF", 5): 5.0,
                   ("FF", 3): 1.8, ("FF", 4): 3.2, ("FF", 5): 5.0,
                   ("OM", 3): 1.8, ("OM", 4): 3.2, ("OM", 5): 5.0,
                   ("PR", 2): 0.8, ("PR", 3): 1.8, ("PR", 4): 3.2, ("PR", 5): 5.0}}


def load(path):
    def clean(df):
        df = df[df.index.notna()].copy()
        df.index = df.index.astype(str).str.strip()
        return df[~df.index.isin(["", "nan"]) & ~df.index.str.startswith("<-")]
    xl = pd.ExcelFile(path)
    counts = clean(pd.read_excel(xl, "counts", index_col=0))
    counts.columns = [str(c).strip() for c in counts.columns]
    taxa = clean(pd.read_excel(xl, "taxa", index_col=0))
    taxa.columns = [c.strip().lower() for c in taxa.columns]
    samples = None
    if "samples" in xl.sheet_names:
        samples = clean(pd.read_excel(xl, "samples", index_col=0))
        samples.columns = [c.strip().lower() for c in samples.columns]
    return counts, taxa, samples


def validate(counts, taxa, samples=None):
    p = []
    if counts.empty:
        return ["'counts' sheet is empty."]
    bad = [c for c in counts.columns if not pd.api.types.is_numeric_dtype(counts[c])]
    if bad:
        p.append(f"Non-numeric sample columns: {bad}")
    if (counts.fillna(0) < 0).any().any():
        p.append("Negative counts found.")
    miss = sorted(set(counts.index) - set(taxa.index))
    if miss:
        p.append(f"Taxa in 'counts' missing from 'taxa' sheet: {miss}")
    bt = taxa.loc[taxa.index.isin(counts.index) & ~taxa["trophic"].isin(TROPHIC)]
    if len(bt):
        p.append(f"Bad trophic codes for {list(bt.index)}; use {sorted(TROPHIC)}")
    cp = pd.to_numeric(taxa["cp"], errors="coerce")
    bc = taxa.loc[taxa.index.isin(counts.index) & ~cp.isin([1, 2, 3, 4, 5])]
    if len(bc):
        p.append(f"c-p must be 1-5. Check: {list(bc.index)}")
    if samples is not None:
        un = sorted(set(counts.columns) - set(samples.index))
        if un:
            p.append(f"Sample columns with no row in 'samples': {un}")
    return p


def norton(counts, pv_divisor=1.0):
    """Norton (1978) community descriptors. PV = D*sqrt(F) / pv_divisor."""
    c = counts.fillna(0).astype(float)
    n = c.shape[1]
    pres = (c > 0).sum(axis=1)
    o = pd.DataFrame({"samples_positive": pres.astype(int),
                      "total_individuals": c.sum(axis=1),
                      "absolute_frequency_pct": pres / n * 100,
                      "mean_density": c.sum(axis=1) / n})
    fs, ds = o["absolute_frequency_pct"].sum(), o["mean_density"].sum()
    o["relative_frequency_pct"] = o["absolute_frequency_pct"] / fs * 100 if fs else 0.0
    o["relative_density_pct"] = o["mean_density"] / ds * 100 if ds else 0.0
    o["prominence_value"] = (o["mean_density"]
                             * np.sqrt(o["absolute_frequency_pct"]) / pv_divisor)
    o["importance_value"] = o["relative_frequency_pct"] + o["relative_density_pct"]
    return o.sort_values("prominence_value", ascending=False)


def diversity(counts):
    c = counts.fillna(0).astype(float)
    r = {}
    for s in c.columns:
        x = c[s].values
        x = x[x > 0]
        N, S = x.sum(), len(x)
        if N == 0:
            r[s] = dict(richness_S=0, total_N=0.0, shannon_H=np.nan,
                        simpson_D=np.nan, simpson_1minusD=np.nan,
                        pielou_J=np.nan, margalef_d=np.nan)
            continue
        p = x / N
        H = float(-(p * np.log(p)).sum())
        D = float((p ** 2).sum())
        r[s] = dict(richness_S=int(S), total_N=float(N), shannon_H=H,
                    simpson_D=D, simpson_1minusD=1 - D,
                    pielou_J=H / np.log(S) if S > 1 else np.nan,
                    margalef_d=(S - 1) / np.log(N) if N > 1 else np.nan)
    return pd.DataFrame(r).T


def faunal(counts, taxa):
    meta = taxa.reindex(counts.index)
    key = pd.MultiIndex.from_arrays(
        [meta["trophic"].values,
         pd.to_numeric(meta["cp"], errors="coerce").astype("Int64").values],
        names=["trophic", "cp"])
    g = counts.fillna(0).astype(float).copy()
    g.index = key
    g = g.groupby(level=["trophic", "cp"]).sum()

    def mi(sr):
        t = sr.sum()
        return float(sum(cp * v for (_, cp), v in sr.items()) / t) if t else np.nan

    r = {}
    for s in g.columns:
        x = g[s]
        x = x[x > 0]
        if x.sum() == 0:
            r[s] = {}
            continue
        free = x[[t != "PP" for t, _ in x.index]]
        pp = x[[t == "PP" for t, _ in x.index]]
        f25 = free[[cp >= 2 for _, cp in free.index]]
        e = sum(W["enrichment"].get(k, 0) * v for k, v in x.items())
        b = sum(W["basal"].get(k, 0) * v for k, v in x.items())
        st = sum(W["structure"].get(k, 0) * v for k, v in x.items())
        fu2 = float(x.get(("FF", 2), 0.0))
        cid = 3.2 * float(x.get(("BF", 1), 0.0)) + 0.8 * fu2
        bt = x.groupby(level="trophic").sum()
        tot = float(x.sum())
        ba, fu = float(bt.get("BF", 0)), float(bt.get("FF", 0))
        r[s] = {"MI": mi(free), "MI_2_5": mi(f25), "sigma_MI": mi(x), "PPI": mi(pp),
                "PPI_MI_ratio": mi(pp) / mi(free) if mi(free) else np.nan,
                "EI": 100 * e / (e + b) if e + b else np.nan,
                "SI": 100 * st / (st + b) if st + b else np.nan,
                "CI": 100 * (0.8 * fu2) / cid if cid else np.nan,
                "BI": 100 * b / (e + st + b) if e + st + b else np.nan,
                "NCR": ba / (ba + fu) if ba + fu else np.nan,
                **{f"pct_{t}": 100 * float(bt.get(t, 0)) / tot for t in TROPHIC}}
    return pd.DataFrame(r).T


def extra_indices(counts, taxa):
    c = counts.fillna(0).astype(float)
    tr = taxa["trophic"].reindex(c.index)
    cp = pd.to_numeric(taxa["cp"], errors="coerce").reindex(c.index)
    r = {}
    for s in c.columns:
        x = c[s]
        bt = x.groupby(tr).sum()
        pp = float(bt.get("PP", 0))
        fb = float(bt.get("BF", 0)) + float(bt.get("FF", 0))
        tot = float(bt.sum())
        p = bt[bt > 0] / tot if tot else bt
        r[s] = {"wasilewska_index": fb / pp if pp else np.nan,
                "trophic_diversity_H": float(-(p * np.log(p)).sum()) if tot else np.nan,
                "functional_guild_richness": len({(t, k) for t, k, v
                                                  in zip(tr, cp, x) if v > 0})}
    return pd.DataFrame(r).T


def bray_curtis(counts):
    from scipy.spatial.distance import pdist, squareform
    d = squareform(pdist(counts.fillna(0).astype(float).T.values, metric="braycurtis"))
    return pd.DataFrame(d, index=counts.columns, columns=counts.columns)


def pca(counts, k=2):
    X = counts.fillna(0).astype(float).T
    H = np.sqrt(X.div(X.sum(axis=1).replace(0, np.nan), axis=0)).fillna(0)
    Xc = H - H.mean(axis=0)
    U, S, _ = np.linalg.svd(Xc.values, full_matrices=False)
    k = min(k, len(S))
    sc = pd.DataFrame(U[:, :k] * S[:k], index=H.index,
                      columns=[f"PC{i+1}" for i in range(k)])
    return sc, (S ** 2) / (S ** 2).sum() * 100


def compare(tbl, groups):
    from scipy import stats
    g = groups.reindex(tbl.index)
    out = []
    for col in tbl.columns:
        arr = [tbl.loc[g == l, col].dropna().values for l in g.dropna().unique()]
        arr = [a for a in arr if len(a) > 1]
        if len(arr) < 2:
            out.append(dict(index=col, anova_F=np.nan, anova_p=np.nan,
                            kruskal_H=np.nan, kruskal_p=np.nan))
            continue
        try:
            F, pa = stats.f_oneway(*arr)
        except Exception:
            F, pa = np.nan, np.nan
        try:
            Hs, pk = stats.kruskal(*arr)
        except Exception:
            Hs, pk = np.nan, np.nan
        out.append(dict(index=col, anova_F=F, anova_p=pa, kruskal_H=Hs, kruskal_p=pk))
    return pd.DataFrame(out).set_index("index")


# ---------------------------------------------------------------- NSH index
# Ghaderi, Hayden, Jayaramaiah, Hu & He (2025) An Innovative Framework Fosters
# Practical Application of Nematode-Based Indices in Soil Health Assessment.
# European Journal of Soil Science 76: e70149.  doi:10.1111/ejss.70149
#
# Bands and scores transcribed from the 'NSH' sheet of Supporting Information 1
# and cross-checked against the live IF() formulas in the Dataset-1..5 sheets.
# Seven subindicators, each 1-5 except PPI (max 4) and CI (max 3), so the total
# runs 8-32 exactly as the paper states.
#
# Interpretation given by the authors: NSH < 15 degraded soil health; 15-24
# moderate; > 25 well-functioning. The authors state the index still requires
# calibration on more datasets, and their abundance bands come from global biome
# ranges (van den Hoogen et al. 2019) that under-represent tropical cropland.
NSH_BANDS = {
    #            (lower, upper, score) ; upper is exclusive except the last band
    "MI":  [(1.0, 1.5, 1), (1.5, 2.0, 2), (2.0, 2.5, 3), (2.5, 3.0, 4), (3.0, 5.0, 5)],
    "PPI": [(2.0, 2.5, 4), (2.5, 3.0, 2), (3.0, 3.5, 1), (3.5, 5.0, 3)],
    "EI":  [(0, 20, 1), (20, 40, 2), (40, 60, 3), (60, 80, 4), (80, 100, 5)],
    "SI":  [(0, 20, 1), (20, 40, 2), (40, 60, 3), (60, 80, 4), (80, 100, 5)],
    "BI":  [(0, 20, 5), (20, 40, 4), (40, 60, 3), (60, 80, 2), (80, 100, 1)],
    "CI":  [(0, 20, 2), (20, 80, 3), (80, 100, 2)],
    # MF, total abundance, or total biomass of the whole assemblage per 100 g dry soil
    "MF":  [(-np.inf, 200, 1), (200, 500, 2), (500, 1000, 3),
            (1000, 2000, 4), (2000, np.inf, 5)],
}
# Alternative MF bands if you supply total biomass (ug per 100 g dry soil)
NSH_BIOMASS_BANDS = [(-np.inf, 140, 1), (140, 350, 2), (350, 700, 3),
                     (700, 1400, 4), (1400, np.inf, 5)]


def _band_score(value, bands):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    # Guard the band edges against floating-point residue. Values such as
    # 2.4999999999999996 and 19.999999999999996 occur in the authors' own
    # supplementary workbook and would otherwise fall into the band below.
    value = round(float(value), 6)
    for i, (lo, hi, sc) in enumerate(bands):
        last = i == len(bands) - 1
        if (lo <= value < hi) or (last and lo <= value <= hi):
            return sc
    return np.nan          # outside every band - reported, never silently zeroed


def nsh_index(faunal_table, abundance, mf_basis="abundance"):
    """Nematode Soil Health index, Ghaderi et al. (2025).

    faunal_table : per-sample MI, PPI, EI, SI, BI, CI (output of faunal()).
    abundance    : per-sample total nematodes (or MF, or biomass) per 100 g dry soil.
    mf_basis     : 'abundance'/'mf' uses the 200/500/1000/2000 bands;
                   'biomass' uses the 140/350/700/1400 ug bands.
    """
    bands = dict(NSH_BANDS)
    if mf_basis == "biomass":
        bands["MF"] = NSH_BIOMASS_BANDS
    out = pd.DataFrame(index=faunal_table.index)
    for k in ("MI", "PPI", "EI", "SI", "BI", "CI"):
        out[f"{k}_score"] = [_band_score(v, bands[k]) for v in faunal_table[k]]
    ab = pd.Series(abundance).reindex(faunal_table.index)
    out["MF_score"] = [_band_score(v, bands["MF"]) for v in ab]
    out["NSH"] = out.sum(axis=1, skipna=False)
    out["interpretation"] = pd.cut(out["NSH"], [-np.inf, 15, 25, np.inf],
                                   labels=["degraded", "moderate", "well-functioning"],
                                   right=False)
    return out


# ---------------------------------------------------------------- biomass & footprints
# VERIFIED against the primary sources on 14 Aug 2026:
#   Andrassy (1956): W = (L * D^2) / (1.6e6), W in ug, L and D in um.
#   Ferris (2010) Eur J Soil Biol 46:97-104, section 2.1.2.3:
#       F = SUM Nt * [ 0.1*(Wt/mt) + 0.273*Wt^0.75 ]
#       where Nt = abundance, Wt = fresh weight (ug), mt = cp class of taxon t.
#       0.1  = C as fraction of fresh weight (20% dry weight x 52% C)
#       0.273 = 12/44, C from mass of CO2 evolved
#       0.75  = allometric exponent for respiration
#   Enrichment / structure guilds per Ferris et al. (2001):
#       enrichment = Ba1, Fu2         structure = Ba3-5, Fu3-5, Om3-5, Ca2-5
ANDRASSY_DIVISOR = 1.6e6
MF = {"production_c": 0.1, "respiration_c": 0.273, "resp_exp": 0.75}


def biomass(counts, taxa):
    """Andrassy (1956) individual fresh weight (ug) and per-sample biomass.

    Needs length_um and diameter_um in the taxa sheet; alternatively supply
    length_um and de Man's ratio a_ratio (a = L / max body diameter)."""
    L = pd.to_numeric(taxa.get("length_um"), errors="coerce").reindex(counts.index)
    if "diameter_um" in taxa.columns:
        D = pd.to_numeric(taxa["diameter_um"], errors="coerce").reindex(counts.index)
    elif "a_ratio" in taxa.columns:
        D = L / pd.to_numeric(taxa["a_ratio"], errors="coerce").reindex(counts.index)
    else:
        raise ValueError("taxa sheet needs diameter_um (or a_ratio) for biomass.")
    if L is None or L.isna().all():
        raise ValueError("taxa sheet needs length_um for biomass.")
    w = (L * D ** 2) / ANDRASSY_DIVISOR
    ok = w.notna() & (w > 0)
    c = counts.fillna(0).astype(float)
    bm = c.loc[ok].mul(w[ok], axis=0)
    out = pd.DataFrame({"individual_weight_ug": w,
                        "total_individuals": c.sum(axis=1),
                        "total_biomass_ug": bm.sum(axis=1).reindex(counts.index)})
    tot = out["total_biomass_ug"].sum()
    out["relative_biomass_pct"] = out["total_biomass_ug"] / tot * 100 if tot else 0.0
    out.attrs["missing_measurements"] = sorted(counts.index[~ok])
    out.attrs["per_sample_biomass"] = bm
    return out


def footprints(counts, taxa, include_fu2=False):
    """Metabolic footprints after Ferris (2010). Returns per-sample footprints."""
    bm = biomass(counts, taxa)
    w = bm["individual_weight_ug"]
    cp = pd.to_numeric(taxa["cp"], errors="coerce").reindex(counts.index)
    tr = taxa["trophic"].reindex(counts.index)
    per_ind = MF["production_c"] * w / cp + MF["respiration_c"] * w ** MF["resp_exp"]
    fp = counts.fillna(0).astype(float).mul(per_ind, axis=0)

    # NOTE: the enrichment INDEX (EI) uses Ba1 + Fu2, but reproducing Table 4 of
    # Ferris (2010) requires the enrichment FOOTPRINT to be Ba1 only:
    #   sample A -> Ba1 only gives 3.66 (paper: 4) and functional 143.5 (paper: 144)
    #               Ba1+Fu2 gives 5.29 and functional 207.6  -- does not match
    #   sample B -> Ba1 only gives 91.5 (paper: 91) and functional 370.1 (paper: 371)
    # Default follows the paper's own worked example. Set include_fu2=True to
    # use the EI guild set instead, and say which you used in your methods.
    enrich = (tr == "BF") & (cp == 1)
    if include_fu2:
        enrich = enrich | ((tr == "FF") & (cp == 2))
    struct = (((tr.isin(["BF", "FF", "OM"])) & (cp >= 3))
              | ((tr == "PR") & (cp >= 2)))

    Fe, Fs = fp.loc[enrich].sum(axis=0), fp.loc[struct].sum(axis=0)
    out = pd.DataFrame({
        "composite_footprint": fp.sum(axis=0),
        "enrichment_footprint": Fe,
        "structure_footprint": Fs,
        "functional_footprint": Fe * Fs / 2,     # area of the rhomboid, Ferris (2010)
        "herbivore_footprint": fp.loc[tr == "PP"].sum(axis=0),
        "bacterial_footprint": fp.loc[tr == "BF"].sum(axis=0),
        "fungal_footprint": fp.loc[tr == "FF"].sum(axis=0),
        "predator_footprint": fp.loc[tr == "PR"].sum(axis=0),
        "omnivore_footprint": fp.loc[tr == "OM"].sum(axis=0),
    })
    out.attrs["missing_measurements"] = bm.attrs["missing_measurements"]
    return out


# ---------------------------------------------------------------- templates
def write_template(path="NEMO-NEMA_input_template.xlsx"):
    readme = pd.DataFrame({"NEMO-NEMA input template": [
        "Fill the sheets below. Do not rename sheets or first-row headers.",
        "",
        "counts : taxa in rows, samples in columns, individuals recovered.",
        "         Use ONE consistent soil basis for every sample (e.g. per 200 cc).",
        "         Blank or 0 = absent. Numbers only.",
        "",
        "taxa   : one row per taxon in 'counts'; names must match exactly.",
        "         trophic = PP/BF/FF/OM/PR (Yeates et al. 1993)",
        "         cp      = 1-5 (Bongers 1990; Bongers & Bongers 1998)",
        "         source  = WRITE THE REFERENCE. Nothing is pre-filled on purpose:",
        "                   an unverified c-p table corrupts every index silently.",
        "         length_um / diameter_um = OPTIONAL, for biomass (Andrassy 1956).",
        "",
        "samples: one row per sample column; 'group' = the factor you compare.",
        "",
        "Delete the example rows before entering your own data."]})
    counts = pd.DataFrame({"Taxon": ["Helicotylenchus dihystera",
                                     "Tylenchorhynchus sp.", "Acrobeloides sp."],
                           "S01": [41, 7, 120], "S02": [12, 0, 96],
                           "S03": [0, 15, 143], "S04": [8, 3, 88]})
    taxa = pd.DataFrame({"Taxon": ["Helicotylenchus dihystera",
                                   "Tylenchorhynchus sp.", "Acrobeloides sp."],
                         "trophic": ["PP", "PP", "BF"], "cp": [3, 3, 2],
                         "length_um": [None] * 3, "diameter_um": [None] * 3,
                         "source": ["<enter your reference>"] * 3})
    samples = pd.DataFrame({"Sample": ["S01", "S02", "S03", "S04"],
                            "group": ["Untreated"] * 2 + ["Amended"] * 2,
                            "locality": ["Site A"] * 4, "crop": ["Rice"] * 4,
                            "replicate": [1, 2, 1, 2]})
    with pd.ExcelWriter(path, engine="openpyxl") as x:
        readme.to_excel(x, sheet_name="README", index=False)
        counts.to_excel(x, sheet_name="counts", index=False)
        taxa.to_excel(x, sheet_name="taxa", index=False)
        samples.to_excel(x, sheet_name="samples", index=False)
    return path


def write_demo(path="NEMO-NEMA_demo_data.xlsx"):
    rng = np.random.default_rng(7)
    tx = [("Meloidogyne graminicola", "PP", 3), ("Hirschmanniella oryzae", "PP", 3),
          ("Helicotylenchus dihystera", "PP", 3), ("Tylenchorhynchus mashhoodi", "PP", 3),
          ("Hoplolaimus indicus", "PP", 3), ("Xiphinema basiri", "PP", 5),
          ("Acrobeloides sp.", "BF", 2), ("Rhabditis sp.", "BF", 1),
          ("Cephalobus sp.", "BF", 2), ("Aphelenchus avenae", "FF", 2),
          ("Aporcelaimellus sp.", "OM", 5), ("Mononchus sp.", "PR", 4)]
    grp = ["Untreated"] * 4 + ["FYM amended"] * 4 + ["Nematicide"] * 4
    smp = [f"S{i:02d}" for i in range(1, 13)]
    base = {"Untreated": [180, 90, 60, 40, 25, 6, 90, 30, 70, 45, 12, 5],
            "FYM amended": [70, 40, 55, 35, 20, 9, 260, 140, 190, 120, 30, 14],
            "Nematicide": [25, 12, 18, 10, 6, 1, 40, 12, 30, 20, 3, 1]}
    data = {s: rng.poisson(np.array(base[g], float) * rng.uniform(.75, 1.25))
            for s, g in zip(smp, grp)}
    counts = pd.DataFrame(data, index=[t[0] for t in tx])
    taxa = pd.DataFrame({"Taxon": [t[0] for t in tx], "trophic": [t[1] for t in tx],
                         "cp": [t[2] for t in tx],
                         "source": ["SYNTHETIC DEMO — not a real assignment"] * len(tx)})
    samples = pd.DataFrame({"Sample": smp, "group": grp,
                            "replicate": [1, 2, 3, 4] * 3})
    with pd.ExcelWriter(path, engine="openpyxl") as x:
        counts.to_excel(x, sheet_name="counts", index_label="Taxon")
        taxa.to_excel(x, sheet_name="taxa", index=False)
        samples.to_excel(x, sheet_name="samples", index=False)
    return path


# ---------------------------------------------------------------- CLI
def run_cli(path, pv_divisor=1.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False})
    COL = {"PP": "#C44E52", "BF": "#4C72B0", "FF": "#DD8452",
           "OM": "#55A868", "PR": "#8172B3"}

    counts, taxa, samples = load(path)
    probs = validate(counts, taxa, samples)
    if probs:
        print("INPUT PROBLEMS — nothing computed:")
        for p in probs:
            print("  •", p)
        return
    n = norton(counts, pv_divisor).join(taxa[["trophic", "cp"]])
    d, f, ex = diversity(counts), faunal(counts, taxa), extra_indices(counts, taxa)
    nsh = None
    ab = counts.fillna(0).astype(float).sum(axis=0)
    try:
        nsh = nsh_index(f[["MI","PPI","EI","SI","BI","CI"]], ab)
    except Exception as exc:
        print("NSH skipped:", exc)
    fps = None
    if "length_um" in taxa.columns and taxa["length_um"].notna().any():
        try:
            fps = footprints(counts, taxa)
        except Exception as exc:
            print("Footprints skipped:", exc)
    bc = bray_curtis(counts)
    sc, var = pca(counts)
    g = samples["group"] if samples is not None and "group" in samples else None

    with pd.ExcelWriter("NEMO-NEMA_results.xlsx", engine="openpyxl") as x:
        n.to_excel(x, sheet_name="Norton_community")
        d.to_excel(x, sheet_name="Diversity")
        f.to_excel(x, sheet_name="Faunal_indices")
        ex.to_excel(x, sheet_name="Extra_indices")
        bc.to_excel(x, sheet_name="BrayCurtis")
        sc.to_excel(x, sheet_name="PCA_scores")
        if nsh is not None:
            nsh.to_excel(x, sheet_name="NSH_index")
        if fps is not None:
            fps.to_excel(x, sheet_name="Metabolic_footprints")
            biomass(counts, taxa).to_excel(x, sheet_name="Biomass")
        if g is not None:
            compare(d.join(f), g).to_excel(x, sheet_name="Group_tests")

    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    a = ax[0, 0]
    nn = n.sort_values("prominence_value")
    a.barh(range(len(nn)), nn["prominence_value"],
           color=[COL[t] for t in nn["trophic"]])
    a.set_yticks(range(len(nn)))
    a.set_yticklabels(nn.index, fontsize=7)
    a.set_xlabel("Prominence value")
    a.set_title("A. Ranked prominence value", loc="left", fontweight="bold")

    a = ax[0, 1]
    a.scatter(n["absolute_frequency_pct"], n["mean_density"],
              s=n["prominence_value"] / n["prominence_value"].max() * 400 + 30,
              c=[COL[t] for t in n["trophic"]], alpha=.75, edgecolor="w")
    a.set_xlabel("Absolute frequency (%)")
    a.set_ylabel("Mean density")
    a.set_title("B. Frequency vs density (bubble = PV)", loc="left", fontweight="bold")

    a = ax[1, 0]
    comp = f[[c for c in f.columns if c.startswith("pct_")]]
    bot = np.zeros(len(comp))
    for c in comp.columns:
        k = c.split("_")[1]
        a.bar(range(len(comp)), comp[c], bottom=bot, color=COL[k], label=TROPHIC[k])
        bot += comp[c].values
    a.set_xticks(range(len(comp)))
    a.set_xticklabels(comp.index, rotation=90, fontsize=7)
    a.set_ylabel("% of community")
    a.legend(fontsize=7, frameon=False, ncol=2)
    a.set_title("C. Trophic composition", loc="left", fontweight="bold")

    a = ax[1, 1]
    a.scatter(f["EI"], f["SI"], s=80, color="#4C72B0", edgecolor="w")
    a.axhline(50, ls=":", c="grey")
    a.axvline(50, ls=":", c="grey")
    a.set_xlim(0, 100)
    a.set_ylim(0, 100)
    a.set_xlabel("Enrichment index (EI)")
    a.set_ylabel("Structure index (SI)")
    for xx, yy, l in [(10, 93, "Maturing"), (72, 93, "Structured"),
                      (10, 4, "Degraded"), (72, 4, "Disturbed")]:
        a.text(xx, yy, l, color="grey", fontsize=8)
    a.set_title("D. Faunal profile (Ferris et al. 2001)", loc="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig("NEMO-NEMA_figures.png", dpi=170)

    print(f"{counts.shape[0]} taxa x {counts.shape[1]} samples")
    print("\n--- Norton community analysis ---")
    print(n.round(2).to_string())
    print("\n--- Faunal indices ---")
    print(f.round(2).to_string())
    if nsh is not None:
        print("\n--- Nematode Soil Health index (Ghaderi et al. 2025) ---")
        print(nsh.to_string())
        print("  NSH <15 degraded | 15-24 moderate | >25 well-functioning")
        print("  Abundance must be per 100 g dry soil for the MF band to be valid.")
    print("\nWrote NEMO-NEMA_results.xlsx and NEMO-NEMA_figures.png")


# ---------------------------------------------------------------- Streamlit
def _nn_footer(st):
    st.divider()
    st.markdown(
        "<div style='font-size:0.87rem; line-height:1.65; color:#3E4A44'>"
        "<b>Developed by</b><br>"
        "Ashish Kumar Singh &middot; Kavita Jain &middot; "
        "Vishal Singh Somvanshi &middot; Rashid Pervez &middot; "
        "Anil Sirohi &middot; Pankaj<br>"
        "<span style='color:#5B6B62'>Division of Nematology, ICAR-Indian "
        "Agricultural Research Institute, New Delhi 110012</span></div>",
        unsafe_allow_html=True)


def run_app():
    import streamlit as st
    from nemonema_manual import MANUAL
    import nemonema_report as rpt
    import matplotlib.pyplot as plt
    import nemonema_plots as npl
    st.set_page_config(page_title="NEMO-NEMA", page_icon="\U0001F52C", layout="wide")

    st.markdown("""<style>
      .block-container {padding-top: 2.2rem; max-width: 1500px;}
      .stTabs [data-baseweb="tab-list"] {gap: 2px; border-bottom: 1px solid #DED6C6;}
      .stTabs [data-baseweb="tab"] {height: 42px; padding: 0 16px;
        background: transparent; border-radius: 6px 6px 0 0; font-size: 0.93rem;}
      .stTabs [aria-selected="true"] {background: #F1ECE1;
        border-bottom: 2px solid #2E7D5B; color: #17453A;}
      .nn-rule {height:3px; width:74px;
        background:linear-gradient(90deg,#17453A,#8CBB72);
        border-radius:2px; margin:2px 0 14px 0;}
      .nn-sub {color:#5B6B62; font-size:0.92rem; line-height:1.5;}
      .nn-cp {display:inline-block; width:9px; height:9px; border-radius:50%;
        margin-right:5px; vertical-align:middle;}
    </style>""", unsafe_allow_html=True)
    import os as _os
    _h = st.columns([1, 9])
    with _h[0]:
        if _os.path.exists("icar_logo.png"):
            st.image("icar_logo.png", width=92)
    with _h[1]:
        st.markdown("<h1 style='margin-bottom:0'>NEMO-NEMA</h1>"
                    "<div class='nn-rule'></div>"
                    "<div class='nn-sub'><b>N</b>ematode <b>E</b>cological "
                    "<b>M</b>etrics and <b>O</b>rdination for <b>N</b>ematode "
                    "<b>E</b>cosystem <b>M</b>onitoring and <b>A</b>ssessment"
                    "<br><span style='color:#8A9690'>"
                    + "".join(f"<span class='nn-cp' style='background:{c}'></span>"
                              for c in ["#D55E00","#E69F00","#F0E442","#56B4E9","#17453A"])
                    + "c-p 1 to 5 &nbsp;&middot;&nbsp; version 1.0</span></div>",
                    unsafe_allow_html=True)
    st.write("")

    with st.sidebar:
        st.header("1. Your data")
        up = st.file_uploader("Excel workbook (.xlsx)", type=["xlsx"])
        st.download_button("Download blank template",
                           data=_template_bytes(), file_name="NEMO-NEMA_template.xlsx",
                           mime="application/vnd.openxmlformats-officedocument."
                                "spreadsheetml.sheet")
        st.download_button("Download demo dataset",
                           data=_demo_bytes(), file_name="NEMO-NEMA_demo.xlsx",
                           mime="application/vnd.openxmlformats-officedocument."
                                "spreadsheetml.sheet")
        st.header("2. Figure style")
        pal = st.selectbox("Colour palette", list(npl.PALETTES.keys()))
        st.markdown(" ".join(
            f'<span style="display:inline-block;width:22px;height:14px;'
            f'background:{c};border:1px solid #999;margin-right:2px"></span>'
            for c in npl.palette(pal, 8)), unsafe_allow_html=True)
        errtype = st.radio("Error bars", ["SE", "SD", "CI95"], horizontal=True)
        fsize = st.slider("Font size", 8, 16, 10)
        dpi = st.select_slider("Figure resolution (dpi)", [100, 150, 300, 600], 300)
        st.header("3. Options")
        pv = st.radio("Prominence value", ["D x sqrt(F)", "D x sqrt(F) / 10"],
                      help="Both are cited as Norton (1978). State which you used.")
        basis = st.selectbox("Counts are expressed per",
                             ["100 g dry soil", "200 cc soil", "250 g soil", "other"],
                             help="NSH abundance bands assume per 100 g dry soil.")

    if up is None:
        import os as _os
        if _os.path.exists("micrograph.jpg"):
            _c = st.columns([3, 2])
            with _c[1]:
                st.image("micrograph.jpg",
                         caption="Soil nematode, 4x. Division of Nematology, "
                                 "ICAR-IARI.")
            _tgt = _c[0]
        else:
            _tgt = st.container()
        with _tgt:
            st.info("Upload a workbook to begin, or download the demo dataset from the "
                    "sidebar to see what the output looks like.")
            st.markdown("""
    **Your workbook needs three sheets**

    | Sheet | Contents |
    |---|---|
    | `counts` | taxa in rows, samples in columns, individuals recovered |
    | `taxa` | one row per taxon: `trophic` (PP/BF/FF/OM/PR), `cp` (1-5), `source` |
    | `samples` | one row per sample: `group` = the factor being compared |

    Optional `length_um` and `diameter_um` columns in `taxa` unlock biomass and
    metabolic footprints.
    """)
        _nn_footer(st)
        return

    try:
        counts, taxa, samples = load(up)
    except Exception as exc:
        st.error(f"Could not read that workbook: {exc}")
        return

    problems = validate(counts, taxa, samples)
    if problems:
        st.error("Input problems - nothing has been computed:")
        for p in problems:
            st.write("-", p)
        return

    st.success(f"{counts.shape[0]} taxa x {counts.shape[1]} samples "
               f"({int(counts.fillna(0).values.sum()):,} individuals)")

    divisor = 10.0 if pv.endswith("/ 10") else 1.0
    nrt = norton(counts, divisor).join(taxa[["trophic", "cp"]])
    dv, fn, ex = diversity(counts), faunal(counts, taxa), extra_indices(counts, taxa)
    grp = samples["group"] if samples is not None and "group" in samples else None

    ab = counts.fillna(0).astype(float).sum(axis=0)
    try:
        nsh = nsh_index(fn[["MI", "PPI", "EI", "SI", "BI", "CI"]], ab)
    except Exception:
        nsh = None
    fps = None
    if "length_um" in taxa.columns and taxa["length_um"].notna().any():
        try:
            fps = footprints(counts, taxa)
        except Exception as exc:
            st.warning(f"Footprints unavailable: {exc}")

    tabs = st.tabs(["Community", "Diversity", "Faunal analysis",
                    "Soil health (NSH)", "Footprints", "Statistics", "Download", "Report", "Manual"])

    with tabs[0]:
        st.subheader("Norton (1978) community descriptors")
        st.dataframe(nrt.round(3), use_container_width=True)
        st.bar_chart(nrt["prominence_value"])

    with tabs[1]:
        st.subheader("Per-sample diversity")
        st.dataframe(dv.round(3), use_container_width=True)
        m = st.selectbox("Metric to plot", ["shannon_H", "simpson_1minusD",
                                            "pielou_J", "margalef_d",
                                            "richness_S", "total_N"])
        if grp is None:
            st.bar_chart(dv[m])
        else:
            npl.apply_style(fsize)
            tk = npl.tukey_posthoc(dv[m], grp)
            letters = npl.compact_letters(tk) if not tk.empty else None
            fig, ax = plt.subplots(figsize=(7, 4.5))
            npl.bar_with_error(dv[m], grp, err=errtype, pal=pal, ylabel=m,
                               title=m, letters=letters, ax=ax)
            st.pyplot(fig, use_container_width=False)
            st.download_button("Download figure (PNG)", _fig_bytes(fig, dpi),
                               f"{m}.png", "image/png")
            plt.close(fig)

    with tabs[2]:
        st.subheader("Maturity and food-web indices")
        st.caption("MI, MI(2-5), sigmaMI, PPI: Bongers (1990), Bongers & Bongers "
                   "(1998). EI, SI, CI, BI: Ferris et al. (2001). Trophic groups: "
                   "Yeates et al. (1993).")
        st.dataframe(fn.round(3), use_container_width=True)
        st.scatter_chart(fn.reset_index(), x="EI", y="SI")
        st.caption("Faunal profile: upper right = structured, lower left = degraded.")
        st.dataframe(ex.round(3), use_container_width=True)

    with tabs[3]:
        st.subheader("Nematode Soil Health index")
        st.caption("Ghaderi et al. (2025) Eur J Soil Sci 76:e70149. "
                   "Range 8-32. Below 15 degraded, 15-24 moderate, above 25 "
                   "well-functioning.")
        if nsh is None:
            st.info("NSH needs MI, PPI, EI, SI, BI and CI.")
        else:
            if basis != "100 g dry soil":
                st.warning(f"Your counts are per {basis}. The NSH abundance band "
                           "assumes per 100 g dry soil, so the MF score - and the "
                           "total - will be wrong until you convert.")
            st.dataframe(nsh, use_container_width=True)
            if grp is not None:
                npl.apply_style(fsize)
                tk = npl.tukey_posthoc(nsh["NSH"], grp)
                fig, ax = plt.subplots(figsize=(7, 4.5))
                npl.bar_with_error(nsh["NSH"], grp, err=errtype, pal=pal,
                                   ylabel="NSH index", title="Soil health",
                                   letters=npl.compact_letters(tk)
                                   if not tk.empty else None, ax=ax)
                ax.axhline(15, ls="--", c="#C44E52", lw=1)
                ax.axhline(25, ls="--", c="#55A868", lw=1)
                st.pyplot(fig, use_container_width=False)
                plt.close(fig)
            else:
                st.bar_chart(nsh["NSH"])
            st.divider()
            st.subheader("Calibration check")
            st.dataframe(npl.subscore_contributions(nsh).round(2),
                         use_container_width=True)
            if grp is not None:
                disc = npl.nsh_discrimination(nsh, grp)
                c1, c2, c3 = st.columns(3)
                c1.metric("Scale used", f"{disc.get('pct_of_scale_used', 0):.0f}%")
                c2.metric("eta squared",
                          f"{disc.get('eta_squared', float('nan')):.2f}")
                c3.metric("Categories",
                          f"{disc.get('n_categories_occupied', 0)} of 3")
            st.caption("Thresholds were calibrated on Australian and European "
                       "datasets; the authors state further calibration is needed. "
                       "Interpret with care outside those systems.")

    with tabs[4]:
        st.subheader("Biomass and metabolic footprints")
        st.caption("Andrassy (1956) biomass; Ferris (2010) footprints.")
        if fps is None:
            st.info("Add length_um and diameter_um columns to your taxa sheet to "
                    "compute biomass and metabolic footprints.")
        else:
            miss = fps.attrs.get("missing_measurements") or []
            if miss:
                st.warning(f"Excluded for lack of measurements: {', '.join(miss)}")
            st.dataframe(fps.round(2), use_container_width=True)

    with tabs[5]:
        st.subheader("Group comparisons")
        if grp is None:
            st.info("Add a `group` column to the samples sheet to enable this.")
        else:
            combined = dv.join(fn)
            if nsh is not None:
                combined = combined.join(nsh[["NSH"]])
            st.dataframe(compare(combined, grp).round(4), use_container_width=True)
            st.caption("ANOVA and Kruskal-Wallis are both shown because nematode "
                       "counts are often non-normal. Choose one a priori and "
                       "correct for multiple testing.")

    with tabs[6]:
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as x:
            nrt.to_excel(x, sheet_name="Norton_community")
            dv.to_excel(x, sheet_name="Diversity")
            fn.to_excel(x, sheet_name="Faunal_indices")
            ex.to_excel(x, sheet_name="Extra_indices")
            if nsh is not None:
                nsh.to_excel(x, sheet_name="NSH_index")
            if fps is not None:
                fps.to_excel(x, sheet_name="Metabolic_footprints")
            bray_curtis(counts).to_excel(x, sheet_name="BrayCurtis")
            pca(counts)[0].to_excel(x, sheet_name="PCA_scores")
            if grp is not None:
                compare(dv.join(fn), grp).to_excel(x, sheet_name="Group_tests")
            pd.DataFrame({"setting": ["prominence_value", "soil_basis",
                                      "n_samples", "n_taxa"],
                          "value": [pv, basis, counts.shape[1], counts.shape[0]]
                          }).to_excel(x, sheet_name="Methods_log", index=False)
        st.download_button("Download all results (.xlsx)", buf.getvalue(),
                           file_name="NEMO-NEMA_results.xlsx",
                           mime="application/vnd.openxmlformats-officedocument."
                                "spreadsheetml.sheet")

    with tabs[7]:
        st.subheader("Full report")
        st.write("One PDF with every table, the figures, the interpretation, "
                 "the methods log and the references.")
        try:
            _s = compare(dv.join(fn), grp) if grp is not None else None
            _p = rpt.build_report(counts, taxa, samples, nrt, dv, fn, ex, nsh,
                                  fps, _s, plots=npl, palette=pal, err=errtype,
                                  settings={"prominence_value": pv, "soil_basis": basis})
            st.download_button("Download full report (PDF)", _p,
                               "NEMO-NEMA_report.pdf", "application/pdf")
        except Exception as e:
            st.warning(f"Report could not be built: {e}")

    with tabs[8]:
        st.subheader("User manual")
        st.markdown(MANUAL)

    _nn_footer(st)

    st.divider()
    st.caption("Verify trophic and c-p assignments against primary sources before "
               "publishing: every maturity and food-web index depends on them.")


def _fig_bytes(fig, dpi=300):
    b = io.BytesIO()
    fig.savefig(b, format="png", dpi=dpi, bbox_inches="tight")
    return b.getvalue()


def _template_bytes():
    import tempfile, os
    p = os.path.join(tempfile.mkdtemp(), "t.xlsx")
    write_template(p)
    return open(p, "rb").read()


def _demo_bytes():
    import tempfile, os
    p = os.path.join(tempfile.mkdtemp(), "d.xlsx")
    write_demo(p)
    return open(p, "rb").read()


def _under_streamlit():
    try:
        from streamlit.runtime.scriptrunner import get_script_run_ctx
        return get_script_run_ctx() is not None
    except Exception:
        return False


if _under_streamlit():
    run_app()
elif __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "help"
    if cmd == "template":
        print("written:", write_template())
    elif cmd == "demo":
        print("written:", write_demo())
    elif cmd == "run":
        run_cli(sys.argv[2], 10.0 if "--pv10" in sys.argv else 1.0)
    else:
        print(HELP)
