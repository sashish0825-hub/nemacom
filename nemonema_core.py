"""
nemonema_core.py — calculation engine for NEMO-NEMA.

Nematode Ecological Metrics and Ordination for
Nematode Ecosystem Monitoring and Assessment

Deliberately UI-free so every function can be tested against hand-worked values.

VERIFIED AGAINST PRIMARY SOURCES
  Ferris, Bongers & de Goede (2001) Appl Soil Ecol 18:13-29, Fig. 1
      guild weights and the b/e/s components — checked value by value
  Ferris (2010) Eur J Soil Biol 46:97-104, Table 4
      metabolic footprints — reproduces all three worked samples
  Ghaderi et al. (2025) Eur J Soil Sci 76:e70149, Supporting Information 1
      NSH scoring — reproduces all 294 published samples exactly

NOT VERIFIED
  Norton (1978): two prominence-value variants circulate, both attributed to
  that book. The variant is a user choice and is recorded in the methods log.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

VERSION = "NEMO-NEMA 1.1"

TROPHIC = {"PP": "plant parasite", "BF": "bacterivore", "FF": "fungivore",
           "OM": "omnivore", "PR": "predator"}

# Ferris et al. (2001) Fig. 1. Keys are (trophic group, c-p class); values are
# weighting coefficients. Guild abundances are supplied separately.
W = {
    "enrichment": {("BF", 1): 3.2, ("FF", 2): 0.8},
    "basal":      {("BF", 2): 0.8, ("FF", 2): 0.8},
    "structure":  {("BF", 3): 1.8, ("BF", 4): 3.2, ("BF", 5): 5.0,
                   ("FF", 3): 1.8, ("FF", 4): 3.2, ("FF", 5): 5.0,
                   ("OM", 3): 1.8, ("OM", 4): 3.2, ("OM", 5): 5.0,
                   ("PR", 2): 0.8, ("PR", 3): 1.8, ("PR", 4): 3.2, ("PR", 5): 5.0},
}

# Ferris et al. (2001) Fig. 1 and Table 5. x = SI, y = EI, split at 50/50.
QUADRANTS = {
    "A": ("Disturbed", "high EI, low SI — enriched but structurally simple; "
                       "bacterial decomposition channel"),
    "B": ("Structured", "high EI, high SI — enriched and structured; low to "
                        "moderate disturbance"),
    "C": ("Stable", "low EI, high SI — undisturbed and mature; fungal "
                    "decomposition channel"),
    "D": ("Degraded", "low EI, low SI — depleted and structurally simple; "
                      "stressed"),
}


# --------------------------------------------------------------------------
# Loading and validation
# --------------------------------------------------------------------------
def load(path_or_buffer):
    """Read a NEMO-NEMA workbook. Returns (counts, taxa, samples)."""
    def clean(df):
        df = df[df.index.notna()].copy()
        df.index = df.index.astype(str).str.strip()
        return df[~df.index.isin(["", "nan", "None"])
                  & ~df.index.str.startswith(("^", "<-"))]

    xl = pd.ExcelFile(path_or_buffer)
    if "counts" not in xl.sheet_names:
        raise ValueError(
            "No sheet named 'counts'. The workbook needs three sheets: "
            "counts, taxa and samples. Sheets found: "
            + ", ".join(xl.sheet_names))
    counts = clean(pd.read_excel(xl, "counts", index_col=0))
    counts.columns = [str(c).strip() for c in counts.columns]

    taxa = clean(pd.read_excel(xl, "taxa", index_col=0)) \
        if "taxa" in xl.sheet_names else pd.DataFrame(index=counts.index)
    taxa.columns = [str(c).strip().lower() for c in taxa.columns]

    samples = None
    if "samples" in xl.sheet_names:
        samples = clean(pd.read_excel(xl, "samples", index_col=0))
        samples.columns = [str(c).strip().lower() for c in samples.columns]
    return counts, taxa, samples


def validate(counts, taxa, samples=None) -> list:
    """Problems that must be fixed before any result can be trusted."""
    p = []
    if counts.empty:
        return ["The 'counts' sheet is empty."]
    bad = [c for c in counts.columns if not pd.api.types.is_numeric_dtype(counts[c])]
    if bad:
        p.append(f"Non-numeric sample columns in 'counts': {bad}")
    if (counts.fillna(0) < 0).any().any():
        p.append("Negative counts found in 'counts'.")
    if counts.index.duplicated().any():
        p.append(f"Duplicate taxon names: "
                 f"{sorted(counts.index[counts.index.duplicated()].unique())}")
    if "trophic" not in taxa.columns or "cp" not in taxa.columns:
        p.append("The 'taxa' sheet needs 'trophic' and 'cp' columns. "
                 "Use the Auto-assign tab to fill them from taxon names.")
        return p
    missing = sorted(set(counts.index) - set(taxa.index))
    if missing:
        p.append(f"Taxa in 'counts' with no row in 'taxa': {missing}")
    sub = taxa.loc[taxa.index.isin(counts.index)]
    bt = sub.loc[~sub["trophic"].isin(TROPHIC)]
    if len(bt):
        p.append(f"Unrecognised trophic codes for {list(bt.index)}. "
                 f"Use one of {sorted(TROPHIC)}.")
    cp = pd.to_numeric(sub["cp"], errors="coerce")
    bc = sub.loc[~cp.isin([1, 2, 3, 4, 5])]
    if len(bc):
        p.append(f"c-p must be an integer 1-5. Check: {list(bc.index)}")
    if samples is not None:
        un = sorted(set(counts.columns) - set(samples.index))
        if un:
            p.append(f"Sample columns with no row in 'samples': {un}")
    return p


# --------------------------------------------------------------------------
# 1. Norton (1978) community descriptors
# --------------------------------------------------------------------------
def norton(counts: pd.DataFrame, pv_divisor: float = 1.0) -> pd.DataFrame:
    c = counts.fillna(0).astype(float)
    m = c.shape[1]
    if m == 0:
        raise ValueError("No sample columns.")
    k = (c > 0).sum(axis=1)
    n = c.sum(axis=1)
    out = pd.DataFrame({
        "samples_positive": k.astype(int),
        "total_individuals": n,
        "absolute_frequency_pct": k / m * 100,
        "mean_density": n / m,
    })
    fs, ds = out["absolute_frequency_pct"].sum(), out["mean_density"].sum()
    out["relative_frequency_pct"] = (out["absolute_frequency_pct"] / fs * 100
                                     if fs else np.nan)
    out["relative_density_pct"] = (out["mean_density"] / ds * 100 if ds else np.nan)
    out["prominence_value"] = (out["mean_density"]
                               * np.sqrt(out["absolute_frequency_pct"]) / pv_divisor)
    out["importance_value"] = (out["relative_frequency_pct"]
                               + out["relative_density_pct"])
    # flag the case where frequency saturates and PV collapses onto density
    out.attrs["pv_uninformative"] = bool((out["absolute_frequency_pct"] == 100).all())
    return out.sort_values("prominence_value", ascending=False)


# --------------------------------------------------------------------------
# 2. Diversity
# --------------------------------------------------------------------------
def diversity(counts: pd.DataFrame) -> pd.DataFrame:
    """Per-sample diversity. Degenerate cases return NaN, never a misleading 0."""
    c = counts.fillna(0).astype(float)
    rows = {}
    for s in c.columns:
        x = c[s].values
        x = x[x > 0]
        N, S = x.sum(), len(x)
        if N == 0:
            rows[s] = dict(richness_S=0, total_N=0.0, shannon_H=np.nan,
                           simpson_D=np.nan, simpson_1minusD=np.nan,
                           pielou_J=np.nan, margalef_d=np.nan)
            continue
        p = x / N
        H = float(-(p * np.log(p)).sum())
        D = float((p ** 2).sum())
        rows[s] = dict(
            richness_S=int(S), total_N=float(N), shannon_H=H,
            simpson_D=D, simpson_1minusD=1 - D,
            # ln(1) = 0, so evenness is undefined at S = 1, not equal to 1
            pielou_J=(H / np.log(S)) if S > 1 else np.nan,
            margalef_d=((S - 1) / np.log(N)) if N > 1 else np.nan)
    return pd.DataFrame(rows).T


# --------------------------------------------------------------------------
# 3. Maturity and food-web indices
# --------------------------------------------------------------------------
def _guilds(counts, taxa):
    meta = taxa.reindex(counts.index)
    key = pd.MultiIndex.from_arrays(
        [meta["trophic"].values,
         pd.to_numeric(meta["cp"], errors="coerce").astype("Int64").values],
        names=["trophic", "cp"])
    g = counts.fillna(0).astype(float).copy()
    g.index = key
    return g.groupby(level=["trophic", "cp"]).sum()


def faunal(counts: pd.DataFrame, taxa: pd.DataFrame) -> pd.DataFrame:
    g = _guilds(counts, taxa)

    def wmean(sr):
        t = sr.sum()
        return float(sum(cp * v for (_, cp), v in sr.items()) / t) if t else np.nan

    rows = {}
    for s in g.columns:
        x = g[s]
        x = x[x > 0]
        if x.sum() == 0:
            rows[s] = {}
            continue
        free = x[[t != "PP" for t, _ in x.index]]
        pp = x[[t == "PP" for t, _ in x.index]]
        free25 = free[[cp >= 2 for _, cp in free.index]]

        b = sum(W["basal"].get(k, 0) * v for k, v in x.items())
        e = sum(W["enrichment"].get(k, 0) * v for k, v in x.items())
        st = sum(W["structure"].get(k, 0) * v for k, v in x.items())

        Ba1 = float(x.get(("BF", 1), 0.0))
        Fu2 = float(x.get(("FF", 2), 0.0))
        ci_den = 3.2 * Ba1 + 0.8 * Fu2

        bt = x.groupby(level="trophic").sum()
        tot = float(x.sum())
        Ba, Fu = float(bt.get("BF", 0)), float(bt.get("FF", 0))

        rows[s] = {
            "MI": wmean(free), "MI_2_5": wmean(free25), "sigma_MI": wmean(x),
            "PPI": wmean(pp),
            "PPI_MI_ratio": (wmean(pp) / wmean(free)) if wmean(free) else np.nan,
            "EI": 100 * e / (e + b) if (e + b) else np.nan,
            "SI": 100 * st / (st + b) if (st + b) else np.nan,
            # CI is undefined, not zero, when neither enrichment guild is present
            "CI": 100 * (0.8 * Fu2) / ci_den if ci_den else np.nan,
            "BI": 100 * b / (e + st + b) if (e + st + b) else np.nan,
            "NCR": Ba / (Ba + Fu) if (Ba + Fu) else np.nan,
            **{f"pct_{t}": 100 * float(bt.get(t, 0)) / tot for t in TROPHIC},
        }
    return pd.DataFrame(rows).T


def quadrant(ei, si) -> str:
    """Faunal profile quadrat after Ferris et al. (2001) Fig. 1 and Table 5."""
    if pd.isna(ei) or pd.isna(si):
        return ""
    if ei >= 50:
        return "B" if si >= 50 else "A"
    return "C" if si >= 50 else "D"


def extra_indices(counts: pd.DataFrame, taxa: pd.DataFrame) -> pd.DataFrame:
    c = counts.fillna(0).astype(float)
    tr = taxa["trophic"].reindex(c.index)
    cp = pd.to_numeric(taxa["cp"], errors="coerce").reindex(c.index)
    rows = {}
    for s in c.columns:
        x = c[s]
        bt = x.groupby(tr).sum()
        pp = float(bt.get("PP", 0))
        fb = float(bt.get("BF", 0)) + float(bt.get("FF", 0))
        tot = float(bt.sum())
        p = bt[bt > 0] / tot if tot else bt
        rows[s] = {
            "wasilewska_index": fb / pp if pp else np.nan,
            "trophic_diversity_H": float(-(p * np.log(p)).sum()) if tot else np.nan,
            "functional_guild_richness": len({(t, k) for t, k, v
                                              in zip(tr, cp, x) if v > 0}),
        }
    return pd.DataFrame(rows).T


# --------------------------------------------------------------------------
# 4. Nematode Soil Health index
# --------------------------------------------------------------------------
NSH_BANDS = {
    "MI":  [(1.0, 1.5, 1), (1.5, 2.0, 2), (2.0, 2.5, 3), (2.5, 3.0, 4), (3.0, 5.0, 5)],
    "PPI": [(2.0, 2.5, 4), (2.5, 3.0, 2), (3.0, 3.5, 1), (3.5, 5.0, 3)],
    "EI":  [(0, 20, 1), (20, 40, 2), (40, 60, 3), (60, 80, 4), (80, 100, 5)],
    "SI":  [(0, 20, 1), (20, 40, 2), (40, 60, 3), (60, 80, 4), (80, 100, 5)],
    "BI":  [(0, 20, 5), (20, 40, 4), (40, 60, 3), (60, 80, 2), (80, 100, 1)],
    "CI":  [(0, 20, 2), (20, 80, 3), (80, 100, 2)],
    "MF":  [(-np.inf, 200, 1), (200, 500, 2), (500, 1000, 3),
            (1000, 2000, 4), (2000, np.inf, 5)],
}
NSH_BIOMASS_BANDS = [(-np.inf, 140, 1), (140, 350, 2), (350, 700, 3),
                     (700, 1400, 4), (1400, np.inf, 5)]


def _band(value, bands):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    # guard band edges against floating-point residue such as 2.4999999999999996,
    # which occurs in the authors' own supplementary workbook
    value = round(float(value), 6)
    for i, (lo, hi, sc) in enumerate(bands):
        last = i == len(bands) - 1
        if (lo <= value < hi) or (last and lo <= value <= hi):
            return sc
    return np.nan


def nsh_index(faunal_tbl, abundance, mf_basis="abundance") -> pd.DataFrame:
    bands = dict(NSH_BANDS)
    if mf_basis == "biomass":
        bands["MF"] = NSH_BIOMASS_BANDS
    out = pd.DataFrame(index=faunal_tbl.index)
    for k in ("MI", "PPI", "EI", "SI", "BI", "CI"):
        out[f"{k}_score"] = [_band(v, bands[k]) for v in faunal_tbl[k]]
    ab = pd.Series(abundance).reindex(faunal_tbl.index)
    out["MF_score"] = [_band(v, bands["MF"]) for v in ab]
    out["NSH"] = out.sum(axis=1, skipna=False)
    out["interpretation"] = pd.cut(out["NSH"], [-np.inf, 15, 25, np.inf],
                                   labels=["degraded", "moderate",
                                           "well-functioning"], right=False)
    return out


# --------------------------------------------------------------------------
# 5. Biomass and metabolic footprints
# --------------------------------------------------------------------------
ANDRASSY_DIVISOR = 1.6e6
MF_COEF = {"production": 0.1, "respiration": 0.273, "exponent": 0.75}


def biomass(counts: pd.DataFrame, taxa: pd.DataFrame) -> pd.DataFrame:
    """Andrassy (1956): W = (L x D^2) / 1.6e6, W in ug, L and D in um."""
    L = pd.to_numeric(taxa.get("length_um"), errors="coerce").reindex(counts.index)
    if L is None or L.isna().all():
        raise ValueError("The taxa sheet needs 'length_um'.")
    if "diameter_um" in taxa.columns and taxa["diameter_um"].notna().any():
        D = pd.to_numeric(taxa["diameter_um"], errors="coerce").reindex(counts.index)
    elif "a_ratio" in taxa.columns:
        D = L / pd.to_numeric(taxa["a_ratio"], errors="coerce").reindex(counts.index)
    else:
        raise ValueError("The taxa sheet needs 'diameter_um' or 'a_ratio'.")
    w = (L * D ** 2) / ANDRASSY_DIVISOR
    ok = w.notna() & (w > 0)
    c = counts.fillna(0).astype(float)
    bm = c.loc[ok].mul(w[ok], axis=0)
    out = pd.DataFrame({"individual_weight_ug": w,
                        "total_individuals": c.sum(axis=1),
                        "total_biomass_ug": bm.sum(axis=1).reindex(counts.index)})
    tot = out["total_biomass_ug"].sum()
    out["relative_biomass_pct"] = out["total_biomass_ug"] / tot * 100 if tot else np.nan
    out.attrs["missing"] = sorted(counts.index[~ok])
    return out


def footprints(counts: pd.DataFrame, taxa: pd.DataFrame,
               include_fu2: bool = False) -> pd.DataFrame:
    """Ferris (2010): F = SUM Nt [0.1 (Wt/mt) + 0.273 Wt^0.75].

    Reproducing Table 4 of that paper requires the ENRICHMENT footprint to use
    Ba1 only, although the enrichment INDEX uses Ba1 + Fu2:
        sample A -> Ba1 only 3.66 (paper 4); Ba1+Fu2 gives 5.29
        sample B -> Ba1 only 91.5 (paper 91); functional 370.1 (paper 371)
    The default follows the paper's own worked example.
    """
    bm = biomass(counts, taxa)
    w = bm["individual_weight_ug"]
    cp = pd.to_numeric(taxa["cp"], errors="coerce").reindex(counts.index)
    tr = taxa["trophic"].reindex(counts.index)
    per_ind = (MF_COEF["production"] * w / cp
               + MF_COEF["respiration"] * w ** MF_COEF["exponent"])
    fp = counts.fillna(0).astype(float).mul(per_ind, axis=0)

    enrich = (tr == "BF") & (cp == 1)
    if include_fu2:
        enrich = enrich | ((tr == "FF") & (cp == 2))
    struct = ((tr.isin(["BF", "FF", "OM"])) & (cp >= 3)) | ((tr == "PR") & (cp >= 2))

    Fe, Fs = fp.loc[enrich].sum(axis=0), fp.loc[struct].sum(axis=0)
    out = pd.DataFrame({
        "composite_footprint": fp.sum(axis=0),
        "enrichment_footprint": Fe,
        "structure_footprint": Fs,
        "functional_footprint": Fe * Fs / 2,
        "herbivore_footprint": fp.loc[tr == "PP"].sum(axis=0),
        "bacterial_footprint": fp.loc[tr == "BF"].sum(axis=0),
        "fungal_footprint": fp.loc[tr == "FF"].sum(axis=0),
        "predator_footprint": fp.loc[tr == "PR"].sum(axis=0),
        "omnivore_footprint": fp.loc[tr == "OM"].sum(axis=0),
    })
    out.attrs["missing"] = bm.attrs["missing"]
    return out


# --------------------------------------------------------------------------
# 6. Multivariate
# --------------------------------------------------------------------------
def bray_curtis(counts: pd.DataFrame) -> pd.DataFrame:
    from scipy.spatial.distance import pdist, squareform
    d = squareform(pdist(counts.fillna(0).astype(float).T.values,
                         metric="braycurtis"))
    return pd.DataFrame(d, index=counts.columns, columns=counts.columns)


def pca(counts: pd.DataFrame, k: int = 2):
    """PCA on Hellinger-transformed relative abundances."""
    X = counts.fillna(0).astype(float).T
    H = np.sqrt(X.div(X.sum(axis=1).replace(0, np.nan), axis=0)).fillna(0)
    Xc = H - H.mean(axis=0)
    U, S, _ = np.linalg.svd(Xc.values, full_matrices=False)
    k = min(k, len(S))
    scores = pd.DataFrame(U[:, :k] * S[:k], index=H.index,
                          columns=[f"PC{i+1}" for i in range(k)])
    return scores, (S ** 2) / (S ** 2).sum() * 100
