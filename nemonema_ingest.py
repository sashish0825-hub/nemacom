"""
nemonema_ingest.py — Automatic Dataset Import & Standardisation Layer for NEMO-NEMA.

Drop-in module. Adds nothing to and removes nothing from the existing analysis
code: it only converts an arbitrary user upload into the schema the rest of the
application already consumes.

Governing principle: FORMAT CAN CHANGE. DATA CANNOT.

Public API
----------
    set_reference_taxa(iterable)        inject the app's own genus reference
    detect_dataset_structure(...)  ->   DatasetStructure
    map_columns(...)               ->   MappingResult
    convert_to_nemonema_schema(...)->   Standardized
    validate_conversion(...)       ->   ValidationReport
    generate_conversion_report(...)->   pandas.DataFrame
    build_standardized_workbook(...)->  bytes (.xlsx)
    render_import_ui(st)           ->   Standardized | None   (Streamlit screen)

Nothing here imports streamlit at module level, so the layer is testable headless.
"""

from __future__ import annotations

import io
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

INGEST_VERSION = "1.0"

# --------------------------------------------------------------------------
# 1. Reference vocabularies
# --------------------------------------------------------------------------

# Fallback genus list, used only when the host application has not injected its
# own reference via set_reference_taxa(). NEMO-NEMA ships a curated 188-genus
# reference; that one should always take precedence, because it is the list the
# trophic/c-p assignments are keyed to.
_FALLBACK_GENERA = """
Meloidogyne Heterodera Globodera Cactodera Punctodera Meloidodera Pratylenchus
Radopholus Hirschmanniella Zygotylenchus Pratylenchoides Nacobbus Ditylenchus
Anguina Subanguina Nothotylenchus Tylenchorhynchus Merlinius Quinisulcius
Bitylenchus Amplimerlinius Sauertylenchus Trophurus Geocenamus Helicotylenchus
Rotylenchus Pararotylenchus Scutellonema Hoplolaimus Aorolaimus Peltamigratus
Rotylenchulus Tylenchulus Belonolaimus Dolichodorus Xiphinema Paralongidorus
Longidorus Xiphidorus Trichodorus Paratrichodorus Criconemoides Mesocriconema
Criconemella Hemicriconemoides Discocriconemella Ogma Xenocriconemella Bakernema
Hemicycliophora Caloosia Paratylenchus Gracilacus Cacopaurus Tylenchus Filenchus
Malenchus Coslenchus Aglenchus Basiria Boleodorus Psilenchus Neopsilenchus
Ecphyadophora Ottolenchus Miculenchus Tylodorus Aphelenchus Aphelenchoides
Paraphelenchus Seinura Bursaphelenchus Rhabditis Rhabditella Mesorhabditis
Protorhabditis Pelodera Diploscapter Oscheius Caenorhabditis Steinernema
Heterorhabditis Panagrolaimus Panagrellus Cephalobus Eucephalobus
Heterocephalobus Acrobeles Acrobeloides Cervidellus Chiloplacus Zeldia Seleborca
Stegelletina Teratocephalus Metateratocephalus Euteratocephalus Plectus
Anaplectus Wilsonema Chronogaster Rhabdolaimus Prismatolaimus Bastiania Alaimus
Amphidelus Paramphidelus Monhystera Eumonhystera Geomonhystera Achromadora
Prodesmodora Ethmolaimus Aulolaimus Bunonema Diplogaster Diplogasteroides
Mononchoides Butlerius Fictor Mononchus Clarkus Coomansus Mylonchulus
Prionchulus Anatonchus Iotonchus Miconchus Tripyla Trischistoma Tobrilus Ironus
Dorylaimus Mesodorylaimus Prodorylaimus Eudorylaimus Epidorylaimus Laimydorus
Aporcelaimus Aporcelaimellus Discolaimus Discolaimoides Thornia Labronema
Nygolaimus Aquatides Tylencholaimus Tylencholaimellus Dorylaimoides Leptonchus
Diphtherophora Belondira Axonchium Dorylaimellus Oxydirus Longidorella
Xiphinemella Nordia Pungentus Enchodelus Crassolabium Microdorylaimus
""".split()

_REFERENCE_TAXA: set = {g.lower() for g in _FALLBACK_GENERA}
_REFERENCE_IS_INJECTED = False


def set_reference_taxa(taxa: Iterable[str]) -> None:
    """Replace the fallback genus list with the host application's reference.

    Call once at start-up, e.g. from the 188-genus table in the Reference
    database tab. Genus or full binomial entries are both accepted.
    """
    global _REFERENCE_TAXA, _REFERENCE_IS_INJECTED
    cleaned = set()
    for t in taxa:
        if not isinstance(t, str):
            continue
        t = t.strip()
        if not t:
            continue
        cleaned.add(t.lower())
        cleaned.add(t.split()[0].lower())  # genus of a binomial
    if cleaned:
        _REFERENCE_TAXA = cleaned
        _REFERENCE_IS_INJECTED = True


def reference_status() -> str:
    n = len(_REFERENCE_TAXA)
    src = "host application" if _REFERENCE_IS_INJECTED else "built-in fallback"
    return f"{n} reference names loaded ({src})"


# Negative evidence. These are the sources of the Slovakia-as-a-taxon failure:
# place names, crops, treatments and soil variables are text, and text alone
# must never be sufficient evidence of a taxon.
_GAZETTEER = {
    # countries seen in nematode survey data
    "india", "slovakia", "czechia", "czech republic", "slovak republic", "poland",
    "germany", "france", "spain", "italy", "netherlands", "belgium", "austria",
    "hungary", "romania", "bulgaria", "greece", "turkey", "iran", "iraq", "egypt",
    "kenya", "nigeria", "ghana", "south africa", "brazil", "argentina", "chile",
    "mexico", "usa", "united states", "canada", "china", "japan", "korea",
    "vietnam", "thailand", "philippines", "indonesia", "malaysia", "bangladesh",
    "pakistan", "nepal", "sri lanka", "bhutan", "myanmar", "australia",
    "new zealand", "uk", "united kingdom", "ireland", "portugal", "switzerland",
    "sweden", "norway", "denmark", "finland", "russia", "ukraine", "israel",
    # Indian states and UTs
    "andhra pradesh", "arunachal pradesh", "assam", "bihar", "chhattisgarh",
    "goa", "gujarat", "haryana", "himachal pradesh", "jharkhand", "karnataka",
    "kerala", "madhya pradesh", "maharashtra", "manipur", "meghalaya", "mizoram",
    "nagaland", "odisha", "orissa", "punjab", "rajasthan", "sikkim", "tamil nadu",
    "telangana", "tripura", "uttar pradesh", "uttarakhand", "west bengal",
    "delhi", "new delhi", "jammu", "kashmir", "ladakh", "puducherry", "chandigarh",
    "andaman", "nicobar", "lakshadweep", "dadra", "nagar haveli", "daman", "diu",
    # crops
    "rice", "wheat", "maize", "barley", "sorghum", "pearl millet", "finger millet",
    "chickpea", "pigeonpea", "lentil", "mungbean", "urdbean", "soybean",
    "groundnut", "mustard", "sunflower", "sesame", "cotton", "sugarcane", "jute",
    "potato", "tomato", "brinjal", "eggplant", "okra", "chilli", "capsicum",
    "onion", "garlic", "cabbage", "cauliflower", "carrot", "radish", "spinach",
    "cucumber", "pumpkin", "bottle gourd", "bitter gourd", "banana", "mango",
    "guava", "citrus", "grape", "apple", "papaya", "pomegranate", "kiwi",
    "kiwifruit", "tea", "coffee", "coconut", "arecanut", "rubber", "black pepper",
    "turmeric", "ginger", "cardamom", "fallow", "grassland", "forest", "orchard",
    "polyhouse", "greenhouse", "nursery",
    # treatments / management
    "control", "untreated", "check", "organic", "inorganic", "conventional",
    "integrated", "treated", "carbofuran", "fluopyram", "neem cake", "fym",
    "vermicompost", "compost", "biochar", "npk", "urea", "farmyard manure",
    "trichoderma", "pseudomonas", "purpureocillium", "paecilomyces", "bacillus",
    "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10",
    # soil / measurement variables
    "ph", "ec", "oc", "organic carbon", "moisture", "clay", "sand", "silt",
    "loam", "sandy loam", "clay loam", "silty clay", "texture", "bulk density",
    "nitrogen", "phosphorus", "potassium", "sulphur", "zinc", "iron", "manganese",
    "copper", "boron", "temperature", "rainfall", "depth", "yield",
}

# Column-name synonyms. Matching is done on a normalised form, so
# "Sample ID", "sample_id", "SampleID" and "SAMPLE-ID" all collapse together.
ROLE_SYNONYMS: Dict[str, List[str]] = {
    "sample_id": ["sample", "sampleid", "sample no", "sample number", "sample name",
                  "samplename", "sampleno", "id", "plot", "plot id", "plotid",
                  "site", "site id", "siteid", "code", "sample code", "samplecode",
                  "field id", "unit", "observation", "obs"],
    "group": ["group", "grouping", "category", "class", "system", "management",
              "agroecosystem", "ecosystem", "land use", "landuse", "cropping system",
              "habitat", "site type", "sitetype"],
    "treatment": ["treatment", "treat", "trt", "amendment", "management practice",
                  "practice", "regime", "dose", "application"],
    "replicate": ["replicate", "rep", "reps", "replication", "block", "r"],
    "taxon": ["taxon", "taxa", "genus", "genera", "species", "nematode",
              "nematode taxon", "nematode genus", "organism", "name",
              "scientific name", "identification", "identified as"],
    "count": ["count", "counts", "number", "no", "n", "abundance", "individuals",
              "individual", "population", "density", "value", "reading",
              "no of individuals", "number of individuals", "nematode count",
              "population density", "pop density", "juveniles", "j2"],
    "trophic_group": ["trophic group", "trophicgroup", "trophic", "feeding group",
                      "feeding type", "feeding habit", "guild", "functional group"],
    "cp": ["cp", "c p", "cp value", "cpvalue", "c p value", "colonizer persister",
           "colonizer persister value", "cp class", "maturity class"],
    "source": ["source", "reference", "citation", "authority", "assigned by",
               "basis"],
    "location": ["location", "place", "village", "district", "state", "region",
                 "country", "province", "zone", "gps", "latitude", "longitude",
                 "lat", "long", "lon", "altitude", "elevation", "locality"],
    "crop": ["crop", "host", "host plant", "hostplant", "plant", "cultivar host",
             "crop species"],
    "variety": ["variety", "var", "cultivar", "cv", "genotype", "accession",
                "line", "hybrid"],
    "soil_type": ["soil type", "soiltype", "soil", "texture", "soil texture",
                  "soil class"],
    "organic_inorganic": ["organic inorganic", "organic or inorganic", "input type",
                          "inputtype", "farming system", "farming type"],
    "date": ["date", "sampling date", "date of sampling", "collection date",
             "month", "year", "season", "time", "timepoint", "dap", "das"],
    "length": ["length", "body length", "l", "total length"],
    "diameter": ["diameter", "body diameter", "width", "body width", "greatest width"],
}

MAIN_ROLES = ("sample_id", "taxon", "count")
METADATA_ROLES = ("group", "treatment", "replicate", "location", "crop", "variety",
                  "soil_type", "organic_inorganic", "date")
TAXON_ATTR_ROLES = ("trophic_group", "cp", "source")
MEASUREMENT_ROLES = ("length", "diameter")

# Taxonomic morphology of a name. Weak evidence on its own — deliberately so.
_TAXON_SUFFIXES = ("tylenchus", "laimus", "oides", "dorus", "nema", "rhabditis",
                   "onchus", "chus", "ellus", "ella", "phelenchus", "cephalobus",
                   "idae", "inae", "olaimus", "orhynchus")
_TAXON_TOKENS = ("sp.", "spp.", "sp", "spp", "juveniles", "cyst")


def normalise_header(name: Any) -> str:
    """Lower-case, strip punctuation, collapse whitespace."""
    s = str(name) if name is not None else ""
    s = s.replace("_", " ").replace("-", " ").replace(".", " ")
    s = re.sub(r"[^\w\s]", " ", s, flags=re.UNICODE)
    s = re.sub(r"\s+", " ", s).strip().lower()
    return s


# --------------------------------------------------------------------------
# 2. Result containers
# --------------------------------------------------------------------------

@dataclass
class SheetProfile:
    name: str
    n_rows: int
    n_cols: int
    columns: List[str]
    dtypes: Dict[str, str]
    missing_pct: Dict[str, float]
    unique_counts: Dict[str, int]
    orientation: str            # long | wide_samples_rows | wide_taxa_rows | unknown
    orientation_score: float    # 0-1 heuristic match score, NOT a probability
    content_role: str           # data | taxa_reference | sample_metadata | unknown
    content_score: float
    notes: List[str] = field(default_factory=list)


@dataclass
class DatasetStructure:
    filename: str
    filetype: str
    sheets: List[SheetProfile]
    primary_sheet: Optional[str]
    ambiguous: bool
    notes: List[str] = field(default_factory=list)

    def profile(self, sheet: str) -> Optional[SheetProfile]:
        for s in self.sheets:
            if s.name == sheet:
                return s
        return None


@dataclass
class ColumnMapping:
    column: str
    role: str            # one of the roles above, or "metadata" / "ignore"
    score: float         # 0-1 heuristic
    evidence: str


@dataclass
class MappingResult:
    orientation: str
    mappings: List[ColumnMapping]
    taxon_columns: List[str]      # wide layouts only
    warnings: List[str] = field(default_factory=list)
    needs_confirmation: bool = False

    def by_role(self, role: str) -> List[str]:
        return [m.column for m in self.mappings if m.role == role]

    def first(self, role: str) -> Optional[str]:
        cols = self.by_role(role)
        return cols[0] if cols else None

    def as_frame(self) -> pd.DataFrame:
        return pd.DataFrame([
            {"Original column": m.column, "Detected role": m.role,
             "Match score": round(m.score, 2), "Evidence": m.evidence}
            for m in self.mappings
        ])


@dataclass
class Standardized:
    counts_long: pd.DataFrame     # sample_id, taxon, count
    counts_wide: pd.DataFrame     # NEMO-NEMA native: taxa rows x sample columns
    samples: pd.DataFrame         # sample_id + metadata
    taxa: pd.DataFrame            # taxon, trophic_group, cp, source, provenance
    log: List[Dict[str, Any]] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)
    blocking: List[str] = field(default_factory=list)
    original: Optional[pd.DataFrame] = None
    source_info: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationCheck:
    name: str
    passed: bool
    original: Any
    converted: Any
    detail: str = ""


@dataclass
class ValidationReport:
    checks: List[ValidationCheck]

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def as_frame(self) -> pd.DataFrame:
        return pd.DataFrame([
            {"Check": c.name, "Result": "PASS" if c.passed else "FAIL",
             "Original": c.original, "Converted": c.converted, "Detail": c.detail}
            for c in self.checks
        ])


# --------------------------------------------------------------------------
# 3. Cell-level helpers — numeric handling never guesses
# --------------------------------------------------------------------------

_NUMERIC_RE = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")


def parse_count_cell(v: Any) -> Tuple[Optional[float], Optional[str]]:
    """Return (value, problem). Whitespace is stripped; nothing else is altered.

    Anything that is not unambiguously numeric is reported as a problem rather
    than coerced, so that '<1', 'ND' or '1,200' can never silently become a
    number the user did not write.
    """
    if v is None:
        return None, None
    if isinstance(v, (int, np.integer, float, np.floating)):
        if isinstance(v, float) and np.isnan(v):
            return None, None
        return float(v), None
    s = str(v).strip()
    if s == "":
        return None, None
    if _NUMERIC_RE.match(s):
        return float(s), None
    return None, s


def profile_series(s: pd.Series) -> Dict[str, Any]:
    non_null = s.dropna()
    numeric_like = 0
    for v in non_null:
        val, prob = parse_count_cell(v)
        if val is not None:
            numeric_like += 1
    n = len(non_null)
    return {
        "n": n,
        "n_missing": int(s.isna().sum()),
        "missing_pct": float(s.isna().mean() * 100) if len(s) else 0.0,
        "n_unique": int(non_null.nunique()),
        "numeric_frac": (numeric_like / n) if n else 0.0,
        "dtype": str(s.dtype),
    }


def looks_like_taxon_string(value: Any) -> float:
    """Score a single string on how much it looks like a nematode taxon name.

    Reference-database membership dominates. Morphology of the name contributes
    only a little, and place/crop/treatment vocabulary subtracts, because the
    original failure mode was text being treated as taxonomy.
    """
    if not isinstance(value, str):
        return 0.0
    raw = value.strip()
    if not raw:
        return 0.0
    low = raw.lower()
    genus = low.split()[0] if low.split() else low
    genus = genus.rstrip(",;:")

    if low in _GAZETTEER or genus in _GAZETTEER:
        return -1.0
    if low in _REFERENCE_TAXA or genus in _REFERENCE_TAXA:
        return 1.0

    score = 0.0
    tokens = low.replace(",", " ").split()
    if any(t in _TAXON_TOKENS for t in tokens[1:]):
        score += 0.35
    if genus.endswith(_TAXON_SUFFIXES):
        score += 0.30
    if len(tokens) == 2 and raw[:1].isupper() and tokens[1].islower() and len(tokens[1]) > 3:
        score += 0.20          # binomial shape
    if re.search(r"\d", raw):
        score -= 0.30
    if len(raw) < 5:
        score -= 0.20
    return max(min(score, 0.85), 0.0)


def taxon_evidence_for_values(values: Sequence[Any]) -> Tuple[float, int, int]:
    """(mean score, n matched against reference, n scored)."""
    vals = [v for v in values if isinstance(v, str) and v.strip()]
    if not vals:
        return 0.0, 0, 0
    sample = vals[:400]
    scores = [looks_like_taxon_string(v) for v in sample]
    matched = sum(1 for s in scores if s >= 1.0)
    return float(np.mean(scores)), matched, len(sample)


def taxon_evidence_for_header(name: str) -> Tuple[float, str]:
    """Is this *column header* itself a taxon name (wide layout)?"""
    low = normalise_header(name)
    genus = low.split()[0] if low.split() else low
    if low in _GAZETTEER or genus in _GAZETTEER:
        return -1.0, "matches place/crop/treatment vocabulary"
    if low in _REFERENCE_TAXA or genus in _REFERENCE_TAXA:
        return 1.0, "in taxon reference database"
    for role, syns in ROLE_SYNONYMS.items():
        if low in syns:
            return -1.0, f"header is a known {role} label"
    s = looks_like_taxon_string(str(name))
    if s > 0:
        return s, "taxonomic name pattern only (not in reference)"
    return 0.0, "no taxonomic evidence"


def role_from_header(name: str) -> Tuple[Optional[str], float, str]:
    low = normalise_header(name)
    if not low:
        return None, 0.0, "blank header"
    for role, syns in ROLE_SYNONYMS.items():
        if low in syns:
            return role, 0.95, "exact header match"
    for role, syns in ROLE_SYNONYMS.items():
        for syn in syns:
            if len(syn) < 3:
                continue
            if low.startswith(syn + " ") or low.endswith(" " + syn) or f" {syn} " in f" {low} ":
                return role, 0.75, f"header contains '{syn}'"
    return None, 0.0, "no header match"


# --------------------------------------------------------------------------
# 4. detect_dataset_structure
# --------------------------------------------------------------------------

def _read_all_sheets(src: Any, filename: str) -> Tuple[Dict[str, pd.DataFrame], str]:
    ext = os.path.splitext(filename)[1].lower()
    if ext in (".csv", ".txt", ".tsv"):
        sep = "\t" if ext == ".tsv" else None
        if hasattr(src, "seek"):
            src.seek(0)
        df = pd.read_csv(src, sep=sep, engine="python")
        return {"csv": df}, "csv"
    if hasattr(src, "seek"):
        src.seek(0)
    book = pd.read_excel(src, sheet_name=None)
    return book, "excel"


def _classify_sheet(name: str, df: pd.DataFrame) -> SheetProfile:
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    cols = list(df.columns)
    prof = {c: profile_series(df[c]) for c in cols}

    header_roles = {c: role_from_header(c) for c in cols}
    header_taxon = {c: taxon_evidence_for_header(c) for c in cols}

    numeric_cols = [c for c in cols if prof[c]["numeric_frac"] >= 0.9 and prof[c]["n"] > 0]
    taxon_headers = [c for c in numeric_cols if header_taxon[c][0] >= 0.6]

    has_taxon_col = any(header_roles[c][0] == "taxon" for c in cols)
    has_count_col = any(header_roles[c][0] == "count" for c in cols)
    has_sample_col = any(header_roles[c][0] == "sample_id" for c in cols)

    # A column of taxon *values* (long layout, or native taxa-in-first-column)
    value_taxon_cols = []
    for c in cols:
        if prof[c]["numeric_frac"] > 0.5:
            continue
        mean_s, matched, n = taxon_evidence_for_values(df[c].tolist())
        if n >= 3 and (matched / n >= 0.5 or mean_s >= 0.55):
            value_taxon_cols.append(c)

    orientation, score, notes = "unknown", 0.0, []

    first_col = cols[0] if cols else None
    rest_numeric = numeric_cols and all(c in numeric_cols for c in cols[1:]) if len(cols) > 1 else False

    if has_count_col and (has_taxon_col or value_taxon_cols) and has_sample_col:
        orientation, score = "long", 0.95
        notes.append("sample, taxon and count columns all present")
    elif has_count_col and (has_taxon_col or value_taxon_cols):
        orientation, score = "long", 0.75
        notes.append("taxon and count columns present; sample column uncertain")
    elif first_col in value_taxon_cols and rest_numeric:
        orientation, score = "wide_taxa_rows", 0.90
        notes.append("taxa in first column, samples across the header row (NEMO-NEMA native)")
    elif taxon_headers and len(taxon_headers) >= 2:
        orientation, score = "wide_samples_rows", 0.90
        notes.append(f"{len(taxon_headers)} column headers matched taxon names")
    elif len(numeric_cols) >= 3 and value_taxon_cols:
        orientation, score = "wide_samples_rows", 0.45
        notes.append("numeric block present but headers not confirmed as taxa")

    # What is this sheet for?
    content_role, content_score = "unknown", 0.3
    n_taxon_attr = sum(1 for c in cols if header_roles[c][0] in TAXON_ATTR_ROLES)
    n_meta = sum(1 for c in cols if header_roles[c][0] in METADATA_ROLES)
    if orientation in ("long", "wide_taxa_rows", "wide_samples_rows") and score >= 0.7:
        content_role, content_score = "data", score
    elif (has_taxon_col or value_taxon_cols) and n_taxon_attr >= 1 and not has_count_col:
        content_role, content_score = "taxa_reference", 0.85
        notes.append("taxon attributes without counts")
    elif has_sample_col and n_meta >= 1 and not has_count_col and not taxon_headers:
        content_role, content_score = "sample_metadata", 0.85
        notes.append("sample identifiers with descriptive columns only")
    elif orientation != "unknown":
        content_role, content_score = "data", max(score, 0.4)

    return SheetProfile(
        name=name,
        n_rows=int(len(df)),
        n_cols=int(len(cols)),
        columns=cols,
        dtypes={c: prof[c]["dtype"] for c in cols},
        missing_pct={c: round(prof[c]["missing_pct"], 1) for c in cols},
        unique_counts={c: prof[c]["n_unique"] for c in cols},
        orientation=orientation,
        orientation_score=score,
        content_role=content_role,
        content_score=content_score,
        notes=notes,
    )


def detect_dataset_structure(src: Any, filename: str = "upload.xlsx") -> Tuple[DatasetStructure, Dict[str, pd.DataFrame]]:
    """Inspect an upload without altering it. Returns (structure, {sheet: df})."""
    book, ftype = _read_all_sheets(src, filename)
    profiles = [_classify_sheet(name, df) for name, df in book.items()]

    data_sheets = [p for p in profiles if p.content_role == "data"]
    data_sheets.sort(key=lambda p: (p.content_score, p.n_rows), reverse=True)

    primary = data_sheets[0].name if data_sheets else (profiles[0].name if profiles else None)
    ambiguous = False
    notes = []
    if len(data_sheets) > 1:
        top, second = data_sheets[0], data_sheets[1]
        if top.content_score - second.content_score < 0.15:
            ambiguous = True
            notes.append("More than one sheet could hold the principal observations — "
                         "please choose rather than let the app assume.")
    if data_sheets and data_sheets[0].content_score < 0.7:
        ambiguous = True
        notes.append("Structure of the best candidate sheet is not clear-cut; confirm the mapping.")

    return DatasetStructure(filename=filename, filetype=ftype, sheets=profiles,
                            primary_sheet=primary, ambiguous=ambiguous, notes=notes), book


# --------------------------------------------------------------------------
# 5. map_columns
# --------------------------------------------------------------------------

def map_columns(df: pd.DataFrame, orientation: str) -> MappingResult:
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    cols = list(df.columns)
    prof = {c: profile_series(df[c]) for c in cols}
    mappings: List[ColumnMapping] = []
    taxon_cols: List[str] = []
    warnings: List[str] = []

    if orientation == "wide_taxa_rows":
        # First column holds taxa; every other column is a sample.
        first = cols[0]
        mappings.append(ColumnMapping(first, "taxon", 0.92,
                                      "first column of a taxa-in-rows matrix"))
        for c in cols[1:]:
            if prof[c]["numeric_frac"] >= 0.9 or prof[c]["n"] == 0:
                mappings.append(ColumnMapping(c, "sample_column", 0.90,
                                              "numeric column in a taxa-in-rows matrix"))
            else:
                role, sc, ev = role_from_header(c)
                mappings.append(ColumnMapping(c, role or "metadata", sc or 0.4,
                                              ev if role else "non-numeric column in a count matrix"))
                warnings.append(f"'{c}' is not numeric and was not treated as a sample column.")
        return MappingResult(orientation, mappings, taxon_cols, warnings,
                             needs_confirmation=bool(warnings))

    for c in cols:
        role, score, evidence = role_from_header(c)
        p = prof[c]

        if orientation == "wide_samples_rows":
            t_score, t_ev = taxon_evidence_for_header(c)
            numeric = p["numeric_frac"] >= 0.9 and p["n"] > 0
            if numeric and t_score >= 0.6 and role not in ("count", "cp", "replicate", "length", "diameter"):
                taxon_cols.append(c)
                mappings.append(ColumnMapping(c, "taxon_column", t_score, t_ev))
                continue
            # A reference-matched header with some unparseable cells is still a
            # taxon column. Keeping it means the offending cells get reported;
            # dropping it would lose real observations without saying so.
            if t_score >= 1.0 and 0.5 <= p["numeric_frac"] < 0.9:
                taxon_cols.append(c)
                mappings.append(ColumnMapping(c, "taxon_column", t_score,
                                              t_ev + "; some cells are not numeric"))
                warnings.append(f"'{c}' is a known taxon but "
                                f"{round((1 - p['numeric_frac']) * 100)}% of its cells are not "
                                f"numeric. Those cells will be reported, not converted.")
                continue
            if numeric and t_score > 0 and role is None:
                taxon_cols.append(c)
                mappings.append(ColumnMapping(c, "taxon_column", t_score, t_ev + " — confirm"))
                warnings.append(f"'{c}' looks taxonomic only by name pattern; it is not in the "
                                f"reference database. Confirm before analysis.")
                continue

        if role is None:
            # Value-level evidence for a taxon column in a long layout
            if orientation == "long" and p["numeric_frac"] < 0.5:
                mean_s, matched, n = taxon_evidence_for_values(df[c].tolist())
                if n >= 3 and (matched / n >= 0.5 or mean_s >= 0.6):
                    frac = matched / n if n else 0
                    mappings.append(ColumnMapping(
                        c, "taxon", min(0.55 + frac * 0.4, 0.95),
                        f"{matched}/{n} values matched the taxon reference"))
                    continue
            role, score, evidence = "metadata", 0.35, "unrecognised column, preserved as metadata"

        # Guard: a role claim that contradicts the data type is downgraded.
        if role == "count" and p["numeric_frac"] < 0.9 and p["n"] > 0:
            warnings.append(f"'{c}' is named like a count column but {round((1-p['numeric_frac'])*100)}% "
                            f"of its values are not numeric.")
            score = min(score, 0.45)
        if role == "taxon" and p["numeric_frac"] > 0.5:
            warnings.append(f"'{c}' is named like a taxon column but holds numbers.")
            score = min(score, 0.35)
        if role == "sample_id" and p["n"] and p["n_unique"] == 1:
            warnings.append(f"'{c}' is named like a sample column but has a single repeated value.")
            score = min(score, 0.4)

        mappings.append(ColumnMapping(c, role, score, evidence))

    if orientation == "long":
        if not any(m.role == "taxon" for m in mappings):
            warnings.append("No taxon column identified.")
        if not any(m.role == "count" for m in mappings):
            warnings.append("No count column identified.")
        if not any(m.role == "sample_id" for m in mappings):
            warnings.append("No sample column identified.")
    if orientation == "wide_samples_rows" and len(taxon_cols) < 2:
        warnings.append("Fewer than two taxon columns detected in a wide layout.")

    low_conf = [m.column for m in mappings if m.role not in ("metadata", "ignore") and m.score < 0.6]
    needs = bool(warnings) or bool(low_conf)
    return MappingResult(orientation, mappings, taxon_cols, warnings, needs_confirmation=needs)


# --------------------------------------------------------------------------
# 6. convert_to_nemonema_schema
# --------------------------------------------------------------------------

def _clean_label(v: Any) -> Any:
    """Whitespace normalisation only. Case and spelling are left alone."""
    if isinstance(v, str):
        return re.sub(r"\s+", " ", v).strip()
    return v


def convert_to_nemonema_schema(df: pd.DataFrame,
                               mapping: MappingResult,
                               blank_policy: str = "absent",
                               duplicate_policy: str = "flag",
                               source_info: Optional[Dict[str, Any]] = None) -> Standardized:
    """Restructure to the internal schema. No observation is altered.

    blank_policy    'absent' (default) omits empty cells from the long table —
                    arithmetically identical to a zero for every index NEMO-NEMA
                    computes, and it invents nothing. 'zero' writes explicit
                    zeros and records that it did.
    duplicate_policy 'flag' refuses to proceed when the same sample x taxon
                    appears twice. Summing duplicates silently would change the
                    data; the user must decide.
    """
    src = dict(source_info or {})
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    log: List[Dict[str, Any]] = []
    issues: List[str] = []
    blocking: List[str] = []
    orientation = mapping.orientation

    def record(step: str, detail: str, before: Any = "", after: Any = ""):
        log.append({"Step": step, "Detail": detail, "Rows before": before, "Rows after": after})

    record("Input", f"{orientation} layout, {len(df)} rows x {len(df.columns)} columns", len(df), len(df))

    # ---------------- wide: taxa in rows, samples in columns ----------------
    if orientation == "wide_taxa_rows":
        taxon_col = mapping.first("taxon") or df.columns[0]
        sample_cols = [m.column for m in mapping.mappings if m.role == "sample_column"]
        meta_cols = [m.column for m in mapping.mappings
                     if m.role not in ("taxon", "sample_column", "ignore")]
        records, bad_cells = [], []
        for _, row in df.iterrows():
            taxon = _clean_label(row[taxon_col])
            if taxon is None or (isinstance(taxon, float) and np.isnan(taxon)) or taxon == "":
                continue
            for sc in sample_cols:
                val, prob = parse_count_cell(row[sc])
                if prob is not None:
                    bad_cells.append((str(taxon), sc, prob))
                    continue
                if val is None:
                    if blank_policy == "zero":
                        records.append({"sample_id": _clean_label(sc), "taxon": taxon, "count": 0.0})
                    continue
                records.append({"sample_id": _clean_label(sc), "taxon": taxon, "count": val})
        counts_long = pd.DataFrame(records, columns=["sample_id", "taxon", "count"])
        taxa_meta = df[[taxon_col] + meta_cols].copy() if meta_cols else df[[taxon_col]].copy()
        taxa_meta.columns = ["taxon"] + meta_cols
        taxa_meta["taxon"] = taxa_meta["taxon"].map(_clean_label)
        sample_meta = pd.DataFrame({"sample_id": [_clean_label(c) for c in sample_cols]})
        record("Transpose", f"matrix of {len(df)} taxa x {len(sample_cols)} samples unstacked to long form",
               len(df), len(counts_long))

    # ---------------- wide: samples in rows, taxa in columns ----------------
    elif orientation == "wide_samples_rows":
        sample_col = mapping.first("sample_id")
        taxon_cols = mapping.taxon_columns
        if sample_col is None:
            df = df.reset_index(drop=True)
            df["__row_id__"] = [f"row{i+1}" for i in range(len(df))]
            sample_col = "__row_id__"
            issues.append("No sample column was found; row numbers were used as sample "
                          "identifiers. Check this before relying on the results.")
        meta_cols = [m.column for m in mapping.mappings
                     if m.role not in ("taxon_column", "ignore") and m.column != sample_col]
        records, bad_cells = [], []
        for _, row in df.iterrows():
            sid = _clean_label(row[sample_col])
            for tc in taxon_cols:
                val, prob = parse_count_cell(row[tc])
                if prob is not None:
                    bad_cells.append((str(sid), tc, prob))
                    continue
                if val is None:
                    if blank_policy == "zero":
                        records.append({"sample_id": sid, "taxon": _clean_label(tc), "count": 0.0})
                    continue
                records.append({"sample_id": sid, "taxon": _clean_label(tc), "count": val})
        counts_long = pd.DataFrame(records, columns=["sample_id", "taxon", "count"])
        sample_meta = df[[sample_col] + meta_cols].copy()
        sample_meta.columns = ["sample_id"] + meta_cols
        sample_meta["sample_id"] = sample_meta["sample_id"].map(_clean_label)
        taxa_meta = pd.DataFrame({"taxon": [_clean_label(c) for c in taxon_cols]})
        record("Unpivot", f"{len(taxon_cols)} taxon columns melted against '{sample_col}'",
               len(df), len(counts_long))

    # ---------------- long ----------------
    elif orientation == "long":
        sample_col = mapping.first("sample_id")
        taxon_col = mapping.first("taxon")
        count_col = mapping.first("count")
        for label, col in (("sample", sample_col), ("taxon", taxon_col), ("count", count_col)):
            if col is None:
                blocking.append(f"Long-format data needs a {label} column; none was identified.")
        if blocking:
            return Standardized(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(),
                                log, issues, blocking, original=df, source_info=src)
        bad_cells = []
        rows = []
        for _, row in df.iterrows():
            val, prob = parse_count_cell(row[count_col])
            sid, txn = _clean_label(row[sample_col]), _clean_label(row[taxon_col])
            if prob is not None:
                bad_cells.append((str(sid), str(txn), prob))
                continue
            if val is None:
                continue
            rows.append({"sample_id": sid, "taxon": txn, "count": val})
        counts_long = pd.DataFrame(rows, columns=["sample_id", "taxon", "count"])
        meta_cols = [m.column for m in mapping.mappings
                     if m.column not in (sample_col, taxon_col, count_col)
                     and m.role in METADATA_ROLES + ("metadata",) + MEASUREMENT_ROLES]
        if meta_cols:
            sample_meta = df[[sample_col] + meta_cols].copy()
            sample_meta.columns = ["sample_id"] + meta_cols
            sample_meta["sample_id"] = sample_meta["sample_id"].map(_clean_label)
            sample_meta = sample_meta.drop_duplicates(subset=["sample_id"], keep="first")
        else:
            sample_meta = pd.DataFrame({"sample_id": counts_long["sample_id"].unique()})
        attr_cols = [m.column for m in mapping.mappings if m.role in TAXON_ATTR_ROLES]
        if attr_cols:
            taxa_meta = df[[taxon_col] + attr_cols].copy()
            taxa_meta.columns = ["taxon"] + [m.role for m in mapping.mappings if m.role in TAXON_ATTR_ROLES]
            taxa_meta["taxon"] = taxa_meta["taxon"].map(_clean_label)
            taxa_meta = taxa_meta.drop_duplicates(subset=["taxon"], keep="first")
        else:
            taxa_meta = pd.DataFrame({"taxon": counts_long["taxon"].unique()})
        record("Pass-through", "already in long form; column names mapped, values untouched",
               len(df), len(counts_long))
    else:
        blocking.append("Dataset layout could not be determined. Set the mapping manually.")
        return Standardized(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(),
                            log, issues, blocking, original=df, source_info=src)

    # ---------------- shared post-processing ----------------
    if bad_cells:
        shown = "; ".join(f"{a}/{b} = '{c}'" for a, b, c in bad_cells[:8])
        blocking.append(
            f"{len(bad_cells)} cell(s) in the count field are not numeric and were NOT converted: "
            f"{shown}{' …' if len(bad_cells) > 8 else ''}. Correct them in the source file — the "
            f"importer will not guess what they mean.")

    if not counts_long.empty:
        dup = counts_long.duplicated(subset=["sample_id", "taxon"], keep=False)
        if dup.any():
            n_pairs = counts_long.loc[dup, ["sample_id", "taxon"]].drop_duplicates().shape[0]
            if duplicate_policy == "flag":
                blocking.append(
                    f"{n_pairs} sample x taxon combination(s) appear more than once. These may be "
                    f"replicates that need their own sample IDs, or genuine duplicates. Nothing was "
                    f"summed or dropped — decide in the source file, or choose 'sum duplicates'.")
            elif duplicate_policy == "sum":
                before = len(counts_long)
                counts_long = (counts_long.groupby(["sample_id", "taxon"], as_index=False)["count"]
                               .sum())
                record("Duplicates", f"{n_pairs} repeated sample x taxon pairs summed on user instruction",
                       before, len(counts_long))
                issues.append(f"{n_pairs} duplicate sample x taxon pairs were summed at your request.")

    if counts_long.empty:
        blocking.append(
            "The conversion produced no observations at all. Either no count field was "
            "identified, or every value in it was unreadable. Nothing was analysed.")

    if blank_policy == "absent":
        record("Blanks", "empty cells omitted from the long table (absence), not written as zeros")
    else:
        record("Blanks", "empty cells written as explicit zeros on user instruction")

    # Native NEMO-NEMA matrix: taxa down rows, samples across columns.
    if counts_long.empty:
        counts_wide = pd.DataFrame()
    else:
        # Order of first appearance, not alphabetical: the user's row and column
        # order carries meaning (survey sequence, treatment layout) and is part
        # of what "structural only" has to preserve.
        taxon_order = list(dict.fromkeys(counts_long["taxon"]))
        sample_order = list(dict.fromkeys(counts_long["sample_id"]))
        counts_wide = (counts_long.pivot_table(index="taxon", columns="sample_id",
                                               values="count", aggfunc="sum")
                       .reindex(index=taxon_order, columns=sample_order)
                       .fillna(0))
        counts_wide.index.name = "Taxon"
        counts_wide = counts_wide.reset_index()

    taxa_meta = taxa_meta.drop_duplicates(subset=["taxon"]).reset_index(drop=True)
    for c in ("trophic_group", "cp", "source"):
        if c not in taxa_meta.columns:
            taxa_meta[c] = np.nan
    taxa_meta["provenance"] = np.where(taxa_meta["trophic_group"].notna(),
                                       "User-provided", "Not supplied")
    sample_meta = sample_meta.drop_duplicates(subset=["sample_id"]).reset_index(drop=True)

    record("Output", f"{counts_long['sample_id'].nunique() if not counts_long.empty else 0} samples, "
                     f"{counts_long['taxon'].nunique() if not counts_long.empty else 0} taxa",
           len(df), len(counts_long))

    return Standardized(counts_long=counts_long, counts_wide=counts_wide,
                        samples=sample_meta, taxa=taxa_meta, log=log,
                        issues=issues, blocking=blocking, original=df, source_info=src)


def mark_inferred_taxa(std: Standardized, assign_fn) -> Standardized:
    """Fill trophic_group / cp from the host app's auto-assign function.

    `assign_fn(taxon) -> dict(trophic_group=..., cp=..., source=...)` or None.
    Anything filled this way is labelled Inferred/Reference-derived and is never
    presented as a user observation.
    """
    if std.taxa.empty:
        return std
    t = std.taxa.copy()
    for i, row in t.iterrows():
        if pd.notna(row.get("trophic_group")):
            continue
        got = assign_fn(row["taxon"])
        if not got:
            continue
        t.at[i, "trophic_group"] = got.get("trophic_group")
        t.at[i, "cp"] = got.get("cp")
        t.at[i, "source"] = got.get("source")
        t.at[i, "provenance"] = "Inferred/Reference-derived"
    std.taxa = t
    return std


# --------------------------------------------------------------------------
# 7. validate_conversion
# --------------------------------------------------------------------------

def _approx(a: float, b: float, tol: float = 1e-6) -> bool:
    if a is None or b is None:
        return False
    if np.isnan(a) and np.isnan(b):
        return True
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def validate_conversion(original: pd.DataFrame,
                        mapping: MappingResult,
                        std: Standardized) -> ValidationReport:
    """Recompute every summary statistic from the untouched original and compare."""
    checks: List[ValidationCheck] = []
    orientation = mapping.orientation
    orig = original.copy()
    orig.columns = [str(c) for c in orig.columns]

    # Rebuild the original observation set independently of the conversion code.
    vals: List[float] = []
    o_samples, o_taxa = set(), set()
    if orientation == "wide_taxa_rows":
        taxon_col = mapping.first("taxon") or orig.columns[0]
        sample_cols = [m.column for m in mapping.mappings if m.role == "sample_column"]
        o_samples = {_clean_label(c) for c in sample_cols}
        for _, row in orig.iterrows():
            t = _clean_label(row[taxon_col])
            if t in (None, "") or (isinstance(t, float) and np.isnan(t)):
                continue
            o_taxa.add(t)
            for sc in sample_cols:
                v, p = parse_count_cell(row[sc])
                if v is not None:
                    vals.append(v)
    elif orientation == "wide_samples_rows":
        sample_col = mapping.first("sample_id")
        tcols = mapping.taxon_columns
        o_taxa = {_clean_label(c) for c in tcols}
        for idx, row in orig.reset_index(drop=True).iterrows():
            sid = _clean_label(row[sample_col]) if sample_col else f"row{idx+1}"
            o_samples.add(sid)
            for tc in tcols:
                v, p = parse_count_cell(row[tc])
                if v is not None:
                    vals.append(v)
    elif orientation == "long":
        sc, tc, cc = mapping.first("sample_id"), mapping.first("taxon"), mapping.first("count")
        for _, row in orig.iterrows():
            v, p = parse_count_cell(row[cc])
            if v is None:
                continue
            vals.append(v)
            o_samples.add(_clean_label(row[sc]))
            o_taxa.add(_clean_label(row[tc]))

    c_vals = std.counts_long["count"].tolist() if not std.counts_long.empty else []
    c_samples = set(std.counts_long["sample_id"].unique()) if not std.counts_long.empty else set()
    c_taxa = set(std.counts_long["taxon"].unique()) if not std.counts_long.empty else set()

    # If the user asked for explicit zeros, the observation count legitimately grows.
    zero_padded = any(l["Detail"].startswith("empty cells written") for l in std.log)
    summed = any(l["Step"] == "Duplicates" for l in std.log)

    n_obs_ok = (len(c_vals) >= len(vals)) if zero_padded else (
        len(c_vals) == len(vals) if not summed else len(c_vals) <= len(vals))
    checks.append(ValidationCheck("Observations retained", n_obs_ok, len(vals), len(c_vals),
                                  "explicit zeros added" if zero_padded else
                                  ("duplicates summed" if summed else "one-to-one")))

    checks.append(ValidationCheck("Total abundance", _approx(float(np.sum(vals)) if vals else 0.0,
                                                             float(np.sum(c_vals)) if c_vals else 0.0),
                                  round(float(np.sum(vals)), 6) if vals else 0.0,
                                  round(float(np.sum(c_vals)), 6) if c_vals else 0.0,
                                  "sum(original) must equal sum(converted)"))

    nz_orig = [v for v in vals if v != 0]
    nz_conv = [v for v in c_vals if v != 0]
    checks.append(ValidationCheck("Minimum count",
                                  _approx(min(nz_orig) if nz_orig else 0.0,
                                          min(nz_conv) if nz_conv else 0.0),
                                  min(nz_orig) if nz_orig else "-",
                                  min(nz_conv) if nz_conv else "-", "non-zero values"))
    checks.append(ValidationCheck("Maximum count",
                                  _approx(max(vals) if vals else 0.0, max(c_vals) if c_vals else 0.0),
                                  max(vals) if vals else "-", max(c_vals) if c_vals else "-", ""))

    checks.append(ValidationCheck("Sample identifiers", o_samples == c_samples,
                                  len(o_samples), len(c_samples),
                                  "" if o_samples == c_samples
                                  else f"differs: {sorted(o_samples ^ c_samples)[:5]}"))
    checks.append(ValidationCheck("Taxon names", o_taxa == c_taxa, len(o_taxa), len(c_taxa),
                                  "" if o_taxa == c_taxa
                                  else f"differs: {sorted(o_taxa ^ c_taxa)[:5]}"))

    # Matrix round-trip: the native wide table must reproduce the long table.
    if not std.counts_wide.empty:
        wsum = float(std.counts_wide.drop(columns=["Taxon"]).to_numpy(dtype=float).sum())
        checks.append(ValidationCheck("Matrix round-trip",
                                      _approx(wsum, float(np.sum(c_vals)) if c_vals else 0.0),
                                      round(float(np.sum(c_vals)), 6) if c_vals else 0.0,
                                      round(wsum, 6),
                                      "long table and native matrix must agree"))

    # Metadata must still point at the right samples.
    if not std.samples.empty:
        orphan = c_samples - set(std.samples["sample_id"])
        checks.append(ValidationCheck("Metadata alignment", not orphan,
                                      len(c_samples), len(set(std.samples["sample_id"])),
                                      "" if not orphan else f"{len(orphan)} samples without a metadata row"))

    for b in std.blocking:
        checks.append(ValidationCheck("Blocking issue", False, "-", "-", b))

    return ValidationReport(checks)


# --------------------------------------------------------------------------
# 8. Reporting and export
# --------------------------------------------------------------------------

def generate_conversion_report(structure: DatasetStructure,
                               sheet: str,
                               mapping: MappingResult,
                               std: Standardized,
                               report: ValidationReport) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    def add(section, item, value):
        rows.append({"Section": section, "Item": item, "Value": value})

    add("Provenance", "Importer version", f"nemonema_ingest {INGEST_VERSION}")
    add("Provenance", "Converted at", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    add("Provenance", "Original filename", structure.filename)
    add("Provenance", "Original sheet", sheet)
    add("Provenance", "Taxon reference", reference_status())
    prof = structure.profile(sheet)
    if prof:
        add("Original", "Rows", prof.n_rows)
        add("Original", "Columns", prof.n_cols)
        add("Original", "Column names", " | ".join(prof.columns))
        add("Original", "Detected layout", f"{prof.orientation} (score {prof.orientation_score:.2f})")
    for m in mapping.mappings:
        add("Column roles", m.column, f"{m.role} (score {m.score:.2f}; {m.evidence})")
    for w in mapping.warnings:
        add("Mapping warnings", "warning", w)
    for entry in std.log:
        add("Transformations", entry["Step"],
            f"{entry['Detail']}" + (f" [{entry['Rows before']} -> {entry['Rows after']}]"
                                    if entry["Rows before"] != "" else ""))
    for i in std.issues:
        add("Notes", "note", i)
    for b in std.blocking:
        add("Blocking", "blocked", b)
    for c in report.checks:
        add("Validation", c.name,
            f"{'PASS' if c.passed else 'FAIL'} — original {c.original}, converted {c.converted}"
            + (f" ({c.detail})" if c.detail else ""))
    return pd.DataFrame(rows)


# nemonema_core.load() reads every sheet with index_col=0 and lower-cases the
# column names of `taxa` and `samples`. The rest of the app expects the trophic
# column to be called `trophic` and to hold PP/BF/FF/OM/PR, and the grouping
# factor to be called `group`. The export below is written to that contract, not
# to a schema of my own choosing.

_TROPHIC_VOCAB = {
    "pp": "PP", "plant parasite": "PP", "plant-parasite": "PP",
    "plant parasitic": "PP", "plant-parasitic": "PP", "herbivore": "PP",
    "plant feeder": "PP",
    "bf": "BF", "bacterivore": "BF", "bacteriovore": "BF",
    "bacterial feeder": "BF", "bacteria feeder": "BF", "bacterivorous": "BF",
    "ff": "FF", "fungivore": "FF", "fungal feeder": "FF", "fungivorous": "FF",
    "mycophagous": "FF",
    "om": "OM", "omnivore": "OM", "omnivorous": "OM",
    "pr": "PR", "predator": "PR", "predatory": "PR", "carnivore": "PR",
}


def normalise_trophic(value: Any) -> Tuple[Any, Optional[str]]:
    """Map a spelt-out feeding group onto the app's five codes.

    Returns (value, note). Anything not in the table is passed through
    untouched so the Data check tab can flag it, rather than being guessed at.
    """
    if not isinstance(value, str):
        return value, None
    key = re.sub(r"\s+", " ", value).strip().lower().rstrip("s")
    if key in _TROPHIC_VOCAB:
        code = _TROPHIC_VOCAB[key]
        return code, (None if value.strip() == code else f"'{value.strip()}' -> {code}")
    if value.strip().upper() in ("PP", "BF", "FF", "OM", "PR"):
        return value.strip().upper(), None
    return value, f"'{value}' is not one of PP/BF/FF/OM/PR — left as written"


def build_standardized_workbook(std: Standardized,
                                conversion_log: pd.DataFrame,
                                group_col: Optional[str] = None) -> Tuple[bytes, List[str]]:
    """Write a workbook that nemonema_core.load() can read directly.

    `taxa` and `samples` are omitted entirely when there is nothing to put in
    them — that is the counts-only path the app already supports, and it is
    safer than writing a sheet full of blanks. Returns (bytes, notes).
    """
    notes: List[str] = []
    buf = io.BytesIO()

    counts = std.counts_wide.copy()
    if not counts.empty:
        counts = counts.set_index("Taxon")
        counts = counts.apply(pd.to_numeric, errors="coerce")
    else:
        counts = pd.DataFrame(index=pd.Index([], name="Taxon"))

    # ---- taxa: only if the user actually supplied trophic or c-p ----
    taxa = std.taxa.copy()
    has_traits = (taxa.get("trophic_group") is not None
                  and taxa["trophic_group"].notna().any()) or \
                 (taxa.get("cp") is not None and taxa["cp"].notna().any())
    taxa_out = None
    if has_traits:
        t = pd.DataFrame({"Taxon": taxa["taxon"]})
        mapped = [normalise_trophic(v) for v in taxa["trophic_group"]]
        t["trophic"] = [m[0] for m in mapped]
        for m in mapped:
            if m[1]:
                notes.append(m[1])
        t["cp"] = pd.to_numeric(taxa["cp"], errors="coerce")
        bad_cp = taxa["cp"].notna() & t["cp"].isna()
        if bad_cp.any():
            notes.append(f"{int(bad_cp.sum())} c-p value(s) were not numeric and were left blank "
                         f"for you to set in the Auto-assign tab.")
        t["source"] = taxa.get("source")
        t["provenance"] = taxa.get("provenance")
        taxa_out = t.set_index("Taxon")
    else:
        notes.append("No trophic group or c-p was supplied, so no `taxa` sheet was written. "
                     "Use the Auto-assign tab to propose them.")

    # ---- samples: only if there is metadata beyond the identifier ----
    samples = std.samples.copy()
    meta_cols = [c for c in samples.columns if c != "sample_id"]
    samples_out = None
    if meta_cols:
        s = samples.copy()
        s.columns = ["Sample"] + [re.sub(r"[^0-9a-zA-Z]+", "_", str(c)).strip("_").lower()
                                  for c in meta_cols]
        if group_col:
            src = re.sub(r"[^0-9a-zA-Z]+", "_", str(group_col)).strip("_").lower()
            if src in s.columns and src != "group":
                s["group"] = s[src]
                notes.append(f"'{group_col}' was copied to a column named `group`, which is the "
                             f"factor the comparison tabs use. The original column is kept too.")
            elif src not in s.columns:
                notes.append(f"'{group_col}' was not found among the metadata columns.")
        if "group" not in s.columns:
            notes.append("No column named `group` was produced, so the tabs that compare groups "
                         "will treat all samples as one set.")
        samples_out = s.set_index("Sample")
    else:
        notes.append("No sample metadata was found, so no `samples` sheet was written.")

    with pd.ExcelWriter(buf, engine="openpyxl") as xl:
        counts.to_excel(xl, sheet_name="counts")
        if taxa_out is not None:
            taxa_out.to_excel(xl, sheet_name="taxa")
        if samples_out is not None:
            samples_out.to_excel(xl, sheet_name="samples")
        std.counts_long.to_excel(xl, sheet_name="counts_long", index=False)
        conversion_log.to_excel(xl, sheet_name="conversion_log", index=False)
    return buf.getvalue(), notes


# --------------------------------------------------------------------------
# 9. Streamlit confirmation screen
# --------------------------------------------------------------------------

_ROLE_CHOICES = (["sample_id", "taxon", "count", "taxon_column", "sample_column"]
                 + list(METADATA_ROLES) + list(TAXON_ATTR_ROLES) + list(MEASUREMENT_ROLES)
                 + ["metadata", "ignore"])


def render_import_ui(st, key_prefix: str = "ingest"):
    """Upload -> inspect -> map -> confirm -> convert -> validate.

    Returns (workbook_bytes, label) once a conversion has passed validation, or
    (None, None). The result is held in session state, so it survives the reruns
    Streamlit performs on every later widget interaction.
    """
    done_key = f"{key_prefix}_done"

    if done_key in st.session_state:
        wb, label, summary = st.session_state[done_key]
        st.success(summary)
        c1, c2 = st.columns([3, 1])
        c1.download_button("Download standardized NemaNema dataset", wb,
                           "NemaNema_standardized_dataset.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                           key=f"{key_prefix}_dl2")
        if c2.button("Import a different file", key=f"{key_prefix}_reset"):
            del st.session_state[done_key]
            st.rerun()
        return wb, label

    st.subheader("Import a dataset in any layout")
    st.caption("Upload the spreadsheet you already have. The importer restructures it "
               "into a NEMO-NEMA workbook; it never edits your observations.")

    up = st.file_uploader("Excel or CSV", type=["xlsx", "xls", "csv", "tsv"],
                          key=f"{key_prefix}_file")
    if up is None:
        return None, None

    try:
        structure, book = detect_dataset_structure(up, up.name)
    except Exception as e:
        st.error(f"Could not read that file: {e}")
        return None, None

    st.markdown("#### Detected dataset")
    st.dataframe(pd.DataFrame([{
        "Sheet": p.name, "Rows": p.n_rows, "Columns": p.n_cols,
        "Likely layout": p.orientation, "Layout score": round(p.orientation_score, 2),
        "Sheet holds": p.content_role, "Sheet score": round(p.content_score, 2),
    } for p in structure.sheets]), use_container_width=True, hide_index=True)

    for n in structure.notes:
        st.warning(n)

    names = [p.name for p in structure.sheets]
    idx = names.index(structure.primary_sheet) if structure.primary_sheet in names else 0
    sheet = st.selectbox("Sheet holding the observations", names, index=idx,
                         key=f"{key_prefix}_sheet")
    df = book[sheet]
    prof = structure.profile(sheet)

    st.markdown("#### Preview")
    st.dataframe(df.head(8), use_container_width=True)

    layouts = ["long", "wide_samples_rows", "wide_taxa_rows"]
    li = layouts.index(prof.orientation) if prof.orientation in layouts else 0
    orientation = st.selectbox(
        "Layout", layouts, index=li, key=f"{key_prefix}_layout",
        format_func=lambda x: {"long": "Long (one row per sample x taxon)",
                               "wide_samples_rows": "Wide (samples in rows, taxa in columns)",
                               "wide_taxa_rows": "Matrix (taxa in rows, samples in columns)"}[x])

    mapping = map_columns(df, orientation)

    st.markdown("#### Dataset interpretation")
    st.caption("Scores are heuristic match strengths from name and value patterns — "
               "not statistical confidence. Anything below 0.60 deserves a look.")
    base = mapping.as_frame()
    edited = st.data_editor(
        base,
        column_config={"Detected role": st.column_config.SelectboxColumn(
            "Detected role", options=_ROLE_CHOICES, required=True)},
        disabled=["Original column", "Match score", "Evidence"],
        hide_index=True, use_container_width=True, key=f"{key_prefix}_map")

    if not edited.equals(base):
        lookup = {r["Original column"]: r["Detected role"] for _, r in edited.iterrows()}
        new_maps, new_taxa = [], []
        for m in mapping.mappings:
            role = lookup.get(m.column, m.role)
            changed = role != m.role
            new_maps.append(ColumnMapping(m.column, role,
                                          1.0 if changed else m.score,
                                          "set manually" if changed else m.evidence))
            if role == "taxon_column":
                new_taxa.append(m.column)
        mapping = MappingResult(orientation, new_maps, new_taxa, mapping.warnings, False)

    for w in mapping.warnings:
        st.warning(w)

    meta_choices = [m.column for m in mapping.mappings
                    if m.role in METADATA_ROLES + ("metadata",)]
    default_group = mapping.first("group") or (mapping.first("treatment") or "")
    group_col = st.selectbox(
        "Which column is the grouping factor for comparisons?",
        ["(none)"] + meta_choices,
        index=(meta_choices.index(default_group) + 1) if default_group in meta_choices else 0,
        key=f"{key_prefix}_group",
        help="This becomes the `group` column the Statistics, Multivariate and "
             "Summary tabs compare across. Leave as (none) if the samples are not "
             "split into groups.")

    c1, c2 = st.columns(2)
    blank_policy = c1.selectbox(
        "Empty count cells", ["absent", "zero"], key=f"{key_prefix}_blank",
        format_func=lambda x: {"absent": "Treat as absent (recommended)",
                               "zero": "Write explicit zeros"}[x])
    duplicate_policy = c2.selectbox(
        "Repeated sample x taxon rows", ["flag", "sum"], key=f"{key_prefix}_dup",
        format_func=lambda x: {"flag": "Stop and report (recommended)",
                               "sum": "Sum them"}[x])

    if not st.button("Accept detected structure and convert", key=f"{key_prefix}_go"):
        return None, None

    std = convert_to_nemonema_schema(df, mapping, blank_policy, duplicate_policy,
                                     {"filename": structure.filename, "sheet": sheet})
    report = validate_conversion(df, mapping, std)

    st.markdown("#### Data integrity check")
    for c in report.checks:
        line = f"{c.name}: original {c.original}, converted {c.converted}"
        if c.detail:
            line += f" — {c.detail}"
        (st.success if c.passed else st.error)(("PASS — " if c.passed else "FAIL — ") + line)

    log_df = generate_conversion_report(structure, sheet, mapping, std, report)

    if not report.passed:
        st.error("Conversion failed validation. Nothing has been analysed — correct the "
                 "source file or the mapping above and convert again.")
        st.download_button("Download conversion log", log_df.to_csv(index=False).encode(),
                           "NemaNema_conversion_log.csv", "text/csv",
                           key=f"{key_prefix}_faillog")
        return None, None

    wb, notes = build_standardized_workbook(
        std, log_df, None if group_col == "(none)" else group_col)

    for n in notes:
        st.info(n)
    with st.expander("Conversion log"):
        st.dataframe(log_df, use_container_width=True, hide_index=True)

    summary = (f"Imported {structure.filename}: {std.counts_long['sample_id'].nunique()} samples, "
               f"{std.counts_long['taxon'].nunique()} taxa, {len(std.counts_long)} observations. "
               f"Your uploaded file is untouched.")
    st.session_state[done_key] = (wb, f"imported::{structure.filename}", summary)
    st.rerun()
