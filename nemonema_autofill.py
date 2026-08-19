"""
nemonema_autofill.py — propose trophic group and c-p value from a taxon name.

WHY THIS IS SAFE TO BUNDLE, AND WHAT IS NOT
Trophic group and c-p value are facts published in Yeates et al. (1993) J
Nematol 25:315-331 and Bongers & Bongers (1998) Appl Soil Ecol 10:239-251, and
restated in dozens of later papers including Ferris (2010) Table 1. Reproducing
them with attribution is reproducing standard science.

Body weights are NOT bundled. Those in Ferris (2010) Table 1 are his own
computation from 1368 species measurements — original work, not a restatement.
Users who need metabolic footprints transcribe what they need themselves, or
measure their own specimens, which is better anyway.

WHY EVERY PROPOSAL STILL NEEDS CONFIRMATION
Assignments here are at FAMILY level. Within-family variation is real, and
Bongers & Bongers (1998) themselves flag disputed cases. A silent lookup would
turn a genus misspelling or a family-level approximation into a wrong maturity
index that looks entirely normal. So the module proposes; the user confirms.
"""

from __future__ import annotations

import difflib
import numpy as np
import pandas as pd

SRC_STD = ("Yeates et al. (1993) J Nematol 25:315-331 (feeding habit); "
           "Bongers & Bongers (1998) Appl Soil Ecol 10:239-251 (c-p); "
           "as tabulated in Ferris (2010) Eur J Soil Biol 46:97-104 Table 1")
SRC_PROV = "PROVISIONAL — family not tabulated in Ferris (2010) Table 1. VERIFY."
SRC_CONT = ("CONTESTED — sources disagree. Ferris (2010) Table 1 gives feeding "
            "habit 1 (plant-feeding) for Aphelenchoididae; most soil "
            "Aphelenchoides are fungal feeders and the published classification "
            "of Cerevkova et al. (2024) treats them as such. Decide which "
            "applies to your species and state it in your methods.")

# genus -> (family, trophic, cp, confidence)
REFERENCE = {
    # plant parasites
    "meloidogyne": ("Meloidogynidae", "PP", 3, "VERIFY"),
    "heterodera": ("Heteroderidae", "PP", 3, "OK"),
    "globodera": ("Heteroderidae", "PP", 3, "OK"),
    "rotylenchulus": ("Rotylenchulidae", "PP", 3, "OK"),
    "tylenchulus": ("Tylenchulidae", "PP", 3, "VERIFY"),
    "hoplolaimus": ("Hoplolaimidae", "PP", 3, "OK"),
    "helicotylenchus": ("Hoplolaimidae", "PP", 3, "OK"),
    "rotylenchus": ("Hoplolaimidae", "PP", 3, "OK"),
    "scutellonema": ("Hoplolaimidae", "PP", 3, "OK"),
    "pratylenchus": ("Pratylenchidae", "PP", 3, "OK"),
    "radopholus": ("Pratylenchidae", "PP", 3, "OK"),
    "hirschmanniella": ("Pratylenchidae", "PP", 3, "OK"),
    "tylenchorhynchus": ("Telotylenchidae", "PP", 3, "OK"),
    "merlinius": ("Telotylenchidae", "PP", 3, "OK"),
    "quinisulcius": ("Telotylenchidae", "PP", 3, "OK"),
    "dolichodorus": ("Dolichodoridae", "PP", 3, "OK"),
    "xiphinema": ("Longidoridae", "PP", 5, "OK"),
    "longidorus": ("Longidoridae", "PP", 5, "OK"),
    "trichodorus": ("Trichodoridae", "PP", 4, "OK"),
    "paratrichodorus": ("Trichodoridae", "PP", 4, "OK"),
    "criconemoides": ("Criconematidae", "PP", 3, "OK"),
    "criconemella": ("Criconematidae", "PP", 3, "OK"),
    "hemicriconemoides": ("Criconematidae", "PP", 3, "OK"),
    "hemicycliophora": ("Hemicycliophoridae", "PP", 3, "OK"),
    "paratylenchus": ("Paratylenchidae", "PP", 2, "OK"),
    "tylenchus": ("Tylenchidae", "PP", 2, "OK"),
    "filenchus": ("Tylenchidae", "PP", 2, "OK"),
    "basiria": ("Tylenchidae", "PP", 2, "OK"),
    "psilenchus": ("Psilenchidae", "PP", 2, "OK"),
    "ditylenchus": ("Anguinidae", "PP", 2, "OK"),
    "anguina": ("Anguinidae", "PP", 2, "OK"),
    # CONTESTED. Ferris (2010) Table 1 records Aphelenchoididae with feeding
    # habit 1 (plant-feeding). Most soil Aphelenchoides are fungal feeders, and
    # validating NEMO-NEMA against Cerevkova et al. (2024) as published by
    # Ghaderi et al. (2025) showed the original authors classified it FF:
    # treating it as PP gave CI r = 0.931 with a bias of -4.76 units, while FF
    # gave r = 0.999 with a bias of -0.21. Defaulting to FF and flagging the
    # choice, because the genus contains both foliar plant parasites and fungal
    # feeders and the right call depends on which species are present.
    "aphelenchoides": ("Aphelenchoididae", "FF", 2, "CONTESTED"),
    # fungivores
    "aphelenchus": ("Aphelenchidae", "FF", 2, "OK"),
    "tylencholaimus": ("Leptonchidae", "FF", 4, "OK"),
    "dorylaimoides": ("Leptonchidae", "FF", 4, "OK"),
    "diphtherophora": ("Diphtherophoridae", "FF", 3, "OK"),
    # bacterivores
    "rhabditis": ("Rhabditidae", "BF", 1, "OK"),
    "rhabditidae": ("Rhabditidae", "BF", 1, "OK"),
    "mesorhabditis": ("Mesorhabditidae", "BF", 1, "OK"),
    "panagrolaimus": ("Panagrolaimidae", "BF", 1, "OK"),
    "diplogaster": ("Diplogasteridae", "BF", 1, "OK"),
    "cephalobus": ("Cephalobidae", "BF", 2, "OK"),
    "eucephalobus": ("Cephalobidae", "BF", 2, "OK"),
    "acrobeloides": ("Cephalobidae", "BF", 2, "OK"),
    "acrobeles": ("Cephalobidae", "BF", 2, "OK"),
    "chiloplacus": ("Cephalobidae", "BF", 2, "OK"),
    "plectus": ("Plectidae", "BF", 2, "OK"),
    "plectidae": ("Plectidae", "BF", 2, "OK"),
    "monhystera": ("Monhysteridae", "BF", 2, "OK"),
    "prismatolaimus": ("Prismatolaimidae", "BF", 3, "OK"),
    "teratocephalus": ("Teratocephalidae", "BF", 3, "OK"),
    "wilsonema": ("Wilsonematidae", "BF", 2, "VERIFY"),
    "alaimus": ("Alaimidae", "BF", 4, "OK"),
    # omnivores
    "dorylaimus": ("Dorylaimidae", "OM", 4, "OK"),
    "mesodorylaimus": ("Dorylaimidae", "OM", 4, "OK"),
    "eudorylaimus": ("Qudsianematidae", "OM", 4, "OK"),
    "labronema": ("Qudsianematidae", "OM", 4, "OK"),
    "aporcelaimellus": ("Aporcelaimidae", "OM", 5, "OK"),
    "aporcelaimus": ("Aporcelaimidae", "OM", 5, "OK"),
    "thornenema": ("Thornenematidae", "OM", 5, "OK"),
    "chrysonema": ("Chrysonematidae", "OM", 5, "OK"),
    # predators
    "mononchus": ("Mononchidae", "PR", 4, "OK"),
    "clarkus": ("Mononchidae", "PR", 4, "OK"),
    "mylonchulus": ("Mylonchulidae", "PR", 4, "OK"),
    "iotonchus": ("Iotonchidae", "PR", 4, "OK"),
    "discolaimus": ("Discolaimidae", "PR", 5, "OK"),
    "nygolaimus": ("Nygolaimidae", "PR", 5, "OK"),
    "actinolaimus": ("Actinolaimidae", "PR", 5, "OK"),
    "tripyla": ("Tripylidae", "PR", 3, "OK"),

    # ---- expansion: additional genera reported from soils worldwide.
    # Family-level cp and feeding habit as tabulated in Ferris (2010) Table 1,
    # which attributes cp to Bongers & Bongers (1998) and feeding habit to
    # Yeates et al. (1993). Families absent from that table are flagged VERIFY.
    # plant parasites
    "meloidoderita": ("Meloidoderitidae", "PP", 3, "VERIFY"),
    "punctodera": ("Heteroderidae", "PP", 3, "OK"),
    "cactodera": ("Heteroderidae", "PP", 3, "OK"),
    "afenestrata": ("Heteroderidae", "PP", 3, "OK"),
    "aphasmatylenchus": ("Hoplolaimidae", "PP", 3, "OK"),
    "aorolaimus": ("Hoplolaimidae", "PP", 3, "OK"),
    "peltamigratus": ("Hoplolaimidae", "PP", 3, "OK"),
    "zygotylenchus": ("Pratylenchidae", "PP", 3, "OK"),
    "pratylenchoides": ("Pratylenchidae", "PP", 3, "OK"),
    "nacobbus": ("Pratylenchidae", "PP", 3, "OK"),
    "hoplotylus": ("Pratylenchidae", "PP", 3, "OK"),
    "bitylenchus": ("Telotylenchidae", "PP", 3, "OK"),
    "trophurus": ("Telotylenchidae", "PP", 3, "OK"),
    "telotylenchus": ("Telotylenchidae", "PP", 3, "OK"),
    "neodolichodorus": ("Dolichodoridae", "PP", 3, "OK"),
    "belonolaimus": ("Dolichodoridae", "PP", 3, "OK"),
    "paralongidorus": ("Longidoridae", "PP", 5, "OK"),
    "xiphidorus": ("Longidoridae", "PP", 5, "OK"),
    "monotrichodorus": ("Trichodoridae", "PP", 4, "OK"),
    "nanidorus": ("Trichodoridae", "PP", 4, "OK"),
    "mesocriconema": ("Criconematidae", "PP", 3, "OK"),
    "discocriconemella": ("Criconematidae", "PP", 3, "OK"),
    "ogma": ("Criconematidae", "PP", 3, "OK"),
    "criconema": ("Criconematidae", "PP", 3, "OK"),
    "caloosia": ("Hemicycliophoridae", "PP", 3, "OK"),
    "gracilacus": ("Paratylenchidae", "PP", 2, "OK"),
    "cacopaurus": ("Paratylenchidae", "PP", 2, "OK"),
    "coslenchus": ("Tylenchidae", "PP", 2, "OK"),
    "aglenchus": ("Tylenchidae", "PP", 2, "OK"),
    "malenchus": ("Tylenchidae", "PP", 2, "OK"),
    "boleodorus": ("Tylenchidae", "PP", 2, "OK"),
    "neopsilenchus": ("Psilenchidae", "PP", 2, "OK"),
    "atylenchus": ("Atylenchidae", "PP", 2, "OK"),
    "ecphyadophora": ("Ecphyadophoridae", "PP", 2, "OK"),
    "subanguina": ("Anguinidae", "PP", 2, "OK"),
    "nothanguina": ("Anguinidae", "PP", 2, "OK"),
    "hemicaloosia": ("Hemicycliophoridae", "PP", 3, "OK"),
    "rotylenchoides": ("Hoplolaimidae", "PP", 3, "OK"),
    "helicotylenchoides": ("Hoplolaimidae", "PP", 3, "OK"),
    # fungivores
    "aphelenchoides_group": ("Aphelenchoididae", "FF", 2, "OK"),
    "paraphelenchus": ("Aphelenchidae", "FF", 2, "OK"),
    "seinura": ("Aphelenchoididae", "PR", 2, "OK"),
    "tylenchodorus": ("Leptonchidae", "FF", 4, "OK"),
    "botalium": ("Leptonchidae", "FF", 4, "OK"),
    "doryllium": ("Leptonchidae", "FF", 4, "OK"),
    "diphtherophoroides": ("Diphtherophoridae", "FF", 3, "OK"),
    "tylolaimophorus": ("Diphtherophoridae", "FF", 3, "OK"),
    "nothotylenchus": ("Neotylenchidae", "FF", 2, "OK"),
    "deladenus": ("Neotylenchidae", "FF", 2, "OK"),
    # bacterivores
    "caenorhabditis": ("Rhabditidae", "BF", 1, "OK"),
    "pelodera": ("Peloderidae", "BF", 1, "OK"),
    "protorhabditis": ("Protorhabditidae", "BF", 1, "OK"),
    "cruznema": ("Rhabditidae", "BF", 1, "OK"),
    "bursilla": ("Rhabditidae", "BF", 1, "OK"),
    "diploscapter": ("Diploscapteridae", "BF", 1, "OK"),
    "bunonema": ("Bunonematidae", "BF", 1, "OK"),
    "cylindrolaimus": ("Cylindrolaimidae", "BF", 1, "OK"),
    "myolaimus": ("Myolaimidae", "BF", 2, "OK"),
    "acrobeloides_group": ("Cephalobidae", "BF", 2, "OK"),
    "cervidellus": ("Cephalobidae", "BF", 2, "OK"),
    "heterocephalobus": ("Cephalobidae", "BF", 2, "OK"),
    "zeldia": ("Cephalobidae", "BF", 2, "OK"),
    "stegelletina": ("Cephalobidae", "BF", 2, "OK"),
    "seleborca": ("Cephalobidae", "BF", 2, "OK"),
    "nothacrobeles": ("Cephalobidae", "BF", 2, "OK"),
    "panagrellus": ("Panagrolaimidae", "BF", 1, "OK"),
    "turbatrix": ("Panagrolaimidae", "BF", 1, "OK"),
    "anaplectus": ("Plectidae", "BF", 2, "OK"),
    "ceratoplectus": ("Plectidae", "BF", 2, "OK"),
    "wilsonema_group": ("Wilsonematidae", "BF", 2, "VERIFY"),
    "chronogaster": ("Chronogastridae", "BF", 3, "OK"),
    "rhabdolaimus": ("Rhabdolaimidae", "BF", 3, "OK"),
    "aphanolaimus": ("Aphanolaimidae", "BF", 3, "OK"),
    "leptolaimus": ("Leptolaimidae", "BF", 3, "OK"),
    "bastiania": ("Bastianiidae", "BF", 3, "OK"),
    "aulolaimus": ("Aulolaimidae", "BF", 3, "OK"),
    "euteratocephalus": ("Metateratocephalidae", "BF", 3, "OK"),
    "metateratocephalus": ("Metateratocephalidae", "BF", 3, "OK"),
    "geomonhystera": ("Monhysteridae", "BF", 2, "OK"),
    "eumonhystera": ("Monhysteridae", "BF", 2, "OK"),
    "amphidelus": ("Amphidelidae", "BF", 4, "OK"),
    "paramphidelus": ("Amphidelidae", "BF", 4, "OK"),
    # omnivores
    "aporcelaimium": ("Aporcelaimidae", "OM", 5, "OK"),
    "epacrolaimus": ("Aporcelaimidae", "OM", 5, "OK"),
    "makatinus": ("Aporcelaimidae", "OM", 5, "OK"),
    "prodorylaimus": ("Dorylaimidae", "OM", 4, "OK"),
    "laimydorus": ("Dorylaimidae", "OM", 4, "OK"),
    "thonus": ("Qudsianematidae", "OM", 4, "OK"),
    "microdorylaimus": ("Qudsianematidae", "OM", 4, "OK"),
    "allodorylaimus": ("Qudsianematidae", "OM", 4, "OK"),
    "crassolabium": ("Qudsianematidae", "OM", 4, "OK"),
    "ecumenicus": ("Qudsianematidae", "OM", 4, "OK"),
    "pungentus": ("Nordiidae", "OM", 4, "OK"),
    "enchodelus": ("Nordiidae", "OM", 4, "OK"),
    "longidorella": ("Nordiidae", "OM", 4, "OK"),
    "axonchium": ("Belondiridae", "PP", 5, "OK"),
    "dorylaimellus": ("Dorylaimellidae", "PP", 5, "OK"),
    "oxydirus": ("Belondiridae", "PP", 5, "OK"),
    "carcharolaimus": ("Crateronematidae", "OM", 4, "OK"),
    "mydonomus": ("Mydonomidae", "OM", 5, "OK"),
    "thornia": ("Thorniidae", "OM", 4, "OK"),
    "campydora": ("Campydoridae", "OM", 3, "OK"),
    "isolaimium": ("Isolaimidae", "BF", 5, "OK"),
    # predators
    "prionchulus": ("Mononchidae", "PR", 4, "OK"),
    "coomansus": ("Mononchidae", "PR", 4, "OK"),
    "anatonchus": ("Anatonchidae", "PR", 4, "OK"),
    "miconchus": ("Anatonchidae", "PR", 4, "OK"),
    "iotonchus_group": ("Iotonchidae", "PR", 4, "OK"),
    "sporonchulus": ("Mylonchulidae", "PR", 4, "OK"),
    "cobbonchus": ("Cobbonchidae", "PR", 4, "OK"),
    "discolaimoides": ("Discolaimidae", "PR", 5, "OK"),
    "discolaimium": ("Discolaimidae", "PR", 5, "OK"),
    "nygolaimellus": ("Nygolaimidae", "PR", 5, "OK"),
    "aquatides": ("Nygolaimidae", "PR", 5, "OK"),
    "paravulvus": ("Nygolaimidae", "PR", 5, "OK"),
    "paractinolaimus": ("Actinolaimidae", "PR", 5, "OK"),
    "trachypleurosum": ("Actinolaimidae", "PR", 5, "OK"),
    "tripylina": ("Tripylidae", "PR", 3, "OK"),
    "trischistoma": ("Tripylidae", "PR", 3, "OK"),
    "tobrilus": ("Tobrilidae", "PR", 3, "OK"),
    "ironus": ("Ironidae", "PR", 4, "OK"),
    "prismatolaimus_group": ("Prismatolaimidae", "BF", 3, "OK"),
    "onchulus": ("Onchulidae", "PR", 3, "OK"),
    "mononchulus": ("Bathyodontidae", "BF", 4, "OK"),

    # ---- genera present in Cerevkova et al. (2024) Data in Brief 57:111098
    # that were absent from the lookup. Family placements below are given with
    # the confidence I can actually support; VERIFY means confirm before use.
    "acrolobus": ("Cephalobidae", "BF", 2, "OK"),
    "pseudacrobeles": ("Cephalobidae", "BF", 2, "OK"),
    "amplimerlinius": ("Telotylenchidae", "PP", 3, "OK"),
    "geocenamus": ("Telotylenchidae", "PP", 3, "OK"),
    "nagelus": ("Telotylenchidae", "PP", 3, "OK"),
    "lelenchus": ("Tylenchidae", "PP", 2, "OK"),
    "macropostonia": ("Criconematidae", "PP", 3, "OK"),
    "epidorylaimus": ("Qudsianematidae", "OM", 4, "OK"),
    "eudiplogaster": ("Diplogasteridae", "BF", 1, "OK"),
    "tylencholaimellus": ("Leptonchidae", "FF", 4, "VERIFY"),
    # Punctodora is a CHROMADORID, not the cyst nematode Punctodera. Ferris
    # (2010) Table 1 gives Chromadoridae cp 3, feeding habit 6 = unicellular
    # eukaryote feeder, which has no exact equivalent among PP/BF/FF/OM/PR.
    # Recorded here as BF with a VERIFY flag; decide for yourself how to treat
    # unicellular-eukaryote feeders and say so in your methods.
    "punctodora": ("Chromadoridae", "BF", 3, "VERIFY"),
    "ereptonema": ("Rhabditidae", "BF", 1, "VERIFY"),
    "odontolaimus": ("Odontolaimidae", "BF", 3, "VERIFY"),
    "paraxonchium": ("Paraxonchiidae", "OM", 5, "VERIFY"),
}

TROPHIC_NAME = {"PP": "plant parasite", "BF": "bacterivore", "FF": "fungivore",
                "OM": "omnivore", "PR": "predator"}


def _key(name: str) -> str:
    """First word of the taxon name, lowercased. 'Meloidogyne incognita' -> 'meloidogyne'."""
    return str(name).strip().split()[0].lower().rstrip(".,;") if str(name).strip() else ""


def propose(taxon_names, cutoff: float = 0.86) -> pd.DataFrame:
    """Propose trophic and c-p for each name. Nothing is applied automatically.

    match_type
      exact   name matched a reference entry outright
      fuzzy   close but not identical — usually a spelling difference. CONFIRM.
      none    no match. The user must supply the values.

    confidence
      OK       family tabulated in the cited source
      VERIFY   family not in that table; provisional, confirm before use
    """
    rows = []
    keys = list(REFERENCE)
    for name in taxon_names:
        k = _key(name)
        if k in REFERENCE:
            fam, tro, cp, conf = REFERENCE[k]
            rows.append({"taxon": name, "matched_genus": k, "match_type": "exact",
                         "family": fam, "trophic": tro, "cp": cp,
                         "trophic_name": TROPHIC_NAME[tro],
                         "confidence": conf,
                         "source": (SRC_STD if conf == "OK" else
                                    SRC_CONT if conf == "CONTESTED" else SRC_PROV),
                         "action": ("accept" if conf == "OK" else
                                    "CHOOSE — sources disagree"
                                    if conf == "CONTESTED" else "review")})
            continue
        near = difflib.get_close_matches(k, keys, n=1, cutoff=cutoff)
        if near:
            fam, tro, cp, conf = REFERENCE[near[0]]
            rows.append({"taxon": name, "matched_genus": near[0], "match_type": "fuzzy",
                         "family": fam, "trophic": tro, "cp": cp,
                         "trophic_name": TROPHIC_NAME[tro],
                         "confidence": "CHECK SPELLING",
                         "source": SRC_STD if conf == "OK" else SRC_PROV,
                         "action": "CONFIRM — name differs from the reference"})
        else:
            rows.append({"taxon": name, "matched_genus": "", "match_type": "none",
                         "family": "", "trophic": "", "cp": np.nan,
                         "trophic_name": "", "confidence": "NO MATCH",
                         "source": "",
                         "action": "enter trophic and cp yourself, with a source"})
    return pd.DataFrame(rows)


def summarise(proposal: pd.DataFrame) -> dict:
    """Counts by match type, and whether the dataset is ready to analyse."""
    n = len(proposal)
    exact = int((proposal["match_type"] == "exact").sum())
    fuzzy = int((proposal["match_type"] == "fuzzy").sum())
    none = int((proposal["match_type"] == "none").sum())
    verify = int((proposal["confidence"] == "VERIFY").sum())
    contested = int((proposal["confidence"] == "CONTESTED").sum())
    return {"n_taxa": n, "exact": exact, "fuzzy": fuzzy, "unmatched": none,
            "provisional": verify,
            "contested": contested,
            "ready": none == 0,
            "needs_attention": fuzzy + none + verify + contested,
            "message": (f"{exact} of {n} matched exactly."
                        + (f" {fuzzy} matched only approximately — check the spelling."
                           if fuzzy else "")
                        + (f" {none} did not match and must be filled in by hand."
                           if none else "")
                        + (f" {verify} are provisional and need verifying against a "
                           "primary source." if verify else ""))}


def apply_proposal(taxa_df: pd.DataFrame, proposal: pd.DataFrame,
                   overwrite: bool = False) -> pd.DataFrame:
    """Write accepted proposals into a taxa table.

    By default this fills only BLANK cells: a value the user typed is never
    silently replaced by a lookup.
    """
    out = taxa_df.copy()
    p = proposal.set_index("taxon")
    for col in ("trophic", "cp", "family", "source"):
        if col not in out.columns:
            out[col] = pd.Series([np.nan] * len(out), index=out.index, dtype=object)
        elif out[col].isna().all():
            # an all-NaN column is float64; assigning a string into it raises.
            # Cast to object first so text assignments are legal.
            out[col] = out[col].astype(object)
    for tx in out.index:
        if tx not in p.index or p.loc[tx, "match_type"] == "none":
            continue
        for col in ("trophic", "cp", "family", "source"):
            cur = out.loc[tx, col]
            blank = pd.isna(cur) or (isinstance(cur, str) and cur.strip() == "")
            if overwrite or blank:
                out.loc[tx, col] = p.loc[tx, col]
    # cp back to numeric where possible, without breaking blanks
    out["cp"] = pd.to_numeric(out["cp"], errors="coerce")
    return out


def reference_table() -> pd.DataFrame:
    """The full lookup as a browsable table, with morphometric columns blank."""
    rows = []
    for gen, (fam, tro, cp, conf) in sorted(REFERENCE.items()):
        if gen.endswith("_group"):
            continue
        rows.append({
            "genus": gen.capitalize(), "family": fam, "trophic": tro,
            "trophic_name": TROPHIC_NAME[tro], "cp": cp,
            "trait_confidence": conf,
            "trait_source": (SRC_STD if conf == "OK" else
                             SRC_CONT if conf == "CONTESTED"
                             else SRC_PROV.format(fam=fam)),
            "length_um": np.nan, "diameter_um": np.nan, "a_ratio": np.nan,
            "n_measured": np.nan, "stage_measured": "",
            "morphometry_source": "", "family_weight_ug": np.nan,
            "weight_source": "", "notes": "",
        })
    return pd.DataFrame(rows)
