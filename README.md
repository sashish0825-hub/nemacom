# NEMO-NEMA

**N**ematode **E**cological **M**etrics and **O**rdination for
**N**ematode **E**cosystem **M**onitoring and **A**ssessment

An open, browser-based tool for nematode community and soil-health analysis.
Upload one Excel workbook; get community descriptors, diversity indices, soil
food-web indices, metabolic footprints, a soil health score, design-aware
statistics, soil-nutrient correlation and a full PDF report.

Developed by Ashish Kumar Singh, Kavita Jain, Vishal Singh Somvanshi,
Rashid Pervez, Anil Sirohi and Pankaj.
Division of Nematology, ICAR-Indian Agricultural Research Institute,
New Delhi 110012.

## What it computes

| Module | Contents |
|---|---|
| Community | Norton (1978): frequency, density, prominence and importance values |
| Diversity | Shannon, Simpson, Pielou, Margalef, richness |
| Food web | MI, MI(2–5), σMI, PPI, EI, SI, CI, BI, NCR, faunal profile quadrats |
| Soil health | Nematode Soil Health index (Ghaderi et al. 2025) with calibration checks |
| Footprints | Andrássy (1956) biomass; Ferris (2010) metabolic footprints |
| Statistics | ANOVA, Kruskal–Wallis, assumption diagnostics, Tukey, Dunn, FDR correction |
| Multivariate | Bray–Curtis, PCA, PERMANOVA, PERMDISP |
| Soil nutrients | Correlation with soil chemistry, RDA, forward selection, soil PCA |
| Reference | 188 genera with family, trophic group and c-p value |

## Verification

Computations are checked against the primary sources:

- Ferris (2010) Table 4 metabolic footprints — all three worked samples reproduce
- Ghaderi et al. (2025) NSH scoring — all 294 published samples reproduce exactly
- Ferris et al. (2001) faunal profile quadrats — verified against Fig. 1 and Table 5

**Not verified:** two prominence-value variants circulate, both attributed to
Norton (1978). The variant is a user choice and is recorded in the methods log.

## Input

Three sheets: `counts` (taxa × samples), `taxa` (trophic, cp, source),
`samples` (group, and field/plot where available). An optional `nutrients`
sheet enables the soil-correlation module. If `taxa` is incomplete the
Auto-assign tab proposes trophic group and c-p from taxon names, with source
and confidence on every row; nothing is applied without confirmation.

## Run locally

```bash
pip install -r requirements.txt
streamlit run nemonema_app.py
```

## Citing

See the Cite tab in the app, or `CITATION.cff`. Cite the primary sources of
whatever you report as well as the software — it implements published methods
and does not replace them.

## Licence

MIT. Trait assignments are attributed to their published sources; no
third-party trait database is redistributed.
