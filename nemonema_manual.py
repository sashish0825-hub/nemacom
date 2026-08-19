"""nemonema_manual.py — the complete user manual, rendered in the Manual tab."""

MANUAL = r"""
# NEMO-NEMA user manual

**N**ematode **E**cological **M**etrics and **O**rdination for
**N**ematode **E**cosystem **M**onitoring and **A**ssessment · version 1.1

Division of Nematology, ICAR-Indian Agricultural Research Institute, New Delhi

---

## Part 1 · What this does, and what it does not

You supply one Excel workbook of nematode counts. NEMO-NEMA returns community
descriptors, diversity indices, soil food-web indices, metabolic footprints, a
soil health score, design-aware statistics, soil-nutrient correlation, and a
full PDF report.

**It does not** identify your nematodes, decide their trophic group or c-p
value, or recommend management. Those remain yours. The software computes; the
biology is your responsibility.

**Nothing here replaces the primary literature.** Every index implemented is
someone else's published method — cite them, not only this tool. See the Cite
tab.

---

## Part 2 · Getting in

**On the web.** Open the app link in any browser, on a phone or computer.
Nothing to install.

**On your own computer.** Needs Python. In a terminal:

```
cd ~/nemonema
source venv/bin/activate
streamlit run nemonema_app.py
```

Then open `localhost:8501`. Leave that terminal running while you work; press
Ctrl+C when finished. Use a second terminal for anything else.

---

## Part 3 · The data format

### Sheet `counts` — required

Taxa down the first column, samples across the top row. Each cell is the number
of individuals of that taxon recovered from that sample.

| Taxon | S01 | S02 | S03 | S04 |
|---|---|---|---|---|
| Meloidogyne | 180 | 143 | 96 | 70 |
| Helicotylenchus | 41 | 12 | 0 | 55 |
| Aphelenchus | 22 | 31 | 18 | 60 |
| Cephalobus | 120 | 96 | 143 | 210 |

Four rules:

1. **One consistent soil basis for every sample.** Per 100 g dry soil is
   strongly preferred, because the soil health index assumes it. If you counted
   per 200 cc or per 250 g, convert before entering. The software cannot detect
   a mixed basis.
2. **Blank or 0 means absent.** Never enter text, never enter a negative number.
3. **Keep taxonomic resolution consistent.** Do not enter a family *and* a genus
   that sits inside it — `Plectidae` and `Plectus` together double-counts.
4. **Sample names must be unique.** Abbreviating state names is a common way to
   create duplicates: Uttarakhand and Uttar Pradesh both shorten to `UTT`. The
   software checks for this and refuses to run.

### Sheet `taxa` — recommended

One row per taxon in `counts`, spelled **exactly** the same.

| Taxon | trophic | cp | source |
|---|---|---|---|
| Meloidogyne | PP | 3 | Bongers & Bongers (1998) |
| Helicotylenchus | PP | 3 | Ferris (2010) Table 1 |
| Aphelenchus | FF | 2 | Ferris (2010) Table 1 |
| Cephalobus | BF | 2 | Ferris (2010) Table 1 |

- **trophic** — PP plant parasite, BF bacterivore, FF fungivore, OM omnivore,
  PR predator. Source: Yeates et al. (1993).
- **cp** — colonizer–persister value, a whole number 1 to 5. Source: Bongers
  (1990); Bongers & Bongers (1998).
- **source** — where you got each assignment.

That last column matters more than it looks. Every maturity and food-web index
is calculated entirely from `trophic` and `cp`. A wrong assignment produces a
plausible-looking wrong index, and no software can detect it.

**If this sheet is missing or incomplete**, the Auto-assign tab proposes values
from the taxon names. See Part 5.

### Sheet `samples` — recommended

One row per sample column in `counts`, names matching exactly.

| Sample | group | field | replicate |
|---|---|---|---|
| S01 | Organic | Punjab_F1 | 1 |
| S02 | Organic | Punjab_F2 | 2 |
| S03 | Inorganic | Punjab_F3 | 1 |
| S04 | Inorganic | Punjab_F4 | 2 |

- **group** — the factor you are comparing: treatment, land use, site, season.
  Without it, no statistics can run.
- **field**, **plot**, **block** — strongly recommended. Without them the
  software cannot tell an independent replicate from a subsample of the same
  plot. If three rows share a plot they are subsamples, and treating them as
  replicates makes every p-value too small. This is the commonest fatal design
  error in nematode community papers, and it cannot be fixed after the fact.

Aim for **at least 3 independent replicates per group**; 4 or 5 is safer.

### Sheet `nutrients` — optional

Soil chemistry, physics and biology. Enables the Soil nutrients tab.

| Sample | pH | organic_C | NO3_N | Olsen_P | microbial_biomass_C |
|---|---|---|---|---|---|
| S01 | 6.5 | 1.88 | 63 | 49 | 620 |
| S02 | 5.6 | 1.72 | 77 | 40 | 540 |

Any numeric column is treated as a soil variable. You may instead put extra
numeric columns straight into the `samples` sheet.

### Optional columns for biomass

Add to the `taxa` sheet. Only needed for the Footprints tab; every other index
works without them.

- **Route A, best** — `length_um` and `diameter_um` from your own measurements,
  plus `n_measured` and `stage_measured`. Measure 20–30 individuals per taxon
  across the life stages actually present.
- **Route B** — `length_um` and `a_ratio` from a taxonomic description, where
  de Man's *a* = length ÷ greatest diameter. Diameter is derived.
- **Route C, weakest** — `family_weight_ug` and `weight_source` from Ferris
  (2010) Table 1. Note his own caveats: family means hide large variation
  (Dorylaimidae CV near 200%) and derive from adults, so juvenile-rich field
  samples are overestimated.

---

## Part 4 · Running the analysis

**Step 1 — upload.** Sidebar → *1. Your data* → Browse files. A green message
confirms what was read. A red message lists problems; fix them in Excel and
upload again. The software refuses to compute on bad input rather than
producing quiet nonsense.

**Step 2 — set the options.**

*2. Figure style*: colour palette (the first three are safe for colour-vision
deficiency, which several journals now expect), error bar type, font size,
export resolution.

*3. Options*: which prominence-value formula, and the soil basis your counts
are on.

**Step 3 — read Data check first.** Before any result. See Part 6.

**Step 4 — work through the tabs.**

**Step 5 — download.** Report tab for the PDF, or the Excel button at the foot
of the page for every table.

---

## Part 5 · The tabs

### Auto-assign
Proposes trophic group and c-p for each taxon from a 188-genus reference.
Green = exact match, accept. Amber = approximate match, check your spelling.
Red = no match, fill it in yourself. Rows flagged VERIFY belong to families not
in the cited source and are provisional.

Click **Apply proposals to blank cells only** — values you typed are never
overwritten by a lookup.

Confirmation is required by design. A fuzzy match on a misspelling, applied
silently, becomes a wrong maturity index that looks entirely normal.

### Data check
What your results can and cannot support. Read before anything else.

- **Pseudoreplication** — whether your rows are truly independent.
- **Degenerate indices** — indices that compute but mean nothing here. With no
  fungivores, CI evaluates to 0 and NCR to 1.0. Both are arithmetically correct
  and neither carries information; reporting "decomposition is entirely
  bacterial" from them would be an error.
- **Dominance** — samples where one taxon exceeds 70%. Diversity indices then
  largely report that taxon's share.
- **Trait provenance** — your assignments and their sources.

### Community
Norton (1978) descriptors per taxon. Warns when every taxon occurs in every
sample, because prominence value then collapses onto density and adds nothing.

### Diversity
Shannon, Simpson, Pielou, Margalef, richness. Bars carry Tukey letters; bars
sharing a letter do not differ.

### Faunal analysis
MI, MI(2–5), σMI, PPI, EI, SI, CI, BI, NCR, trophic composition, and the
faunal profile with its four quadrats.

### Soil health
The NSH index, plus a calibration panel showing how much of the 8–32 scale your
samples occupy and whether any subscore is constant. A constant subscore adds a
fixed amount to every sample and discriminates nothing.

### Footprints
Biomass and metabolic footprints, if measurements were supplied.

### Summary
Group means with SD, SE and 95% CI, in publication layout: `2.21 ± 0.03 (n=12)`.

### Statistics
ANOVA and Kruskal–Wallis with n, degrees of freedom, F, η², assumption
diagnostics, post-hoc tests, and your choice of Benjamini–Hochberg or Holm
correction.

### Multivariate
Bray–Curtis, PCA, PERMANOVA and PERMDISP.

### Soil nutrients
Correlation of indices against soil variables, redundancy analysis, forward
selection, soil PCA.

### Reference
188 genera with family, trophic group and c-p. Searchable, filterable,
downloadable.

### Validation
Thirty internal consistency checks, each recomputing from the raw counts.
Verdict: PASS, PASS WITH WARNINGS, or FAIL.

### Report
The full PDF.

### Cite
How to cite the software, and which paper to cite for each index.

---

## Part 6 · Reading the indices

| Index | Low means | High means |
|---|---|---|
| MI | disturbed, coloniser-dominated | mature, undisturbed |
| PPI | ectoparasites dominant | longidorids and trichodorids dominant |
| EI | few resources available | recently enriched |
| SI | simple food web | structured and connected |
| BI | little stress | stressed and depleted |
| CI | bacterial decomposition | fungal decomposition |

**Faunal profile quadrats** (Ferris et al. 2001), x = SI, y = EI, split at 50/50:

| | Low SI | High SI |
|---|---|---|
| **High EI** | **A Disturbed** — enriched, structurally simple, bacterial channel | **B Structured** — enriched and structured |
| **Low EI** | **D Degraded** — depleted and simple, stressed | **C Stable** — undisturbed, mature, fungal channel |

**NSH**: below 15 degraded, 15–24 moderate, above 25 well-functioning. Treat as
comparative, not absolute — the index is new and its authors state it needs
further calibration.

---

## Part 7 · Choosing a statistical test

Decide from the *type of variable*, before you see any p-value.

| Variable | Examples | Default model |
|---|---|---|
| Counts | total N, taxon abundance | negative binomial GLM |
| Bounded 0–100 | EI, SI, BI, CI, trophic % | beta regression or permutation test |
| Bounded continuous | MI, PPI, NCR | permutation test |
| Unbounded continuous | Shannon, Margalef | ANOVA if assumptions hold |
| Whole community | composition | PERMANOVA with PERMDISP |

ANOVA and Kruskal–Wallis are both displayed so that **disagreement between them
is visible**. They are diagnostic, not a menu. Choosing whichever gives the
smaller p after seeing the data inflates the false-positive rate and is not
defensible.

**If your design is paired or blocked** — the same field appears in every
treatment — one-way ANOVA is the wrong model. Use a paired test or fit the
block as a random effect. The Data check tab detects this automatically.

**Correct for multiple testing.** Twenty indices at α 0.05 gives roughly a 64%
chance of at least one false positive. Report the corrected column.

---

## Part 8 · Common problems

**"No sheet named 'counts'"** — your sheets are called Sheet1, Sheet2. Rename
them exactly `counts`, `taxa`, `samples`.

**"Taxa in counts with no row in taxa"** — a name differs between the two
sheets. Check for trailing spaces.

**"Duplicate sample names"** — two columns share a name. Common when
abbreviating.

**Footprints tab empty** — no `length_um` and `diameter_um`. This is a correct
result, not a fault.

**CI = 0 and NCR = 1.0 everywhere** — you have no fungivores. Do not report
these as findings.

**Statistics tab empty** — no `group` column in `samples`.

**Everything reads "moderate" under NSH** — your samples occupy a narrow part of
the scale. Check the calibration panel; this is a known limitation of the index.

**A summary differs from the table above it** — check whether the table is
showing every row. A partial table and a full-dataset mean are both correct and
will disagree.

---

## Part 9 · Before you publish

1. Verify every trophic and c-p assignment against a primary source.
2. Cross-check your faunal indices against NINJA (shiny.wur.nl/ninja) on one
   dataset, and report the comparison.
3. State which prominence-value formula you used — two variants circulate, both
   attributed to Norton (1978).
4. State your soil basis and taxonomic resolution.
5. Declare your experimental hierarchy: what was a replicate, what a subsample.
6. Report multiple-testing-corrected p-values and say which correction.
7. Do not report degenerate indices as values.
8. Cite the primary sources, not only this software.

---

## Part 10 · What this software cannot do for you

It cannot tell whether your identifications are right. It cannot tell whether
your sampling design answers your question. It cannot tell whether a c-p value
you entered is correct for your population. And it cannot turn a nematode
community into a management recommendation — no such inference is justified
from community data alone.

The arithmetic is verified. The biology is yours.

---

## References

Andrássy, I. (1956) *Acta Zoologica* 2: 1–15 ·
Anderson, M.J. (2001) *Austral Ecology* 26: 32–46 ·
Anderson, M.J. (2006) *Biometrics* 62: 245–253 ·
Benjamini, Y. & Hochberg, Y. (1995) *J R Stat Soc B* 57: 289–300 ·
Bongers, T. (1990) *Oecologia* 83: 14–19 ·
Bongers, T. & Bongers, M. (1998) *Appl Soil Ecol* 10: 239–251 ·
Ferris, H., Bongers, T. & de Goede, R.G.M. (2001) *Appl Soil Ecol* 18: 13–29 ·
Ferris, H. (2010) *Eur J Soil Biol* 46: 97–104 ·
Ghaderi, R. et al. (2025) *Eur J Soil Sci* 76: e70149 ·
Legendre, P. & Gallagher, E. (2001) *Oecologia* 129: 271–280 ·
Margalef, R. (1958) *Gen Syst* 3: 36–71 ·
Norton, D.C. (1978) *Ecology of Plant-Parasitic Nematodes*, Wiley ·
Pielou, E.C. (1966) *J Theor Biol* 13: 131–144 ·
Shannon, C.E. & Weaver, W. (1949) *The Mathematical Theory of Communication* ·
Simpson, E.H. (1949) *Nature* 163: 688 ·
Yeates, G.W. et al. (1993) *J Nematol* 25: 315–331
"""
