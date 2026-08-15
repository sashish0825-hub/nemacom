MANUAL = """
### 1. Prepare three sheets

**counts** - taxa in rows, samples in columns, individuals recovered.
One consistent soil basis for every sample; per 100 g dry soil preferred.
Blank or 0 means absent. Numbers only.

**taxa** - one row per taxon, names matching counts exactly.
`trophic` = PP / BF / FF / OM / PR (Yeates et al. 1993).
`cp` = 1 to 5 (Bongers 1990; Bongers & Bongers 1998).
`source` = where you got each assignment. Every food-web index is built from
trophic and cp, so a wrong assignment gives a wrong index no software can catch.

**samples** - one row per sample. `group` = the factor you compare.
Add `field` and `plot` so replicates can be told from subsamples.
At least 3 independent replicates per group.

### 2. Optional - biomass

Fill `length_um` and `diameter_um` from your own measurements (best), or
`length_um` and `a_ratio` from a description, or `family_weight_ug` from
Ferris (2010) Table 1 (weakest - family means vary hugely and come from adults).

### 3. Reading the indices

| Index | Low | High |
|---|---|---|
| MI | disturbed | mature |
| PPI | ectoparasites | longidorids, trichodorids |
| EI | few resources | recently enriched |
| SI | simple food web | structured |
| BI | little stress | stressed |
| CI | bacterial channel | fungal channel |

Faunal profile (Ferris et al. 2001): high EI + high SI = structured;
high EI + low SI = disturbed; low EI + high SI = stable; both low = degraded.

### 4. Common problems

Missing taxa - a name differs between sheets; check trailing spaces.
CI = 0 and NCR = 1.0 everywhere - no fungivores; these are undefined, not zero.
Empty Statistics tab - no `group` column.
Everything "moderate" under NSH - narrow scale span, a known index limitation.

### 5. Before publishing

Verify every trophic and cp assignment. Cross-check against NINJA. State your
prominence formula, soil basis and taxonomic resolution. Report corrected
p-values. Declare what was a replicate and what a subsample.
"""
