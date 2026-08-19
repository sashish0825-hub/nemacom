"""
nemonema_summary.py — replicate handling: group means, dispersion, and
publication-ready summary tables.

Two distinct operations, often confused:

  aggregate_subsamples()   averages TECHNICAL subsamples (several aliquots or
                           cores from the same plot) BEFORE any analysis. This
                           is a correction for pseudoreplication: those rows
                           were never independent, so they should never have
                           entered the analysis as separate observations.

  group_summary()          summarises BIOLOGICAL replicates (independent plots
                           within a treatment) AFTER analysis, into mean and a
                           dispersion measure.

Averaging biological replicates would throw away the very variation the
statistics need. Averaging technical subsamples is mandatory. Getting these the
wrong way round is the commonest design error in nematode community papers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
def aggregate_subsamples(counts: pd.DataFrame, samples: pd.DataFrame,
                         by: str, how: str = "mean"):
    """Collapse technical subsamples to one observation per unit of `by`.

    Use ONLY when several rows are aliquots or cores from the same plot. If the
    rows are independent plots, do not use this - you would be discarding real
    replication.
    """
    if by not in samples.columns:
        raise ValueError(f"column '{by}' not in the samples sheet")
    key = samples[by].astype(str)
    grouped = counts.T.groupby(key.reindex(counts.columns).values)
    new_counts = (grouped.mean() if how == "mean" else grouped.sum()).T
    new_counts = new_counts.round(2) if how == "mean" else new_counts

    meta = samples.copy()
    meta[by] = key.values
    agg = {c: "first" for c in meta.columns if c != by}
    new_samples = meta.groupby(by).agg(agg)
    new_samples.index.name = "Sample"

    n_before, n_after = counts.shape[1], new_counts.shape[1]
    new_counts.attrs["note"] = (
        f"{n_before} rows collapsed to {n_after} by '{by}' ({how}). "
        "Only correct if the collapsed rows were technical subsamples.")
    return new_counts, new_samples


# --------------------------------------------------------------------------
def group_summary(table: pd.DataFrame, groups: pd.Series,
                  err: str = "SE", decimals: int = 2) -> pd.DataFrame:
    """Mean, SD, SE, 95% CI and n per group for every numeric column.

    The 95% CI uses the t distribution, not 1.96 — at the n of 3 to 5 typical
    of field trials the normal approximation is noticeably too narrow.
    """
    from scipy import stats as st
    g = groups.reindex(table.index)
    num = table.select_dtypes(include=[np.number])
    rows = []
    for col in num.columns:
        for lev in pd.unique(g.dropna()):
            v = num.loc[g == lev, col].dropna().values
            if len(v) == 0:
                continue
            n = len(v)
            sd = v.std(ddof=1) if n > 1 else np.nan
            se = sd / np.sqrt(n) if n > 1 else np.nan
            tcrit = st.t.ppf(0.975, n - 1) if n > 1 else np.nan
            rows.append({
                "index": col, "group": str(lev), "n": n,
                "mean": round(float(v.mean()), decimals),
                "SD": round(float(sd), decimals) if sd == sd else np.nan,
                "SE": round(float(se), decimals) if se == se else np.nan,
                "CI95_lower": round(float(v.mean() - tcrit * se), decimals)
                              if se == se else np.nan,
                "CI95_upper": round(float(v.mean() + tcrit * se), decimals)
                              if se == se else np.nan,
                "min": round(float(v.min()), decimals),
                "max": round(float(v.max()), decimals),
                "CV_pct": round(float(sd / v.mean() * 100), 1)
                          if sd == sd and v.mean() else np.nan,
            })
    out = pd.DataFrame(rows)
    out.attrs["err"] = err
    return out


def summary_wide(table: pd.DataFrame, groups: pd.Series, err: str = "SE",
                 decimals: int = 2) -> pd.DataFrame:
    """Indices as rows, groups as columns, each cell 'mean ± err (n)'.

    This is the layout journals expect for a results table, and it can be
    pasted straight into a manuscript.
    """
    s = group_summary(table, groups, err, decimals)
    if s.empty:
        return pd.DataFrame()
    s["cell"] = s.apply(
        lambda r: (f"{r['mean']:.{decimals}f} ± {r[err]:.{decimals}f} (n={r['n']})"
                   if r[err] == r[err] else f"{r['mean']:.{decimals}f} (n={r['n']})"),
        axis=1)
    wide = s.pivot(index="index", columns="group", values="cell")
    wide.columns.name = None
    wide.index.name = "Index"
    wide.attrs["caption"] = (
        f"Values are group means ± {err}; n is the number of independent "
        "replicates. Verify that replicates are biologically independent — "
        "technical subsamples must be averaged before this table is built.")
    return wide


# --------------------------------------------------------------------------
def replicate_structure(samples: pd.DataFrame, group_col: str = "group") -> dict:
    """Describe the replication actually present, and whether it supports
    inference."""
    cols = {c.lower(): c for c in samples.columns}
    gc = cols.get(group_col)
    if gc is None:
        return {"ok": False, "message": "No `group` column in the samples sheet."}
    counts = samples[gc].value_counts().sort_index()
    smallest = int(counts.min())
    out = {"replicates_per_group": counts.to_dict(),
           "n_groups": int(len(counts)), "smallest": smallest,
           "balanced": bool(counts.nunique() == 1)}
    if smallest < 2:
        out["ok"] = False
        out["message"] = (f"A group has n={smallest}. No dispersion can be "
                          "estimated and no test can run. Add replicates, or "
                          "redefine `group` at a level that has them.")
    elif smallest < 3:
        out["ok"] = True
        out["message"] = (f"Smallest group n={smallest}. Tests will run but are "
                          "very weakly powered; three is the usual minimum for "
                          "any inferential claim.")
    else:
        out["ok"] = True
        out["message"] = (f"n={smallest} or more per group"
                          + ("; balanced design." if out["balanced"]
                             else "; unbalanced design — state this in methods."))
    return out
