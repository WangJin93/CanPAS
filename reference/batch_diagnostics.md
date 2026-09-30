# Per-cohort expression diagnostics for a multi-cohort analysis

Summarises the expression distribution of one gene (or a gene vector) in
each cohort, together with the cut-point values the pooled analyses use,
so a batch effect can be inspected before anything is pooled.

**Multi-cohort pooling uses within-cohort standardised effect sizes.**
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
is a two-stage meta-analysis: stage 1 fits the Cox model inside each
cohort, standardises the marker within that cohort and returns a log-HR
and its SE; stage 2 pools those per-cohort estimates by inverse
variance. A cohort-wide additive shift of the expression scale (the
usual platform batch effect) cancels inside the cohort and therefore
*does not enter the pooled estimate*. The pooled **Kaplan-Meier** of
[`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md)
splits each cohort at a **within-cohort percentile** (the cohort's own
median, or its own \\(100-p)\\th percentile), so again no absolute
cross-platform threshold is used. The summaries this function returns
make both facts checkable: `median` and the two cut-point columns are
the values used, and the standardised columns are what actually enters
the pool.

## Usage

``` r
batch_diagnostics(
  datasets,
  genes,
  type = "OS",
  top_pct = 25,
  merged = NULL,
  process_duplicates = "max",
  max_try = 3
)
```

## Arguments

- datasets:

  Character vector of dataset accessions or TCGA projects.

- genes:

  One or more gene symbols.

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or a concrete token. It only
  selects which survival columns are joined and which token is reported;
  the expression summaries are independent of it.

- top_pct:

  Percentile used by the `cut = "top_pct"` rule of
  [`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md):
  the cut-point reported is the cohort's \\(100 - top\\pct)\\th
  percentile, i.e. `top_pct = 25` means the top 25% of the cohort is
  "High" (default 25).

- merged:

  Optional named list (cohort -\> merged data.frame) that skips all
  retrieval; useful offline and in tests.

- process_duplicates:

  Probe collapsing rule passed to
  [`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)
  (default `"max"`).

- max_try:

  Retrieval attempts per cohort.

## Value

A list of class `cpas_batch_diagnostics`:

- `table`::

  one row per cohort and gene: `dataset`, `gene`, `endpoint` (the token
  resolved for that cohort), `n` (non-missing expression values),
  `n_missing`, `median`, `q1`, `q3`, `iqr`, `min`, `max`, `mean`, `sd`,
  `cut_median` (the median-split cut-point used), `cut_top_pct` (the
  percentile cut-point used), `n_high_median`, `n_low_median`,
  `n_high_top_pct`, `n_low_top_pct`, and `median_z` (the cohort median
  on the within-cohort standardised scale, i.e. what the pool compares)

- `cohorts`::

  one row per cohort: how many genes were found, the token resolved and
  the sample count

- `values`::

  the long per-sample values the plot and the summaries are built from:
  `dataset`, `gene`, `value`

- `batch`::

  the cross-cohort spread: per gene, `min_median`, `max_median`,
  `median_shift` (max minus min, the absolute cross-platform difference
  in expression units), and `median_z_shift` (the same shift after
  within-cohort standardisation)

- `errors`, `input`, `manifest`::

  cohorts that could not be read, the echoed inputs, and the analysis
  manifest

## Details

The cut-point columns are computed exactly as
[`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md)
computes them, on the same rows: the median-split value is
`stats::median(value)` and the percentile value is
`stats::quantile(value, 1 - top_pct/100)` with the default type-7
interpolation. They are reported per cohort, never pooled, because
pooling an absolute expression threshold across platforms is exactly the
mistake the within-cohort rules avoid.

Cohorts are read through
[`cohort_merged`](https://wangjin93.github.io/CanPAS/reference/cohort_merged.md),
which goes through the same cached accessors as the analyses
([`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md),
[`merge_surv_expr`](https://wangjin93.github.io/CanPAS/reference/merge_surv_expr.md),
[`get_data`](https://wangjin93.github.io/CanPAS/reference/get_data.md),
[`tcga_surv_table`](https://wangjin93.github.io/CanPAS/reference/tcga_surv_table.md)).
A cohort whose tables are already in the local cache is therefore
analysed without a network call; passing `merged` skips retrieval
entirely.

## See also

[`plot_batch_diagnostics`](https://wangjin93.github.io/CanPAS/reference/plot_batch_diagnostics.md),
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md),
[`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md),
[`cohort_merged`](https://wangjin93.github.io/CanPAS/reference/cohort_merged.md)

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Two real lung cohorts; the second run reads the cache offline.
   b <- batch_diagnostics(c("GSE31210", "GSE37745"), genes = "GAPDH",
                          type = "DFS", top_pct = 25)
   b$table[, c("dataset", "gene", "n", "median", "cut_median", "cut_top_pct")]
   b$batch                      # how far apart the absolute scales are
   plot_batch_diagnostics(b)    # the distributions side by side
} # }
```
