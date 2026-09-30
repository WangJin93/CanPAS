# Pooled Kaplan-Meier (marker split inside every cohort)

Splits each cohort into High/Low at the marker and then pools by
`method`: `"ipd"` pools the patients directly into one `survfit` (the
log-rank test is reported both unstratified and stratified by cohort);
`"meta"` two-stage: at every time point the cohorts' KM survival
probabilities S(t) are transformed with log(-log S) and pooled by
inverse variance (RE/FE), giving a pooled survival curve and the
1/3/5-year table; `"both"` returns both routes.

## Usage

``` r
cpas_km_pooled(
  merged,
  marker,
  type = "OS",
  method = c("both", "ipd", "meta"),
  landmarks = c(1, 3, 5),
  meta_method = "RE",
  cut = c("median", "top_pct", "custom"),
  top_pct = NULL,
  cut_value = NULL
)
```

## Arguments

- merged:

  Named list of the cohorts' `merge_surv_expr` merged data.frames. The
  names must be catalog accessions so the concrete endpoint token of
  each cohort can be resolved. The result carries `dataset_endpoints`,
  the token every dataset actually used.

- marker:

  Gene column name.

- type:

  Endpoint family or token (OS/RFS/...).

- method:

  `"ipd"`, `"meta"` or `"both"`.

- landmarks:

  Time points, in years, at which the pooled survival probabilities are
  tabulated and drawn (default 1, 3 and 5 years).

- meta_method:

  Pooling method for the time-point survival probabilities:
  `"RE"`/`"DL"` (DerSimonian-Laird random effects, the default kept for
  continuity with the published curves), `"REML"`, `"HK"` or `"FE"`
  (fixed effect).

- cut:

  High/low split rule, one of three: `"median"` (default) takes the top
  50% inside every cohort, i.e. High = marker above **that cohort's
  own** median; `"top_pct"` sorts each cohort by expression and takes
  the highest `top_pct`% as High (threshold = that cohort's \\100 -
  top_pct\\ percentile, so `top_pct = 25` means the top 25% are high
  expressors); `"custom"` uses the absolute threshold `cut_value`, with
  High above it and Low at or below it. An absolute threshold is
  meaningful here because every matrix served by the mirror is on the
  log2 scale (Section 2.3 of the accompanying paper), so one value means
  the same thing in every cohort; percentile splits are nevertheless the
  conventional reading of "high versus low expression", and the
  application offers only median and top_pct. A per-cohort search for
  the best cut point is deliberately **not** offered: it would repeat,
  cohort by cohort, the cut-point search quantified in Section 3.4 and
  compound its inflation when pooled.

- top_pct:

  Percentage for `cut = "top_pct"`: a single finite value between 1 and
  99; the patients with the highest `top_pct`% of expression are High.
  The threshold is computed with the default
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html)
  interpolation (type 7); when several patients sit exactly on the
  threshold the High group can be slightly larger than `top_pct`%, and
  the actual threshold and counts are returned (`cohort_thresholds`,
  `cutpoint`). A cohort left one-sided by the rule is recorded in
  `empty_cohorts` with its reason and does not enter the pool.

- cut_value:

  Threshold for `cut = "custom"`: a single finite value on the log2
  expression scale (or on the signature score). A cohort with no patient
  on one side of the threshold is recorded in `empty_cohorts` with its
  reason and does not enter the pool.

## Value

An ordinary `list` (no class attribute; not an S3 object, so there is no
method dispatch):

- `df`::

  the pooled analysis data (time/status/marker/dataset/group)

- `datasets`, `n_high`, `n_low`, `method`, `cutpoint`::

  the cohorts included and the size of each group

- `dataset_endpoints`::

  the endpoint token every dataset actually used (named vector)

- `cut`, `top_pct`, `cut_value`, `cutpoint`::

  the split rule actually used, its parameter or threshold, and a
  printable description of the rule; `cohort_thresholds` gives the
  per-cohort threshold when `cut = "top_pct"`

- `n_dropped`, `empty_cohorts`, `empty_reasons`::

  the split rule is exhaustive for the cohorts that entered (`n_dropped`
  is always 0); when a threshold or percentile rule leaves one side
  empty in a cohort, that cohort is recorded in `empty_cohorts` with its
  reason and does not enter the pool

- `skipped_cohorts`, `skipped_reasons`::

  cohorts excluded for insufficient data and why - the endpoint could
  not be resolved, a column is missing, fewer than 10 complete (time,
  status, marker) rows, or the marker takes a single value in that
  cohort. These cohorts do not enter the pool but are **never silently
  dropped**

- `fit`, `logrank_p`, `logrank_p_stratified`::

  the `survfit` and the log-rank p (pooled / stratified by cohort) when
  `method` is ipd or both

- `meta_landmarks`, `meta_curve`::

  the time-point pooled survival table and the fine-grid curve when
  `method` is meta or both

- `manifest`::

  the analysis manifest
  ([`cpas_manifest`](https://wangjin93.github.io/CanPAS/reference/cpas_manifest.md)):
  the cut-point rule, the number of cut-points searched (0: no search),
  the cohorts skipped and why, the versions and the timestamp

## Details

Each cohort is split on its own marker median (High = marker \> that
cohort's median), so what is compared is the risk above versus below
each cohort's own median and not an absolute threshold across cohorts;
when an absolute threshold is wanted, use
[`plot_km()`](https://wangjin93.github.io/CanPAS/reference/plot_km.md)
cohort by cohort. The IPD route reports the log-rank p both unstratified
and stratified by cohort; reporting the stratified one is recommended.
In the time-point route, a time point beyond a cohort's longest
follow-up reuses that cohort's last observed S(t) (`extend = TRUE`), and
the number of cohorts actually contributing at each time point is
recorded in `meta_landmarks$k`.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Pooled Kaplan-Meier from an explicitly built list of merged cohorts
   ## (names must be catalog accessions so the endpoint is resolved per cohort).
   merged <- list(GSE31210 = cohort_merged("GSE31210", "GAPDH", type = "RFS"),
                   GSE37745 = cohort_merged("GSE37745", "GAPDH", type = "RFS"))
   km <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS", method = "both")
   km$datasets; km$dataset_endpoints; km$logrank_p

   plot_cpas_km(km)                       # pooled curve + landmark table
   plot_cpas_km_perdataset(km, ncol = 2)  # one panel per cohort
} # }
```
