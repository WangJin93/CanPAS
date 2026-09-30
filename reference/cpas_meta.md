# Two-stage cross-cohort integrative (meta) analysis

Runs a two-stage integrative analysis over several cohorts of one cancer
type: stage 1 fits `coxph(Surv(time, status) ~ marker + confounders)` in
every cohort, giving log-HR +/- SE per cohort; stage 2 pools them by
inverse variance and reports the pooled HR (95% CI), Z/p, the
heterogeneity statistics Q/df/I\\^2\\/tau\\^2\\ and two prediction
intervals.

The pooling method is chosen by `method`: `"REML"` (default, iterative
restricted maximum likelihood for tau\\^2\\), `"DL"` (DerSimonian-Laird,
the previous default, kept for continuity), `"HK"`
(Hartung-Knapp-Sidik-Jonkman adjusted variance of the pooled estimate,
with its characteristic wider interval) or `"FE"` (fixed effect). `"RE"`
is still accepted as a synonym of `"DL"`.

## Usage

``` r
cpas_meta(
  datasets,
  marker,
  type = "OS",
  method = c("REML", "DL", "HK", "FE"),
  pi_method = c("t", "normal", "HK"),
  confounders = NULL,
  min_events = 5,
  max_try = 3,
  merged = NULL,
  pooling = c("family", "exact"),
  overlap = c("warn", "refuse", "dedupe"),
  auto_repair = FALSE,
  class_table = NULL,
  overlap_table = NULL
)
```

## Arguments

- datasets:

  Character vector of dataset accessions (matching
  `dataset_info$Accession`).

- marker:

  One gene symbol (e.g. `"TP53"`) or a signature formula string (e.g.
  `"0.5*TP53+0.3*GAPDH"`).

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or a raw endpoint token
  (OS/DSS/DFS/RFS/PFS/MFS/DFI/PFI/DRFS/EFS, ...). A family is resolved
  per cohort to the concrete token that cohort provides (see
  [`endpoint_resolve`](https://wangjin93.github.io/CanPAS/reference/endpoint_resolve.md))
  and recorded in `per_dataset$endpoint` and in
  `per_dataset$pooling_class`.

- method:

  Pooling method: `"REML"` (default), `"DL"`, `"HK"` or `"FE"`; `"RE"`
  is accepted as `"DL"`.

- pi_method:

  Construction rule of the 95% prediction interval for a new cohort
  (spec B6): `"t"` (default) is the documented \\\hat\mu \pm
  t\_{0.975,k-2}\sqrt{se\_{pooled}^2+\tau^2}\\ interval (needs at least
  3 cohorts); `"normal"` is the normal approximation \\\hat\mu \pm
  1.96\sqrt{se\_{pooled}^2+\tau^2}\\ (needs at least 2); `"HK"` uses the
  Hartung-Knapp-Sidik-Jonkman adjusted standard error of the pooled
  estimate with \\t\_{0.975,k-1}\\ (needs at least 2). The rule actually
  used is recorded in `$pooled$pi_method` and `$pooled$pi_rule`, and the
  other construction is returned beside it, clearly labelled, in
  `$pooled$pi_alt_lower`/`pi_alt_upper` with `$pooled$pi_alt_rule`. Only
  the prediction interval is affected: `pi_method` does not change the
  pooled HR, its confidence interval, or any heterogeneity statistic.

- confounders:

  Optional character vector of covariate names to adjust for (they must
  exist in the merged data; the covariates available in each cohort are
  intersected). With the default `auto_repair = FALSE` a cohort in which
  a requested covariate is constant is reported as `not estimable`
  instead of being adjusted for a reduced covariate set.

- min_events:

  Minimum number of events for a cohort to enter (default 5).

- max_try:

  Retries per cohort when fetching fails (API jitter).

- merged:

  Optional named list (cohort name -\> already merged data.frame). When
  supplied nothing is fetched from the network.

- pooling:

  `"family"` (default) pools the Exact-equivalent and the
  Clinically-related rows of the endpoint family, i.e. the cohorts that
  contribute RFS/EFS/DFI to DFS, CSS/BCSS to DSS, PFI to PFS and DRFS to
  MFS as well as those whose token is the family itself. `"exact"` pools
  only rows whose pooling class is `Exact-equivalent`. Both modes are
  available because `"family"` is what makes the cross-token claim (the
  DFS pool would fall from 98 cohorts to the ~20 whose token is
  literally DFS); the mode actually used is recorded in
  `$pooled$pooling`, in `$per_dataset$pooling_class` and in `$manifest`.

- overlap:

  What to do about cohort pairs known to share patients: `"warn"`
  (default) proceeds but warns and records the pairs in
  `$pooled$overlap_pairs`, `"refuse"` stops with an actionable error,
  `"dedupe"` keeps one member of each overlapping group (the larger
  cohort, ties by accession sort) and records what it dropped in
  `$excluded` and `$manifest$overlap_dropped`.

- auto_repair:

  `FALSE` (default) is the fail-safe: a cohort whose model would need
  covariates dropped, levels merged or rows excluded is reported as
  `not estimable` with its reasons and offending terms in `$excluded`,
  and no model is returned for it. `TRUE` restores the previous
  behaviour (repair the model) and records every modification in
  `$modifications`.

- class_table:

  Optional endpoint-semantics table (the frozen schema of
  [`endpoint_semantics`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md))
  used instead of the shipped copy.

- overlap_table:

  Optional shared-patient register (the frozen schema of
  [`cohort_overlap`](https://wangjin93.github.io/CanPAS/reference/cohort_overlap.md))
  used instead of the shipped copy.

## Value

A list of class `cpas_meta`:

- `per_dataset`::

  one row per pooled cohort with `dataset`, `endpoint` (the token used),
  `pooling_class`, `n`, `events`, `HR`, `lower`, `upper`, `logHR`, `se`,
  `p`

- `pooled`::

  one row: the method used, `pi_method` (the prediction-interval rule
  used), `k`, `HR`, `lower`, `upper`, `p`, `Q`, `df`, `p_heterogeneity`,
  `I2`, `tau2`, the primary prediction interval `pi_lower`/`pi_upper`
  (by default t with k-2 df) and the alternative interval
  `pi_alt_lower`/`pi_alt_upper` (by default the normal approximation),
  plus `pooling`, `pooling_classes`, `overlap_pairs` and `overlap_mode`

- `excluded`::

  cohorts that did not enter and why (`cohort`, `stage`, `reason`);
  `stage = "model"` rows are the fail-safe refusals

- `modifications`::

  every modification `auto_repair = TRUE` made (`cohort`, `variable`,
  `detail`, `reason`)

- `errors`, `input`, `manifest`::

  fetch/endpoint errors per cohort, the echoed inputs, and the analysis
  manifest
  ([`cpas_manifest`](https://wangjin93.github.io/CanPAS/reference/cpas_manifest.md))

## See also

[`endpoint_semantics`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md),
[`cohort_overlap`](https://wangjin93.github.io/CanPAS/reference/cohort_overlap.md),
[`cpas_manifest`](https://wangjin93.github.io/CanPAS/reference/cpas_manifest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Real cohorts, real endpoint: each cohort contributes the DFS-family
   ## token it actually has, which the function resolves and reports.
   m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
   m$per_dataset[, c("dataset", "endpoint", "pooling_class", "n", "events", "HR")]
   m$pooled
   plot_meta_forest(m)
   loo_meta(m)               # leave-one-out sensitivity
   cpas_manifest(m)          # what was pooled, and under which defaults
   print(m)

   ## Strict mode: only cohorts whose token IS the family definition.
   m2 <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS",
                   pooling = "exact")
   m2$per_dataset$pooling_class

   ## Prediction-interval rule (spec B6): the default t(k-2) interval, the
   ## normal approximation, or the HK-based interval. The interval used is
   ## recorded in $pooled$pi_method.
   m_normal <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53",
                         type = "DFS", pi_method = "normal")
   m_normal$pooled[, c("pi_method", "pi_lower", "pi_upper", "pi_alt_lower")]
} # }
```
