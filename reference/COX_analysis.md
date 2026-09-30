# Univariate or Multivariate Cox Regression

Fits Cox proportional-hazards models on a data.frame of survival plus
continuous and/or categorical covariates, and formats a tidy results
table (hazard ratios, 95% CIs, Wald p-values, sample sizes per
variable).

## Usage

``` r
COX_analysis(
  df,
  type = NULL,
  cont_Variates = NULL,
  cate_Variates = NULL,
  method = c("uni", "multi"),
  precision = 3,
  min_level_n = 5L,
  min_covariates = 2L,
  drop_nonestimable = NULL,
  auto_repair = FALSE
)
```

## Arguments

- df:

  Data.frame. If `type` is given it must contain `ID`, `<type>_time` and
  `<type>_status`; if `type = NULL` the columns must already be named
  `time` and `status`.

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token; the family is
  resolved to the token available in `df`. `NULL` when `df` already
  contains `time`/`status`. Default `NULL`.

- cont_Variates:

  Character vector of continuous covariate columns.

- cate_Variates:

  Character vector of categorical covariate columns (the first factor
  level is the reference level).

- method:

  Either `"uni"` (one model per covariate) or `"multi"` (single
  multivariable model with all covariates).

- precision:

  Decimal places used to format the display columns (`HR_text`,
  `P_text`); numeric columns keep full precision.

- min_level_n:

  Integer. Minimum number of patients required for a categorical level
  to keep its own coefficient in the multivariable model. Rarer levels
  are merged into an `"Other"` level when that leaves the covariate
  usable; if merging is not possible the patients belonging to those
  levels are excluded from the multivariable model and the covariate is
  kept. Default `5`.

- min_covariates:

  Integer. Number of covariates below which the multivariable model is
  no longer regarded as an adjusted model. When the automatic repair can
  only keep fewer than `min_covariates` covariates, the model is still
  fitted and returned (with a warning) so that the analysis can be
  inspected, but `metadata$reduced` is `TRUE` and
  `metadata$reduced_note` states which covariate(s) survived and why the
  others could not be co-estimated. Default `2`. Only a model with no
  estimable covariate at all is refused.

- drop_nonestimable:

  Deprecated. The old name of `auto_repair`; when supplied it sets
  `auto_repair` (and warns). Use `auto_repair` instead.

- auto_repair:

  Logical. `FALSE` (the default, and the fail-safe) is the new
  behaviour: a multivariable model that would need covariates dropped,
  levels merged or patients excluded is **not** fitted. The function
  returns an explicit `not estimable` result instead - class
  `cpas_not_estimable` (also carrying `cpas_COX`), with an empty
  `models` list, `status = "not estimable"`, the reason(s) in `$reasons`
  and the offending term(s) in `$offending_terms`, plus what the repair
  *would* have done in `$metadata$would_have_dropped`. `TRUE` restores
  the previous behaviour: constant or nearly constant covariates and
  perfectly separated or very rare levels are screened out first,
  covariates that still keep the joint model from converging are then
  removed one at a time (or, when the offender cannot be isolated, by
  trying each removal in turn, down to a single covariate), and every
  modification is recorded in `metadata$dropped_covariates`,
  `metadata$drop_notes` and the analysis manifest.

## Value

List of class `cpas_COX`:

- `input_params`::

  echoed inputs

- `processed_data`::

  analysis data after numeric/factor coercion

- `models`::

  list of fitted `coxph` models

- `summaries`::

  corresponding `summary.coxph` objects

- `results_table`::

  tidy results with columns `Var1`, `Variates`, `Level`, `Type` (`cont`,
  `cate_header`, `cate_level`), `N`, `HR`, `HR95L`, `HR95H`, `Pvalue`,
  `HR_text`, `P_text`

- `metadata`::

  sample sizes, events, and, for multivariable models,
  `dropped_covariates` and `drop_notes` describing every covariate that
  the automatic repair removed and why, plus `final_covariates`,
  `reduced` and `reduced_note` which flag a model that had to be reduced
  below `min_covariates`; `estimable` and `status` state whether a model
  was produced at all

- `manifest`::

  the analysis manifest
  ([`cpas_manifest`](https://wangjin93.github.io/CanPAS/reference/cpas_manifest.md)):
  covariates, rows and covariates dropped with reasons, the
  proportional-hazards check, the versions and the timestamp

- `reasons`, `offending_terms`, `estimable`, `status`::

  present on every result; for a fail-safe refusal
  (`auto_repair = FALSE`) they carry the reason(s) and the offending
  term(s) and `estimable` is `FALSE`, which is how a caller detects the
  refusal without a
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html)

## Details

Status is treated as `1 = event, 0 = censored`; rows with missing
time/status are dropped. For categorical covariates the first factor
level is the reference (no row is printed for it); header rows carry
`HR = NA` and are intended for plot annotations. Univariate models are
estimated on the complete cases of `(time, status, variable)`;
multivariate models on the complete cases of all covariates together.

The default is the fail-safe: `auto_repair = FALSE`. A multivariable
model that is not estimable exactly as requested comes back as an object
that says so (`estimable = FALSE`, `status = "not estimable"`, with the
reasons and the offending terms) rather than as a quietly reduced model.
This is a deliberate change of default; `auto_repair = TRUE` restores
the automatic repair described next and records every modification.

A multivariable model can fail not because the data are unusable but
because one covariate is: a constant or nearly constant continuous
covariate, a categorical level with very few patients, a level in which
all (or no) patients have an event, or a covariate collinear with the
others. With `drop_nonestimable = TRUE` such covariates are located and
removed, and the remaining model is fitted and reported together with
the reason for each removal; coefficients are never reported for a model
that did not converge, because those estimates are meaningless rather
than merely imprecise. The search continues down to a single covariate:
if no pair of covariates can be co-estimated, the largest estimable
model is returned with `metadata$reduced = TRUE`, a warning, and
`metadata$reduced_note` explaining that the result is not adjusted for
confounding. Only when nothing at all can be estimated does the function
stop, and it then reports both the covariates it had already removed and
the reason the final fit failed.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A real cohort from the CanPAS mirror: Lung Cancer, OS + DSS, n = 133.
   di <- dataset_info[dataset_info$Accession == "GSE14814", ]
   di[, c("Accession", "Type", "N", "EndpointFamilies")]

   ## Prerequisite chain: expression (probes collapsed) then survival, and
   ## with clin = TRUE the clinical covariates this cohort really carries.
   d <- cohort_merged("GSE13507", c("GAPDH", "ACTB"), type = "OS", clin = TRUE)
   colnames(d)

   ## One model per covariate.
   r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                      cate_Variates = c("grade", "N"), method = "uni")
   head(r$results_table)

   ## Joint model: covariates that cannot be co-estimated are located and
   ## removed, and every decision is recorded together with its reason.
   r2 <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                       cate_Variates = c("grade", "N"), method = "multi")
   r2$metadata$dropped_covariates
   r2$metadata$reduced
} # }
```
