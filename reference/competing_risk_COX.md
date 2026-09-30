# Cause-specific and subdistribution (Fine-Gray) Cox models

Fits, for the same covariates, (1) the cause-specific Cox model
(`Surv(time, status == etype)`) and (2) the Fine-Gray subdistribution
hazard model obtained with
[`survival::finegray()`](https://rdrr.io/pkg/survival/man/finegray.html)
and a weighted `coxph()`. The two answer different questions and are
reported side by side.

## Usage

``` r
competing_risk_COX(
  df,
  time,
  status,
  covariates = NULL,
  etype = 1,
  conf_level = 0.95,
  max_levels = 20L
)
```

## Arguments

- df:

  Data.frame with a time column and a cause-coded status column.

- time:

  Name of the numeric follow-up time column (years).

- status:

  Name of the status column: `0` = censored, `1, 2, ...` = event types.

- covariates:

  Character vector of covariates (the first one is usually the marker of
  interest). A factor stays categorical; a character column whose values
  all parse as numbers is treated as continuous (mirrored clinical
  columns arrive as text, e.g. `"55.4"`), any other character column
  becomes a factor. The classification actually used is returned in
  `covariate_types`.

- etype:

  Event type of interest (default 1). Competing events are all other
  non-zero codes.

- conf_level:

  Confidence level for the reported intervals.

- max_levels:

  Maximum number of levels accepted for a categorical covariate (default
  20). A covariate with one coefficient per level is usually a
  continuous variable that was supplied as text, so the call stops with
  an explanatory error instead of reporting meaningless hazard ratios;
  raise the limit for a deliberate fine stratification.

## Value

Object of class `cpas_competing`: list with

- `cause_specific`:

  tidy data.frame (`Variates`, `Level`, `HR`, `HR95L`, `HR95H`,
  `Pvalue`, `N`)

- `subdistribution`:

  the same columns for the Fine-Gray model (hazard ratios are
  subdistribution hazard ratios, sHR)

- `n`, `events_interest`, `events_competing`, `n_censored`, `etype`:

  sample and event accounting

- `covariate_types`:

  named vector: how each covariate was used (`"continuous"`,
  `"categorical"` or the text-coercion cases)

- `diagnostics`:

  per model, `ok` and, when `FALSE`, the reason (non-convergence,
  separation, infinite coefficient): the affected coefficient is then
  omitted from the table and a warning is raised, so a meaningless
  hazard ratio is never reported as a result

- `fits`:

  the two fitted `coxph` objects

## Details

The cause-specific hazard is the rate of the event among patients still
event-free; it answers "does the marker act on the disease process". The
subdistribution hazard keeps patients with a competing event in the risk
set and is directly linked to the cumulative incidence; it answers "does
the marker change the probability of dying of the disease". Reporting
only the cause-specific hazard can overstate the clinical effect when
competing events are common (and vice versa).

A numeric covariate contributes one coefficient. If a continuous
variable is supplied as text (as the mirrored clinical columns are), it
would otherwise be read as a factor with one level per value;
`competing_risk_COX()` restores the evident type, reports it in
`covariate_types`, and refuses (rather than silently reporting) a
categorical covariate with more than `max_levels` levels.
Non-convergence, separation and infinite coefficients are reported
through `diagnostics` and excluded from the tables.

The Fine-Gray standard errors are clustered by patient, because
[`survival::finegray()`](https://rdrr.io/pkg/survival/man/finegray.html)
expands every patient who has a competing event into several rows;
without the cluster term the reported standard error is 17-23% too
small. The point estimates and the clustered standard errors agree with
[`cmprsk::crr()`](https://rdrr.io/pkg/cmprsk/man/crr.html) to within
0.1%.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A real cohort that carries both OS and DSS, through the prerequisite
   ## chain: cause 1 = cancer death (DSS event), cause 2 = other death.
   d  <- cohort_merged("GSE14814", "GAPDH", type = "OS", clin = TRUE)
   cr <- data.frame(
     time   = d$OS_time,
     status = ifelse(d$OS_status == 0, 0, ifelse(d$DSS_status == 1, 1, 2)),
     age    = d$age,     # numeric: one coefficient, as for a Cox model
     sex    = factor(d$sex))
   table(cr$status)   # 0 censored, 1 cancer death, 2 other death

   cc <- competing_risk_COX(cr, time = "time", status = "status",
                            covariates = c("age", "sex"), etype = 1)
   cc$covariate_types    # how each covariate was used (continuous / categorical)
   cc$cause_specific     # cause-specific hazard ratios
   cc$subdistribution    # Fine-Gray subdistribution hazard ratios
   cc$diagnostics        # ok + reason per model; unusable coefficients are omitted
} # }
```
