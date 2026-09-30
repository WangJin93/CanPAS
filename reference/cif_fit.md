# Cumulative incidence functions for competing risks

Estimates the Aalen-Johansen cumulative incidence function (CIF) of each
event type with
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
on a multi-state response `Surv(time, factor(status))`, optionally by
group.

## Usage

``` r
cif_fit(
  df,
  time,
  status,
  group = NULL,
  times = NULL,
  conf_level = 0.95,
  group_cut = c("none", "median", "tertile", "quartile"),
  max_levels = 20L
)
```

## Arguments

- df:

  Data.frame with a time column and a cause-of-event status column.

- time:

  Name of the numeric follow-up time column (years).

- status:

  Name of the status column: `0` = censored, `1, 2, ...` = event types
  (cause of the event).

- group:

  Optional name of a **categorical** column used to split the CIF (e.g.
  a marker group or sex). A continuous column is refused: it would
  become one stratum per distinct value.

- times:

  Optional numeric vector of time points at which the CIF is tabulated
  (default: deciles of the observed event times).

- conf_level:

  Confidence level for the pointwise interval.

- group_cut:

  What to do when `group` names a continuous column: `"none"` (default)
  stops with an explanatory error, `"median"` splits at the median into
  low/high, `"tertile"` into low/middle/high and `"quartile"` into
  Q1-Q4. The cut points are reported in the result (`cut_points`), so
  the split is reproducible.

- max_levels:

  More distinct values than this make a group column count as continuous
  (default 20).

## Value

Object of class `cpas_cif`: list with `table` (data.frame: `group`,
`cause`, `time`, `n_risk`, `cif`, `lower`, `upper`), `fit` (the
`survfit` object), `causes`, `events`, `n`, `group_levels` and
`group_sizes` (how the strata were formed) and the echoed inputs.

## Details

With competing events, the complement of the Kaplan-Meier estimate is
**not** the incidence of the event of interest: patients who die of
another cause are removed from the risk set in the Kaplan-Meier
estimator and their future event probability is redistributed over the
remaining causes. The CIF does not redistribute that probability, which
is why it is the appropriate descriptive curve whenever competing events
exist.

`group` names a grouping variable, not a covariate: the column is
factored as supplied, so a continuous variable such as age would produce
one stratum per distinct value (and one curve per stratum). Such a
column is refused unless `group_cut` asks for a documented quantile
split, in which case the threshold(s) used are returned with the result.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A real cohort that carries both OS and DSS, through the prerequisite
   ## chain: cause 1 = cancer death (DSS event), cause 2 = other death.
   d  <- cohort_merged("GSE14814", "GAPDH", type = "OS", clin = TRUE)
   cr <- data.frame(
     time   = d$OS_time,
     status = ifelse(d$OS_status == 0, 0, ifelse(d$DSS_status == 1, 1, 2)),
     age    = d$age,                        # continuous
     sex    = factor(d$sex),                # categorical
     age_group = factor(ifelse(d$age > stats::median(d$age, na.rm = TRUE),
                                 "older", "younger")))
   table(cr$status)   # 0 censored, 1 cancer death, 2 other death

   ci <- cif_fit(cr, time = "time", status = "status", times = c(1, 3, 5))
   ci$table

   ## Stratified by a real categorical column (sex) ...
   ci_sex <- cif_fit(cr, time = "time", status = "status", group = "sex")
   ci_sex$group_levels; ci_sex$group_sizes
   plot_cif(ci_sex)

   ## ... or by a documented split of the CONTINUOUS age column. Passing age
   ## directly as group is refused: it would create one stratum per distinct
   ## age. group_cut says where to cut, and returns the cut point used.
   ci_age <- cif_fit(cr, time = "time", status = "status",
                     group = "age", group_cut = "median")
   ci_age$cut_points; ci_age$group_levels
} # }
```
