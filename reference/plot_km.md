# Kaplan-Meier Survival Plot by Marker

Splits samples into `High` / `Low` groups around a marker cut-point
(default: median) and draws a Kaplan-Meier curve with
[`survminer::ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html).

## Usage

``` r
plot_km(df, type = "OS", marker, cutpoint = "median", ...)
```

## Arguments

- df:

  Data.frame with columns `ID`, `<type>_time`, `<type>_status` and
  `marker`.

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token (default
  `"OS"`).

- marker:

  Name of the numeric marker column used for stratification.

- cutpoint:

  Cut-point rule: a number (fixed threshold), `"median"` (default) or
  `"mean"`.

- ...:

  Further arguments passed to
  [`survminer::ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html)
  (e.g. `pval = TRUE`, `risk.table = TRUE`, `legend = "bottom"`).

## Value

A `ggsurvplot` object (access the ggplot via `$plot`).

## Details

Rows with missing time/status/marker are removed; groups are
`High = marker > cutpoint` and `Low = marker <= cutpoint`.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## The prerequisite chain returns one analysis-ready table per cohort.
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")

   plot_km(d, type = "OS", marker = "GAPDH", pval = TRUE)

   ## A signature or a specific cutpoint works on the same table.
   plot_km(d, type = "OS", marker = "GAPDH", cutpoint = "median")
} # }
```
