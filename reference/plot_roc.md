# Time-dependent ROC Curve for a Marker

Computes a time-dependent ROC curve at a given prediction horizon with
[`survivalROC::survivalROC`](https://rdrr.io/pkg/survivalROC/man/survivalROC.html)
and draws it with ggplot2.

## Usage

``` r
plot_roc(
  df,
  type = "OS",
  marker,
  predict.time = 1,
  method = c("KM", "NNE"),
  span = NULL
)
```

## Arguments

- df:

  Data.frame with columns `<type>_time`, `<type>_status` and `marker`.

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token (default
  `"OS"`).

- marker:

  Name of the numeric marker column.

- predict.time:

  Prediction horizon (same unit as `<type>_time`, typically years).
  Default 1.

- method:

  ROC method passed to `survivalROC`: `"KM"` (default) or `"NNE"`.

- span:

  Span for `method = "NNE"` (ignored for KM).

## Value

A ggplot object; the estimated AUC is attached as attribute `"AUC"` of
the returned object.

## Details

Rows with missing time/status/marker are removed first. The
`<type>_time` values are interpreted as the same unit as `predict.time`.

## Examples

``` r
if (FALSE) { # \dontrun{
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")

   g <- plot_roc(d, type = "OS", marker = "GAPDH", predict.time = 3)
   g$auc
} # }
```
