# Sensitivity analysis: leave-one-out

Refits the pooled estimate once per cohort, each time leaving one cohort
out, so the influence of a single cohort on the overall result can be
assessed.

## Usage

``` r
loo_meta(x)
```

## Arguments

- x:

  An object of class `cpas_meta`.

## Value

data.frame: `left_out` names the cohort that was removed and the
remaining columns are the re-pooled summary for that omission, computed
with the same pooling method (`x$input$method`) as the original
analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
   m <- cpas_meta(c("GSE31210", "GSE37745", "GSE42127"), marker = "TP53", type = "DFS")
   loo_meta(m)
} # }
```
