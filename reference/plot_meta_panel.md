# Plot a per-gene meta-analysis panel

Forest-style plot of the pooled hazard ratios of a
[`cpas_meta_panel`](https://wangjin93.github.io/CanPAS/reference/cpas_meta_panel.md)
result: one row per gene, point size scaled by the number of pooled
datasets, colour by FDR significance, with the heterogeneity (I2) shown
next to each estimate.

## Usage

``` r
plot_meta_panel(x, annotate = TRUE, digits = 4)
```

## Arguments

- x:

  An object returned by
  [`cpas_meta_panel`](https://wangjin93.github.io/CanPAS/reference/cpas_meta_panel.md).

- annotate:

  Add the per-gene I2 (and the FDR-adjusted p for the significant genes)
  to the axis labels.

- digits:

  Decimal places for the printed estimates (default 4).

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
   p <- cpas_meta_panel(c("GSE31210", "GSE37745"), genes = c("TP53", "GAPDH"),
                         type = "DFS")
   plot_meta_panel(p, annotate = TRUE, digits = 4)
} # }
```
