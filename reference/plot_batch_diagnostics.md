# Plot per-cohort expression distributions

Draws the per-cohort expression distribution of every gene in a
[`batch_diagnostics`](https://wangjin93.github.io/CanPAS/reference/batch_diagnostics.md)
result as side-by-side box plots on one shared axis, with the
within-cohort median-split and percentile cut-points marked, so the
cross-platform shift and the cut-points that absorb it are visible in
the same panel.

## Usage

``` r
plot_batch_diagnostics(
  x,
  genes = NULL,
  show_cut_points = TRUE,
  free_y = FALSE,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`batch_diagnostics`](https://wangjin93.github.io/CanPAS/reference/batch_diagnostics.md).

- genes:

  Optional subset of genes to draw (default: all).

- show_cut_points:

  Draw the per-cohort median (solid segment) and percentile (dashed
  segment) cut-points (default `TRUE`).

- free_y:

  Let each gene panel keep its own y scale (default `FALSE`: one shared
  axis, which is what shows the cross-cohort differences).

- ...:

  Reserved, currently unused.

## Value

A `ggplot` object.

## See also

[`batch_diagnostics`](https://wangjin93.github.io/CanPAS/reference/batch_diagnostics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
   b <- batch_diagnostics(c("GSE31210", "GSE37745"), genes = "GAPDH", type = "DFS")
   plot_batch_diagnostics(b)
} # }
```
