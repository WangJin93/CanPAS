# Quick Kaplan-Meier Plots for a Dataset

One-stop helper: retrieves expression of the requested genes for a
dataset, merges it with the survival table, and draws a grid of
Kaplan-Meier plots (median split) for every gene that maps on the
platform.

## Usage

``` r
quick_km(
  dataset,
  genes,
  process_duplicates = "max",
  type = "OS",
  ncol = 3,
  cutpoint = "median",
  ...
)
```

## Arguments

- dataset:

  Dataset accession in `dataset_info` (GEO), or `"TCGA-<PROJECT>"` (e.g.
  `"TCGA-LUAD"`).

- genes:

  Character vector of gene symbols.

- process_duplicates:

  Collapsing rule for multiple probes per gene (see
  [`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)).
  Default `"max"`.

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token. A family is
  resolved to the token available in the cohort (e.g. DFS -\>
  RFS/DFI/EFS).

- ncol:

  Number of columns in the plot grid.

- cutpoint:

  Cut-point rule passed to
  [`plot_km`](https://wangjin93.github.io/CanPAS/reference/plot_km.md).

- ...:

  Further arguments passed to
  [`plot_km`](https://wangjin93.github.io/CanPAS/reference/plot_km.md) /
  [`survminer::ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html).

## Value

Invisibly, the list of `ggsurvplot` objects; the combined grid is drawn
to the active graphics device.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## quick_km() runs the same prerequisite chain internally, for one cohort
   ## and one or more genes, and draws the curves.
   quick_km("GSE14814", c("GAPDH", "ACTB"), pval = TRUE, legend = "bottom")

   ## The same data through the explicit chain, for customisation.
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
   plot_km(d, type = "OS", marker = "GAPDH", cutpoint = 12)   # absolute cut-off
} # }
```
