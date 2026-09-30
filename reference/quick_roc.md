# Quick Time-dependent ROC Plots for a Dataset

One-stop helper: retrieves expression of the requested genes for a
dataset, merges it with the survival table, and draws a grid of
time-dependent ROC curves for every gene that maps on the platform.

## Usage

``` r
quick_roc(
  dataset,
  genes,
  process_duplicates = "max",
  type = "OS",
  predict.time = 1,
  method = c("KM", "NNE"),
  ncol = 3,
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

- predict.time:

  Prediction horizon (same time unit as the survival column, typically
  years).

- method:

  ROC method passed to
  [`plot_roc`](https://wangjin93.github.io/CanPAS/reference/plot_roc.md).

- ncol:

  Number of columns in the plot grid.

- ...:

  Further arguments passed to
  [`plot_roc`](https://wangjin93.github.io/CanPAS/reference/plot_roc.md).

## Value

Invisibly, the list of ggplot objects; the combined grid is drawn to the
active graphics device.

## Examples

``` r
if (FALSE) { # \dontrun{
   quick_roc("GSE14814", c("GAPDH", "ACTB"), predict.time = 3)

   ## Equivalent, using the prerequisite chain for a single cohort.
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
   plot_roc(d, type = "OS", marker = "GAPDH", predict.time = 3)$auc
} # }
```
