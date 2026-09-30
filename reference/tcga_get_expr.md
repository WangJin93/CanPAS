# Multi-gene TCGA expression (wide format)

Fetches each requested gene for a TCGA project from UCSC Xena and merges
them into one wide data.frame (rows = samples).

## Usage

``` r
tcga_get_expr(dataset = "LUAD", genes = "TP53")
```

## Arguments

- dataset:

  TCGA dataset id (`"LUAD"` or `"TCGA-LUAD"`).

- genes:

  Character vector of gene symbols.

## Value

data.frame with first column `sample` and one numeric column per
requested gene (columns follow the order of `genes`).

## Examples

``` r
if (FALSE) { # \dontrun{
   ## The expression step alone, for a TCGA project (same schema as
   ## get_expr_data() for GEO cohorts).
   e <- tcga_get_expr("LUAD", c("TP53", "GAPDH"))
   colnames(e$expr_data)
   head(e$expr_data)
} # }
```
