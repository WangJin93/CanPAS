# TCGA cohort merged with expression (CanPAS format)

Combines on-demand Xena expression
([`tcga_get_expr`](https://wangjin93.github.io/CanPAS/reference/tcga_get_expr.md))
with the local clinical/survival table
([`tcga_surv_table`](https://wangjin93.github.io/CanPAS/reference/tcga_surv_table.md))
into a single data.frame whose schema matches GEO merged data: first
column `ID`, then `<type>_time` / `<type>_status` and clinical columns,
then one numeric column per gene.

## Usage

``` r
tcga_merged(dataset = "LUAD", genes = c("TP53"), type = "OS")
```

## Arguments

- dataset:

  TCGA dataset id (`"LUAD"` or `"TCGA-LUAD"`).

- genes:

  Gene symbols.

- type:

  Endpoint used downstream (must exist for the project).

## Value

data.frame.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## TCGA keeps its clinical columns, so covariates are usable directly.
   d <- tcga_merged("LUAD", c("TP53", "GAPDH"), type = "OS")
   colnames(d)

   r <- COX_analysis(d, type = "OS", cont_Variates = c("TP53", "age"),
                      cate_Variates = c("sex", "stage"), method = "uni")
   head(r$results_table)
} # }
```
