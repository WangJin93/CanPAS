# TCGA projects available through the TCGA helpers

Named vector mapping TCGA project abbreviations to the cancer type
labels used in CanPAS. It is derived from the shipped catalog (every
[`dataset_info`](https://wangjin93.github.io/CanPAS/reference/dataset_info.md)
row whose `Accession` starts with `"TCGA-"`; the `Type` column supplies
the label), so it always agrees with the catalog and every listed
project is accepted by
[`tcga_project_dataset`](https://wangjin93.github.io/CanPAS/reference/tcga_project_dataset.md),
[`tcga_surv_table`](https://wangjin93.github.io/CanPAS/reference/tcga_surv_table.md),
[`tcga_merged`](https://wangjin93.github.io/CanPAS/reference/tcga_merged.md),
[`tcga_get_expr`](https://wangjin93.github.io/CanPAS/reference/tcga_get_expr.md)
and
[`canonical_type`](https://wangjin93.github.io/CanPAS/reference/canonical_type.md).

## Usage

``` r
tcga_retained
```

## Format

Named character vector, one entry per TCGA cohort in the catalog.
