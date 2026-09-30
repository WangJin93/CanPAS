# Canonical CanPAS cancer type of a TCGA project

Canonical CanPAS cancer type of a TCGA project

## Usage

``` r
canonical_type(dataset)
```

## Arguments

- dataset:

  TCGA dataset id (`"LUAD"` or `"TCGA-LUAD"`). A GEO or CGGA accession
  is not a TCGA id and returns `NA`.

## Value

The canonical cancer-type label or `NA_character_` when the project is
unknown.

## Examples

``` r
  ## Real TCGA accessions from the catalog, and the type label CanPAS reports.
  canonical_type("TCGA-LUAD")
#> [1] "Lung Cancer"
  canonical_type("TCGA-ACC")
#> [1] "Adrenocortical Cancer"
  canonical_type("TCGA-XXXX")   # unknown project -> NA
#> [1] NA
```
