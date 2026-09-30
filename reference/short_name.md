# Short dataset name of a TCGA project

Prefixes a TCGA project id with `"TCGA-"`. Only the first element of
`dataset` is used, and the input must be a TCGA id (`"LUAD"` or
`"TCGA-LUAD"`) — a GEO accession such as `"GSE13507"` would come back as
`"TCGA-GSE13507"`, which is why GEO cohorts are labelled by their
catalog accessions instead.

## Usage

``` r
short_name(dataset)
```

## Arguments

- dataset:

  TCGA dataset id (`"LUAD"` or `"TCGA-LUAD"`).

## Value

A single string, `"TCGA-<PROJECT>"` (only the first element of `dataset`
is used).

## Examples

``` r
  ## Short labels used in figures and tables (TCGA ids only).
  short_name("LUAD")
#> [1] "TCGA-LUAD"
  short_name("TCGA-LUAD")
#> [1] "TCGA-LUAD"
```
