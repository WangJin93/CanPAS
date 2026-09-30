# Endpoint semantics of the curated cohorts

Returns the companion table that records, for every cohort and pooling
family, what the endpoint actually measures and which pooling class it
belongs to. The table is shipped with the package
(`inst/extdata/endpoint_semantics.csv`) and is a byte-identical copy of
the curation pipeline's `data/endpoint_semantics.csv`.

[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
uses it to decide what `pooling = "exact"` keeps;
[`endpoint_pooling_class`](https://wangjin93.github.io/CanPAS/reference/endpoint_pooling_class.md)
is the one-cell accessor.

## Usage

``` r
endpoint_semantics(accession = NULL, family = NULL, path = NULL)
```

## Arguments

- accession:

  Optional accession filter (e.g. `"GSE31210"`).

- family:

  Optional pooling-family filter (`OS`, `DSS`, `DFS`, `PFS`, `MFS`).

- path:

  Optional path of the companion CSV; defaults to the shipped copy.

## Value

data.frame with the frozen columns `Accession`, `Type`, `Family`,
`Token`, `TokenRole`, `SourceField`, `EventDefinition`, `TimeOrigin`,
`CensoringRule`, `CompetingEvents`, `Derived`, `PoolingClass`,
`Evidence`, `Note`. Errors when the shipped table is absent (it is
produced by the curation pipeline, not by the package).

## See also

[`endpoint_pooling_class`](https://wangjin93.github.io/CanPAS/reference/endpoint_pooling_class.md),
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  ## One row per cohort x endpoint family, with the pooling class.
  es <- endpoint_semantics(family = "DFS")
  table(es$PoolingClass)

  ## The class of one cell, the value cpas_meta(pooling = "exact") uses.
  endpoint_pooling_class("GSE31210", "DFS")   # "Clinically-related" (RFS)
  endpoint_pooling_class("GSE13507", "OS")    # "Exact-equivalent"
  endpoint_pooling_class("GSE13507", "DFS")   # "Absent" (no DFS endpoint)
} # }
```
