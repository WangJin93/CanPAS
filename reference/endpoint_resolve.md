# Resolve an endpoint family for one cohort

Returns the concrete endpoint token available for a cohort and endpoint
(family name or raw token), using the columns of the packaged
`dataset_info` catalog (e.g. `EP_DFS`).

## Usage

``` r
endpoint_resolve(dataset, type)
```

## Arguments

- dataset:

  Dataset accession (e.g. `"GSE31210"`, `"TCGA-LGG"`).

- type:

  Family name (`"OS"`, `"DSS"`, `"DFS"`, `"PFS"`, `"MFS"`) or a raw
  endpoint token.

## Value

Endpoint token (e.g. `"RFS"`). A family is resolved to the token the
cohort actually provides; a raw token is returned only when it is a
known endpoint token, otherwise `NA_character_` (which is also returned
when the cohort does not provide the family).

## Examples

``` r
  ## The token a real cohort actually contributes for a family.
  endpoint_resolve("GSE31210", "DFS")   # "RFS": that cohort is RFS-based
#> [1] "RFS"
  endpoint_resolve("GSE13507", "OS")    # "OS"
#> [1] "OS"
  endpoint_resolve("TCGA-LUAD", "PFS")  # "PFI"
#> [1] "PFI"

  ## NA is returned when the cohort has no endpoint of that family.
  endpoint_resolve("GSE31210", "DSS")
#> [1] NA
```
