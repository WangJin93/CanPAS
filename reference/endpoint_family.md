# Pooling family of a survival endpoint

Maps a raw endpoint token (e.g. `"RFS"`) to the pooling family used by
meta-analysis and pooled Kaplan-Meier (`OS`, `DSS`, `DFS`, `PFS`,
`MFS`). Family names are returned unchanged, unknown tokens give `NA`.

## Usage

``` r
endpoint_family(type)
```

## Arguments

- type:

  Endpoint token or family name.

## Value

Character scalar (family) or `NA_character_`.

## Examples

``` r
  ## Families group the tokens that different cohorts use for the same
  ## endpoint, so one analysis can span cohorts.
  endpoint_family("RFS")   # "DFS"
#> [1] "DFS"
  endpoint_family("PFI")   # "PFS"
#> [1] "PFS"
  endpoint_family("OS")    # "OS"
#> [1] "OS"

  ## The catalog records which families each real cohort can be analysed on.
  di <- dataset_info[dataset_info$Accession == "GSE13507", ]
  di[, c("Accession", "Type", "N", "EndpointFamilies")]
#>    Accession           Type   N EndpointFamilies
#> 81  GSE13507 Bladder Cancer 165               OS
```
