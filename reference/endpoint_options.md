# Endpoints available for a cohort, grouped by family

Returns one row per endpoint family available for a cohort, with the
concrete endpoint token, whether the token is a TCGA-derived endpoint,
the analysis family and a display label such as `"DFS (RFS)"` or
`"DFS (DFI) [derived]"`.

## Usage

``` r
endpoint_options(dataset, family = NULL)
```

## Arguments

- dataset:

  Dataset accession (e.g. `"GSE31210"`, `"TCGA-LGG"`).

- family:

  Optional family filter (one or more of `"OS"`, `"DSS"`, `"DFS"`,
  `"PFS"`, `"MFS"`).

## Value

data.frame with columns `family`, `token`, `derived`, `label`.

## Examples

``` r
  ## What a real cohort can be analysed on, with the token per family.
  endpoint_options("TCGA-LGG")
#>   family token derived               label
#> 1     OS    OS   FALSE                  OS
#> 2    DSS   DSS   FALSE                 DSS
#> 3    DFS   DFI    TRUE DFS (DFI) [derived]
#> 4    PFS   PFI    TRUE PFS (PFI) [derived]
  endpoint_options("GSE13507")
#>   family token derived label
#> 1     OS    OS   FALSE    OS

  ## Only one family at a time, when that is all you need.
  endpoint_options("GSE31210", family = "DFS")
#>   family token derived     label
#> 1    DFS   RFS   FALSE DFS (RFS)
```
