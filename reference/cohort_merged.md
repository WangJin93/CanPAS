# Unified CanPAS cohort reader (GEO or TCGA)

Returns one merged data.frame for a cohort regardless of the source:
datasets prefixed `"TCGA-"` are handled through the TCGA helpers (Xena +
local clinical tables); anything else is treated as a GEO accession
served by the CanPAS API. The output schema is identical for both
sources, so downstream functions of the package work unchanged.

## Usage

``` r
cohort_merged(
  dataset,
  genes = "TP53",
  type = "OS",
  process_duplicates = "max",
  clin = FALSE
)
```

## Arguments

- dataset:

  Dataset identifier, e.g. `"GSE14814"` or `"TCGA-LUAD"`.

- genes:

  Gene symbols to retrieve.

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token; families are
  resolved to the token available in the cohort and recorded as the
  `endpoint` attribute of the returned data.frame.

- process_duplicates:

  Probe collapsing rule used for GEO queries.

- clin:

  If `TRUE`, clinical covariates supplied for that cohort are appended:
  `age`, `sex`, `stage`, `T`, `N`, `M`, `grade`, `histology`, whichever
  the cohort carries. GEO cohorts take them from the mirrored clinical
  table that
  [`merge_surv_expr()`](https://wangjin93.github.io/CanPAS/reference/merge_surv_expr.md)
  returns in `raw_surv_data`; TCGA cohorts already include them. Columns
  missing for a cohort are simply absent. The mirror stores these
  columns as text, so they are converted: a column whose values all
  parse as numbers (`age`, often `T`, `N`, `M`) becomes numeric, every
  other one becomes a factor. Default `FALSE` (survival + genes only).

## Value

data.frame: first column `ID` followed by survival columns, one numeric
column per gene, and (with `clin = TRUE`) the clinical covariates
available for that cohort.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## One cohort, one call: the endpoint token is resolved from the catalog,
   ## expression is fetched (mirror for GEO, Xena for TCGA), survival merged.
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
   attr(d, "family"); attr(d, "endpoint")
   colnames(d)

   ## clin = TRUE adds the clinical covariates that cohort really carries, so
   ## a multivariable model can use them directly.
   d2 <- cohort_merged("GSE13507", "GAPDH", type = "OS", clin = TRUE)
   colnames(d2)

   ## A TCGA project goes through the same function and has the same schema.
   tcga <- cohort_merged("TCGA-LUAD", c("TP53", "GAPDH"), type = "OS")
   colnames(tcga)[1:6]
} # }
```
