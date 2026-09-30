# Local TCGA clinical/survival tables used by the TCGA helpers

Reports which local clinical and survival tables the TCGA helpers
([`tcga_surv_table`](https://wangjin93.github.io/CanPAS/reference/tcga_surv_table.md),
[`tcga_merged`](https://wangjin93.github.io/CanPAS/reference/tcga_merged.md),
[`cohort_merged`](https://wangjin93.github.io/CanPAS/reference/cohort_merged.md))
read, and where each one comes from. The two tables are small (196 KB
together) and are **bundled with the package**, so TCGA cohorts can be
analysed as soon as CanPAS is installed, without setting any environment
variable. `CPAS_DATA_ROOT` is an optional override that points the
helpers at a project checkout instead.

Each file is resolved by taking the first candidate that exists, in this
order: (a) an explicit override - `options(CanPAS.tcga_clinical_rda = )`
/ `options(CanPAS.tcga_survival_rda = )`, or a `TCGA_CLI_RDA` /
`TCGA_SURV_RDA` object assigned in the global environment; (b)
`CPAS_DATA_ROOT` when set and `<root>/data/tcga/<file>` exists; (c) the
copy bundled inside the installed package
(`system.file("extdata/tcga", package = "CanPAS")`); (d) otherwise the
file is unresolvable and the helpers stop with a message that names both
ways to supply it and lists the paths that were tried.

## Usage

``` r
tcga_local_tables()
```

## Value

A list with `clinical` and `survival` (the resolved paths,
`NA_character_` when a file could not be resolved), `source` (named
character vector, one of `"override"`, `"CPAS_DATA_ROOT"` or `"package"`
per file), `found` (single logical: both files exist) and `tried` (named
list of the candidate paths per file, named by the rule that produced
them).

## Examples

``` r
  ## Which tables the TCGA helpers will read, and where they come from.
  t <- tcga_local_tables()
  t$found
#> [1] TRUE
  t$source
#>  clinical  survival 
#> "package" "package" 
  basename(t$clinical)
#> [1] "tcga_clinical.rda"
```
