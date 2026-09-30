# TCGA clinical + survival table

Loads the local TCGA clinical and survival tables (bundled with the
package; see
[`tcga_local_tables`](https://wangjin93.github.io/CanPAS/reference/tcga_local_tables.md)
for the resolution order, and set `CPAS_DATA_ROOT` to override them with
a project checkout) and returns one row per tumour sample of a project,
with CanPAS-compatible columns: `ID`, `<endpoint>_status` /
`<endpoint>_time` (years), `age`, `sex`, `stage` and, when present,
`histology` and `grade`.

## Usage

``` r
tcga_surv_table(dataset = "LUAD", endpoints = c("OS", "DSS", "DFI", "PFI"))
```

## Arguments

- dataset:

  TCGA dataset id (`"LUAD"` or `"TCGA-LUAD"`).

- endpoints:

  Endpoint prefixes to export.

## Value

data.frame (see description).

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Local clinical table of a project: all endpoint families at once.
   s <- tcga_surv_table("LUAD")
   colnames(s)
   head(s)

   ## Only the endpoints the project really carries.
   s2 <- tcga_surv_table("LUAD", endpoints = c("OS", "DSS"))
   colnames(s2)
} # }
```
