# Shared-patient overlap register

Returns the register of cohort pairs that share patients, with the
number of shared patients, the basis of the finding and the evidence it
was read from. The table is shipped with the package
(`inst/extdata/cohort_overlap.csv`) and is a byte-identical copy of the
curation pipeline's `data/suppl/cohort_overlap.csv`.

[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
consults it through its `overlap` argument.

## Usage

``` r
cohort_overlap(path = NULL)
```

## Arguments

- path:

  Optional path of the register CSV; defaults to the shipped copy.

## Value

data.frame with the frozen columns `AccessionA`, `AccessionB`,
`SharedPatients`, `Basis`, `Evidence`. Errors when the shipped table is
absent (it is produced by the curation pipeline, not by the package).

## See also

[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  ov <- cohort_overlap()
  head(ov)
  ## Which of my cohorts overlap at all?
  ov[ov$AccessionA %in% c("GSE11969", "GSE13213") |
     ov$AccessionB %in% c("GSE11969", "GSE13213"), ]
} # }
```
