# Analysis manifest of a CanPAS result

Returns the analysis manifest of a result: the cohorts, the endpoint
family and the token each cohort actually contributed, the pooling
classes that were pooled, how the cohorts were selected, which rows and
covariates were dropped and why, the cut-point rule and how many
cut-points were searched, the search-adjusted p-value when one was
computed, the proportional-hazards (`cox.zph`) result, the meta-analysis
method, tau\\^2\\, \\I^2\\, the prediction interval(s), the R and
package versions, and the timestamp.

Every analysis result carries the same object as its `$manifest`
element; `cpas_manifest(result)` reads it. The manifest prints as a
readable block and converts with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) into a
two-column `field`/`value` table; the structured pieces (per-cohort
pooling classes, dropped rows/covariates, overlapping pairs) are
data.frames inside the manifest itself.

## Usage

``` r
cpas_manifest(x = NULL, ...)
```

## Arguments

- x:

  A result object returned by an analysis function
  ([`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md),
  [`COX_analysis`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md),
  [`COX_screen_adjust`](https://wangjin93.github.io/CanPAS/reference/COX_screen_adjust.md),
  [`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md),
  ...). When omitted, an empty manifest skeleton with the same field
  list is returned.

- ...:

  Unused; accepted so that `cpas_manifest(result)` and `cpas_manifest()`
  take the same shape as the other accessors.

## Value

An object of class `cpas_manifest` (a named list).

## See also

[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
  cpas_manifest(m)                 # readable block
  as.data.frame(cpas_manifest(m))  # field / value table
  cpas_manifest(m)$pooling_classes # what was pooled
} # }
```
