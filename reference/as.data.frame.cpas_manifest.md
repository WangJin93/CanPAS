# Convert an analysis manifest to a data.frame

Flattens the manifest into a two-column `field`/`value` table so it can
be written to CSV, compared between runs, or embedded in a report.
Table-valued fields are collapsed into one row each; the structured
tables remain in the manifest itself (`$pooling_table`, `$dropped_rows`,
`$dropped_covariates`, `$overlap_pairs`).

## Usage

``` r
# S3 method for class 'cpas_manifest'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `cpas_manifest`.

- row.names, optional, ...:

  Passed to `as.data.frame` for compatibility.

## Value

data.frame with columns `field` and `value`.
