# Print Method for cpas_COX_by_genes Class

Custom print method for the result object returned by COX_by_genes().

## Usage

``` r
# S3 method for class 'cpas_COX_by_genes'
print(x, ...)
```

## Arguments

- x:

  An object of class cpas_COX_by_genes.

- ...:

  Additional arguments passed to print.

## Value

No return value, prints a summary of the univariable Cox regression
results.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A real per-gene screen on a real cohort.
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB", "TP53"), type = "OS")
   r <- COX_by_genes(d, type = "OS", genes = c("GAPDH", "ACTB", "TP53"))
   print(r)
} # }
```
