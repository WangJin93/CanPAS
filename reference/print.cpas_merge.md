# Print Method for cpas_merge Class

Custom print method for the result object returned by merge_surv_expr
function.

## Usage

``` r
# S3 method for class 'cpas_merge'
print(x, ...)
```

## Arguments

- x:

  An object of class cpas_merge.

- ...:

  Additional arguments passed to print.

## Value

No return value, prints a summary of the merged data.

## Examples

``` r
if (FALSE) { # \dontrun{
   e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
   m <- merge_surv_expr("GSE14814", e)
   print(m)
} # }
```
