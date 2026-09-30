# Print Method for cpas_get_expr Class

Custom print method for the result object returned by get_expr_data
function.

## Usage

``` r
# S3 method for class 'cpas_get_expr'
print(x, ...)
```

## Arguments

- x:

  An object of class cpas_get_expr.

- ...:

  Additional arguments passed to print.

## Value

No return value, prints a summary of the gene expression data retrieval.

## Examples

``` r
if (FALSE) { # \dontrun{
   e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
   print(e)
} # }
```
