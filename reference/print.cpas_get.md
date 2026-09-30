# Print Method for cpas_get Class

Custom print method for the result object returned by get_data function.

## Usage

``` r
# S3 method for class 'cpas_get'
print(x, ...)
```

## Arguments

- x:

  An object of class cpas_get.

- ...:

  Additional arguments passed to print.

## Value

No return value, prints a summary of the API response.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## The raw API accessor, printed.
   r <- get_data("GSE14814", "surv_data")
   print(r)
} # }
```
