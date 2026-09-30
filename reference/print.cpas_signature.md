# Print Method for cpas_signature Class

Custom print method for the result object returned by
get_signature_value function.

## Usage

``` r
# S3 method for class 'cpas_signature'
print(x, ...)
```

## Arguments

- x:

  An object of class cpas_signature.

- ...:

  Additional arguments passed to print.

## Value

No return value, prints a summary of the gene signature calculation.

## Examples

``` r
if (FALSE) { # \dontrun{
   s <- get_signature_value("0.5*GAPDH + 0.5*ACTB", "GSE14814")
   print(s)
} # }
```
