# Print Method for cpas_COX_by_datasets Class

Custom print method for the result object returned by COX_by_datasets().

## Usage

``` r
# S3 method for class 'cpas_COX_by_datasets'
print(x, ...)
```

## Arguments

- x:

  An object of class cpas_COX_by_datasets.

- ...:

  Additional arguments passed to print.

## Value

No return value, prints a summary of the multi-dataset COX analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
   r <- COX_by_datasets(c("GSE14814", "GSE31210"), gene = "GAPDH", type = "OS")
   print(r)
} # }
```
