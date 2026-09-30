# Print Method for cpas_COX Class

Custom print method for the result object returned by COX_analysis
function.

## Usage

``` r
# S3 method for class 'cpas_COX'
print(x, ...)
```

## Arguments

- x:

  An object of class cpas_COX.

- ...:

  Additional arguments passed to print.

## Value

No return value, prints a summary of the COX regression results.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A real Cox fit: expression + survival + clinical covariates.
   d <- cohort_merged("GSE13507", "GAPDH", type = "OS", clin = TRUE)
   r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                      cate_Variates = c("grade", "N"), method = "uni")
   print(r)
} # }
```
