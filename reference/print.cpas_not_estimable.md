# Print Method for cpas_not_estimable Class

Prints a fail-safe refusal: the multivariable model was not estimable
exactly as requested and was therefore not fitted. The reasons and the
offending terms are printed, together with what the automatic repair
would have changed (`auto_repair = TRUE`).

## Usage

``` r
# S3 method for class 'cpas_not_estimable'
print(x, ...)
```

## Arguments

- x:

  An object of class `cpas_not_estimable` (also a `cpas_COX` object).

- ...:

  Additional arguments passed to print.

## Value

Invisibly returns `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A covariate that cannot be co-estimated: the model is refused, with
   ## reasons, instead of being repaired (the package default).
   d <- cohort_merged("GSE13507", "GAPDH", type = "OS", clin = TRUE)
   d$flat <- 1
   r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "flat"),
                     method = "multi")
   r$estimable            # FALSE
   r$offending_terms      # "flat"
   print(r)
} # }
```
