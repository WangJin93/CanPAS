# Plot cumulative incidence curves

Draws the cumulative incidence functions returned by
[`cif_fit`](https://wangjin93.github.io/CanPAS/reference/cif_fit.md) as
step curves with pointwise confidence ribbons.

## Usage

``` r
plot_cif(x, cause = NULL, ...)
```

## Arguments

- x:

  An object returned by
  [`cif_fit`](https://wangjin93.github.io/CanPAS/reference/cif_fit.md).

- cause:

  Event type (cause code) to plot; default: all.

- ...:

  Unused, kept for S3 compatibility.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A real cohort that carries both OS and DSS, through the prerequisite
   ## chain: cause 1 = cancer death (DSS event), cause 2 = other death.
   d  <- cohort_merged("GSE14814", "GAPDH", type = "OS", clin = TRUE)
   cr <- data.frame(
     time   = d$OS_time,
     status = ifelse(d$OS_status == 0, 0, ifelse(d$DSS_status == 1, 1, 2)),
     age    = d$age,                        # continuous
     sex    = factor(d$sex),                # categorical
     age_group = factor(ifelse(d$age > stats::median(d$age, na.rm = TRUE),
                                 "older", "younger")))
   table(cr$status)   # 0 censored, 1 cancer death, 2 other death

   x <- cif_fit(cr, time = "time", status = "status", times = c(1, 3, 5))
   plot_cif(x)             # cumulative incidence of cause 1
   plot_cif(x, cause = 2)

   ## The same figure per group, and per tertile of a continuous variable.
   plot_cif(cif_fit(cr, time = "time", status = "status", group = "sex"))
   plot_cif(suppressMessages(
     cif_fit(cr, time = "time", status = "status", group = "age",
             group_cut = "tertile")))
} # }
```
