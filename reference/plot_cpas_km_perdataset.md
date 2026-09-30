# Grid of per-cohort Kaplan-Meier panels

Splits a
[`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md)
result back into one Kaplan-Meier panel per dataset (each cohort split
at its own median) and lays them out as a grid.

## Usage

``` r
plot_cpas_km_perdataset(km, ncol = 4, draw = TRUE)
```

## Arguments

- km:

  An object returned by
  [`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md).

- ncol:

  Number of columns in the grid.

- draw:

  `TRUE` (default) draws the grid on the current device with
  [`gridExtra::grid.arrange`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html);
  `FALSE` returns a composable `patchwork` object instead, so the caller
  can combine it with other figures (passing the return value of
  `grid.arrange` to patchwork does not error but is silently dropped,
  which is why the composed case needs `draw = FALSE`).

## Value

Invisibly `NULL` when `draw = TRUE` (the figure has been drawn); a
`patchwork` object when `draw = FALSE`.
