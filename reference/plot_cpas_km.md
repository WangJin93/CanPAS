# Plot pooled / meta Kaplan-Meier result

Visualizes a
[`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md)
result: when method `"ipd"`/`"both"` was used a classical stratified
Kaplan-Meier plot is drawn; otherwise (method `"meta"`) the merged meta
survival curve with confidence band is plotted.

## Usage

``` r
plot_cpas_km(km, which = c("auto", "ipd", "meta"))
```

## Arguments

- km:

  An object returned by
  [`cpas_km_pooled`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md).

- which:

  Which route to draw: `"auto"` (default) draws the
  individual-patient-data curve when that route was run and the
  time-point meta curve otherwise, exactly as before; `"ipd"` or
  `"meta"` selects a route explicitly. A `method = "both"` result holds
  both, so a caller that wants to show them together (as the Shiny
  application does) can draw each in turn instead of only ever getting
  the IPD curve.

## Value

A `ggsurvplot` object for the IPD route (with its risk table) or a
ggplot object for the meta route, carrying the requested landmark
estimates and their confidence intervals as points on the curve.
