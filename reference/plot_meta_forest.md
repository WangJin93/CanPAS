# Forest plot of an integrative (meta) analysis

Draws the per-cohort HR (95% CI) together with the pooled HR of a
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
result as a forest plot.

## Usage

``` r
plot_meta_forest(
  x,
  digits = 4,
  show_stars = TRUE,
  label_size = 4.2,
  star_size = NULL,
  y_headroom = 3.6,
  x_frac = 0.02,
  label_lift = 1,
  label_where = c("inside", "subtitle"),
  ...
)
```

## Arguments

- x:

  An object of class `cpas_meta`.

- digits:

  Decimal places used for the annotations in the figure (default 4).

- show_stars:

  Whether to draw the per-cohort significance stars at the left edge of
  the panel (`*` p\<0.05, `**` p\<0.01, `***` p\<0.001); default `TRUE`.

- label_size:

  Font size of the pooled HR / 95% PI annotation (default 4.2).

- star_size:

  Font size of the stars and of the `Overall` label; defaults to
  `0.85 * label_size`.

- y_headroom:

  Head-room at the top of the panel so the `Overall` diamond and its
  confidence interval are not clipped (default 3.6).

- x_frac:

  Horizontal position, as a fraction from the left edge of the panel,
  used to anchor the stars and the pooled annotation (default 0.02).

- label_lift:

  Vertical offset of the pooled annotation relative to the `Overall` row
  (default 1).

- label_where:

  Where to place the pooled annotation: `"inside"` puts it in the
  top-left corner of the panel and `"subtitle"` above the panel (default
  `"inside"`).

- ...:

  Reserved, currently unused.

## Value

A `ggplot` object.

## Details

The stars and the pooled annotation are anchored to a **finite** data
value derived from the fixed axis expansion, not to `x = -Inf`: an
infinite value becomes `NaN` through the log10 scale, and the text layer
is then silently dropped with a warning. The per-cohort significance
stars sit at the left edge inside the panel, and the pooled HR / 95% PI
annotation is placed in the top-left corner of the panel by default.

## Examples

``` r
if (FALSE) { # \dontrun{
   m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
   plot_meta_forest(m, digits = 4)
} # }
```
