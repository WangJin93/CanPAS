# Forest Plot of Cox Regression Results

Draws a forest plot (point = HR, error bar = 95% CI) from a tidy Cox
results table as produced by
[`COX_analysis`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md)
or
[`COX_screen_adjust`](https://wangjin93.github.io/CanPAS/reference/COX_screen_adjust.md).
Rows are coloured by significance and risk direction; categorical
variable names are shown as group headers.

## Usage

``` r
forest_plot(
  COX_out,
  which = c("uni", "multi"),
  HR_threshold = 1,
  p_threshold = 0.05,
  colors = c("#00AFBB", "black", "#FC4E07"),
  log_x = TRUE,
  HR_order = c("decrease", "increase", "none"),
  group_levels = c("label", "block"),
  point_size = 2.2,
  bar_height = 0.3,
  digits = 4
)
```

## Arguments

- COX_out:

  A data.frame with (at least) columns `Variates` (or `gene`, as
  returned by
  [`COX_by_genes`](https://wangjin93.github.io/CanPAS/reference/COX_by_genes.md)),
  `HR`, `HR95L`, `HR95H`, `Pvalue`; optionally columns `Type` (`cont` /
  `cate_header` / `cate_level`) and `Level`. If the columns `HR_uni` ...
  `Pvalue_uni` from
  [`COX_screen_adjust`](https://wangjin93.github.io/CanPAS/reference/COX_screen_adjust.md)
  are present instead, use `which = "uni"` or `"multi"`.

- which:

  Which set of columns to plot when `COX_out` contains both univariate
  and multivariate estimates (`"uni"` or `"multi"`).

- HR_threshold:

  Reference value for risk direction (dashed vertical line). Default 1.

- p_threshold:

  Significance threshold for colouring. Default 0.05.

- colors:

  Vector of three colours: \[low-risk, non-significant, high-risk\].
  Defaults `c("#00AFBB", "black", "#FC4E07")`.

- log_x:

  Whether to plot the x-axis on a log2 scale. Default TRUE.

- HR_order:

  Ordering of estimator rows: `"none"` (input order), `"decrease"`
  (default) or `"increase"` by HR.

- group_levels:

  How the rows of a categorical covariate stay identifiable once the
  rows are reordered (only used when `HR_order` is not `"none"`).
  `"label"` (default) sorts strictly by HR: a categorical level is then
  no longer adjacent to its covariate name, so every level is labelled
  `"covariate: level"` and the separate header rows are dropped.
  `"block"` keeps the levels of each covariate together, sorts them
  inside their block and orders the blocks by HR (by the most extreme
  level in the requested direction), so the bold covariate header still
  sits directly above its own levels.

- point_size:

  Point size.

- bar_height:

  Height of the error bars.

- digits:

  Decimal places for the printed estimates (default 4).

## Value

A ggplot object. The plotted rows and their order are in the returned
object's `data`, so a caller can verify that each label and its estimate
belong together.

## Details

A categorical covariate contributes one row per level, and a level name
alone ("G3") is only meaningful next to its covariate name. Sorting rows
by HR therefore interleaves covariates; the previous implementation
re-attached the covariate headers by moving whole position ranges, which
could place a level under another covariate's header. Now either the
level labels carry the covariate name (`group_levels = "label"`,
default, strictly sorted) or the covariate blocks stay contiguous
(`group_levels = "block"`). `HR_order = "none"` leaves the input order
untouched, header rows included.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Fit a real model first (expression + survival + clinical covariates
   ## come from the CanPAS mirror through cohort_merged()).
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS", clin = TRUE)
   r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                      cate_Variates = c("sex", "stage"), method = "uni")

   forest_plot(r$results_table)

   ## The same figure for the multivariable model, on a linear scale.
   r2 <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                       cate_Variates = c("sex", "stage"), method = "multi")
   forest_plot(r2$results_table, log_x = FALSE, p_threshold = 0.01)
} # }
```
