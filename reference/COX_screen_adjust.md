# Univariate -\> Multivariate Cox Analysis Workflow

Runs univariate Cox models for all supplied covariates, selects
covariates significant at `p.threshold` (main effect of a categorical
variable = any of its levels), then - if at least two covariates
qualify - runs one multivariable Cox model and returns both univariate
and multivariate estimates side by side, optionally as a formatted
`flextable`.

## Usage

``` r
COX_screen_adjust(
  df,
  type = "OS",
  cont_Variates = NULL,
  cate_Variates = NULL,
  precision = 3,
  p.threshold = 0.05,
  auto_repair = FALSE
)
```

## Arguments

- df:

  Data.frame (see
  [`COX_analysis`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md)
  for the required layout, including `type`).

- type:

  Survival endpoint prefix, default `"OS"`.

- cont_Variates:

  Continuous covariate column names.

- cate_Variates:

  Categorical covariate column names.

- precision:

  Decimal places for formatted values.

- p.threshold:

  Significance threshold for entering the multivariate step.

- auto_repair:

  Passed to
  [`COX_analysis`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md)
  for the multivariate step. `FALSE` (default) is the fail-safe: when
  the multivariate model is not estimable exactly as requested, it is
  refused with its reasons and offending terms
  (`multi_metadata$estimable` is `FALSE`) instead of being silently
  reduced. `TRUE` restores the automatic repair.

## Value

List with:

- `result`::

  long data.frame with univariate and (if available) multivariate
  HR/CI/p per covariate and level

- `print_result`::

  a `flextable` of `result` (a three-line table ready for
  [`flextable::save_as_docx()`](https://rdrr.io/pkg/flextable/man/save_as_docx.html))
  when the suggested flextable package is installed, otherwise the plain
  `data.frame` with the same numbers, formatted as a three-line table
  (top rule, rule under the header, bottom rule) with the univariable
  and multivariable columns grouped under their own header row, ready
  for
  [`flextable::save_as_docx()`](https://rdrr.io/pkg/flextable/man/save_as_docx.html)
  or printing

- `sig_variates`::

  covariates that entered the multivariate step

- `uni_table`, `multi_table`::

  raw tables of each step

- `multi_metadata`::

  metadata of the multivariate step, including `reduced` and
  `reduced_note` when the model had to be reduced below `min_covariates`
  (see
  [`COX_analysis`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md));
  a message states the same

- `cont_Variates`, `cate_Variates`::

  echoed inputs

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Real cohort with real covariates: expression + survival + clinical
   ## columns are fetched and merged by the prerequisite functions.
   d <- cohort_merged("GSE13507", "GAPDH", type = "OS", clin = TRUE)

   q <- COX_screen_adjust(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                           cate_Variates = c("grade", "N"), p.threshold = 0.05)
   q$sig_variates
   head(q$result)
   ## multi_metadata$reduced is TRUE when the joint model had to be reduced
   q$multi_metadata$reduced_note
} # }
```
