# Cox Regression for Many Genes (Panel Screening)

Fits one univariable Cox proportional-hazards model per gene
(`Surv(time, status) ~ gene`) and returns hazard ratios, 95% confidence
intervals, Wald p-values and Benjamini-Hochberg FDRs for every gene,
together with the fitted model objects. This is a panel screen: each
gene is estimated separately and the genes are never adjusted for one
another. Use
[`COX_analysis`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md)`(method = "multi")`
for a joint model with adjusted hazard ratios, and
[`COX_by_datasets`](https://wangjin93.github.io/CanPAS/reference/COX_by_datasets.md)
for the same analysis across several datasets.

## Usage

``` r
COX_by_genes(df, type = "OS", genes)
```

## Arguments

- df:

  Data.frame whose first column is `ID` and that contains the columns
  `<type>_time` (years), `<type>_status` (0/1) and one numeric column
  per marker.

- type:

  Endpoint family (`"OS"`, `"DSS"`, `"DFS"`, `"PFS"`, `"MFS"`) or a
  concrete token (`"RFS"`, `"DFI"`, `"PFI"` ...). Families are resolved
  to the token available in `df`.

- genes:

  Character vector of gene (expression) columns to test.

## Value

List of class `cpas_COX_by_genes`:

- `input_params`::

  echoed inputs and analysis time

- `processed_data`::

  analysis data frame after cleaning

- `individual_models`::

  named list of `coxph` fits

- `individual_summaries`::

  named list of model summaries

- `results_table`::

  data.frame with columns `gene`, `HR`, `HR95L`, `HR95H`, `Pvalue` and
  `P_adj` (Benjamini-Hochberg FDR across all genes tested) plus
  `P_adj_text`

- `metadata`::

  sample size, events, gene count, and `failed_genes` /
  `failure_reasons` for genes that could not be estimated (separation,
  non-convergence, constant gene)

## Details

HRs are per one raw unit of the marker. For cross-platform comparability
consider per-SD standardization (see
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)).
Every gene is tested separately, so the table also carries the
Benjamini-Hochberg FDR (`P_adj`): report it when the genes are a
screening panel rather than a pre-specified hypothesis. Rows with
missing time/status or missing gene values are dropped; genes that
produce a degenerate model are skipped with a warning and reported as
`NA` rows in `results_table` plus a note in the metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## One real cohort, several genes: the prerequisite chain gives the table.
   d <- cohort_merged("GSE14814", c("GAPDH", "ACTB", "TP53"), type = "OS")

   r <- COX_by_genes(d, type = "OS", genes = c("GAPDH", "ACTB", "TP53"))
   head(r$results_table)
} # }
```
