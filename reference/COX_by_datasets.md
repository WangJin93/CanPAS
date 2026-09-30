# Cox Regression for Many Datasets (Single Gene)

Runs the same univariable Cox analysis (one gene) on several datasets
(GEO accessions or `TCGA-<PROJECT>` projects) and returns one combined
tidy table plus the per-dataset objects. The estimates are reported per
dataset and are NOT pooled: use
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
for an inverse-variance meta-analysis with heterogeneity statistics.
Expected hazard ratios are per raw expression unit, so they are not
comparable across platforms; `cpas_meta` standardises the gene within
each dataset (per SD) before pooling.

## Usage

``` r
COX_by_datasets(
  datasets,
  gene,
  type = "OS",
  precision = 3,
  process_duplicates = "max"
)
```

## Arguments

- datasets:

  Character vector of dataset accessions or TCGA projects (e.g.
  `c("GSE14814", "GSE31210")` or `"TCGA-LUAD"`).

- gene:

  Single gene symbol to test.

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token; families are
  resolved per cohort and the token used is reported in `endpoint`.

- precision:

  Decimal places for formatted display columns.

- process_duplicates:

  Probe collapsing rule (see
  [`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)).
  Default `"max"`.

## Value

List of class `cpas_multi_cox`:

- `input_params`::

  echoed inputs

- `individual_results`::

  per successful dataset, a list with `expr_data`, `merged_data` and
  `cox_analysis`

- `combined_results`::

  data.frame with columns `dataset`, `endpoint` (the token actually
  used), `Variates` (gene), `Level`, `N`, `HR`, `HR95L`, `HR95H`,
  `Pvalue` and a Benjamini-Hochberg `P_adj` across all tested dataset x
  gene combinations

- `metadata`, `errors`::

  bookkeeping

## Examples

``` r
if (FALSE) { # \dontrun{
   ## The same gene in real cohorts of one cancer type (lung, OS).
   di <- dataset_info[dataset_info$Type == "Lung Cancer", "Accession"]
   head(di)

   r <- COX_by_datasets(c("GSE14814", "GSE31210"), gene = "GAPDH", type = "OS")
   head(r$results_table)
} # }
```
