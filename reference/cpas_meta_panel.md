# Per-gene meta-analysis panel across datasets

Runs the two-stage meta-analysis of
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
for every gene in `genes` over the same set of `datasets` and collects
the pooled estimates, heterogeneity statistics, prediction intervals and
coverage into one table, with Benjamini-Hochberg FDR across the genes.
Expression and survival tables are fetched \*\*once per dataset for all
genes\*\*, so a panel of G genes over D datasets costs D retrievals
rather than G x D.

This is the pooled counterpart of the per-dataset screening table: for
each gene the per-dataset estimates are pooled (inverse variance, random
or fixed effects) instead of merely being listed, and the gene-level
p-values are adjusted for multiple testing.

**Defaults follow
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md).**
The pooling method defaults to `"REML"` - the same default as
`cpas_meta` (spec B10a: the panel previously defaulted to `"DL"`, which
made two calls with the same arguments disagree; `"DL"` remains
available explicitly, and passing `method = "DL"` reproduces the earlier
panel tables exactly). The prediction-interval rule defaults to
`pi_method = "t"`, matching `cpas_meta`.

## Usage

``` r
cpas_meta_panel(
  datasets,
  genes,
  type = "OS",
  method = c("REML", "DL", "HK", "FE", "RE"),
  pi_method = c("t", "normal", "HK"),
  confounders = NULL,
  min_events = 5,
  process_duplicates = "max",
  fdr_method = "BH",
  merged = NULL,
  max_try = 3,
  progress = TRUE
)
```

## Arguments

- datasets:

  Character vector of dataset accessions or TCGA projects.

- genes:

  Character vector of gene symbols (the panel).

- type:

  Endpoint family (OS/DSS/DFS/PFS/MFS) or a concrete token; the family
  is resolved per dataset and the token used is reported in
  `per_dataset$endpoint` and `endpoints`.

- method:

  Pooling method passed to
  [`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md):
  `"REML"` (the default, matching `cpas_meta`), `"DL"`, `"HK"`, `"FE"`;
  `"RE"` is accepted as `"DL"`.

- pi_method:

  Prediction-interval rule passed to
  [`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md):
  `"t"` (default, \\t\_{k-2}\\), `"normal"` or `"HK"`. The rule used is
  recorded in `input$pi_method` and in the manifest; the per-gene table
  carries the primary interval in `pi_lower`/`pi_upper` and the clearly
  labelled alternative in `pi_alt_lower`/`pi_alt_upper`.

- confounders:

  Optional covariates adjusted for inside every dataset.

- min_events:

  Minimum number of events a dataset must contribute (default 5);
  datasets below it are excluded with a recorded reason.

- process_duplicates:

  Probe collapsing rule used when retrieving expression (see
  [`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)).

- fdr_method:

  Method passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)
  for the gene-level FDR (default `"BH"`).

- merged:

  Optional named list of already merged data frames (dataset -\>
  data.frame), which skips all retrieval (used by the test suite and by
  callers who cache their own data).

- max_try:

  Number of retrieval attempts per dataset.

- progress:

  Print a progress message per gene.

## Value

List of class `cpas_meta_panel`:

- `table`::

  one row per gene: `gene`, `k` (datasets pooled), `datasets_used`,
  `total_n`, `total_events`, `HR`, `lower`, `upper`, `p`, `P_adj`,
  `P_adj_text`, `I2`, `tau2`, `pi_lower`, `pi_upper`, `pi_alt_lower`,
  `pi_alt_upper`, `p_heterogeneity`

- `per_dataset`::

  long table: `gene`, `dataset`, `endpoint`, `n`, `events`, `HR`,
  `lower`, `upper`, `p`

- `endpoints`::

  named list (gene -\> named vector of tokens used)

- `errors`::

  reasons why a gene or a dataset was not pooled

- `fetch_errors`::

  datasets that could not be retrieved

- `input`::

  echoed inputs

## Details

Each gene is standardised within each dataset (per SD) before pooling,
so the hazard ratios are comparable across platforms, and the endpoint
token each dataset contributes is reported per gene (a dataset reporting
RFS and another reporting DFI are both usable under the DFS family; the
mixed tokens are listed). A gene with fewer than two usable datasets is
reported with `k = 1` and no pooled p-value, and the reason is recorded
in `errors`.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Two real lung cohorts with a DFS-family endpoint (RFS where available).
   p <- cpas_meta_panel(c("GSE31210", "GSE37745"), genes = c("TP53", "GAPDH"),
                         type = "DFS")
   p$table
   plot_meta_panel(p)
   ## the same panel with the earlier DerSimonian-Laird default, and with the
   ## normal-approximation prediction interval
   p_dl <- cpas_meta_panel(c("GSE31210", "GSE37745"), genes = "TP53",
                           type = "DFS", method = "DL", pi_method = "normal")
   p_dl$input$method; p_dl$input$pi_method
} # }
```
