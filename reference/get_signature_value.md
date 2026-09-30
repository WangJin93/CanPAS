# Weighted Gene-Signature Score per Sample

Parses a linear signature string such as
`"0.3*GAPDH + 0.7*ACTB + 0.5*RPN1"` and computes the weighted sum of
expression values for every sample of a dataset.

## Usage

``` r
get_signature_value(
  signature,
  dataset,
  process_duplicates = "max",
  allow_missing = FALSE
)
```

## Arguments

- signature:

  Character string of the form `<weight>*<GENE> +/- <weight>*<GENE> ...`
  (weights may be negative).

- dataset:

  Dataset accession (GEO or `"TCGA-<PROJECT>"`).

- process_duplicates:

  Probe collapsing rule passed to
  [`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)
  (default `"max"`).

- allow_missing:

  If `FALSE` (default), all signature genes must be measurable on the
  platform; otherwise the score is computed on the subset of genes that
  are available (with a message).

## Value

List of class `cpas_signature`:

- `input_params`::

  echoed inputs

- `signature_info`::

  `original_signature`, `genes`, `weights` and `weight_gene_pairs`
  (aligned to available genes)

- `expr_data`::

  full `cpas_get_expr` object

- `results_table`::

  data.frame `ID`, `signature`

- `metadata`::

  bookkeeping

## Details

Weights are applied in the order of the signature string. Genes that
cannot be measured are reported; with `allow_missing = TRUE` the score
uses only measured genes (weights of missing genes dropped).

## Examples

``` r
if (FALSE) { # \dontrun{
   ## A signature is scored inside a real cohort: expression is fetched,
   ## probes are collapsed, and missing genes stop the call unless allowed.
   s <- get_signature_value("0.5*GAPDH + 0.5*ACTB", "GSE14814")
   s$metadata$genes_used; s$metadata$sample_count
   s$signature_info$weights         # the parsed weights, in order
   head(s$results_table)            # ID + one score per patient

   ## Genes that are not measurable in that cohort: stop, or score the subset.
   s2 <- get_signature_value("0.5*GAPDH + 0.5*ACTB", "GSE13507",
                              allow_missing = TRUE)
   s2$metadata$genes_used; s2$metadata$genes_missing
} # }
```
