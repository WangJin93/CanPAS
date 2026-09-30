# Retrieve Expression Data for Genes of a Dataset

Resolves gene symbols to probes on the dataset's platform (via the
ID_map table and the platform's probe-to-gene mapping) and downloads the
matching expression values from the CanPAS API. Multiple probes mapping
to the same gene can be collapsed with a user-selected rule.

Datasets whose accession starts with `"TCGA-"` are handled differently:
the mirror holds no TCGA expression table, so values are fetched per
gene from UCSC Xena on demand
([`tcga_gene_expr_df`](https://wangjin93.github.io/CanPAS/reference/tcga_gene_expr_df.md),
the same source the App uses). Xena returns gene-level log2 values, so
there is no probe resolution or collapsing and `process_duplicates` is
ignored.

## Usage

``` r
get_expr_data(
  dataset,
  gene = NULL,
  process_duplicates = c("max", "mean", "median", "min", "no"),
  genes = NULL
)
```

## Arguments

- dataset:

  Dataset accession in `dataset_info`, e.g. `"GSE14814"` or
  `"TCGA-LUAD"`.

- gene:

  Character vector of gene symbols.

- process_duplicates:

  How to collapse several probes that map to the same gene: `"max"`
  (default), `"mean"`, `"median"`, `"min"` - or `"no"` to keep
  probe-level rows (columns then carry probe ids instead of symbols).
  Ignored for TCGA datasets.

- genes:

  One or more gene symbols (character); equivalent to `gene` (supply
  only one of the two).

## Value

List of class `cpas_get_expr`:

- `input_params`::

  echoed inputs and analysis time

- `raw_ids`::

  rows of ID_map matching the requested symbols

- `platform_info`::

  dataset_info row of the dataset

- `ref_ids`::

  probe (row_names) to gene mapping actually used; for TCGA one row per
  gene, since there are no probes

- `expr_data`::

  data.frame, first column `ID` (samples), remaining columns are genes
  (or probes) with numeric expression

- `metadata`::

  genes requested / found / matched / sample count, plus `source`
  (`"mirror"` or `"xena"`)

## Details

The internal tables `ID_map` and `dataset_info` ship with the package
and are loaded lazily; they never need to be attached manually. Genes
that exist in `ID_map` but have no probe on the dataset's platform are
skipped with a message. For TCGA datasets the symbol is passed to Xena
directly, so a gene missing from `ID_map` is still attempted, and genes
absent from Xena are skipped with a message.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## GEO: probe ids are resolved through the platform map and collapsed to
   ## one value per gene.
   e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
   e$metadata
   head(e$expr_data)

   ## TCGA: the mirror holds no expression table, so values are fetched per
   ## gene from UCSC Xena (gene-level log2, no probe collapsing).
   tcga <- get_expr_data("TCGA-LUAD", "TP53")
   head(tcga$expr_data)
} # }
```
