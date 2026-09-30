# Fetch one gene's TCGA expression from UCSC Xena

Downloads expression of a single gene for a TCGA project on demand
(per-gene query, no full-matrix download) using the internal
`UCSCXenaShiny:::get_data` helper, and returns a long data.frame of
tumour samples only.

## Usage

``` r
tcga_gene_expr_df(dataset, gene, max_try = 3, use_cache = NULL)
```

## Arguments

- dataset:

  TCGA dataset id (`"LUAD"` or `"TCGA-LUAD"`).

- gene:

  Single gene symbol.

- max_try:

  Number of attempts before giving up (Xena may be flaky).

- use_cache:

  Reuse the local file cache. `NULL` (default) follows
  `getOption("CanPAS.cache", TRUE)` and `CANPAS_CACHE` / `CPAS_CACHE`;
  `TRUE` or `FALSE` forces the choice. Xena is queried one gene at a
  time and is intermittently unreachable, so a cached gene table is
  reused after `CanPAS.cache_ttl` seconds only when the live request
  fails.

## Value

data.frame with columns `sample`, `value`.

## Examples

``` r
if (FALSE) { # \dontrun{
   x <- tcga_gene_expr_df("LUAD", "TP53")   # long table: sample, value
   head(x)

   ## The project id and the Xena dataset id it maps to, plus which projects
   ## CanPAS retains.
   tcga_project_dataset("TCGA-LUAD")
   names(tcga_retained)
} # }
```
