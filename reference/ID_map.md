# Gene-symbol to Entrez mapping (internal reference)

Mapping of gene symbols to Entrez gene ids (and Ensembl gene ids when
available) used to resolve probes across GEO platforms. Loaded lazily by
the package; there is no need to attach it manually.

## Usage

``` r
ID_map
```

## Format

A data.frame with columns `Symbol`, `Ensembl_gene` and `gene_id`.
