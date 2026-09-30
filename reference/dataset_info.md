# CanPAS dataset catalog

One row per cohort in the CanPAS catalog - 197 rows: 147 GEO series, 14
EMBL-EBI cohorts (ArrayExpress/BioStudies), 31 TCGA projects, 3 CGGA
glioma cohorts and 2 cBioPortal-hosted studies (`A5-PCPG`,
`IMmotion150`) - with the accession, platform (GPL), cancer type and the
endpoint / sample-size columns described below. Used by
[`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)
to find the platform of a dataset.

## Usage

``` r
dataset_info
```

## Format

A data.frame with 32 columns. `Accession`, `Type`, `GPL`, `N`,
`SurvivalTypes` (raw endpoint tokens) and the endpoint-family mapping
added by the analysis plan B: `EndpointFamilies` (browsing families),
`EP_OS`, `EP_DSS`, `EP_DFS`, `EP_PFS`, `EP_MFS` (the concrete token
available for each pooling family), `EndpointPrimary` and
`EndpointDerived` (TCGA-derived DFI/PFI).

Sample-size columns. The delivered local artefacts under `data/` are the
authoritative layer; the mirror and this packaged object are downstream
copies of it, checked against it by
`pipeline/R/16_verify_catalog_mirror.R`. `n_surv` is the row count of
the delivered survival table and `n_events` the primary endpoint's
events among those `N` analysable samples, so it is not the event total
of the survival table:

- `N`:

  Analysable sample count: samples with both expression data in the
  mirror and non-missing `time` and `status` for the dataset's
  `EndpointPrimary`. `NA` when the dataset has no expression table or no
  primary endpoint annotation, i.e. gene-level analysis is not possible.

- `n_expr`, `n_surv`:

  Samples in the mirror expression table and rows in the mirror survival
  table (`NA` when absent).

- `n_events`:

  Events (`status = 1`) for `EndpointPrimary`, counted inside the
  analysable set.

- `n_OS`, `n_DSS`, `n_DFS`, `n_PFS`, `n_MFS`:

  Analysable sample count per pooling family, resolved to the token each
  dataset actually provides (GSE31210 contributes to `n_DFS` through
  RFS). `0` when the family is unavailable.

- `expr_in_mirror`:

  Whether an expression table exists in the mirror.

Two columns describe how a cohort may be used:

- `CohortGroup`:

  Label of the group of cohorts that share patients (same study on
  another platform, or the same series deposited twice); 33 of the 197
  rows carry a label and the rest are `NA`. Computed by
  `pipeline/R/26_cohort_overlap.R` from sample titles.

- `Note`:

  Free-text provenance and use flag for the cohort, written by the
  curation pipeline. 120 of the 197 rows carry no note and are stored as
  the literal `NA`; the other 77 are ASCII prose. A note records
  whichever of the following applies, in the pipeline's own wording: the
  platform and the probe-to-gene coverage of the expression table
  (`GPL...`; `probes->genes: ... n/m non-control = p%`); the sample and
  event accounting behind `N` and `n_events`, with paired or replicate
  designs marked as such because they double-count patients; how the
  endpoint was derived (time unit, status coding, audit or patient-level
  check); an overlap group (`OVERLAPS ...`,
  `SAME SERIES on another platform - do not pool together ...`, or
  `TITLE COLLISION (not a shared-patient group)`); the registration
  batch date (`... batch 2026-09-24`); and `[!]` where the note raises a
  caveat the user should read. The two rows whose `n_convention` is
  `clinical-record` (`GSE325123`, `GSE31312`) additionally state the
  expression-join-restricted `N`/`n_events` pair in a
  `[n_convention=clinical-record: ...]` clause, and
  `pipeline/R/16_verify_catalog_mirror.R` treats those two rows as
  documented convention exceptions, not rule violations. The App shows
  the note on the Datasets page and warns on the multi-dataset pages
  when an overlapping pair is selected.

Two bookkeeping columns are not used by the analysis functions:

- `X`:

  Row index carried over from the pre-removal catalog (values 1-213, 197
  distinct). Harmless residue; kept so the packaged table stays
  cell-identical to `data/dataset_info.csv`.

- `method`:

  Assay / data type recorded for the cohort: `"RNA"` (136),
  `"TCGA-RNAseq"` (33), `"RNA-seq"` (9), `"array"` (1), `"SRA"` (1), and
  `NA` for the 17 lung-cancer GEO cohorts whose method was not recorded.

Earlier releases recorded in `N` the planned/expression cohort size,
which overstated the analysable sample count for many datasets (GSE31210
was listed as 133 while 226 tumours are analysable).

## Examples

``` r
  ## The catalog is the starting point of every real analysis: it records the
  ## cohorts mirrored in CanPAS and, per cohort, which endpoint families can be
  ## analysed and how many samples are usable.
  head(dataset_info[, c("Accession", "Type", "N", "EndpointFamilies")])
#>   Accession        Type   N EndpointFamilies
#> 1  GSE14814 Lung Cancer 133           OS,DSS
#> 2   GSE8894 Lung Cancer 138              DFS
#> 3  GSE31210 Lung Cancer 226           OS,DFS
#> 4  GSE13213 Lung Cancer 117               OS
#> 5  GSE37745 Lung Cancer 196           OS,DFS
#> 6  GSE17710 Lung Cancer  56           OS,DFS

  ## Cohorts of one cancer type with a DFS-family endpoint, largest first.
  lung <- dataset_info[dataset_info$Type == "Lung Cancer" &
                        grepl("DFS", dataset_info$EndpointFamilies), ]
  lung[order(-lung$N), c("Accession", "N", "EndpointFamilies")]
#>     Accession   N EndpointFamilies
#> 130 TCGA-LUAD 565   OS,DSS,DFS,PFS
#> 131 TCGA-LUSC 542   OS,DSS,DFS,PFS
#> 8    GSE30219 293           OS,DFS
#> 10   GSE41271 274           OS,DFS
#> 3    GSE31210 226           OS,DFS
#> 5    GSE37745 196           OS,DFS
#> 14   GSE50081 181           OS,DFS
#> 2     GSE8894 138              DFS
#> 13   GSE74777 107           OS,DFS
#> 6    GSE17710  56           OS,DFS

  ## The catalog drives the analysis helpers, e.g. the endpoint token.
  endpoint_resolve(lung$Accession[1], "DFS")
#> [1] "RFS"
```
