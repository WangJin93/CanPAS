# Xena dataset id of a TCGA project

Builds the UCSC Xena sampleMap/HiSeqV2 dataset id for a TCGA project,
e.g. `"TCGA.LUAD.sampleMap/HiSeqV2"`.

## Usage

``` r
tcga_project_dataset(dataset)
```

## Arguments

- dataset:

  TCGA dataset id, either the project abbreviation (e.g. `"LUAD"`) or
  the accession (`"TCGA-LUAD"`).

## Value

A single string.
