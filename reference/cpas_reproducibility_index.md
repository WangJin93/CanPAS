# Index of the reproducibility material

A tidy data frame indexing the material behind the reported numbers: the
67 curation steps of the shipped pipeline (with the consolidated
`inst/pipeline/` file that provides each one), the exclusion register
(109 records), the repair and defect records, the frozen-state record,
and the Additional-file layout.

The index is assembled from two small CSVs that ship with the package
(`inst/reproducibility/curation_steps.csv` and
`inst/reproducibility/materials.csv`) and resolved with
[`system.file`](https://rdrr.io/r/base/system.file.html), so it works
offline from an installed package. It is an index, not a data dump: it
carries counts and locations, never the registers' rows, and it never
reads the catalog or the mirror.

The registers, defect reports and additional files themselves live in
the project tree (the reproducibility bundle), not in the package; their
`location` is reported relative to the project root (`CPAS_DATA_ROOT`,
default `/home/Jingle/data/Project/CPAS`) and is resolved to an absolute
existing path in `resolved_source` when that root is present on the
machine.

## Usage

``` r
cpas_reproducibility_index(package = "CanPAS", root = NULL)
```

## Arguments

- package:

  Package to resolve the shipped files in (default `"CanPAS"`); mainly
  useful for tests.

- root:

  Optional project root used to resolve the external material (default
  `CPAS_DATA_ROOT`, else the documented project path).

## Value

A data frame, one row per indexed item, with columns:

- `section`:

  `"curation_step"`, `"register"`, `"repair_record"`, `"frozen_state"`
  or `"additional_file"`

- `number`:

  1..67 for the curation steps, `NA` otherwise

- `item`:

  step id, or the item name for the other sections

- `script`:

  the original script of a step (`NA` otherwise)

- `what_it_does`:

  what the step or item is

- `consolidated_file`:

  the `inst/pipeline/` file providing a step (`NA` otherwise)

- `records`:

  the number of records the item carries where that number is part of
  the record (109 for the exclusion register, 197 for Additional file 1,
  ...); `NA` when the item is not a record list

- `location`:

  where the item lives, relative to the project root

- `package_path`:

  the path inside the installed package, when the item ships with the
  package (`NA` otherwise)

- `installed_path`, `available`:

  the resolved
  [`system.file()`](https://rdrr.io/r/base/system.file.html) path and
  whether it exists

- `resolved_source`:

  the absolute path of `location` under the project root when it exists
  there, else `NA`

## See also

[`cpas_pipeline_steps`](https://wangjin93.github.io/CanPAS/reference/cpas_pipeline_steps.md),
[`endpoint_semantics`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md),
[`cohort_overlap`](https://wangjin93.github.io/CanPAS/reference/cohort_overlap.md),
[`cpas_manifest`](https://wangjin93.github.io/CanPAS/reference/cpas_manifest.md)

## Examples

``` r
idx <- cpas_reproducibility_index()
table(idx$section)
#> 
#> additional_file   curation_step    frozen_state        register   repair_record 
#>               4              67               2               4               4 
subset(idx, section == "register", select = c("item", "records", "location"))
#>                        item records
#> 68       exclusion_register     109
#> 69 cohort_curation_register      43
#> 70    executed_run_register      NA
#> 71      rebaseline_register      NA
#>                                                                                                                        location
#> 68                                          run/CanPAS-tool-paper/final/additional-files/Additional_file_2_excluded_cohorts.csv
#> 69    run/CanPAS-tool-paper/final/additional-files/Additional_file_3_reproducibility_bundle/COHORT_CURATION_REGISTER_D01-D43.md
#> 70 run/CanPAS-tool-paper/final/additional-files/Additional_file_3_reproducibility_bundle/EXECUTED_RUN_REGISTER_original_pass.md
#> 71                 run/CanPAS-tool-paper/final/additional-files/Additional_file_3_reproducibility_bundle/REBASELINE_REGISTER.md
subset(idx, section == "curation_step")[1:3, c("number", "script", "consolidated_file")]
#>   number          script    consolidated_file
#> 1      1  01_parse_gse.R    01_ingest_parse.R
#> 2      2    02_gpl_map.R   02_platform_maps.R
#> 3      3 03_surv_table.R 03_survival_tables.R
```
