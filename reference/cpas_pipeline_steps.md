# The 67 curation steps of the shipped pipeline

Returns the step inventory of the curation pipeline that built the
CanPAS mirror and catalog: one row per step with its canonical run
order, its original script name, what the step does, and which
consolidated file under `inst/pipeline/` provides it.

The inventory is a small CSV shipped with the package
(`inst/reproducibility/curation_steps.csv`); the consolidated file of
every step is resolved with
[`system.file`](https://rdrr.io/r/base/system.file.html), so the result
is usable from an installed package. When the shipped entry point
`inst/pipeline/run_pipeline.R` is present, each row is cross-checked
against its step table as well, so drift between the inventory and the
entry point is visible rather than silent.

## Usage

``` r
cpas_pipeline_steps(package = "CanPAS")
```

## Arguments

- package:

  Package to resolve the shipped files in (default `"CanPAS"`); mainly
  useful for tests.

## Value

A data frame with one row per step (67 rows):

- `number`:

  canonical run order, 1..67

- `id`:

  the step id used by `run_pipeline.R` (e.g. `"07b"`)

- `script`:

  the original script name (e.g. `"07b_split_tnm.R"`)

- `what_it_does`:

  one-line description of the step

- `kind`:

  `"numbered construction / repair step"` (58 rows) or
  `"helper / demo / validation script"` (9 rows)

- `consolidated_file`:

  the `inst/pipeline/` file that provides the step

- `installed_path`:

  the resolved absolute path of that file (`""` when the package does
  not ship it)

- `available`:

  whether `installed_path` exists

- `matches_entry_point`:

  `TRUE`/`FALSE` when the entry point could be read and the row agrees
  with its step table; `NA` when it could not be read

## See also

[`cpas_reproducibility_index`](https://wangjin93.github.io/CanPAS/reference/cpas_reproducibility_index.md)

## Examples

``` r
steps <- cpas_pipeline_steps()
nrow(steps)                                  # 67
#> [1] 67
table(steps$kind)
#> 
#>   helper / demo / validation script numbered construction / repair step 
#>                                   9                                  58 
head(steps[, c("number", "script", "consolidated_file")])
#>   number            script         consolidated_file
#> 1      1    01_parse_gse.R         01_ingest_parse.R
#> 2      2      02_gpl_map.R        02_platform_maps.R
#> 3      3   03_surv_table.R      03_survival_tables.R
#> 4      4    04_qc_report.R 08_qc_screen_and_verify.R
#> 5      5 05_dataset_plan.R 06_catalog_and_registry.R
#> 6      6    06_upload_db.R        07_mirror_upload.R
```
