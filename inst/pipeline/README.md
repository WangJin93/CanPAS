# CanPAS curation pipeline (scripts only)

This directory ships the **data-curation pipeline** that built the CanPAS
mirror and the catalog that the package serves.  It is distributed **inside the
package** so that the curation is part of the released artefact, not a private
side repository.

**Scripts only — no data.** Nothing in this directory is a dataset: the scripts
read and write a *data root* that lives outside the package (see below).  The
only data files that ship with CanPAS are the small curated tables already under
`inst/extdata/tcga/` and `data/`, which are unrelated to this directory.

## What the pipeline does

A dataset enters the mirror through the same six stages the original pipeline
documented:

| stage | what happens | local products (under the data root) |
|---|---|---|
| 1 download + parse | GEO series matrix (or CGGA archive, or a publication table) is parsed | `data/raw/<ACC>.gz`, `data/expr/<ACC>.rds`, `data/pheno/<ACC>.rds`, and a row in `data/dataset_info.csv` |
| 2 platform annotation | platform SOFT -> probe -> ENTREZ map | `data/processed/gpl/<GPL>.rds` |
| 3 survival/clinical table | endpoint columns are chosen per cohort, units become years, status becomes 0/1, clinical columns are standardised | `data/processed/surv/<ACC>_surv.rds` |
| 4 QC | sample/event counts, time range, annotation coverage | `pipeline/out/QC_<date>.csv` |
| 5 clinical / scale harmonisation | clinical normalisation (stage/T/N/M/sex/age/grade) and log2(x+1) expression-scale harmonisation | in-place edits of the tables above, with backups under `pipeline/backup/` |
| 6 upload | expression, survival and platform tables are written to the MySQL mirror | remote tables, plus a log under `pipeline/out/` |

Beyond that main line the pipeline also performs cohort repair, catalog
assembly, mirror verification, candidate screening and the analysis/demo
scripts that exercised the finished package.

## Required data root

Every step reads its own root, exactly as before: either `CPAS_DATA_ROOT`,
or the path hard-coded in that step, or an explicit first argument.  The
project tree is expected to contain:

```
all datasets.csv          master cohort list
data/raw/                 downloaded series matrices
data/expr/  data/pheno/   parsed expression and clinical tables
data/processed/surv/      curated survival/clinical tables
data/processed/gpl/       platform annotation maps
data/dataset_info.csv     the catalog
CanPAS/data/dataset_info.rda, CanPAS/data/ID_map.rda
pipeline/out/             reports and logs written by the pipeline
```

Credentials are read from the environment (`CPAS_DB_PASSWORD`, optional
`CPAS_DB_USER` / `CPAS_DB_HOST`); see `Renviron.example`.

## Layout

The original pipeline was 67 standalone `pipeline/R/*.R` scripts that were run
one process per script.  They are consolidated here into nine files.  Each
original script is embedded **byte-for-byte** inside a runner function
`run_<script>()`, so sourcing a file defines functions only and has no side
effects, and each step still runs in its own process.

| consolidated file | original scripts merged in (in this order) |
|---|---|
| `01_ingest_parse.R` | `01_parse_gse.R`, `10_parse_cgga.R`, `90_build_gse108474_suppl.R`, `91_complete_gse14520_surv.R`, `92_build_geo_expansion_expr_pheno.R`, `93_build_geo_expansion_bespoke.R` |
| `02_platform_maps.R` | `02_gpl_map.R`, `12_gpl_map_symbol.R`, `38_geo_expansion_annotation_check.R`, `102_extend_gpl4133_agilent_name.R`, `103_build_embl_gpl_maps.R` |
| `03_survival_tables.R` | `03_surv_table.R`, `14_split_gse40272_surv.R`, `15_split_gse40272_expr.R`, `17_recod_gse86166_status.R`, `23_fix_surv_ids_and_endpoints.R`, `25_fix_endpoint_annotations.R`, `35_add_tcga_cohorts.R`, `39_build_gse1379_surv.R`, `40_add_tcga_chol_dlbc.R`, `110_build_gse205209.R` |
| `04_cohort_builds.R` | `95_build_suppl_expansion.R`, `96_build_suppl_expansion.R`, `100_build_relaxed_gate.R`, `104_build_embl_cohorts.R`, `106_build_embl_step2.R` |
| `05_clinical_and_scale.R` | `07_standardize_clinical.R`, `07b_split_tnm.R`, `07c_T_usage_table.R`, `08_verify_normalization.R`, `11_endpoint_families.R`, `18_log2_transform.R`, `20_log2_transform_db.R`, `21_reupload_log2_tables.R`, `22_report_log2_impact.R` |
| `06_catalog_and_registry.R` | `05_dataset_plan.R`, `19_fix_catalog_N.R`, `26_cohort_overlap.R`, `27_catalog_notes.R`, `41_register_final3_rows.R`, `94_update_catalog_geo_expansion.R`, `98_update_catalog_suppl_expansion.R`, `101_update_catalog_relaxed_gate.R`, `105_update_catalog_embl.R`, `107_register_embl_step2.R`, `111_register_gse205209.R` |
| `07_mirror_upload.R` | `06_upload_db.R`, `09_update_db_surv.R`, `13_upload_small_cohorts.R`, `24_make_tables_writable.R`, `42_upload_gpl1223.R`, `99_extend_gpl_db.R`, `112_upload_gpl27956.R` |
| `08_qc_screen_and_verify.R` | `04_qc_report.R`, `16_verify_catalog_mirror.R`, `36_screen_geo_candidates.R`, `37_geo_expansion_screen.R`, `113_sweep_gpl_mirror_divergence.R` |
| `09_analysis_and_figures.R` | `batch_integrate.R`, `batch_integrate2.R`, `demo_meta_lung.R`, `demo_tcga_integration.R`, `demo_tcga_ondemand.R`, `demo_unified_reader.R`, `test_cpas_dataset.R`, `test_cpas_GSE44001.R`, `validate_cpas.R` |

`CPAS_full_test_report.Rmd` is the pipeline's own end-to-end test report; it is
carried verbatim and is not merged into any of the nine files.

## Step order

`run_pipeline.R` holds the canonical order (the same order as the original
entry point, extended to the later steps) and runs one step per process:

```sh
export CPAS_DATA_ROOT=/path/to/CPAS
Rscript run_pipeline.R --list                  # every step, with the group file
Rscript run_pipeline.R --check                 # consistency check only (step 16)
Rscript run_pipeline.R --from 07 --to 09       # a range
Rscript run_pipeline.R --only 11,16            # selected steps
Rscript run_pipeline.R --dry-run               # print the commands, run nothing
```

Without `--only`/`--from`/`--to` **every** listed step is selected; use
`--dry-run` first.  Steps that write to the MySQL mirror refuse to run unless
`CPAS_DB_PASSWORD` is set.

An individual step can also be run directly from R:

```r
source(system.file("pipeline", "05_clinical_and_scale.R", package = "CanPAS"))
run_07b_split_tnm()          # equivalent to: Rscript pipeline/R/07b_split_tnm.R
```

## Notes on the consolidation

* Every function definition of the original 67 scripts is present with an
  **identical deparsed body** (270 definitions, 0 altered), and the set of
  original function names is a subset of the merged set.  The whole text of 66
  of the 67 scripts appears byte-for-byte in the consolidated files.
* Names that several scripts defined differently (`note`, `num`, `%||%`,
  `parseGSEMatrix`, ...) are no longer ambiguous: each embedded script keeps its
  own scope, exactly as when it ran as its own `Rscript`.
* Two small pieces of glue are added, both documented in the files themselves:
  `05_clinical_and_scale.R` exposes the clinical normalizers
  (`MISS_TOKENS`, `is_missing_token`, `strip_missing`, `clean_text`,
  `norm_tnm`) at file level, because the original `07b_split_tnm.R` re-read and
  `eval()`-ed them out of `07_standardize_clinical.R` on disk; and
  `04_cohort_builds.R` exposes `scale_verdict()` at file level, because
  `95_build_suppl_expansion.R` used the copy defined by
  `96_build_suppl_expansion.R`.  Both copies are verbatim.
* `data/`, the mirror, the catalog and the shipped TCGA tables are untouched by
  this directory.
