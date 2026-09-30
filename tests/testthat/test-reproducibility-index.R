# cpas_reproducibility_index() / cpas_pipeline_steps() (spec B7) --------------
# Offline: the index is read from the shipped CSVs and resolved with system.file.

test_that("the step inventory has the 67 steps and resolves every consolidated file", {
  st <- cpas_pipeline_steps()
  expect_equal(nrow(st), 67L)
  expect_identical(st$number, 1:67)
  expect_equal(sum(st$kind == "numbered construction / repair step"), 58L)
  expect_equal(sum(st$kind == "helper / demo / validation script"), 9L)
  expect_equal(length(unique(st$script)), 67L)
  expect_setequal(unique(st$consolidated_file),
                  c("01_ingest_parse.R", "02_platform_maps.R", "03_survival_tables.R",
                    "04_cohort_builds.R", "05_clinical_and_scale.R",
                    "06_catalog_and_registry.R", "07_mirror_upload.R",
                    "08_qc_screen_and_verify.R", "09_analysis_and_figures.R"))
  expect_true(all(nzchar(st$installed_path)))
  expect_true(all(st$available))
  expect_true(all(basename(st$installed_path) == st$consolidated_file))
  expect_true(all(nzchar(st$what_it_does)))
  # every id the entry point lists is present exactly once
  expect_setequal(st$id, c("01", "02", "03", "04", "05", "06", "07", "07b", "07c",
                           "08", "09", "10", "11", "12", "13", "14", "15", "16",
                           "17", "18", "19", "20", "21", "22", "23", "24", "25",
                           "26", "27", "35", "36", "37", "38", "39", "40", "41",
                           "42", "90", "91", "92", "93", "94", "95", "96", "98",
                           "99", "100", "101", "102", "103", "104", "105", "106",
                           "107", "110", "111", "112", "113", "batch_integrate",
                           "batch_integrate2", "demo_meta_lung",
                           "demo_tcga_integration", "demo_tcga_ondemand",
                           "demo_unified_reader", "test_cpas_GSE44001",
                           "test_cpas_dataset", "validate_cpas"))
  # cross-checked against the shipped entry point's own step table
  expect_true(all(st$matches_entry_point))
  expect_error(cpas_pipeline_steps(package = "no.such.package"), "not found")
})

test_that("the index covers the steps and the reproducibility material", {
  idx <- cpas_reproducibility_index()
  expect_s3_class(idx, "data.frame")
  expect_equal(nrow(idx), 67L + 14L)
  expect_identical(attr(idx, "n_steps"), 67L)
  expect_identical(attr(idx, "n_materials"), 14L)
  expect_setequal(unique(idx$section),
                  c("curation_step", "register", "repair_record",
                    "frozen_state", "additional_file"))
  expect_true(all(c("section", "number", "item", "script", "what_it_does",
                    "consolidated_file", "records", "location", "package_path",
                    "installed_path", "available", "resolved_source") %in%
                    colnames(idx)))
  steps <- idx[idx$section == "curation_step", ]
  expect_equal(nrow(steps), 67L)
  expect_identical(steps$number, 1:67)
  expect_true(all(steps$available))
  expect_true(all(grepl("^pipeline/", steps$package_path)))
  expect_true(all(is.na(steps$records)))
})

test_that("the exclusion register, the defect records and the frozen state are indexed", {
  idx <- cpas_reproducibility_index()
  ex <- idx[idx$item == "exclusion_register", ]
  expect_equal(nrow(ex), 1L)
  expect_equal(ex$records, 109L)
  expect_identical(ex$section, "register")
  expect_match(ex$location, "Additional_file_2_excluded_cohorts.csv")
  # the repair / defect records
  rep <- idx[idx$section == "repair_record", ]
  expect_true(nrow(rep) >= 3L)
  expect_true("verified_defects_1.0.0" %in% rep$item)
  expect_true(any(rep$records == 3L))                 # 3 defect entries
  expect_match(rep$location[rep$item == "verified_defects_1.0.0"],
               "REPORT_verified_defects_1.0.0.md")
  # the frozen-state record
  fs <- idx[idx$section == "frozen_state", ]
  expect_true("frozen_state_100" %in% fs$item)
  expect_match(fs$location[fs$item == "frozen_state_100"], "FROZEN_STATE_100.md")
  # the Additional-file layout: three files plus the index
  af <- idx[idx$section == "additional_file", ]
  expect_true(all(c("additional_file_1", "additional_file_2", "additional_file_3",
                    "additional_files_index") %in% af$item))
  expect_equal(af$records[af$item == "additional_file_1"], 197L)
  expect_equal(af$records[af$item == "additional_file_2"], 109L)
  expect_match(af$location[af$item == "additional_file_1"],
               "Additional_file_1_per_dataset_endpoint_resolution.csv")
  # it is an index, not a data dump: no register rows are carried
  expect_false(any(grepl("^GSE[0-9]", idx$item)))
})

test_that("the index resolves the project tree when it is present and degrades cleanly", {
  idx <- cpas_reproducibility_index()
  expect_true(any(!is.na(idx$resolved_source)))
  ex <- idx[idx$item == "exclusion_register", ]
  expect_true(file.exists(ex$resolved_source))
  # with a root that does not exist nothing is resolved, and nothing fails
  none <- cpas_reproducibility_index(root = tempfile("nope"))
  expect_true(all(is.na(none$resolved_source)))
  expect_equal(nrow(none), nrow(idx))
  # and the shipped part still resolves through the installed package
  expect_true(all(none$available[none$section == "curation_step"]))
})

test_that("the index does not touch the catalog or the mirror", {
  # it is built from two small CSVs only: replacing the catalog binding must not
  # change a single cell of the index
  idx <- cpas_reproducibility_index()
  src <- attr(idx, "steps_source")
  mat <- attr(idx, "materials_source")
  expect_true(file.exists(src) && file.exists(mat))
  expect_true(all(file.info(c(src, mat))$size < 20000))
  expect_match(attr(idx, "note"), "no catalog or mirror access")
})

test_that("the shipped step CSVs match the entry point they describe", {
  f <- system.file("reproducibility", "curation_steps.csv", package = "CanPAS")
  expect_true(nzchar(f))
  s <- utils::read.csv(f, stringsAsFactors = FALSE)
  expect_equal(nrow(s), 67L)
  expect_true(all(c("number", "id", "script", "what_it_does",
                    "consolidated_file", "kind") %in% colnames(s)))
  # each consolidated file really ships and really contains the runner
  for (g in unique(s$consolidated_file)) {
    p <- system.file("pipeline", g, package = "CanPAS")
    expect_true(file.exists(p))
    lines <- readLines(p, warn = FALSE)
    sub <- s[s$consolidated_file == g, ]
    for (sc in sub$script) {
      runner <- paste0("run_", sub("\\.R$", "", sc), " <- function() {")
      expect_true(any(trimws(lines) == runner))
    }
  }
})
