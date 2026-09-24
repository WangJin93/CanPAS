# Catalog <-> TCGA helper agreement -------------------------------------------
# Regression test for B1 of the 1.0.0 audit (2026-09-24): `tcga_retained` was a
# hard-coded 15-project vector while the catalog carried 31 TCGA cohorts, so 16
# of them failed in every TCGA helper (tcga_project_dataset / tcga_surv_table /
# tcga_merged / tcga_get_expr / get_expr_data) and canonical_type() returned NA
# for them. The project table is now derived from the shipped catalog, and this
# test runs under R CMD check - unlike the network tests in test-tcga.R, which
# are skip_on_cran() - so the next catalog expansion cannot silently break the
# helpers again. It needs no network, no credentials and no local data root.

test_that("every catalogued TCGA cohort is supported by the TCGA helpers", {
  di   <- CanPAS::dataset_info
  acc  <- as.character(di$Accession)
  tcga <- acc[grepl("^TCGA-", acc)]
  expect_gt(length(tcga), 0L)
  proj <- toupper(sub("^TCGA-", "", tcga))

  # the exported project table agrees with the catalog, cell by cell
  expect_setequal(names(tcga_retained), proj)
  expect_identical(unname(tcga_retained[proj]),
                   as.character(di$Type[match(tcga, acc)]))

  # and every one of them is accepted by the helpers (no network, no data root)
  for (i in seq_along(proj)) {
    p <- proj[i]
    expect_identical(tcga_project_dataset(p),
                     paste0("TCGA.", p, ".sampleMap/HiSeqV2"))
    expect_identical(tcga_project_dataset(tcga[i]),
                     paste0("TCGA.", p, ".sampleMap/HiSeqV2"))
    expect_false(is.na(canonical_type(tcga[i])))
    expect_identical(canonical_type(tcga[i]), as.character(di$Type[match(tcga[i], acc)]))
  }
})

test_that("the built-in fallback TCGA list matches the catalog", {
  # tcga_retained is refreshed from the catalog when the package is loaded; the
  # built-in vector is only the fallback for an unreadable catalog, so it must
  # be kept in step too. This is the assertion that fails on catalog drift.
  di <- CanPAS::dataset_info
  tc <- di[grepl("^TCGA-", as.character(di$Accession)), , drop = FALSE]
  fb <- CanPAS:::.cpas_tcga_retained_builtin
  expect_identical(names(fb), toupper(sub("^TCGA-", "", as.character(tc$Accession))))
  expect_identical(unname(fb), as.character(tc$Type))
})

test_that("the TCGA helpers still reject unknown projects", {
  expect_error(tcga_project_dataset("XXXX"), "Unsupported TCGA project")
  expect_true(is.na(canonical_type("TCGA-XXXX")))
  expect_true(is.na(canonical_type("GSE13507")))   # GEO accession, not a TCGA id
})
