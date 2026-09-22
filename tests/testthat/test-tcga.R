# TCGA accessions in the GEO-oriented readers ---------------------------------
# The mirror carries no TCGA expression table, so get_expr_data() fetches from
# UCSC Xena and merge_surv_expr() reads the local clinical table. These tests
# need the network (Xena) and the optional UCSCXenaShiny package; they skip
# otherwise, like test-signature.R.

xena_ready <- function() {
  if (!requireNamespace("UCSCXenaShiny", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({ CanPAS::tcga_gene_expr_df("LUAD", "TP53"); TRUE },
                  error = function(e) FALSE))
}

test_that("get_expr_data() routes TCGA accessions to Xena", {
  skip_on_cran()
  skip_if_not(xena_ready(), "UCSC Xena / UCSCXenaShiny unavailable")
  e <- get_expr_data("TCGA-LUAD", "TP53")
  expect_s3_class(e, "cpas_get_expr")
  expect_equal(e$metadata$source, "xena")
  expect_equal(e$metadata$dataset_accession, "TCGA-LUAD")
  expect_true(all(c("ID", "TP53") %in% colnames(e$expr_data)))
  expect_gt(nrow(e$expr_data), 100L)
  expect_true(all(is.finite(e$expr_data$TP53)))
  # Xena HiSeqV2 is log2(x + 1): values are in the single/low double digits
  expect_true(max(e$expr_data$TP53) < 25)
  # no probe collapsing for TCGA: one ref_ids row per gene
  expect_equal(nrow(e$ref_ids), 1L)
  expect_equal(e$ref_ids$Symbol, "TP53")
})

test_that("get_expr_data() skips genes absent from Xena and keeps the rest", {
  skip_on_cran()
  skip_if_not(xena_ready(), "UCSC Xena / UCSCXenaShiny unavailable")
  e <- suppressMessages(get_expr_data("TCGA-LUAD", c("TP53", "NOTAGENE123")))
  expect_equal(colnames(e$expr_data), c("ID", "TP53"))
})

test_that("merge_surv_expr() uses the local TCGA clinical table", {
  skip_on_cran()
  skip_if_not(xena_ready(), "UCSC Xena / UCSCXenaShiny unavailable")
  e <- get_expr_data("TCGA-LUAD", "TP53")
  m <- merge_surv_expr("TCGA-LUAD", e)
  expect_s3_class(m, "cpas_merge")
  expect_true(all(c("ID", "OS_status", "OS_time", "TP53") %in% colnames(m$merged_data)))
  expect_gt(nrow(m$merged_data), 100L)
  expect_true(all(m$merged_data$OS_status %in% c(0, 1)))
  # status columns must be 0/1 where present (NA = endpoint not available for
  # that sample, which is legitimate: TCGA DSS/DFI cover fewer patients)
  st <- grep("_status$", colnames(m$merged_data), value = TRUE)
  expect_true(all(vapply(st, function(cc) {
    v <- m$merged_data[[cc]]
    all(v[!is.na(v)] %in% c(0, 1))
  }, logical(1))))
})

test_that("TCGA and GEO cohorts can be analysed side by side", {
  skip_on_cran()
  skip_if_not(xena_ready(), "UCSC Xena / UCSCXenaShiny unavailable")
  online <- isTRUE(tryCatch({ CanPAS::get_data("GSE31210", "surv_data"); TRUE },
                            error = function(e) FALSE))
  skip_if_not(online, "CanPAS API unreachable")
  r <- COX_by_datasets(datasets = c("TCGA-LUAD", "GSE31210"), gene = "TP53", type = "OS")
  expect_s3_class(r, "cpas_COX_by_datasets")
  expect_equal(nrow(r$combined_results), 2L)
  expect_length(r$errors, 0L)
  expect_true(all(is.finite(r$combined_results$HR)))
  expect_equal(r$combined_results$endpoint, c("OS", "OS"))
})

test_that("unknown TCGA accessions still fail fast", {
  expect_error(get_expr_data("TCGA-XXXX", "TP53"), "not found in dataset_info")
})
