# cpas_meta_panel(): per-gene meta-analysis panel ---------------------------

mk_panel_data <- function(n = 140, seed = 1, genes = c("G1", "G2", "G3")) {
  set.seed(seed)
  x1 <- rnorm(n)                       # G1: harmful
  x2 <- rnorm(n)                       # G2: protective
  data.frame(ID = sprintf("S%04d", seq_len(n)),
             RFS_time = round(rexp(n, 0.2) * exp(-0.5 * x1 + 0.4 * x2), 3),
             RFS_status = rbinom(n, 1, 0.6),
             G1 = x1, G2 = x2, G3 = rnorm(n), stringsAsFactors = FALSE)
}

test_that("cpas_meta_panel pools every gene and adjusts the p-values", {
  m <- list(A = mk_panel_data(seed = 1), B = mk_panel_data(seed = 2),
            C = mk_panel_data(seed = 3))
  p <- cpas_meta_panel(datasets = c("A", "B", "C"),
                       genes = c("G1", "G2", "G3"), type = "RFS",
                       merged = m, progress = FALSE)
  expect_s3_class(p, "cpas_meta_panel")
  expect_equal(nrow(p$table), 3L)
  expect_true(all(c("gene", "k", "total_n", "total_events", "HR", "lower",
                    "upper", "p", "P_adj", "I2", "tau2", "pi_lower") %in%
                    colnames(p$table)))
  ok <- !is.na(p$table$p)
  expect_equal(p$table$P_adj[ok], stats::p.adjust(p$table$p[ok], "BH"),
               tolerance = 1e-12)
  expect_true(all(p$table$P_adj[ok] >= p$table$p[ok] - 1e-12))
  expect_true(all(p$table$k[ok] == 3L))          # all three datasets pooled
  expect_equal(sum(p$table$total_events), sum(p$per_dataset$events))
  # the harmful gene should come out above 1, the protective one below
  expect_gt(p$table$HR[p$table$gene == "G1"], 1)
  expect_lt(p$table$HR[p$table$gene == "G2"], 1)
  # per-dataset table and token map
  expect_equal(nrow(p$per_dataset), 9L)
  expect_equal(sort(unique(p$per_dataset$gene)), c("G1", "G2", "G3"))
  expect_true(all(p$per_dataset$endpoint == "RFS"))
  expect_equal(names(p$endpoints), c("G1", "G2", "G3"))
  expect_equal(p$metadata$genes_pooled, 3L)
  expect_output(print(p), "per-gene meta-analysis panel")
  expect_s3_class(plot_meta_panel(p), "ggplot")
})

test_that("cpas_meta_panel reports genes that cannot be pooled", {
  m <- list(A = mk_panel_data(seed = 1), B = mk_panel_data(seed = 2))
  m$A$G4 <- stats::rnorm(nrow(m$A))               # usable in A
  m$B$G4 <- 1                                     # constant in B -> B excluded
  p <- cpas_meta_panel(datasets = c("A", "B"), genes = c("G1", "G4"), type = "RFS",
                       merged = m, progress = FALSE)
  expect_equal(p$table$k[p$table$gene == "G1"], 2L)
  expect_equal(p$table$k[p$table$gene == "G4"], 1L)   # only A contributes
  expect_true(grepl("constant", p$errors[["G4"]], ignore.case = TRUE))
  expect_equal(p$metadata$genes_single_dataset, 1L)
  expect_true(is.na(p$table$P_adj[p$table$gene == "G4"]))
  # a gene that is absent everywhere is reported as failed, not dropped
  p2 <- cpas_meta_panel(datasets = c("A", "B"), genes = c("NOPE"), type = "RFS",
                        merged = m, progress = FALSE)
  expect_equal(p2$table$k, 0L)
  expect_true(grepl("gene missing", p2$errors[["NOPE"]], ignore.case = TRUE))
  expect_equal(p2$metadata$genes_failed, 1L)
})

test_that("cpas_meta_panel validates its inputs and honours fdr_method", {
  m <- list(A = mk_panel_data(seed = 1), B = mk_panel_data(seed = 2))
  expect_error(cpas_meta_panel(genes = "G1", merged = m), "'datasets'")
  expect_error(cpas_meta_panel(datasets = "A", merged = m), "'genes'")
  expect_error(cpas_meta_panel(datasets = "A", genes = "G1", merged = unname(m)),
               "named list")
  p <- cpas_meta_panel(datasets = c("A", "B"), genes = c("G1", "G2", "G3"),
                       type = "RFS", merged = m, fdr_method = "none",
                       progress = FALSE)
  expect_equal(p$table$P_adj, p$table$p)
  expect_error(cpas_meta_panel(datasets = c("A", "B"), genes = "G1",
                               type = "RFS", merged = m, method = "nope",
                               progress = FALSE), "arg")
  expect_error(plot_meta_panel(list(table = data.frame())), "cpas_meta_panel")
})

test_that("plot_meta_forest renders from a cpas_meta object", {
  m <- list(A = mk_panel_data(seed = 11), B = mk_panel_data(seed = 12),
            C = mk_panel_data(seed = 13))
  r <- suppressMessages(cpas_meta(datasets = c("A", "B", "C"), marker = "G1",
                                  type = "RFS", merged = m))
  expect_s3_class(plot_meta_forest(r), "ggplot")          # regression: bare `table`
  expect_s3_class(plot_meta_forest(r, digits = 4), "ggplot")
  p <- plot_meta_forest(r, digits = 4)
  # the per-dataset axis must carry the dataset labels, not the base table() function
  expect_true(all(r$per_dataset$dataset %in% levels(p$data$dataset)))
})
