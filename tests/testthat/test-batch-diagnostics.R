# batch_diagnostics() / plot_batch_diagnostics() (spec B5) -------------------
# Offline: cohorts are supplied through 'merged', so nothing is retrieved.

mk_bd_cohort <- function(n, shift = 0, seed = 1, na_at = integer(0)) {
  set.seed(seed)
  v <- stats::rnorm(n, 10 + shift, 1)
  v[na_at] <- NA_real_
  data.frame(ID = sprintf("S%04d", seq_len(n)),
             OS_time = round(stats::rexp(n, 0.15) + 0.1, 3),
             OS_status = stats::rbinom(n, 1, 0.5),
             GAPDH = v,
             TP53 = stats::rnorm(n, 8 + shift, 1.4),
             stringsAsFactors = FALSE)
}

test_that("batch_diagnostics summarises every cohort and gene", {
  m <- list(A = mk_bd_cohort(60, 0, seed = 1), B = mk_bd_cohort(50, 2.5, seed = 2))
  b <- batch_diagnostics(c("A", "B"), genes = c("GAPDH", "TP53"), type = "OS",
                         merged = m, top_pct = 25)
  expect_s3_class(b, "cpas_batch_diagnostics")
  expect_equal(nrow(b$table), 4L)
  expect_true(all(c("dataset", "gene", "endpoint", "n", "n_missing", "median",
                    "q1", "q3", "iqr", "min", "max", "mean", "sd",
                    "cut_median", "cut_top_pct", "n_high_median", "n_low_median",
                    "n_high_top_pct", "n_low_top_pct", "median_z") %in%
                    colnames(b$table)))
  # the summaries are the ones the values give
  for (i in seq_len(nrow(b$table))) {
    v <- m[[b$table$dataset[i]]][[b$table$gene[i]]]
    expect_equal(b$table$n[i], sum(!is.na(v)))
    expect_equal(b$table$median[i], stats::median(v, na.rm = TRUE))
    expect_equal(b$table$iqr[i], diff(stats::quantile(v, c(.25, .75), na.rm = TRUE,
                                                      names = FALSE)))
    expect_equal(b$table$min[i], min(v, na.rm = TRUE))
    expect_equal(b$table$max[i], max(v, na.rm = TRUE))
    # the two cut-points are exactly the rules cpas_km_pooled() uses
    expect_equal(b$table$cut_median[i], stats::median(v, na.rm = TRUE))
    expect_equal(b$table$cut_top_pct[i],
                 stats::quantile(v, 1 - 25 / 100, na.rm = TRUE, names = FALSE))
    expect_equal(b$table$n_high_median[i], sum(v > b$table$cut_median[i], na.rm = TRUE))
    expect_equal(b$table$n_low_median[i], sum(v <= b$table$cut_median[i], na.rm = TRUE))
  }
  # the cross-cohort spread is reported per gene
  expect_equal(nrow(b$batch), 2L)
  expect_true(all(b$batch$median_shift > 0))
  g <- b$batch[b$batch$gene == "GAPDH", ]
  expect_equal(g$min_median, min(b$table$median[b$table$gene == "GAPDH"]))
  expect_equal(g$max_median, max(b$table$median[b$table$gene == "GAPDH"]))
  expect_equal(g$median_shift, g$max_median - g$min_median)
})

test_that("the reported cut-points agree with cpas_km_pooled()", {
  m <- list(A = mk_bd_cohort(80, 0, seed = 11), B = mk_bd_cohort(70, 1.8, seed = 12))
  b <- batch_diagnostics(c("A", "B"), genes = "GAPDH", type = "OS", merged = m,
                         top_pct = 30)
  km <- cpas_km_pooled(m, marker = "GAPDH", type = "OS", method = "ipd",
                       cut = "median")
  expect_equal(unname(km$cohort_medians[as.character(b$table$dataset)]),
               b$table$cut_median, tolerance = 1e-12)
  km2 <- cpas_km_pooled(m, marker = "GAPDH", type = "OS", method = "ipd",
                        cut = "top_pct", top_pct = 30)
  expect_equal(unname(km2$cohort_thresholds[as.character(b$table$dataset)]),
               b$table$cut_top_pct, tolerance = 1e-12)
})

test_that("missing values are counted, not dropped silently", {
  m <- list(A = mk_bd_cohort(40, 0, seed = 3, na_at = 1:5),
            B = mk_bd_cohort(30, 1, seed = 4))
  b <- batch_diagnostics(c("A", "B"), genes = "GAPDH", merged = m)
  expect_equal(b$table$n_missing[b$table$dataset == "A"], 5L)
  expect_equal(b$table$n[b$table$dataset == "A"], 35L)
  expect_equal(b$table$n_missing[b$table$dataset == "B"], 0L)
  expect_equal(nrow(b$values), 65L)
})

test_that("cohorts and genes that cannot be summarised are reported, not fatal", {
  m <- list(A = mk_bd_cohort(40, 0, seed = 5), B = mk_bd_cohort(30, 1, seed = 6))
  b <- batch_diagnostics(c("A", "B", "MISSING"), genes = c("GAPDH", "NOPE"),
                         merged = m)
  expect_true("MISSING" %in% names(b$errors))
  expect_true("A/NOPE" %in% names(b$errors))
  expect_true(all(b$table$gene == "GAPDH"))
  expect_match(b$manifest$selection_rule, "summarised per cohort")
  expect_error(batch_diagnostics(character(0), "GAPDH", merged = m), "'datasets'")
  expect_error(batch_diagnostics("A", character(0), merged = m), "'genes'")
  expect_error(batch_diagnostics("A", "GAPDH", merged = unname(m)), "named list")
  expect_error(batch_diagnostics("A", "GAPDH", merged = m, top_pct = 0), "'top_pct'")
  expect_error(batch_diagnostics("A", "GAPDH", merged = list()), "named list")
})

test_that("the diagnostic plot draws one axis with the cut-points marked", {
  m <- list(A = mk_bd_cohort(60, 0, seed = 7), B = mk_bd_cohort(50, 2.5, seed = 8))
  b <- batch_diagnostics(c("A", "B"), genes = "GAPDH", merged = m)
  p <- plot_batch_diagnostics(b)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 110L)
  # one panel, one shared axis: the cohort factor keeps the order
  expect_identical(levels(p$data$dataset), c("A", "B"))
  p2 <- plot_batch_diagnostics(b, show_cut_points = FALSE)
  expect_s3_class(p2, "ggplot")
  # a second gene adds a facet
  b2 <- batch_diagnostics(c("A", "B"), genes = c("GAPDH", "TP53"), merged = m)
  p3 <- plot_batch_diagnostics(b2)
  expect_s3_class(p3, "ggplot")
  expect_error(plot_batch_diagnostics(list()), "cpas_batch_diagnostics")
  expect_error(plot_batch_diagnostics(b, genes = "NOPE"), "No value to plot")
  expect_output(print(b), "batch diagnostics")
})

test_that("the manifest states the within-cohort pooling and cut-point rules", {
  m <- list(A = mk_bd_cohort(40, 0, seed = 9))
  b <- batch_diagnostics("A", "GAPDH", merged = m)
  mm <- cpas_manifest(b)
  expect_s3_class(mm, "cpas_manifest")
  expect_match(mm$cut_rule, "within-cohort percentile cut-points")
  expect_match(mm$cut_rule, "no absolute")
  expect_true(any(grepl("within-cohort standardised effect sizes", mm$notes)))
  expect_true(any(grepl("own percentile", mm$notes)))
  expect_identical(mm$analysis, "batch_diagnostics")
})
