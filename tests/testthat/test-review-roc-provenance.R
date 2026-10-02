test_that("ROC plot retains analysis counts and provenance", {
  d <- mk_cohort(120, seed = 813)
  d$GAPDH[1:3] <- NA
  g <- plot_roc(d, marker = "GAPDH", predict.time = 1)
  md <- attr(g, "analysis_metadata")
  expect_equal(md$n_input, 120L)
  expect_equal(md$n_dropped, 3L)
  expect_equal(md$n, 117L)
  expect_s3_class(attr(g, "manifest"), "cpas_manifest")
  expect_true(is.data.frame(attr(g, "roc_data")))
  expect_true(is.finite(attr(g, "AUC")))
})
