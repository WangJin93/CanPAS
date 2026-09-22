test_that("cpas_meta pools per-cohort per-SD log-HRs (RE)", {
  a <- mk_cohort(200, "A", seed = 5)
  b <- mk_cohort(150, "B", seed = 6)
  r <- cpas_meta(datasets = c("A", "B"), marker = "GAPDH", type = "OS",
                 method = "RE", merged = list(A = a, B = b))
  expect_equal(nrow(r$per_dataset), 2)
  expect_true(is.finite(r$pooled$HR) && r$pooled$HR > 0)
  # pooled estimate lies between the two per-cohort estimates
  expect_gte(r$pooled$HR, min(r$per_dataset$HR))
  expect_lte(r$pooled$HR, max(r$per_dataset$HR))
})

test_that("cpas_meta multivariate honours complete cases", {
  a <- mk_cohort(200, "A", seed = 7)
  a$age[1:10] <- NA
  r <- cpas_meta(datasets = "A", marker = "GAPDH", type = "OS", method = "FE",
                 confounders = c("age", "sex", "stage"), merged = list(A = a))
  expect_equal(r$per_dataset$n, sum(complete.cases(a[c("OS_time", "OS_status",
                                                      "GAPDH", "age", "sex", "stage")])))
})

test_that("loo_meta removes one cohort at a time", {
  a <- mk_cohort(180, "A", seed = 8)
  b <- mk_cohort(180, "B", seed = 9)
  c <- mk_cohort(180, "C", seed = 10)
  r <- cpas_meta(c("A", "B", "C"), "GAPDH", "OS", "FE", merged = list(A = a, B = b, C = c))
  loo <- loo_meta(r)
  expect_equal(nrow(loo), 3)
  expect_setequal(loo$left_out, c("A", "B", "C"))
})

test_that("cpas_km_pooled runs ipd + meta landmarks", {
  a <- mk_cohort(200, "A", seed = 11)
  b <- mk_cohort(150, "B", seed = 12)
  km <- cpas_km_pooled(list(A = a, B = b), marker = "GAPDH", type = "OS",
                       method = "both", landmarks = c(1, 3))
  expect_s3_class(km$fit, "survfit")
  expect_true(is.finite(km$logrank_p))
  expect_true(is.finite(km$logrank_p_stratified))
  expect_true(is.data.frame(km$meta_landmarks))
  expect_true(nrow(km$meta_landmarks) >= 4)   # 2 groups x 2 landmarks
  expect_true(all(c("High", "Low") %in% levels(km$df$group)))
})

test_that("forest plot builds from cox results table", {
  d <- mk_cohort(300, seed = 13)
  r <- COX_analysis(d, type = "OS", cont_Variates = "GAPDH",
                    cate_Variates = "stage", method = "uni")
  p <- forest_plot(r$results_table)
  expect_s3_class(p, "ggplot")
})
