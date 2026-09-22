test_that("COX_by_genes returns consistent HRs and S3 class", {
  d <- mk_cohort(300, seed = 1)
  r <- COX_by_genes(d, "OS", c("GAPDH", "TNS1", "PTEN"))
  expect_s3_class(r, "cpas_COX_by_genes")
  expect_true(all(c("gene", "HR", "HR95L", "HR95H", "Pvalue", "P_adj") %in%
                    colnames(r$results_table)))
  # HR equals exp(coef) of the fitted model
  b <- stats::coef(r$individual_models[["GAPDH"]])
  expect_equal(as.numeric(r$results_table$HR[r$results_table$gene == "GAPDH"]),
               as.numeric(exp(b)), tolerance = 1e-9)
})

test_that("COX_analysis univariate handles categorical covariates", {
  d <- mk_cohort(300, seed = 2)
  r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                    cate_Variates = c("sex", "stage"), method = "uni")
  expect_s3_class(r, "cpas_COX")
  expect_true("cate_header" %in% r$results_table$Type)
  expect_true(any(r$results_table$Type == "cate_level"))
})

test_that("COX_analysis multivariate matches survival::coxph", {
  d <- mk_cohort(300, seed = 3)
  r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                    cate_Variates = "sex", method = "multi")
  fit <- survival::coxph(survival::Surv(OS_time, OS_status) ~ GAPDH + age + sex, data = d)
  g <- r$results_table$HR[r$results_table$Type == "cont" &
                            r$results_table$Var1 == "GAPDH"]
  expect_equal(g, unname(exp(stats::coef(fit)["GAPDH"])), tolerance = 1e-9)
})

test_that("COX_screen_adjust degrades gracefully with no signal", {
  d <- mk_cohort(200, seed = 4)
  q <- COX_screen_adjust(d, type = "OS", cont_Variates = c("GAPDH", "age"),
                           cate_Variates = "sex", p.threshold = 1e-9)
  expect_true(is.data.frame(q$result))
  expect_null(q$multi_table)
})
