# Focused offline regression tests for competing-risk model reporting.
review_competing_data <- function(censored = TRUE) {
  set.seed(202603)
  n <- 400L
  x <- stats::rnorm(n)
  interest <- stats::rexp(n, rate = 0.15 * exp(0.4 * x))
  competing <- stats::rexp(n, rate = 0.12 * exp(-0.2 * x))
  censor <- if (censored) stats::rexp(n, rate = 0.06) else rep(Inf, n)
  time <- pmin(interest, competing, censor)
  status <- ifelse(censor == time, 0L, ifelse(interest == time, 2L, 5L))
  data.frame(time = time, status = status, x = x)
}

test_that("separated models withhold all inference and keep debugging fits", {
  d <- review_competing_data()
  d$separated <- as.numeric(d$status == 2L)
  warnings <- character(0)
  result <- withCallingHandlers(
    competing_risk_COX(d, "time", "status", c("x", "separated"), etype = 2),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("cause-specific model is not fully estimable", warnings)))
  expect_true(any(grepl("subdistribution model is not fully estimable", warnings)))
  for (model in c("cause_specific", "subdistribution")) {
    expect_false(result$diagnostics[[model]]$ok)
    expect_true(nzchar(result$diagnostics[[model]]$reason))
    expect_null(result[[model]])
    expect_s3_class(result$fits[[model]], "coxph")
  }
})

test_that("singular models do not publish the remaining finite coefficients", {
  d <- review_competing_data()
  d$duplicate <- d$x
  result <- suppressWarnings(competing_risk_COX(
    d, "time", "status", c("x", "duplicate"), etype = 2
  ))
  for (model in c("cause_specific", "subdistribution")) {
    expect_false(result$diagnostics[[model]]$ok)
    expect_null(result[[model]])
    expect_s3_class(result$fits[[model]], "coxph")
  }
})

test_that("noncontiguous event labels map correctly with and without censoring", {
  for (censored in c(TRUE, FALSE)) {
    d <- review_competing_data(censored)
    recoded <- d
    recoded$status <- ifelse(d$status == 2L, 1L, ifelse(d$status == 5L, 2L, 0L))
    for (etype in c(2L, 5L)) {
      result <- competing_risk_COX(d, "time", "status", "x", etype = etype)
      reference <- competing_risk_COX(
        recoded, "time", "status", "x", etype = if (etype == 2L) 1L else 2L
      )
      expect_equal(result$events_interest, sum(d$status == etype))
      expect_equal(result$events_competing, sum(!d$status %in% c(0L, etype)))
      expect_equal(result$n_censored, sum(d$status == 0L))
      for (model in c("cause_specific", "subdistribution")) {
        expect_true(result$diagnostics[[model]]$ok)
        expect_equal(nrow(result[[model]]), 1L)
        expect_true(all(is.finite(unlist(result[[model]][c("HR", "HR95L", "HR95H", "Pvalue")]))))
        expect_equal(result[[model]], reference[[model]], tolerance = 1e-10)
      }
      expect_equal(sum(result$fits$subdistribution$y[, "status"]), sum(d$status == etype))
    }
  }
})

test_that("Fine-Gray point estimates agree with optional cmprsk references", {
  skip_if_not_installed("cmprsk")
  for (censored in c(TRUE, FALSE)) {
    d <- review_competing_data(censored)
    for (etype in c(2L, 5L)) {
      result <- competing_risk_COX(d, "time", "status", "x", etype = etype)
      reference <- cmprsk::crr(
        ftime = d$time, fstatus = d$status, cov1 = cbind(x = d$x),
        failcode = etype, cencode = 0L, gtol = 1e-7, maxiter = 100L
      )
      expect_true(reference$converged)
      expect_true(result$diagnostics$subdistribution$ok)
      expect_equal(unname(stats::coef(result$fits$subdistribution)),
                   unname(reference$coef), tolerance = 1e-4)
      expect_true(!is.null(result$fits$subdistribution$naive.var))
      expect_equal(result$subdistribution$HR,
                   unname(exp(stats::coef(result$fits$subdistribution))))
    }
  }
})
