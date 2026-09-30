# pi_method: the prediction-interval construction rule (spec B6) --------------
# Offline: meta_pool() is fed fixed effect sizes; cpas_meta()/cpas_meta_panel()
# are fed synthetic cohorts through 'merged'.

fixed_es <- function() {
  list(b = c(0.31, 0.22, 0.42, 0.18, 0.29, 0.35),
       se = c(0.12, 0.15, 0.11, 0.19, 0.14, 0.13))
}

test_that("the default pi_method is 't' and reproduces the documented formula", {
  es <- fixed_es()
  po <- CanPAS:::meta_pool(es$b, es$se, method = "REML")
  expect_identical(po$pi_method, "t")
  tq <- stats::qt(0.975, df = po$k - 2)
  se_pi <- sqrt(po$se^2 + po$tau2)
  expect_equal(po$pi_lower, exp(po$logHR - tq * se_pi), tolerance = 1e-12)
  expect_equal(po$pi_upper, exp(po$logHR + tq * se_pi), tolerance = 1e-12)
  expect_match(po$pi_rule, "k-2")
  expect_match(po$pi_alt_rule, "normal")
  expect_equal(po$pi_alt_lower, exp(po$logHR - 1.96 * se_pi), tolerance = 1e-12)
})

test_that("pi_method = 'normal' swaps primary and alternative", {
  es <- fixed_es()
  t_po <- CanPAS:::meta_pool(es$b, es$se, method = "REML", pi_method = "t")
  n_po <- CanPAS:::meta_pool(es$b, es$se, method = "REML", pi_method = "normal")
  expect_identical(n_po$pi_method, "normal")
  # primary of one is the alternative of the other
  expect_equal(n_po$pi_lower, t_po$pi_alt_lower, tolerance = 1e-12)
  expect_equal(n_po$pi_upper, t_po$pi_alt_upper, tolerance = 1e-12)
  expect_equal(n_po$pi_alt_lower, t_po$pi_lower, tolerance = 1e-12)
  expect_equal(n_po$pi_alt_upper, t_po$pi_upper, tolerance = 1e-12)
  se_pi <- sqrt(n_po$se^2 + n_po$tau2)
  expect_equal(n_po$pi_lower, exp(n_po$logHR - 1.96 * se_pi), tolerance = 1e-12)
  # only the interval changes: the estimate and the heterogeneity do not
  expect_equal(n_po$HR, t_po$HR, tolerance = 1e-15)
  expect_equal(n_po$lower, t_po$lower, tolerance = 1e-15)
  expect_equal(n_po$tau2, t_po$tau2, tolerance = 1e-15)
  expect_equal(n_po$I2, t_po$I2, tolerance = 1e-15)
  expect_equal(n_po$p, t_po$p, tolerance = 1e-15)
})

test_that("pi_method = 'HK' uses the HK-adjusted SE with t(k-1)", {
  es <- fixed_es()
  po <- CanPAS:::meta_pool(es$b, es$se, method = "REML", pi_method = "HK")
  expect_identical(po$pi_method, "HK")
  w2 <- 1 / (es$se^2 + po$tau2)
  bm <- sum(w2 * es$b) / sum(w2)
  hk <- CanPAS:::.cpas_hk(es$b, w2, bm, sqrt(1 / sum(w2)))
  se_pi_hk <- sqrt(hk$se^2 + po$tau2)
  crit <- stats::qt(0.975, df = po$k - 1)
  expect_equal(po$pi_lower, exp(bm - crit * se_pi_hk), tolerance = 1e-12)
  expect_equal(po$pi_upper, exp(bm + crit * se_pi_hk), tolerance = 1e-12)
  expect_match(po$pi_rule, "Hartung-Knapp")
  # the HK interval is not the t(k-2) interval here
  t_po <- CanPAS:::meta_pool(es$b, es$se, method = "REML", pi_method = "t")
  expect_false(isTRUE(all.equal(po$pi_lower, t_po$pi_lower)))
  # the alternative is the clearly labelled t(k-2) construction
  expect_match(po$pi_alt_rule, "k-2")
  expect_equal(po$pi_alt_lower, t_po$pi_lower, tolerance = 1e-12)
})

test_that("two cohorts give a finite normal / HK interval and an NA t interval", {
  b <- c(0.3, 0.1); se <- c(0.1, 0.2)
  tt <- CanPAS:::meta_pool(b, se, method = "REML", pi_method = "t")
  expect_equal(tt$k, 2L)
  expect_true(is.na(tt$pi_lower) && is.na(tt$pi_upper))
  expect_true(is.finite(tt$pi_alt_lower))        # the normal alternative is defined
  nn <- CanPAS:::meta_pool(b, se, method = "REML", pi_method = "normal")
  expect_true(is.finite(nn$pi_lower) && is.finite(nn$pi_upper))
  expect_true(is.na(nn$pi_alt_lower))            # its t alternative needs k >= 3
  hh <- CanPAS:::meta_pool(b, se, method = "REML", pi_method = "HK")
  expect_true(is.finite(hh$pi_lower) && is.finite(hh$pi_upper))
  # a single cohort reports no interval at all
  one <- CanPAS:::meta_pool(0.3, 0.2, method = "REML", pi_method = "normal")
  expect_true(is.na(one$pi_lower) && is.na(one$pi_alt_lower))
  expect_identical(one$pi_method, "normal")
})

test_that("pi_method is validated and case-insensitive", {
  es <- fixed_es()
  expect_error(CanPAS:::meta_pool(es$b, es$se, pi_method = "nope"), "'pi_method' must be one of")
  expect_identical(CanPAS:::.cpas_pi_method("NORMAL"), "normal")
  expect_identical(CanPAS:::.cpas_pi_method(c("t", "HK")), "t")
  expect_identical(CanPAS:::meta_pool(es$b, es$se, pi_method = "hk")$pi_method, "HK")
})

test_that("cpas_meta records and honours pi_method", {
  mk <- function(seed) {
    set.seed(seed)
    n <- 120
    x <- stats::rnorm(n)
    data.frame(ID = sprintf("S%04d", seq_len(n)),
               OS_time = round(stats::rexp(n, 0.2) * exp(-0.4 * x), 3),
               OS_status = stats::rbinom(n, 1, 0.6),
               G1 = x, stringsAsFactors = FALSE)
  }
  m <- list(A = mk(1), B = mk(2), C = mk(3))
  d <- cpas_meta(c("A", "B", "C"), "G1", "OS", merged = m)
  expect_identical(d$pooled$pi_method, "t")
  expect_identical(d$input$pi_method, "t")
  expect_identical(cpas_manifest(d)$pi_method, "t")
  expect_match(cpas_manifest(d)$pi_primary_rule, "k-2")
  n <- cpas_meta(c("A", "B", "C"), "G1", "OS", merged = m, pi_method = "normal")
  expect_identical(n$pooled$pi_method, "normal")
  expect_identical(cpas_manifest(n)$pi_method, "normal")
  expect_equal(n$pooled$HR, d$pooled$HR, tolerance = 1e-15)
  expect_equal(n$pooled$lower, d$pooled$lower, tolerance = 1e-15)
  expect_equal(n$pooled$pi_lower, d$pooled$pi_alt_lower, tolerance = 1e-12)
  expect_true(any(grepl("pi_method = \"normal\"", cpas_manifest(n)$notes)))
  h <- cpas_meta(c("A", "B", "C"), "G1", "OS", merged = m, pi_method = "HK")
  expect_identical(h$pooled$pi_method, "HK")
  expect_equal(h$pooled$HR, d$pooled$HR, tolerance = 1e-15)
  out <- utils::capture.output(print(n))
  expect_true(any(grepl("pi_method: normal", out)))
  expect_true(any(grepl("primary, normal", out)))
  expect_error(cpas_meta(c("A", "B"), "G1", "OS", merged = m, pi_method = "nope"),
               "'pi_method' must be one of")
  # the two-cohort pool keeps the normal interval and reports NA for t
  two <- cpas_meta(c("A", "B"), "G1", "OS", merged = m)
  expect_equal(two$pooled$k, 2L)
  expect_true(is.na(two$pooled$pi_lower))
  expect_true(is.finite(two$pooled$pi_alt_lower))
})

test_that("cpas_meta_panel defaults to REML and passes pi_method through", {
  mkp <- function(seed) {
    set.seed(seed)
    n <- 120
    x <- stats::rnorm(n)
    data.frame(ID = sprintf("S%04d", seq_len(n)),
               OS_time = round(stats::rexp(n, 0.2) * exp(-0.4 * x), 3),
               OS_status = stats::rbinom(n, 1, 0.6),
               G1 = x, G2 = stats::rnorm(n), stringsAsFactors = FALSE)
  }
  m <- list(A = mkp(1), B = mkp(2), C = mkp(3))
  p <- cpas_meta_panel(c("A", "B", "C"), genes = c("G1", "G2"), type = "OS",
                       merged = m, progress = FALSE)
  # spec B10a: the panel default is REML, the cpas_meta default
  expect_identical(p$input$method, "REML")
  expect_identical(cpas_manifest(p)$meta_method, "REML")
  d <- cpas_meta_panel(c("A", "B", "C"), genes = c("G1", "G2"), type = "OS",
                       method = "REML", merged = m, progress = FALSE)
  expect_equal(p$table$HR, d$table$HR, tolerance = 1e-15)
  dl <- cpas_meta_panel(c("A", "B", "C"), genes = c("G1", "G2"), type = "OS",
                        method = "DL", merged = m, progress = FALSE)
  expect_false(isTRUE(all.equal(dl$table$tau2[order(dl$table$gene)],
                                d$table$tau2[order(d$table$gene)])))  # a different estimator
  # pi_method defaults to "t" and is recorded, per gene table included
  expect_identical(p$input$pi_method, "t")
  expect_identical(cpas_manifest(p)$pi_method, "t")
  expect_true(all(c("pi_lower", "pi_upper", "pi_alt_lower", "pi_alt_upper") %in%
                    colnames(p$table)))
  expect_match(cpas_manifest(p)$pi_primary_rule, "k-2")
  n <- cpas_meta_panel(c("A", "B", "C"), genes = "G1", type = "OS",
                       pi_method = "normal", merged = m, progress = FALSE)
  expect_identical(n$input$pi_method, "normal")
  expect_identical(cpas_manifest(n)$pi_method, "normal")
  g1 <- p$table$gene == "G1"
  expect_equal(n$table$pi_lower, p$table$pi_alt_lower[g1], tolerance = 1e-12)
  expect_equal(n$table$pi_alt_lower, p$table$pi_lower[g1], tolerance = 1e-12)
  expect_equal(n$table$pi_upper, p$table$pi_alt_upper[g1], tolerance = 1e-12)
  h <- cpas_meta_panel(c("A", "B", "C"), genes = "G1", type = "OS",
                       pi_method = "HK", merged = m, progress = FALSE)
  expect_match(cpas_manifest(h)$pi_primary_rule, "Hartung-Knapp")
  expect_error(cpas_meta_panel(c("A", "B"), genes = "G1", type = "OS",
                               pi_method = "nope", merged = m, progress = FALSE),
               "'pi_method' must be one of")
})
