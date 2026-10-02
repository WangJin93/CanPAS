test_that("short follow-up cohorts do not interrupt landmark pooling", {
  set.seed(10)
  short <- data.frame(time = runif(50, .01, .3), status = rep(0:1, 25),
                      dataset = "short", group = "Low")
  long <- data.frame(time = runif(50, 1, 5), status = rep(0:1, 25),
                     dataset = "long", group = "Low")
  r <- CanPAS:::km_surv_meta(rbind(short, long), c(1.5, 3))
  expect_true(nrow(r$landmarks) > 0)
  expect_true(all(r$landmarks$k == 1L))
  r0 <- CanPAS:::km_surv_meta(short, c(1, 3))
  expect_null(r0$landmarks)
  expect_null(r0$curve)
})

test_that("one-sided median cohorts are excluded and recorded", {
  a <- mk_cohort(100, seed = 876)
  b <- mk_cohort(100, seed = 877)
  a$GAPDH <- c(rep(1, 70), rep(0, 30))
  r <- cpas_km_pooled(list(A = a, B = b), "GAPDH", method = "ipd")
  expect_true("A" %in% r$empty_cohorts)
  expect_true(all(r$df$dataset == "B"))
  expect_false("A" %in% r$manifest$dropped_rows$cohort)
})

test_that("classification token conflicts cannot enter exact pooling", {
  a <- mk_cohort(120, seed = 878)
  tb <- data.frame(Accession = "SYN", Family = "OS", Token = "CSS",
                   TokenRole = "primary", PoolingClass = "Exact-equivalent")
  expect_error(cpas_meta("SYN", "GAPDH", merged = list(SYN = a),
                        pooling = "exact", class_table = tb), "token conflict")
})

test_that("HK alternative PI uses unadjusted random-effects SE", {
  b <- c(.1, .1001, .1002); se <- rep(.2, 3)
  hk <- CanPAS:::meta_pool(b, se, method = "HK", pi_method = "HK")
  re <- CanPAS:::meta_pool(b, se, method = "REML", pi_method = "t")
  expect_equal(hk$pi_alt_lower, re$pi_lower)
  expect_equal(hk$pi_alt_upper, re$pi_upper)
})
