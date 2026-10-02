test_that("pooled KM counts excluded incomplete rows", {
  d <- mk_cohort(120, seed = 716)
  d$OS_time[1:7] <- NA_real_
  r <- cpas_km_pooled(list(SYN = d), "GAPDH", method = "ipd")
  expect_equal(nrow(r$df), 113L)
  expect_equal(r$n_dropped, 7L)
  expect_equal(sum(r$manifest$dropped_rows$n_dropped), 7L)
})

test_that("HK landmark intervals use the small-sample critical value", {
  set.seed(31)
  d <- do.call(rbind, lapply(1:4, function(j)
    data.frame(time = rexp(100), status = rbinom(100, 1, .65),
               dataset = paste0("S", j), group = "Low")))
  z <- CanPAS:::km_surv_meta(d, c(.5, 1), meta_method = "HK")
  mm <- do.call(rbind, lapply(split(d, d$dataset), function(x) {
    s <- summary(survival::survfit(survival::Surv(time, status) ~ 1, x), times = .5)
    data.frame(S = s$surv, SE = s$std.err)
  }))
  p <- CanPAS:::meta_pool(log(-log(mm$S)), mm$SE / (mm$S * abs(log(mm$S))), method = "HK")
  expect_equal(z$landmarks$lower[1], exp(-p$upper), tolerance = 1e-7)
  expect_equal(z$landmarks$upper[1], exp(-p$lower), tolerance = 1e-7)
})
