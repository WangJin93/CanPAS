# Robustness regressions for the 2026-09-11 audit ---------------------------
# Every test here corresponds to a defect that was reproduced before the fix.

mk <- function(n = 60, hr = 1.4, ev = 0.6, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  data.frame(ID = sprintf("S%04d", seq_len(n)),
             OS_time = round(rexp(n, 0.2) * exp(-log(hr) * x), 3),
             OS_status = rbinom(n, 1, ev),
             RFS_time = round(rexp(n, 0.2), 3),
             RFS_status = rbinom(n, 1, ev),
             marker = x, stringsAsFactors = FALSE)
}

test_that("cpas_meta with a single cohort returns that estimate, not NaN", {
  A <- mk(120, 1.3, 0.6, 1)
  r <- suppressMessages(cpas_meta("A", "marker", type = "OS", merged = list(A = A)))
  expect_equal(r$pooled$k, 1L)
  expect_true(is.finite(r$pooled$HR))
  expect_true(is.finite(r$pooled$lower) && is.finite(r$pooled$upper))
  expect_true(is.finite(r$pooled$p))
  # heterogeneity statistics are undefined for k = 1 -> NA, never NaN
  expect_true(is.na(r$pooled$I2))
  expect_true(is.na(r$pooled$tau2))
  expect_true(is.na(r$pooled$p_heterogeneity))
})

test_that("meta_pool with identical estimates gives I2 = 0 (not NaN)", {
  p <- CanPAS:::meta_pool(c(0.2, 0.2, 0.2), c(0.1, 0.1, 0.1), method = "RE")
  expect_equal(p$I2, 0)
  expect_equal(p$tau2, 0)
  expect_true(is.finite(p$HR))
})

test_that("duplicated cohort ids are rejected", {
  A <- mk(60, 1.3, 0.6, 1); B <- mk(60, 1.8, 0.6, 2)
  expect_error(cpas_meta(c("A", "A"), "marker", type = "OS",
                         merged = list(A = A, A = B)),
               "duplicated accessions")
  expect_error(cpas_meta(c("A", "B"), "marker", type = "OS",
                         merged = list(A = A, A = B)),
               "duplicated names")
})

test_that("mixed endpoint tokens inside one family are reported", {
  a <- mk(60, 1.3, 0.6, 1); b <- mk(60, 1.5, 0.6, 2)
  names(a)[2:3] <- c("RFS_time", "RFS_status")
  names(b)[2:3] <- c("RFS_time", "RFS_status")
  r <- suppressMessages(cpas_meta(c("A", "B"), "marker", type = "RFS",
                                  merged = list(A = a, B = b)))
  expect_equal(unique(r$per_dataset$endpoint), "RFS")
})

test_that("loo_meta works with two cohorts", {
  A <- mk(120, 1.3, 0.6, 1); B <- mk(120, 1.8, 0.6, 2)
  r <- suppressMessages(cpas_meta(c("A", "B"), "marker", type = "OS",
                                  merged = list(A = A, B = B)))
  lo <- loo_meta(r)
  expect_equal(nrow(lo), 2L)
  expect_true(all(is.finite(lo$HR)))
  expect_error(loo_meta(structure(list(per_dataset = r$per_dataset[1, ],
                                       input = r$input), class = "cpas_meta")),
               "at least 2 datasets")
})

test_that("cpas_km_pooled validates status, time and events", {
  A <- mk(60, 1.3, 0.6, 1); B <- mk(60, 1.8, 0.6, 2)
  bad <- A; bad$OS_status <- ifelse(bad$OS_status == 1, 2, 0)
  expect_error(cpas_km_pooled(list(A = bad, B = B), "marker", type = "OS"),
               "must be coded 0")
  neg <- A; neg$OS_time[1] <- -1
  expect_error(cpas_km_pooled(list(A = neg, B = B), "marker", type = "OS"),
               "negative survival time")
  cens <- A; cens$OS_status <- 0; cens2 <- B; cens2$OS_status <- 0
  expect_error(cpas_km_pooled(list(A = cens, B = cens2), "marker", type = "OS"),
               "No event was observed")
})

test_that("landmarks far beyond the follow-up do not return NaN survival", {
  A <- mk(60, 1.3, 0.6, 1); B <- mk(60, 1.8, 0.6, 2)
  km <- suppressWarnings(cpas_km_pooled(list(A = A, B = B), "marker", type = "OS",
                                        method = "meta", landmarks = c(1, 500)))
  expect_true(all(is.finite(km$meta_landmarks$S)))
  expect_true(all(km$meta_landmarks$k >= 1))
})

test_that("plot_km validates time, status and the cut point", {
  d <- mk(40, 1.4, 0.6, 3)
  bad <- d; bad$OS_status[1] <- 3
  expect_error(plot_km(bad, "OS", "marker"), "coded 0")
  neg <- d; neg$OS_time[1] <- -2
  expect_error(plot_km(neg, "OS", "marker"), "negative")
  expect_error(plot_km(d, "OS", "marker", cutpoint = NA_real_), "finite")
  expect_error(plot_km(d, "OS", "marker", cutpoint = c(1, 2)), "single finite number")
})

test_that("plot_km warns when a group has fewer than 5 patients", {
  d <- mk(20, 1.4, 0.6, 4)
  d$marker <- c(rep(0, 19), 100)
  expect_warning(plot_km(d, "OS", "marker", cutpoint = 50), "smallest group")
})

test_that("plot_roc refuses horizons where the AUC is not estimable", {
  d <- mk(200, 1.6, 0.7, 5)
  expect_error(plot_roc(d, "OS", "marker", predict.time = max(d$OS_time) * 2),
               "exceeds the longest follow-up")
  expect_error(plot_roc(d, "OS", "marker", predict.time = 10^-6),
               "No event occurs before")
})

test_that("degenerate markers are reported as failures, not as HRs", {
  d <- mk(60, 1.4, 0.6, 6)
  d$const <- 1
  expect_warning(r <- COX_by_genes(d, "OS", c("marker", "const")), "could not be modelled")
  expect_equal(r$metadata$failed_genes, "const")
  expect_true(is.na(r$results_table$HR[r$results_table$gene == "const"]))
})

test_that("COX_analysis validates time and skips unusable covariates", {
  d <- mk(60, 1.4, 0.6, 7)
  neg <- d; neg$OS_time[1] <- -1
  expect_error(COX_analysis(neg, type = "OS", cont_Variates = "marker"), "negative")
  d$onelev <- factor("A")
  expect_warning(r <- COX_analysis(d, type = "OS", cont_Variates = "marker",
                                   cate_Variates = "onelev", method = "uni"),
                 "skipped")
  expect_true("marker" %in% r$results_table$Variates)
  expect_true(any(grepl("^onelev", r$metadata$skipped_variables)))
})

test_that("COX_analysis keeps non-syntactic covariate names", {
  d <- mk(80, 1.4, 0.7, 8)
  d[["a b"]] <- rnorm(nrow(d))
  r <- COX_analysis(d, type = "OS", cont_Variates = c("marker", "a b"), method = "multi")
  expect_true("a b" %in% r$results_table$Variates)
})

test_that("COX_analysis warns when events per variable is low and reports the PH test", {
  d <- mk(60, 1.4, 0.4, 9)
  d$v1 <- rnorm(nrow(d)); d$v2 <- rnorm(nrow(d)); d$v3 <- rnorm(nrow(d))
  expect_warning(COX_analysis(d, type = "OS", cont_Variates = c("marker", "v1", "v2", "v3"),
                              method = "multi"),
                 "events per variable")
  r <- suppressWarnings(COX_analysis(d, type = "OS", cont_Variates = c("marker", "v1"),
                                     method = "multi"))
  expect_false(is.null(r$metadata$ph_test))
  expect_true(is.finite(r$metadata$complete_cases))
})

test_that("COX_by_genes validates the status coding and negative times", {
  d <- mk(40, 1.4, 0.6, 10)
  bad <- d; bad$OS_status <- ifelse(bad$OS_status == 1, 2, 0)
  expect_error(COX_by_genes(bad, "OS", "marker"), "coded 0")
  neg <- d; neg$OS_time[1] <- -1
  expect_error(COX_by_genes(neg, "OS", "marker"), "negative")
})

test_that("endpoint_resolve validates raw tokens", {
  expect_equal(endpoint_resolve("GSE31210", "DFS"), "RFS")
  expect_true(is.na(endpoint_resolve("GSE31210", "NOPE")))
  expect_true(is.na(endpoint_resolve("NOT_A_COHORT", "OS")))
  expect_equal(endpoint_resolve("GSE31210", "RFS"), "RFS")
})
