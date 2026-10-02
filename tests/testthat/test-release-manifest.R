test_that("ordinary Cox rejects duplicated identifiers and audits missingness", {
  d <- mk_cohort(120, seed = 811)
  bad <- d; bad$ID[2] <- bad$ID[1]
  expect_error(COX_analysis(bad, "OS", "GAPDH"), "duplicated ID")
  d$OS_time[1:4] <- NA
  x <- COX_analysis(d, "OS", c("GAPDH", "age"))
  expect_equal(x$manifest$n_input, 120)
  expect_equal(x$manifest$n_analyzed, 116)
  expect_equal(x$manifest$n_excluded, 4)
  expect_match(x$manifest$dataset_hash, "^[0-9a-f]{32}$")
  expect_equal(x$manifest$model_audit[[1]]$n_analyzed, 116)
  expect_setequal(x$manifest$covariates_requested, c("GAPDH", "age"))
  d$age[5:10] <- NA
  y <- COX_analysis(d, "OS", c("GAPDH", "age"))
  expect_equal(y$manifest$model_audit$GAPDH$n_analyzed, 116)
  expect_equal(y$manifest$model_audit$age$n_analyzed, 110)
  expect_equal(y$manifest$n_analyzed, 116)
  expect_false(identical(x$manifest$dataset_hash, y$manifest$dataset_hash))
})

test_that("meta and pooled KM manifests count requested rows and exclusions", {
  d <- mk_cohort(120, seed = 812); d$OS_time[1:3] <- NA
  z <- cpas_meta(c("A", "B"), "GAPDH", "OS", merged = list(A = d, B = d))
  expect_equal(z$manifest$n_input, 240)
  expect_equal(z$manifest$n_analyzed, 234)
  expect_equal(z$manifest$n_excluded, 6)
  expect_equal(z$manifest$events, sum(z$per_dataset$events))
  expect_match(z$manifest$dataset_hash, "^[0-9a-f]{32}$")
  expect_match(z$manifest$hash_scope, "per-cohort")
  k <- cpas_km_pooled(list(A = d, B = d), "GAPDH", method = "ipd")
  expect_equal(k$manifest$n_input, 240)
  expect_equal(k$manifest$n_analyzed, 234)
  expect_equal(k$manifest$n_excluded, 6)
  expect_equal(k$manifest$events, sum(k$df$status))
})

test_that("manifest schema does not claim observed truth by default", {
  m <- cpas_manifest()
  expect_identical(m$coverage_status, "schema_only")
  expect_true(is.na(m$dataset_hash))
  expect_true(is.na(m$seed))
})

test_that("CIF manifests measure analyzed input and do not fabricate endpoint evidence", {
  d <- data.frame(t = 1:30, s = rep(0:2, 10))
  d$t[1] <- NA_real_
  x <- cif_fit(d, time = "t", status = "s", times = c(5, 10))
  m <- cpas_manifest(x)
  expect_identical(m$primary_path, "cif_fit")
  expect_equal(m$n_input, 30L)
  expect_equal(m$n_analyzed, 29L)
  expect_equal(m$n_excluded, 1L)
  expect_equal(m$events, 20L)
  expect_match(m$dataset_hash, "^[a-f0-9]{32}$")
  expect_true(is.na(m$accession))
  expect_true(is.na(m$raw_token))
  expect_match(m$endpoint_evidence, "not supplied")
  expect_true(is.na(m$seed))
  y <- cif_fit(d, time = "t", status = "s", times = c(5, 10))
  expect_identical(m$dataset_hash, cpas_manifest(y)$dataset_hash)
  d$t[2] <- 2.5
  z <- cif_fit(d, time = "t", status = "s", times = c(5, 10))
  expect_false(identical(m$dataset_hash, cpas_manifest(z)$dataset_hash))
})

test_that("attribute-based manifests are retrieved without changing plot contract", {
  d <- mk_cohort(100, seed = 172)
  p <- plot_roc(d, marker = "GAPDH", predict.time = 1)
  expect_s3_class(cpas_manifest(p), "cpas_manifest")
  k <- plot_km(d, marker = "GAPDH")
  m <- cpas_manifest(k)
  expect_equal(m$n_input, 100L)
  expect_equal(m$n_analyzed, 100L)
  expect_equal(m$n_excluded, 0L)
  expect_equal(m$events, sum(d$OS_status == 1))
  expect_match(m$dataset_hash, "^[a-f0-9]{32}$")
})
