test_that("registered overlaps use offline merged inputs and explicit policy", {
  a <- mk_cohort(120, seed = 914)
  b <- mk_cohort(120, seed = 915)
  ov <- data.frame(AccessionA = "A", AccessionB = "B", SharedPatients = 15,
    Basis = "test", Evidence = "synthetic")
  expect_warning(z <- cpas_meta(c("A", "B"), "GAPDH", "OS", merged = list(A = a, B = b),
    overlap = "warn", overlap_table = ov), "sharing patients")
  expect_equal(nrow(z$overlap_pairs), 1)
  expect_error(cpas_meta(c("A", "B"), "GAPDH", "OS", merged = list(A = a, B = b),
    overlap = "refuse", overlap_table = ov), "Refusing")
})

test_that("primary provenance is independent of attribute versus list storage", {
  d <- mk_cohort(120, seed = 916)
  q <- COX_screen_adjust(d, "OS", c("GAPDH", "age"))
  expect_equal(q$manifest$n_input, 120)
  expect_equal(q$manifest$n_analyzed, 120)
  expect_match(q$manifest$dataset_hash, "^[0-9a-f]{32}$")
  r <- plot_roc(d, "OS", "GAPDH", predict.time = 2)
  expect_equal(cpas_manifest(r)$events, sum(d$OS_status))
  expect_match(cpas_manifest(r)$dataset_hash, "^[0-9a-f]{32}$")
})
