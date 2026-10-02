test_that("cpas_meta rejects non-poolable endpoint classes in family mode", {
  d <- mk_cohort(120, seed = 20260930)
  tb <- data.frame(Accession = "SYN", Family = "OS", Token = "OS",
                   TokenRole = "primary", PoolingClass = "Not-poolable")
  expect_error(cpas_meta("SYN", "GAPDH", merged = list(SYN = d),
                          class_table = tb), "no dataset produced estimable results.*pooling")
})

test_that("cpas_meta rejects nonbinary survival status before fitting", {
  d <- mk_cohort(120, seed = 20260931)
  d$OS_status <- rep(c(1, 2), length.out = nrow(d))
  expect_error(cpas_meta("SYN", "GAPDH", merged = list(SYN = d)),
               "status must be coded 0 .* 1")
})

test_that("missing requested confounders are not silently dropped", {
  d <- mk_cohort(120, seed = 20261001)
  expect_error(cpas_meta("SYN", "GAPDH", confounders = "not_in_data",
                          merged = list(SYN = d)),
               "no dataset produced estimable results.*absent")
  repaired <- cpas_meta("SYN", "GAPDH", confounders = "not_in_data",
                        auto_repair = TRUE, merged = list(SYN = d))
  expect_true(any(repaired$modifications$variable == "not_in_data"))
})
