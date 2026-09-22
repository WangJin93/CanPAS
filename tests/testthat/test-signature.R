test_that("signature parser aligns weights to measured genes", {
  skip_on_cran()
  testthat::skip_if_not(requireNamespace("jsonlite", quietly = TRUE))
  online <- tryCatch({ CanPAS::get_data("GSE14814", "surv_data"); TRUE },
                     error = function(e) FALSE)
  skip_if_not(online, "CanPAS API unreachable")
  s <- get_signature_value("0.3*GAPDH + 0.7*TNS1", "GSE14814")
  expect_s3_class(s, "cpas_signature")
  expect_true(all(c("ID", "signature") %in% colnames(s$results_table)))
  expect_equal(s$signature_info$genes, c("GAPDH", "TNS1"))
})


test_that("the documented signature syntax is parsed as written", {
  # parsing is checked without the network: an unsupported formula must fail
  # on the formula itself, a supported one must not
  ok <- c("0.5*GAPDH + 0.5*ACTB", "0.5*GAPDH+0.5*ACTB", "0.5*GAPDH - 0.3*ACTB",
          "0.5*GAPDH -0.3*ACTB", "1e-2*GAPDH + 2*ACTB", "0.5*GAPDH")
  bad <- c("GAPDH + ACTB", "0.5*GAPDH + ", "0.5*GAPDH & 0.5*ACTB")
  parse_err <- function(x) {
    tryCatch({
      # a fake dataset makes the parser run and then fail at the fetch, which is
      # enough to tell a formula error from an accepted formula
      get_signature_value(x, "NO-SUCH-COHORT")
      NA_character_
    }, error = function(e) conditionMessage(e))
  }
  for (f in ok) expect_false(grepl("Unsupported text|Could not parse|must contain", parse_err(f)),
                             info = f)
  for (f in bad) expect_match(parse_err(f), "Unsupported text|Could not parse|must contain",
                              info = f)
})
