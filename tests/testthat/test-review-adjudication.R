review_adjudication_table <- function(sizes) {
  prefixes <- c(GEO = "GSE", TCGA = "TCGA-", `EMBL-EBI` = "E-",
                CGGA = "CGGA", other = "PUB-")
  accession <- unlist(lapply(names(sizes), function(source) {
    paste0(prefixes[[source]], seq_len(sizes[[source]]))
  }), use.names = FALSE)
  data.frame(Accession = accession, Family = "OS", Token = "OS",
             Evidence = "Synthetic allocation regression fixture",
             stringsAsFactors = FALSE)
}

test_that("Hamilton awards each residual seat to a distinct stratum", {
  tab <- review_adjudication_table(c(GEO = 4L, TCGA = 4L, `EMBL-EBI` = 2L))
  sheet <- endpoint_adjudication(n = 4, seed = 7, table = tab)
  strata <- attr(sheet, "sampling")$strata
  expect_identical(strata$stratum, c("OS | EMBL-EBI", "OS | GEO", "OS | TCGA"))
  # Quotas are 0.8, 1.6, 1.6: EMBL-EBI wins first, then GEO wins the tie.
  expect_identical(strata$drawn, c(1L, 2L, 1L))
  expect_equal(nrow(sheet), 4L)
  expect_identical(as.integer(table(factor(CanPAS:::.cpas_accession_source(sheet$cohort),
                                         levels = c("EMBL-EBI", "GEO", "TCGA")))),
                   c(1L, 2L, 1L))
  expect_identical(sheet, endpoint_adjudication(n = 4, seed = 7, table = tab))
  expect_identical(sheet, endpoint_adjudication(n = 4, seed = 7,
                                               table = tab[nrow(tab):1L, ]))
})

test_that("Hamilton allocation obeys quotas and capacity for every sample size", {
  fixtures <- list(c(GEO = 4L, TCGA = 4L, `EMBL-EBI` = 2L),
                   c(GEO = 2L, TCGA = 1L),
                   c(GEO = 2L, other = 1L))
  for (sizes in fixtures) {
    tab <- review_adjudication_table(sizes)
    for (n in seq_len(nrow(tab) + 2L)) {
      sheet <- endpoint_adjudication(n = n, seed = 11, table = tab)
      strata <- attr(sheet, "sampling")$strata
      target <- min(n, nrow(tab))
      quota <- target * strata$available / sum(strata$available)
      expected <- as.integer(floor(quota))
      remaining <- target - sum(expected)
      if (remaining > 0L) {
        winners <- order(-(quota - floor(quota)), strata$stratum)[seq_len(remaining)]
        expected[winners] <- expected[winners] + 1L
      }
      expect_identical(strata$drawn, expected)
      expect_false(anyNA(strata$drawn))
      expect_true(all(strata$drawn >= 0L & strata$drawn <= strata$available))
      expect_equal(sum(strata$drawn), target)
      expect_equal(nrow(sheet), target)
      expect_equal(length(unique(sheet$record_id)), target)
      observed <- table(factor(paste("OS", CanPAS:::.cpas_accession_source(sheet$cohort),
                                     sep = " | "), levels = strata$stratum))
      expect_identical(as.integer(observed), strata$drawn)
    }
  }
})

test_that("Hamilton samples keep stable record ids across sizes and seeds", {
  tab <- review_adjudication_table(c(GEO = 4L, TCGA = 4L, `EMBL-EBI` = 2L))
  full <- endpoint_adjudication(n = nrow(tab), seed = 99, table = tab)
  for (n in c(1L, 4L, 7L)) {
    sheet <- endpoint_adjudication(n = n, seed = 11, table = tab)
    matched <- match(sheet$record_id, full$record_id)
    expect_false(anyNA(matched))
    expect_identical(sheet$cohort, full$cohort[matched])
    expect_identical(sheet$token, full$token[matched])
    expect_identical(sheet$evidence, full$evidence[matched])
    expect_identical(sheet$record_id, sort(sheet$record_id))
  }
})

test_that("Hamilton sampling preserves existing and absent RNG states", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  tab <- review_adjudication_table(c(GEO = 4L, TCGA = 4L, `EMBL-EBI` = 2L))
  set.seed(123)
  before <- get(".Random.seed", envir = .GlobalEnv)
  sheet <- endpoint_adjudication(n = 4, seed = 42, table = tab)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  rm(".Random.seed", envir = .GlobalEnv)
  expect_identical(sheet, endpoint_adjudication(n = 4, seed = 42, table = tab))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
