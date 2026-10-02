# Hardening round: the new statistical defaults, the fail-safe, the analysis
# manifest, the endpoint-semantics / pooling-class layer and the shared-patient
# overlap handling. Everything here is offline: fixed synthetic cohorts and (for
# the metafor comparison) a fixed small effect-size table.

# --------------------------------------------------------------------------
# B1.1  pooling methods validated against metafor::rma()
# --------------------------------------------------------------------------
# A fixed effect-size table (log-HR and SE), small enough to be a literal.
fixed_es <- function() {
  list(b  = c(0.41, -0.13, 0.62, 0.28, 0.05, 0.55, -0.02),
       se = c(0.18, 0.22, 0.15, 0.26, 0.31, 0.12, 0.20))
}

test_that("method = 'DL' reproduces DerSimonian-Laird exactly, metafor included", {
  es <- fixed_es()
  po <- CanPAS:::meta_pool(es$b, es$se, method = "DL")
  expect_identical(po$method, "DL")
  # closed form, checked independently of the package code
  w <- 1 / es$se^2; bfe <- sum(w * es$b) / sum(w)
  Q <- sum(w * (es$b - bfe)^2); df <- length(es$b) - 1
  tau2 <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
  expect_equal(po$tau2, tau2, tolerance = 1e-12)
  expect_equal(po$logHR, sum(es$b / (es$se^2 + tau2)) / sum(1 / (es$se^2 + tau2)),
               tolerance = 1e-12)
  skip_if_not_installed("metafor")
  mr <- metafor::rma(es$b, es$se^2, method = "DL")
  expect_equal(po$tau2, mr$tau2, tolerance = 1e-12)
  expect_equal(po$logHR, as.numeric(mr$beta[1]), tolerance = 1e-12)
  expect_equal(po$se, mr$se[1], tolerance = 1e-12)
})

test_that("method = 'REML' agrees with metafor::rma(method = 'REML')", {
  es <- fixed_es()
  po <- CanPAS:::meta_pool(es$b, es$se, method = "REML")
  expect_identical(po$method, "REML")
  expect_identical(po$tau2_method, "REML")
  # tau^2 must solve the REML estimating equation of the random-effects model
  g <- function(t2) {
    w <- 1 / (es$se^2 + t2); mu <- sum(w * es$b) / sum(w)
    sum(w^2 * ((es$b - mu)^2 - es$se^2)) / sum(w^2) + 1 / sum(w) - t2
  }
  expect_lt(abs(g(po$tau2)), 1e-8)
  skip_if_not_installed("metafor")
  mr <- metafor::rma(es$b, es$se^2, method = "REML")
  # metafor stops on a very flat REML likelihood, so the two estimates agree to
  # ~1e-5 rather than to machine precision; the likelihood must agree exactly
  expect_equal(po$tau2, mr$tau2, tolerance = 1e-4)
  expect_equal(po$logHR, as.numeric(mr$beta[1]), tolerance = 1e-4)
  expect_equal(po$se, mr$se[1], tolerance = 1e-4)
  nll <- function(t2) {
    w <- 1 / (es$se^2 + t2); mu <- sum(w * es$b) / sum(w)
    0.5 * (sum(log(es$se^2 + t2)) + log(sum(w)) + sum(w * (es$b - mu)^2))
  }
  expect_equal(nll(po$tau2), nll(mr$tau2), tolerance = 1e-8)
  expect_true(po$tau2 >= 0)
})

test_that("method = 'HK' is the Knapp-Hartung variance and widens the interval", {
  es <- fixed_es()
  po <- CanPAS:::meta_pool(es$b, es$se, method = "HK")
  dl <- CanPAS:::meta_pool(es$b, es$se, method = "DL")
  expect_identical(po$method, "HK")
  expect_identical(po$se_method, "Hartung-Knapp-Sidik-Jonkman adjusted")
  # HK shares the REML tau^2 and only adjusts the variance of the pooled estimate
  expect_equal(po$tau2, CanPAS:::meta_pool(es$b, es$se, method = "REML")$tau2,
               tolerance = 1e-12)
  # the t(k-1) interval is wider than the z interval of the random-effects model
  expect_lt(po$lower, dl$lower)
  expect_gt(po$upper, dl$upper)
  skip_if_not_installed("metafor")
  # metafor 5.0.1 has no method = "HK": the Knapp-Hartung adjustment is selected
  # with test = "knha" on top of a tau^2 estimator (REML is its default)
  mr <- metafor::rma(es$b, es$se^2, method = "REML", test = "knha")
  expect_equal(po$tau2, mr$tau2, tolerance = 1e-4)
  expect_equal(po$logHR, as.numeric(mr$beta[1]), tolerance = 1e-4)
  expect_equal(po$se, mr$se[1], tolerance = 1e-4)
  expect_equal(log(po$lower), mr$ci.lb, tolerance = 1e-4)
  expect_equal(log(po$upper), mr$ci.ub, tolerance = 1e-4)
  expect_equal(po$p, mr$pval, tolerance = 1e-4)
})

test_that("method = 'FE' and the legacy alias 'RE' still work", {
  es <- fixed_es()
  fe <- CanPAS:::meta_pool(es$b, es$se, method = "FE")
  w <- 1 / es$se^2
  expect_equal(fe$tau2, 0)
  expect_equal(fe$logHR, sum(w * es$b) / sum(w), tolerance = 1e-12)
  re <- CanPAS:::meta_pool(es$b, es$se, method = "RE")
  expect_identical(re$method, "DL")
  expect_identical(re$method_requested, "RE")
  expect_error(CanPAS:::meta_pool(es$b, es$se, method = "nope"), "'method' must be one of")
})

test_that("k = 1 reports the cohort's own estimate and no heterogeneity", {
  po <- CanPAS:::meta_pool(0.3, 0.2, method = "REML")
  expect_equal(po$k, 1L)
  expect_true(is.na(po$tau2) && is.na(po$I2))
  expect_true(is.na(po$pi_lower))
  expect_true(is.finite(po$HR))
})

# --------------------------------------------------------------------------
# B1.2  the prediction interval is primary + clearly-labelled alternative
# --------------------------------------------------------------------------
test_that("both prediction intervals are returned and labelled", {
  es <- fixed_es()
  po <- CanPAS:::meta_pool(es$b, es$se, method = "REML")
  expect_true(is.finite(po$pi_lower) && is.finite(po$pi_upper))
  expect_true(is.finite(po$pi_alt_lower) && is.finite(po$pi_alt_upper))
  expect_match(po$pi_rule, "k-2")
  expect_match(po$pi_alt_rule, "normal")
  # the primary is the t(k-2) construction and is the wider of the two here
  expect_lt(po$pi_lower, po$pi_alt_lower)
  expect_gt(po$pi_upper, po$pi_alt_upper)
  # and it equals the documented formula
  tq <- stats::qt(0.975, df = po$k - 2)
  se_pi <- sqrt(po$se^2 + po$tau2)
  expect_equal(po$pi_lower, exp(po$logHR - tq * se_pi), tolerance = 1e-12)
  expect_equal(po$pi_alt_upper, exp(po$logHR + 1.96 * se_pi), tolerance = 1e-12)
})

# --------------------------------------------------------------------------
# A1  endpoint semantics / pooling classes
# --------------------------------------------------------------------------
test_that("the shipped companion table has the frozen schema and full coverage", {
  f <- system.file("extdata", "endpoint_semantics.csv", package = "CanPAS")
  skip_if(!nzchar(f), "companion table not shipped in this build")
  es <- endpoint_semantics()
  expect_true(nrow(es) >= 985L)
  expect_true(all(c("Accession", "Type", "Family", "Token", "TokenRole",
                    "SourceField", "EventDefinition", "TimeOrigin",
                    "CensoringRule", "CompetingEvents", "Derived",
                    "PoolingClass", "Evidence", "Note") %in% colnames(es)))
  expect_false(any(is.na(es$PoolingClass) | !nzchar(es$PoolingClass)))
  expect_false(any(duplicated(paste(es$Accession, es$Family))))
  expect_setequal(unique(es$Family), c("OS", "DSS", "DFS", "PFS", "MFS"))
  # the 197 catalog cohorts x 5 families are covered exactly
  di <- dataset_info
  expect_setequal(unique(es$Accession), di$Accession)
  expect_equal(nrow(es), length(di$Accession) * 5L)
})

test_that("endpoint_semantics() filters by cohort and by family", {
  skip_if(!nzchar(system.file("extdata", "endpoint_semantics.csv", package = "CanPAS")),
          "companion table not shipped in this build")
  one <- endpoint_semantics(accession = "GSE14814")
  expect_true(all(one$Accession == "GSE14814"))
  expect_equal(nrow(one), 5L)
  dfs <- endpoint_semantics(accession = "GSE31210", family = "DFS")
  expect_equal(nrow(dfs), 1L)
  expect_identical(as.character(dfs$Token), "RFS")
})

test_that("endpoint_pooling_class() answers every catalog cell with the frozen vocabulary", {
  di <- dataset_info
  fams <- c("OS", "DSS", "DFS", "PFS", "MFS")
  cls <- vapply(di$Accession, function(a)
    vapply(fams, function(f) endpoint_pooling_class(a, f), character(1)),
    character(length(fams)))
  expect_false(any(is.na(cls)))
  expect_true(all(cls %in% c("Exact-equivalent", "Clinically-related",
                             "Not-poolable", "Unknown", "Absent")))
  # the token-level distinction the layer exists for
  expect_identical(endpoint_pooling_class("GSE31210", "DFS"), "Clinically-related")
  expect_identical(endpoint_pooling_class("GSE31210", "OS"), "Exact-equivalent")
  expect_identical(endpoint_pooling_class("GSE4922_GPL96", "DFS"), "Exact-equivalent")
  # an absent family is the absence of an endpoint, not a pooling verdict
  expect_identical(endpoint_pooling_class("GSE14814", "MFS"), "Absent")
  expect_error(endpoint_pooling_class(c("A", "B"), "OS"), "one accession")
})

test_that("endpoint_pooling_class() falls back to the token map without a shipped table", {
  # the fallback is the documented token -> family rule, so it must reproduce the
  # exact/related split from the catalog alone
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  expect_error(endpoint_semantics(path = tmp), "not available")
  p <- file.path(tempdir(), "no_such_dir", "endpoint_semantics.csv")
  expect_error(endpoint_semantics(path = p), "not available")
})

test_that("a malformed companion table is rejected with an actionable message", {
  bad <- tempfile(fileext = ".csv")
  on.exit(unlink(bad), add = TRUE)
  utils::write.csv(data.frame(Accession = "X", Family = "OS", Token = "OS",
                              PoolingClass = "", stringsAsFactors = FALSE),
                   bad, row.names = FALSE)
  expect_error(endpoint_semantics(path = bad), "blank PoolingClass")
  utils::write.csv(data.frame(Accession = "X", Family = "OS", Token = "OS",
                              PoolingClass = "Sort-of", stringsAsFactors = FALSE),
                   bad, row.names = FALSE)
  expect_error(endpoint_semantics(path = bad), "outside the frozen vocabulary")
})

# --------------------------------------------------------------------------
# A1  pooling = "family" (default) vs "exact"
# --------------------------------------------------------------------------
test_that("pooling = 'exact' keeps only Exact-equivalent rows and records the rest", {
  # Both cohorts carry the DFS columns, so the analysis is estimable either way;
  # the fixture marks B's cell as Clinically-related, which is what exact mode
  # has to act on.
  a <- mk_cohort(200, "A", seed = 21)
  b <- mk_cohort(150, "B", seed = 22)
  c3 <- mk_cohort(180, "C", seed = 23)
  add_dfs <- function(d) { d$DFS_time <- d$OS_time; d$DFS_status <- d$OS_status; d }
  a <- add_dfs(a); b <- add_dfs(b); c3 <- add_dfs(c3)
  ct <- data.frame(
    Accession = c("A", "B", "C", "A", "B", "C"),
    Family    = c("OS", "OS", "OS", "DFS", "DFS", "DFS"),
    Token     = c("OS", "OS", "OS", "DFS", "DFS", "DFS"),
    TokenRole = c("primary", "primary", "primary", "primary", "contributing", "primary"),
    PoolingClass = c("Exact-equivalent", "Exact-equivalent", "Exact-equivalent",
                     "Exact-equivalent", "Clinically-related", "Exact-equivalent"),
    stringsAsFactors = FALSE)
  mg <- list(A = a, B = b, C = c3)

  fam <- cpas_meta(c("A", "B", "C"), "GAPDH", "DFS", merged = mg, class_table = ct)
  expect_identical(fam$pooled$pooling, "family")
  expect_equal(fam$pooled$k, 3L)
  expect_true("pooling_class" %in% colnames(fam$per_dataset))
  expect_setequal(strsplit(fam$pooled$pooling_classes, ", ", fixed = TRUE)[[1]],
                  c("Exact-equivalent", "Clinically-related"))
  expect_setequal(fam$manifest$pooling_classes,
                  c("Exact-equivalent", "Clinically-related"))

  exact <- cpas_meta(c("A", "B", "C"), "GAPDH", "DFS", merged = mg,
                     pooling = "exact", class_table = ct)
  expect_identical(exact$pooled$pooling, "exact")
  expect_equal(exact$pooled$k, 2L)
  expect_setequal(exact$per_dataset$dataset, c("A", "C"))
  expect_false("B" %in% exact$per_dataset$dataset)
  expect_true(all(exact$excluded$stage == "pooling"))
  expect_match(exact$excluded$reason[1], 'pooling = "exact"')
  expect_match(exact$excluded$reason[1], "Clinically-related")
  # the classes actually pooled are recorded on the result and in the manifest
  expect_identical(exact$pooled$pooling_classes, "Exact-equivalent")
  expect_identical(exact$manifest$pooling, "exact")
  expect_equal(nrow(exact$manifest$pooling_table), 2L)
  expect_true(all(c("accession", "family", "token", "pooling_class") %in%
                    colnames(exact$manifest$pooling_table)))
  # the per-cohort class is on the long table as well
  expect_identical(fam$per_dataset$pooling_class[fam$per_dataset$dataset == "B"],
                   "Clinically-related")
})

# --------------------------------------------------------------------------
# A6  overlap = warn / refuse / dedupe
# --------------------------------------------------------------------------
ov_fixture <- function() {
  data.frame(AccessionA = c("A", "B", "D"),
             AccessionB = c("B", "C", "E"),
             SharedPatients = c(40L, 12L, 7L),
             Basis = c("title", "platform-pair", "GEO-seed"),
             Evidence = c("titles", "platform ids", "sample ids"),
             stringsAsFactors = FALSE)
}

test_that("the shipped overlap register has the frozen schema", {
  f <- system.file("extdata", "cohort_overlap.csv", package = "CanPAS")
  skip_if(!nzchar(f), "overlap register not shipped in this build")
  ov <- cohort_overlap()
  expect_true(nrow(ov) >= 1L)
  expect_identical(colnames(ov),
                   c("AccessionA", "AccessionB", "SharedPatients", "Basis", "Evidence"))
  expect_true(is.integer(ov$SharedPatients))
  expect_false(any(is.na(ov$AccessionA) | is.na(ov$AccessionB)))
  expect_true(all(c("GSE11969", "GSE13213") %in% c(ov$AccessionA, ov$AccessionB)))
})

test_that("overlap = 'warn' proceeds, warns and records the offending pairs", {
  a <- mk_cohort(200, "A", seed = 31); b <- mk_cohort(150, "B", seed = 32)
  c3 <- mk_cohort(180, "C", seed = 33)
  expect_warning(
    r <- cpas_meta(c("A", "B", "C"), "GAPDH", "OS", merged = list(A = a, B = b, C = c3),
                   overlap_table = ov_fixture()),
    "sharing patients")
  expect_equal(r$pooled$k, 3L)
  expect_identical(r$pooled$overlap_mode, "warn")
  expect_equal(nrow(r$overlap_pairs), 2L)
  expect_equal(nrow(r$pooled$overlap_pairs[[1]]), 2L)
  expect_match(r$pooled$overlap_pairs_text, "A/B")
  expect_equal(nrow(r$manifest$overlap_pairs), 2L)
})

test_that("overlap = 'refuse' stops with an actionable error naming the pairs", {
  a <- mk_cohort(200, "A", seed = 34); b <- mk_cohort(150, "B", seed = 35)
  c3 <- mk_cohort(180, "C", seed = 36)
  err <- tryCatch(cpas_meta(c("A", "B", "C"), "GAPDH", "OS",
                            merged = list(A = a, B = b, C = c3),
                            overlap = "refuse", overlap_table = ov_fixture()),
                  error = function(e) conditionMessage(e))
  expect_match(err, "share patients")
  expect_match(err, "A / B")
  expect_match(err, "dedupe")
})

test_that("overlap = 'dedupe' keeps the larger cohort deterministically and records it", {
  a <- mk_cohort(200, "A", seed = 37); b <- mk_cohort(150, "B", seed = 38)
  c3 <- mk_cohort(180, "C", seed = 39)
  d <- mk_cohort(90, "D", seed = 40); e <- mk_cohort(90, "E", seed = 41)
  r <- cpas_meta(c("A", "B", "C", "D", "E"), "GAPDH", "OS",
                 merged = list(A = a, B = b, C = c3, D = d, E = e),
                 overlap = "dedupe", overlap_table = ov_fixture())
  # A>B>C is one group (A wins); D and E tie on n, so the accession sort decides
  # and keeps the first of the sorted pair, "D"
  expect_setequal(r$per_dataset$dataset, c("A", "D"))
  expect_setequal(r$excluded$cohort, c("B", "C", "E"))
  expect_true(all(r$excluded$stage == "overlap"))
  expect_match(r$excluded$reason[r$excluded$cohort == "B"], "shares patients with A")
  expect_match(r$excluded$reason[r$excluded$cohort == "E"], "shares patients with D")
  expect_equal(nrow(r$manifest$overlap_dropped), 3L)
  # deterministic: the same call in the reverse order gives the same survivors
  r2 <- cpas_meta(rev(c("A", "B", "C", "D", "E")), "GAPDH", "OS",
                  merged = list(A = a, B = b, C = c3, D = d, E = e),
                  overlap = "dedupe", overlap_table = ov_fixture())
  expect_setequal(r2$per_dataset$dataset, r$per_dataset$dataset)
})

test_that("an empty or absent register is handled explicitly", {
  a <- mk_cohort(120, "A", seed = 42); b <- mk_cohort(110, "B", seed = 43)
  empty <- ov_fixture()[0, , drop = FALSE]
  r <- cpas_meta(c("A", "B"), "GAPDH", "OS", merged = list(A = a, B = b),
                 overlap_table = empty)
  expect_equal(nrow(r$overlap_pairs), 0L)
  expect_identical(r$pooled$overlap_pairs_text, "")
  # with no register shipped at all, "refuse" and "dedupe" must refuse rather
  # than pretend the check was made
  missing <- tempfile(fileext = ".csv")
  expect_error(cohort_overlap(path = missing), "not available")
  expect_error(CanPAS:::.cpas_overlap_read(missing), "26_cohort_overlap")
  expect_error(CanPAS:::.cpas_overlap_read(NULL), "not available")
  # cpas_meta consults the shipped path, so with a register present the modes run
  skip_if(!nzchar(system.file("extdata", "cohort_overlap.csv", package = "CanPAS")),
          "register not shipped in this build")
  r2 <- cpas_meta(c("A", "B"), "GAPDH", "OS", merged = list(A = a, B = b),
                  overlap = "dedupe")
  expect_equal(nrow(r2$excluded[r2$excluded$stage == "overlap", ]), 0L)
  expect_gt(nrow(cohort_overlap()), 0L)
})

# --------------------------------------------------------------------------
# B1.3  the fail-safe default in the Cox / multivariable paths
# --------------------------------------------------------------------------
sim_multi_h <- function(n = 400, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(n); age <- stats::rnorm(n, 60, 10)
  grade <- factor(sample(c("G1", "G2", "G3"), n, TRUE))
  hr <- exp(0.4 * x + 0.02 * (age - 60) + 0.5 * (grade == "G3"))
  data.frame(ID = sprintf("S%03d", seq_len(n)), OS_time = stats::rexp(n, hr / 10),
             OS_status = stats::rbinom(n, 1, 0.7), marker = x, age = age, grade = grade)
}

test_that("auto_repair = FALSE returns an explicit not estimable result, with no model", {
  d <- sim_multi_h()
  d$flat <- 1
  r <- suppressMessages(COX_analysis(d, type = "OS",
                                     cont_Variates = c("marker", "age", "flat"),
                                     cate_Variates = "grade", method = "multi"))
  expect_s3_class(r, "cpas_not_estimable")
  expect_s3_class(r, "cpas_COX")
  expect_false(isTRUE(r$estimable))
  expect_identical(r$status, "not estimable")
  expect_length(r$models, 0L)
  expect_equal(nrow(r$results_table), 0L)
  expect_true("flat" %in% r$offending_terms)
  expect_match(r$reasons, "constant")
  # what the repair would have done is still reported
  expect_true("flat" %in% r$metadata$would_have_dropped$variable)
  expect_true(inherits(r$manifest, "cpas_manifest"))
  # and print() says so rather than showing a table
  out <- utils::capture.output(print(r))
  expect_true(any(grepl("NOT ESTIMABLE", out)))
})

test_that("auto_repair = FALSE never silently returns a reduced model", {
  d <- sim_multi_h(seed = 3)
  d$age_copy <- d$age          # perfectly collinear pair
  r <- suppressMessages(COX_analysis(d, type = "OS",
                                     cont_Variates = c("marker", "age", "age_copy"),
                                     method = "multi"))
  expect_false(isTRUE(r$estimable))
  expect_length(r$models, 0L)
  expect_gte(length(r$offending_terms), 1L)
  expect_true(all(r$offending_terms %in% c("age", "age_copy")))
  expect_match(r$status, "not estimable")
})

test_that("auto_repair = TRUE restores the repair and records every modification", {
  d <- sim_multi_h()
  d$flat <- 1
  r <- suppressWarnings(COX_analysis(d, type = "OS",
                                     cont_Variates = c("marker", "age", "flat"),
                                     cate_Variates = "grade", method = "multi",
                                     auto_repair = TRUE))
  expect_true(isTRUE(r$estimable))
  expect_gt(length(r$models), 0L)
  expect_false("flat" %in% r$metadata$final_covariates)
  expect_true("flat" %in% r$metadata$dropped_covariates$variable)
  expect_equal(nrow(r$manifest$dropped_covariates), 1L)
})

test_that("an estimable model is untouched by the fail-safe", {
  d <- sim_multi_h()
  r <- COX_analysis(d, type = "OS", cont_Variates = c("marker", "age"),
                    cate_Variates = "grade", method = "multi")
  expect_true(isTRUE(r$estimable))
  expect_equal(nrow(r$metadata$dropped_covariates), 0L)
  expect_true(all(is.finite(r$results_table$HR[!is.na(r$results_table$HR)])))
})

test_that("the deprecated drop_nonestimable alias still selects the repair", {
  d <- sim_multi_h(); d$flat <- 1
  expect_warning(
    r <- COX_analysis(d, type = "OS", cont_Variates = c("marker", "age", "flat"),
                      method = "multi", drop_nonestimable = TRUE),
    "deprecated")
  expect_true(isTRUE(r$estimable))
})

test_that("cpas_meta reports a cohort whose confounders cannot be co-estimated", {
  a <- mk_cohort(200, "A", seed = 51); b <- mk_cohort(150, "B", seed = 52)
  a$flat <- 1                                   # constant in A only
  b$flat <- seq_len(nrow(b))                    # present and variable in B
  r <- cpas_meta(c("A", "B"), "GAPDH", "OS", confounders = c("age", "flat"),
                 merged = list(A = a, B = b))
  expect_equal(r$pooled$k, 1L)
  expect_true(nrow(r$not_estimable) >= 1L)
  expect_match(r$not_estimable$reason[1], "auto_repair = FALSE")
  expect_match(r$not_estimable$reason[1], "flat")
  # with the repair explicitly requested both cohorts are admissable again
  r2 <- cpas_meta(c("A", "B"), "GAPDH", "OS", confounders = c("age", "flat"),
                  merged = list(A = a, B = b), auto_repair = TRUE)
  expect_equal(r2$pooled$k, 2L)
  expect_equal(nrow(r2$modifications), 1L)
  expect_identical(r2$modifications$variable, "flat")
})

# --------------------------------------------------------------------------
# B1.4  the analysis manifest
# --------------------------------------------------------------------------
test_that("every analysis result carries a manifest with the required fields", {
  a <- mk_cohort(200, "A", seed = 61); b <- mk_cohort(150, "B", seed = 62)
  r <- cpas_meta(c("A", "B"), "GAPDH", "OS", merged = list(A = a, B = b))
  m <- cpas_manifest(r)
  expect_s3_class(m, "cpas_manifest")
  expect_identical(m, r$manifest)
  need <- c("analysis", "cohorts", "n_cohorts", "family", "token", "tokens",
            "token_role", "pooling", "pooling_classes", "pooling_table",
            "overlap_mode", "overlap_source", "overlap_pairs", "overlap_dropped",
            "selection_rule", "dropped_rows", "dropped_covariates", "cut_rule",
            "cut_points_searched", "search_adjusted_p", "ph_test",
            "meta_method", "tau2", "I2", "pi_primary", "pi_primary_rule",
            "pi_alt", "pi_alt_rule", "versions", "timestamp", "notes")
  expect_true(all(need %in% names(m)))
  expect_identical(m$meta_method, "REML")
  expect_identical(m$family, "OS")
  expect_setequal(m$cohorts, c("A", "B"))
  expect_true(is.finite(m$tau2) || !is.na(m$tau2))
  expect_true(all(c("R", "CanPAS", "survival") %in% names(m$versions)))
  expect_s3_class(m$timestamp, "POSIXct")
  expect_match(m$selection_rule, "supplied in 'datasets'")
  expect_match(m$cut_rule, "no cut-point is searched")
  # printable and machine-readable
  out <- utils::capture.output(print(m))
  expect_true(any(grepl("CanPAS analysis manifest", out)))
  expect_true(any(grepl("meta_method", out)))
  df <- as.data.frame(m)
  expect_identical(colnames(df), c("field", "value"))
  expect_true(all(c("meta_method", "I2", "tau2", "timestamp") %in% df$field))
  expect_identical(df$value[df$field == "meta_method"], "REML")
})

test_that("the manifest records the cut-point rule and the number searched", {
  a <- mk_cohort(200, "A", seed = 63); b <- mk_cohort(150, "B", seed = 64)
  km <- cpas_km_pooled(list(A = a, B = b), marker = "GAPDH", type = "OS",
                       method = "both", landmarks = c(1, 3))
  m <- cpas_manifest(km)
  expect_match(m$cut_rule, "50% split")
  expect_identical(m$cut_points_searched, 0L)
  expect_true(is.na(m$search_adjusted_p))
  expect_identical(m$analysis, "cpas_km_pooled")
})

test_that("the manifest is available on a Cox result too", {
  d <- sim_multi_h()
  r <- COX_analysis(d, type = "OS", cont_Variates = "marker", method = "uni")
  m <- cpas_manifest(r)
  expect_identical(m$analysis, "COX_analysis")
  expect_identical(m$cohorts, "d")
  expect_true(is.data.frame(as.data.frame(m)))
  expect_error(cpas_manifest(list(a = 1)), "needs a CanPAS analysis result")
  expect_s3_class(cpas_manifest(), "cpas_manifest")
})

# --------------------------------------------------------------------------
# B1.5  the new defaults are documented
# --------------------------------------------------------------------------
test_that("the documented defaults are the implemented defaults", {
  expect_identical(eval(formals(cpas_meta)$method), c("REML", "DL", "HK", "FE"))
  expect_identical(eval(formals(cpas_meta)$pooling), c("family", "exact"))
  expect_identical(eval(formals(cpas_meta)$overlap), c("warn", "refuse", "dedupe"))
  expect_false(eval(formals(cpas_meta)$auto_repair))
  expect_false(eval(formals(COX_analysis)$auto_repair))
  expect_identical(eval(formals(plot_meta_forest)$digits), 4)
  expect_true(is.function(cpas_manifest))
  expect_true(is.function(cohort_overlap))
  expect_true(is.function(endpoint_semantics))
  expect_true(is.function(endpoint_pooling_class))
  # as.data.frame() is a registered method for the manifest
  expect_true(any(grepl("as.data.frame.cpas_manifest",
                        as.character(utils::methods(class = "cpas_manifest")))))
})

test_that("the man pages of the five previously Chinese topics are in English", {
  rd <- system.file("man", package = "CanPAS")
  skip_if(!nzchar(rd), "installed man/ not available")
  topics <- c("cpas_meta", "plot_meta_forest", "cpas_km_pooled",
              "plot_cpas_km_perdataset", "loo_meta")
  for (tp in topics) {
    p <- file.path(rd, paste0(tp, ".Rd"))
    skip_if(!file.exists(p), paste("no man page for", tp))
    txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
    cjk <- paste0("[", intToUtf8(0x4e00), "-", intToUtf8(0x9fff), "]")
    expect_false(grepl(cjk, txt), info = tp)
    expect_match(txt, tp, fixed = TRUE)
  }
})
