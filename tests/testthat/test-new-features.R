# Tests for the features of this release -------------------------------------
# Competing risks, multiple-testing correction, meta prediction interval and
# the preferred argument names (cohorts / genes).

mk_comp <- function(n = 400, seed = 7) {
  set.seed(seed)
  x <- rnorm(n)
  t1 <- rexp(n, 0.25) * exp(-0.5 * x)      # event of interest
  t2 <- rexp(n, 0.15)                      # competing event
  tt <- pmin(t1, t2)
  cause <- ifelse(t1 < t2, 1, 2)
  cen <- rexp(n, 0.05)
  cause[cen < tt] <- 0
  tt <- pmin(tt, cen)
  data.frame(time = round(tt, 4), status = cause, age = x,
             group = ifelse(x > 0, "high", "low"), stringsAsFactors = FALSE)
}

test_that("cif_fit returns a monotone cumulative incidence per cause", {
  d <- mk_comp()
  ci <- cif_fit(d, "time", "status", times = c(1, 3, 5, 8))
  expect_s3_class(ci, "cpas_cif")
  expect_true(all(c("group", "cause", "time", "n_risk", "cif", "lower", "upper") %in%
                    colnames(ci$table)))
  for (cs in ci$causes) {
    v <- ci$table$cif[ci$table$cause == cs]
    expect_true(all(diff(v) >= -1e-12))          # CIF is non-decreasing
  }
  expect_equal(sort(unique(ci$table$cause)), c(1, 2))
  expect_equal(ci$n, nrow(d))
  expect_equal(sum(ci$events), ci$n - ci$n_censored)
  # at the last time point the two CIFs must stay below the all-cause failure
  expect_true(all(ci$table$cif <= 1))
  # grouped version
  cig <- cif_fit(d, "time", "status", group = "group", times = c(1, 5))
  expect_equal(sort(unique(cig$table$group)), c("high", "low"))
  expect_s3_class(plot_cif(ci), "ggplot")
})

test_that("cif_fit validates its inputs", {
  d <- mk_comp(100)
  expect_error(cif_fit(data.frame(t = 1:20, s = rbinom(20, 1, .5)), "t", "s"),
               "At least two competing event types")
  bad <- d; bad$time[1] <- -1
  expect_error(cif_fit(bad, "time", "status"), "Negative follow-up")
  expect_error(cif_fit(d, "nope", "status"), "not found")
})

test_that("competing_risk_COX matches cmprsk::crr for the subdistribution HR", {
  d <- mk_comp()
  cr <- competing_risk_COX(d, "time", "status", covariates = "age", etype = 1)
  expect_s3_class(cr, "cpas_competing")
  expect_equal(cr$events_interest, sum(d$status == 1))
  expect_equal(cr$events_competing, sum(d$status == 2))
  expect_true(all(is.finite(cr$cause_specific$HR)))
  expect_true(all(is.finite(cr$subdistribution$HR)))
  expect_true(all(cr$subdistribution$HR95L < cr$subdistribution$HR))
  expect_true(all(cr$subdistribution$HR < cr$subdistribution$HR95H))
  expect_output(print(cr), "subdistribution hazard ratios")
  if (requireNamespace("cmprsk", quietly = TRUE)) {
    cc <- cmprsk::crr(d$time, d$status, cov1 = cbind(age = d$age),
                      failcode = 1, cencode = 0)
    ref <- exp(summary(cc)$coef[1])
    # crr uses its own fitting routine (different convergence tolerance); agreement
    # to 1e-4 relative is the level supported by the two implementations
    expect_equal(cr$subdistribution$HR[1], ref, tolerance = 1e-4)
    expect_equal(unname(log(cr$subdistribution$HR[1])), unname(summary(cc)$coef[1]),
                 tolerance = 1e-4)
  }
})

test_that("competing_risk_COX validates the event coding", {
  d <- mk_comp(200)
  expect_error(competing_risk_COX(transform(d, status = as.numeric(status == 1)),
                                  "time", "status", covariates = "age"),
               "at least two event types")
  expect_error(competing_risk_COX(d, "time", "status", covariates = "age", etype = 3),
               "not among the event codes")
  expect_error(competing_risk_COX(d, "time", "status"), "'covariates'")
  small <- d[d$status == 1, ][1:3, ]
  expect_error(competing_risk_COX(rbind(small, d[d$status == 2, ][1:3, ]), "time", "status",
                                  covariates = "age"),
               "Fewer than 5 events of interest")
})

test_that("COX_by_genes reports Benjamini-Hochberg FDR over the panel", {
  set.seed(11)
  n <- 200
  d <- data.frame(ID = sprintf("S%03d", seq_len(n)),
                  OS_time = round(rexp(n, 0.2), 3),
                  OS_status = rbinom(n, 1, 0.7))
  for (g in paste0("G", 1:8)) d[[g]] <- rnorm(n)
  r <- COX_by_genes(d, type = "OS", genes = paste0("G", 1:8))
  tb <- r$results_table
  expect_true(all(c("Pvalue", "P_adj", "P_adj_text") %in% colnames(tb)))
  expect_equal(tb$P_adj, stats::p.adjust(tb$Pvalue, method = "BH"), tolerance = 1e-12)
  ok <- !is.na(tb$Pvalue)
  expect_true(all(tb$P_adj[ok] >= tb$Pvalue[ok] - 1e-12))   # FDR is never smaller
  expect_equal(r$metadata$fdr_method, "BH")
  expect_equal(r$metadata$n_significant_fdr005, sum(tb$P_adj < 0.05, na.rm = TRUE))
})

test_that("cpas_meta reports a prediction interval for k >= 3 and NA below that", {
  mk <- function(n, hr, seed) {
    set.seed(seed); x <- rnorm(n)
    data.frame(ID = sprintf("S%04d", seq_len(n)),
               OS_time = round(rexp(n, 0.2) * exp(-log(hr) * x), 3),
               OS_status = rbinom(n, 1, 0.6), G = x)
  }
  A <- mk(120, 1.5, 1); B <- mk(140, 1.2, 2); C <- mk(100, 2.0, 3)
  r3 <- suppressMessages(cpas_meta(datasets = c("A", "B", "C"), marker = "G",
                                   type = "OS", merged = list(A = A, B = B, C = C)))
  expect_true(is.finite(r3$pooled$pi_lower) && is.finite(r3$pooled$pi_upper))
  expect_true(r3$pooled$pi_lower <= r3$pooled$lower + 1e-12)  # PI is wider than the CI
  expect_true(r3$pooled$pi_upper >= r3$pooled$upper - 1e-12)
  expect_output(print(r3), "prediction interval")

  r2 <- suppressMessages(cpas_meta(datasets = c("A", "B"), marker = "G", type = "OS",
                                   merged = list(A = A, B = B)))
  expect_true(is.na(r2$pooled$pi_lower))
  expect_true(is.na(r2$pooled$pi_upper))
})

test_that("cpas_meta takes its study axis through 'datasets' only", {
  set.seed(5); n <- 120
  mk <- function(k) data.frame(ID = sprintf("S%03d", seq_len(n)),
                               OS_time = round(rexp(n, 0.2), 3),
                               OS_status = rbinom(n, 1, 0.6), G = rnorm(n))
  A <- mk(1); B <- mk(2)
  r <- suppressMessages(cpas_meta(datasets = c("A", "B"), marker = "G", type = "OS",
                                  merged = list(A = A, B = B)))
  expect_equal(r$input$datasets, c("A", "B"))
  expect_true(is.finite(r$pooled$HR))
  expect_equal(nrow(r$per_dataset), 2L)
  expect_error(cpas_meta(marker = "G", type = "OS", merged = list(A = A)),
               "at least one accession")
  # the earlier internal argument names must no longer exist
  expect_error(cpas_meta(tables = c("A", "B"), marker = "G", type = "OS",
                         merged = list(A = A, B = B)), "unused argument")
  expect_error(cpas_meta(cohorts = c("A", "B"), marker = "G", type = "OS",
                         merged = list(A = A, B = B)), "unused argument")
})

test_that("genes is accepted as the preferred name for gene", {
  expect_error(get_expr_data("GSE14814", gene = "A", genes = "B"),
               "both supplied")
  expect_error(get_expr_data("GSE14814"), "Supply the gene symbols")
})


test_that("Fine-Gray standard errors are clustered by patient (match cmprsk::crr)", {
  skip_if_not_installed("cmprsk")
  for (sd_seed in c(7, 23)) {
    d <- mk_comp(400, seed = sd_seed)
    cr <- competing_risk_COX(d, "time", "status", covariates = "age", etype = 1)
    pkg_se <- (log(cr$subdistribution$HR95H) - log(cr$subdistribution$HR95L)) /
      (2 * stats::qnorm(0.975))
    cc <- cmprsk::crr(d$time, d$status, cov1 = cbind(age = d$age),
                      failcode = 1, cencode = 0)
    sm <- summary(cc)
    ref_se <- as.numeric(sm$coef[1, "se(coef)"])
    # the naive (unclustered) variance is 17-23% too small, so require < 2%
    expect_lt(abs(pkg_se / ref_se - 1), 0.02)
    expect_equal(log(cr$subdistribution$HR[1]), as.numeric(sm$coef[1, "coef"]),
                 tolerance = 1e-4)
  }
})

test_that("perfect separation is reported instead of a meaningless HR", {
  n <- 60
  d <- data.frame(ID = sprintf("S%02d", seq_len(n)),
                  OS_time = round(stats::rexp(n, 0.2), 3),
                  OS_status = rep(c(1, 0), each = 30))
  d$sep <- as.numeric(d$OS_status == 1)          # separates events from censored

  # COX_by_genes: NA row + the gene recorded as failed, with a reason
  expect_warning(r1 <- COX_by_genes(d, type = "OS", genes = "sep"),
                 "could not be modelled")
  expect_true(is.na(r1$results_table$HR[r1$results_table$gene == "sep"]))
  expect_equal(r1$metadata$failed_genes, "sep")
  expect_true(grepl("infinite|separation|converged", r1$metadata$failure_reasons[["sep"]],
                    ignore.case = TRUE))

  # COX_analysis, univariable: the covariate is skipped and the reason is kept
  expect_error(COX_analysis(d, type = "OS", cont_Variates = "sep", method = "uni"),
               "No estimable covariate remained")

  # COX_analysis, multivariable: refuse to publish the model
  expect_error(COX_analysis(d, type = "OS", cont_Variates = "sep", method = "multi"),
               "not estimable")

  # a clean covariate in the same data set still works
  set.seed(2); d$ok <- stats::rnorm(n)
  expect_warning(r2 <- COX_analysis(d, type = "OS", cont_Variates = c("ok", "sep"),
                                    method = "uni"), "skipped")
  expect_true("ok" %in% r2$results_table$Variates)
  expect_false("sep" %in% r2$results_table$Variates)
})

test_that("cpas_meta excludes a cohort whose marker effect is not estimable", {
  n <- 120
  set.seed(4)
  good <- data.frame(ID = sprintf("S%03d", seq_len(n)),
                     OS_time = round(stats::rexp(n, 0.2), 3),
                     OS_status = stats::rbinom(n, 1, 0.6), marker = stats::rnorm(n))
  bad <- good
  bad$OS_status <- rep(c(1, 0), each = n / 2)
  bad$marker <- as.numeric(bad$OS_status == 1)   # separated marker
  r <- suppressMessages(cpas_meta(datasets = c("GOOD", "BAD"), marker = "marker",
                                  type = "OS", merged = list(GOOD = good, BAD = bad)))
  expect_equal(nrow(r$per_dataset), 1L)
  expect_equal(r$per_dataset$dataset, "GOOD")
  expect_true(grepl("cox:", r$errors[["BAD"]]))
})

test_that("COX_by_datasets validates its inputs and drops the old argument names", {
  expect_error(COX_by_datasets(datasets = character(0), gene = "TP53"),
               "at least one accession")
  expect_error(COX_by_datasets(gene = "TP53"), "at least one accession")
  expect_error(COX_by_datasets(datasets = c("A", "B"), gene = c("TP53", "PTEN")),
               "single symbol")
  # the earlier internal argument names must no longer exist
  expect_error(COX_by_datasets(tables = "GSE14814", gene = "TP53"), "unused argument")
  expect_error(COX_by_datasets(cohorts = "GSE14814", gene = "TP53"), "unused argument")
})

test_that("get_data keeps the API request parameter 'table' while exposing 'dataset'", {
  # the mirror API only understands ?action=...&table=...; the R argument is
  # 'dataset', so the two must not be conflated (an offline test is the only way
  # to catch a rename that silently changes the external contract)
  err <- tryCatch(get_data(dataset = "GSE14814", action = "surv_data",
                           base_url = "http://127.0.0.1:9/none"),
                  error = function(e) conditionMessage(e))
  expect_true(grepl("table=GSE14814", err, fixed = TRUE))
  expect_false(grepl("dataset=GSE14814", err, fixed = TRUE))
  expect_error(get_data(dataset = "bad id", action = "surv_data"),
               "unsupported characters")
})

# --- covariate typing and safeguards in competing_risk_COX() ---------------
# Mirrored clinical columns arrive as text ("55.4"), and coxph() would read a
# character column as a factor with one level per value: 111 coefficients for
# age. These tests pin the corrected behaviour and the guardrails.

sim_cr <- function(n = 300, seed = 11, text = TRUE) {
  set.seed(seed)
  age_num <- round(stats::rnorm(n, 65, 8), 1)
  sex <- sample(c("female", "male"), n, TRUE)
  t1 <- stats::rexp(n, exp(0.03 * (age_num - 65)) / 8)   # cancer death
  t2 <- stats::rexp(n, 0.02)                             # other death
  status <- ifelse(t1 < t2 & t1 < 20, 1, ifelse(t2 <= t1 & t2 < 20, 2, 0))
  # the mirror returns the clinical columns as text, which is what made coxph()
  # read a continuous variable as a factor with one level per value
  data.frame(time = pmin(t1, t2, 20), status = status,
             age = if (text) as.character(age_num) else age_num,
             sex = sex, stringsAsFactors = FALSE)
}

test_that("a text numeric covariate yields one coefficient, not one per value", {
  d <- sim_cr()
  expect_true(is.character(d$age))
  expect_message(
    r <- competing_risk_COX(d, "time", "status", covariates = c("age", "sex"), etype = 1),
    "inferred from text")
  expect_identical(unname(r$covariate_types["age"]), "continuous (numeric text)")
  expect_true(nrow(r$cause_specific) == 2L)          # age + sex(male), not 300
  expect_true(all(r$cause_specific$Variates %in% c("age", "sex")))
  expect_true(all(is.finite(r$cause_specific$HR)))
  expect_lt(max(abs(log(r$cause_specific$HR)), na.rm = TRUE), 5)
  expect_true(r$diagnostics$cause_specific$ok)
})

test_that("an explicitly numeric covariate stays continuous and silent", {
  d <- sim_cr()
  d$age <- as.numeric(d$age)
  d$sex <- factor(d$sex)                 # already typed: nothing to infer
  expect_silent(r <- competing_risk_COX(d, "time", "status",
                                        covariates = c("age", "sex"), etype = 1))
  expect_identical(unname(r$covariate_types["age"]), "continuous")
  expect_equal(nrow(r$cause_specific), 2L)
})

test_that("a categorical covariate with too many levels is refused, not reported", {
  d <- sim_cr(n = 200)
  d$many <- factor(sprintf("g%03d", seq_len(nrow(d))))    # one level per patient
  expect_error(competing_risk_COX(d, "time", "status", covariates = "many", etype = 1),
               "more than 20 levels")
  # the same call is allowed when the limit is raised deliberately
  r <- suppressWarnings(competing_risk_COX(d, "time", "status", covariates = "many",
                                           etype = 1, max_levels = Inf))
  expect_s3_class(r, "cpas_competing")
})

test_that("a separated covariate is flagged and its coefficient is not published", {
  d <- sim_cr(n = 200)
  d$sep <- factor(ifelse(d$status == 1, "event", "other"))   # perfectly separated
  r <- suppressWarnings(competing_risk_COX(d, "time", "status",
                                           covariates = c("age", "sep"), etype = 1))
  expect_false(r$diagnostics$cause_specific$ok)
  expect_true(nzchar(r$diagnostics$cause_specific$reason))
  cs <- r$cause_specific
  # the unusable level is dropped; nothing absurd is returned
  if (!is.null(cs) && nrow(cs))
    expect_true(all(is.finite(cs$HR)) && max(abs(log(cs$HR)), na.rm = TRUE) < 20)
})

# --- cif_fit(): a group column must be categorical ------------------------
# 'group' names a grouping column and is factored as supplied, so a continuous
# variable would give one stratum per distinct value. These tests pin the
# refusal, the documented quantile splits, and that a real group is unchanged.

test_that("a continuous group column is refused, not turned into 100+ strata", {
  d <- sim_cr(n = 160)
  expect_true(is.character(d$age))
  err <- tryCatch(cif_fit(d, "time", "status", group = "age"), error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "continuous, not categorical")
  expect_match(conditionMessage(err), "group_cut")
  # the same for a genuinely numeric column
  d$age_num <- as.numeric(d$age)
  expect_error(cif_fit(d, "time", "status", group = "age_num"),
               "continuous, not categorical")
})

test_that("group_cut splits a continuous column at documented quantiles", {
  d <- sim_cr(n = 200)
  expect_message(
    r <- cif_fit(d, "time", "status", group = "age", group_cut = "median",
                 times = c(1, 3)),
    "split by median at")
  expect_equal(r$group_cut, "median")
  expect_length(r$cut_points, 1L)
  expect_equal(unname(r$group_levels), c("low", "high"))
  expect_equal(sum(r$group_sizes), nrow(d))
  expect_equal(sort(unique(r$table$group)), c("high", "low"))
  expect_lt(max(r$table$cif, na.rm = TRUE), 1)

  r3 <- suppressMessages(cif_fit(d, "time", "status", group = "age",
                                group_cut = "tertile", times = c(1, 3)))
  expect_equal(unname(r3$group_levels), c("low", "middle", "high"))
  expect_length(r3$cut_points, 2L)
  r4 <- suppressMessages(cif_fit(d, "time", "status", group = "age",
                                 group_cut = "quartile", times = c(1, 3)))
  expect_equal(unname(r4$group_levels), c("Q1", "Q2", "Q3", "Q4"))
  expect_length(r4$cut_points, 3L)
})

test_that("a categorical group column keeps its own levels", {
  d <- sim_cr(n = 200)
  d$age <- as.numeric(d$age)
  r <- cif_fit(d, "time", "status", group = "sex", times = c(1, 3))
  expect_equal(r$group_cut, "none")
  expect_null(r$cut_points)
  expect_equal(sort(unname(r$group_levels)), c("female", "male"))
  expect_equal(sort(unique(r$table$group)), c("female", "male"))
  expect_equal(sum(r$group_sizes), nrow(d))
})
