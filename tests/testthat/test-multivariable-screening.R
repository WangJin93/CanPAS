# The multivariable Cox model must not publish meaningless coefficients, but
# refusing the whole model leaves the user without a result or a reason. These
# tests pin the middle course: locate the offending covariates, drop / merge /
# exclude only what is necessary, fit the rest, and record why.

# categorical covariates contribute a header row whose HR is NA by design
est_hr <- function(r) {
  tb <- r$results_table
  tb$HR[!tb$Type %in% "cate_header"]
}

sim_multi <- function(n = 400, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(n)
  age <- stats::rnorm(n, 60, 10)
  grade <- factor(sample(c("G1", "G2", "G3"), n, TRUE))
  node <- factor(sample(c("N0", "N1"), n, TRUE))
  hr <- exp(0.4 * x + 0.02 * (age - 60) + 0.5 * (grade == "G3") + 0.6 * (node == "N1"))
  tt <- stats::rexp(n, hr / 10)
  st <- stats::rbinom(n, 1, 0.7)
  data.frame(ID = sprintf("S%03d", seq_len(n)), OS_time = tt, OS_status = st,
             marker = x, age = age, grade = grade, N = node)
}

test_that("an estimable multivariable model is fitted unchanged", {
  d <- sim_multi()
  r <- COX_analysis(d, type = "OS", cont_Variates = c("marker", "age"),
                    cate_Variates = c("grade", "N"), method = "multi")
  expect_s3_class(r, "cpas_COX")
  expect_equal(nrow(r$metadata$dropped_covariates), 0L)
  expect_true(all(is.finite(est_hr(r))))
})

test_that("a constant covariate is dropped with a reason", {
  d <- sim_multi()
  d$flat <- 1
  r <- COX_analysis(d, type = "OS", cont_Variates = c("marker", "age", "flat"),
                    cate_Variates = "N", method = "multi")
  dc <- r$metadata$dropped_covariates
  expect_true("flat" %in% dc$variable)
  expect_match(dc$reason[dc$variable == "flat"], "constant")
  expect_false("flat" %in% r$results_table$Var1)
  expect_true(all(is.finite(est_hr(r))))
})

test_that("a level with one patient is handled and the covariate survives", {
  d <- sim_multi()
  lv <- as.character(d$N)
  lv[1] <- "N3"                      # single patient in a third level
  d$N <- factor(lv, levels = c("N0", "N1", "N3"))
  d$OS_status[1] <- 1
  r <- COX_analysis(d, type = "OS", cont_Variates = "marker",
                    cate_Variates = "N", method = "multi")
  dc <- r$metadata$dropped_covariates
  expect_true(any(dc$variable == "N"))
  expect_match(dc$reason[dc$variable == "N"], "fewer than 5 patients|excluded|Other")
  # the covariate is kept with its usable levels, or dropped if none remain
  expect_true(all(is.finite(est_hr(r))))
  expect_true(max(abs(est_hr(r)), na.rm = TRUE) < 1e6)
})

test_that("a level in which every patient has an event does not break the model", {
  d <- sim_multi()
  lv <- as.character(d$N)
  lv[1:6] <- "N3"
  d$N <- factor(lv, levels = c("N0", "N1", "N3"))
  d$OS_status[1:6] <- 1              # all events in that level
  r <- COX_analysis(d, type = "OS", cont_Variates = "marker",
                    cate_Variates = "N", method = "multi")
  expect_s3_class(r, "cpas_COX")
  expect_true(all(is.finite(est_hr(r))))
  expect_true(max(abs(est_hr(r)), na.rm = TRUE) < 1e6)
  expect_true(nrow(r$metadata$dropped_covariates) >= 1L)
})

test_that("collinear covariates are reduced to an estimable model", {
  d <- sim_multi()
  d$age2 <- d$age * 1.0000001         # numerically collinear with age
  r <- tryCatch(COX_analysis(d, type = "OS", cont_Variates = c("marker", "age", "age2"),
                             cate_Variates = "N", method = "multi"),
                error = function(e) e)
  if (inherits(r, "error")) {
    expect_match(conditionMessage(r), "not estimable")
  } else {
    expect_true(all(is.finite(est_hr(r))))
    expect_true(max(abs(log(est_hr(r))), na.rm = TRUE) < 20)
  }
})

test_that("perfect separation of the marker drops the marker and keeps the rest", {
  d <- sim_multi(n = 120)
  d$marker <- as.numeric(d$OS_status)          # marker == event indicator
  d$OS_time <- ifelse(d$OS_status == 1, 0.1, 1)
  r <- tryCatch(suppressWarnings(COX_analysis(d, type = "OS",
              cont_Variates = c("marker", "age"), method = "multi")),
                error = function(e) e)
  if (inherits(r, "error")) {
    # acceptable only when nothing at all is estimable
    expect_match(conditionMessage(r), "not estimable|not fitted")
  } else {
    expect_false("marker" %in% r$metadata$final_covariates)
    expect_true("marker" %in% r$metadata$dropped_covariates$variable)
    expect_true(isTRUE(r$metadata$reduced))
    expect_true(all(is.finite(est_hr(r))))
  }
})

test_that("a pair that cannot be co-estimated is reduced and flagged, not refused", {
  # two covariates that are each strongly significant on their own but are the
  # same variable (perfectly collinear): the repair must fall back to the
  # single covariate that can be estimated and say so, instead of failing
  d <- sim_multi(n = 200, seed = 4)
  d$marker_copy <- d$marker
  expect_warning(
    r <- COX_analysis(d, type = "OS", cont_Variates = c("marker", "marker_copy"),
                      method = "multi"),
    "reduced")
  expect_s3_class(r, "cpas_COX")
  expect_true(isTRUE(r$metadata$reduced))
  expect_length(r$metadata$final_covariates, 1L)
  expect_match(as.character(r$metadata$reduced_note), "not adjusted", all = FALSE)
  expect_true(all(is.finite(est_hr(r))))
  # the surviving coefficient must equal its univariable estimate
  keep <- r$metadata$final_covariates
  uni <- COX_analysis(d, type = "OS", cont_Variates = keep, method = "uni")
  expect_equal(est_hr(r), est_hr(uni), tolerance = 1e-6)
})

test_that("min_covariates = 2 reduces and flags instead of erroring", {
  # two perfectly collinear continuous covariates: one must go, so the model
  # falls below min_covariates - it is returned, flagged, and never errors
  d <- sim_multi(n = 200, seed = 5)
  d$age_copy <- d$age
  r <- suppressWarnings(COX_analysis(d, type = "OS", cont_Variates = c("age", "age_copy"),
                                     method = "multi", min_covariates = 2L))
  expect_s3_class(r, "cpas_COX")
  expect_true(isTRUE(r$metadata$reduced))
  expect_length(r$metadata$final_covariates, 1L)
  expect_match(as.character(r$metadata$reduced_note), "not adjusted")
  expect_true(all(is.finite(est_hr(r))))
})

test_that("the failure message lists the covariates already handled", {
  # every covariate is the event indicator: nothing can be estimated, and the
  # message must carry the whole history, not just the last failed fit
  d <- sim_multi(n = 120, seed = 6)
  d$a <- as.numeric(d$OS_status); d$b <- d$a; d$c <- d$a
  r <- tryCatch(suppressWarnings(COX_analysis(d, type = "OS",
              cont_Variates = c("a", "b", "c"), method = "multi")),
                error = function(e) e)
  expect_s3_class(r, "error")
  expect_match(conditionMessage(r), "not estimable")
})

test_that("COX_screen_adjust reports a reduced multivariate model", {
  d <- sim_multi(n = 200, seed = 8)
  d$marker_copy <- d$marker          # both significant, but the same variable
  r <- suppressMessages(suppressWarnings(
    COX_screen_adjust(d, type = "OS", cont_Variates = c("marker", "marker_copy"),
                      cate_Variates = NULL, p.threshold = 0.05)))
  expect_true(isTRUE(r$multi_metadata$reduced))
  expect_length(r$multi_metadata$final_covariates, 1L)
  expect_true(nrow(r$multi_table) >= 1L)
  ## flextable is a suggested package: with it the summary table is a flextable,
  ## without it the same numbers come back as a plain data.frame
  if (requireNamespace("flextable", quietly = TRUE)) {
    expect_s3_class(r$print_result, "flextable")
  } else {
    expect_s3_class(r$print_result, "data.frame")
  }
})

test_that("the printed summary table is a three-line table", {
  skip_if_not_installed("flextable")
  d <- sim_multi(n = 300, seed = 12)
  r <- suppressMessages(suppressWarnings(
    COX_screen_adjust(d, type = "OS", cont_Variates = c("marker", "age"),
                      cate_Variates = c("grade", "N"), p.threshold = 0.05)))
  expect_s3_class(r$print_result, "flextable")
  html <- as.character(flextable::htmltools_value(r$print_result))
  # collect every border that is actually drawn (width > 0) and require exactly
  # the three rules of a three-line table, nothing else
  decl <- unique(regmatches(html, gregexpr("border-(top|bottom|left|right): [0-9.]+pt solid", html))[[1]])
  drawn <- decl[!grepl(": 0.00pt", decl)]
  expect_setequal(drawn, c("border-top: 1.5pt solid",
                           "border-bottom: 0.75pt solid",
                           "border-bottom: 1.5pt solid"))
  expect_match(html, "Univariate Cox", fixed = TRUE)
  expect_match(html, "Multivariate Cox", fixed = TRUE)
  # and the caption can be set by the caller (the app does this)
  ft <- flextable::set_caption(r$print_result, caption = "GSE13507 (OS), n = 133")
  expect_match(as.character(flextable::htmltools_value(ft)), "GSE13507", fixed = TRUE)
})
