# The COX page fits a univariable model for every covariate and then a
# multivariable model on those that reached the p threshold. When that
# multivariable model is not estimable, the page must still show the
# univariable figure and table, and it must state why - a bare "Analysis
# failed" that hides the plot is what these tests prevent.

cox_app_ready <- function() {
  for (p in c("shiny", "bs4Dash", "DT")) if (!requireNamespace(p, quietly = TRUE)) return(FALSE)
  f <- system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS")
  nzchar(f) && file.exists(f)
}

# synthetic cohort: survival depends on 'marker', and 'age2' is a duplicate of
# 'marker' - both are univariably significant, but they cannot be estimated
# together, which is exactly the situation that used to abort the whole page
cox_stub_df <- function(n = 200, seed = 42) {
  set.seed(seed)
  x <- stats::rnorm(n)
  age <- stats::rnorm(n, 60, 10)
  st <- stats::rbinom(n, 1, 0.7)
  tt <- stats::rexp(n, exp(0.9 * x - 0.02 * (age - 60)) / 12)
  data.frame(ID = sprintf("P%03d", seq_len(n)), OS_time = tt, OS_status = st,
             marker = x, age = age, age2 = x, stringsAsFactors = FALSE)
}

run_cox_app <- function(expr, df = cox_stub_df(), clin = c("age", "age2")) {
  skip_if_not(cox_app_ready(), "shiny app files unavailable")
  # the assertions must be evaluated inside testServer's environment (where
  # input/output/session live), not in this function's frame
  expr_q <- substitute(expr)
  suppressMessages({ library(shiny); library(bs4Dash); library(DT) })
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  env <- new.env(parent = globalenv())
  for (f in c("core.R", "mod_welcome.R", "mod_datasets.R", "mod_km.R", "mod_cox.R",
              "mod_cox_by_genes.R", "mod_cox_by_datasets.R", "mod_pooled_km.R",
              "mod_meta.R", "mod_methods.R"))
    sys.source(file.path(ad, f), envir = env)
  # data access is stubbed: no network, deterministic cohort
  env$.build_marker_df <- function(acc, ep, kind = "gene", gene = NULL, ref = NULL,
                                   formula_txt = NULL, rule = "max", allow_missing = FALSE) df
  env$.loc_clin <- function(acc) data.frame(ID = df$ID, age = df$age, age2 = df$age2)
  rv <- shiny::reactiveValues(sel = list(acc = "TESTCOHORT", type = "Test", family = "OS",
                                         ep = "OS"))
  testServer(get("server_mod_cox", envir = env),
             args = list(id = "cox", rv = rv, dataset_info = CanPAS::dataset_info), {
    session$setInputs(kind = "gene", gene = "MARKER", ref = "__all__", rule = "max",
                      show_all = FALSE, clin = clin, p_thr = 0.05)
    session$setInputs(go = 1)
    session$setInputs(go = 2)
    eval(expr_q, envir = environment())
  })
}

plain <- function(x) gsub("<[^>]*>", " ", paste(as.character(x), collapse = " "))

test_that("a multivariable model reduced to one covariate still plots and explains", {
  run_cox_app({
    st <- plain(output$status)
    expect_match(st, "reduced", ignore.case = TRUE)          # flagged in the status line
    expect_match(st, "NOT adjusted for confounding", fixed = TRUE)
    note <- plain(output$multi_note)
    expect_match(note, "NOT adjusted for confounding", fixed = TRUE)
    expect_match(note, "could not be co-estimated")
    # both figures and both tables are still produced
    expect_false(is.null(output$forest_uni))
    expect_false(is.null(output$forest_multi))
    expect_false(is.null(output$tbl_uni))
    expect_false(is.null(output$tbl_multi))
    # the reason is attached to the covariate that had to go
    expect_match(note, "age2")
  })
})

test_that("a covariate that cannot be estimated leaves the univariable output intact", {
  d <- cox_stub_df(seed = 7)
  d$marker <- as.numeric(d$OS_status)      # marker == event indicator: separation
  d$age2 <- d$marker
  run_cox_app({
    st <- plain(output$status)
    # the univariable table survives even though the marker cannot be estimated
    expect_false(is.null(output$tbl_uni))
    expect_match(st, "skipped|Only one covariate|No covariate reached")
  }, df = d)
})

# --- the Print results tab ------------------------------------------------
# The COX page carries a second tab with the combined univariable +
# multivariable summary that COX_screen_adjust() formats as a three-line table.

test_that("the COX page exposes a Print results tab", {
  skip_if_not(cox_app_ready(), "shiny app files unavailable")
  suppressMessages({ library(shiny); library(bs4Dash) })
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  env <- new.env(parent = globalenv())
  for (f in c("core.R", "mod_km.R", "mod_cox.R", "mod_help.R"))
    sys.source(file.path(ad, f), envir = env)
  h <- paste(as.character(get("ui_mod_cox", envir = env)("cox")), collapse = "")
  expect_match(h, "Print results", fixed = TRUE)
  expect_match(h, "Download table (Word)", fixed = TRUE)
})

test_that("the Print results tab renders a three-line summary table", {
  skip_if_not(cox_app_ready(), "shiny app files unavailable")
  suppressMessages({ library(shiny); library(bs4Dash); library(DT) })
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  env <- new.env(parent = globalenv())
  for (f in c("core.R", "mod_welcome.R", "mod_datasets.R", "mod_km.R", "mod_cox.R",
              "mod_cox_by_genes.R", "mod_cox_by_datasets.R", "mod_pooled_km.R",
              "mod_meta.R", "mod_methods.R", "mod_help.R"))
    sys.source(file.path(ad, f), envir = env)
  set.seed(42); n <- 220
  x <- stats::rnorm(n); age <- stats::rnorm(n, 62, 9)
  st <- stats::rbinom(n, 1, 0.7)
  tt <- stats::rexp(n, exp(0.9 * x + 0.05 * (age - 62)) / 12)   # both significant
  df <- data.frame(ID = sprintf("P%03d", seq_len(n)), OS_time = tt, OS_status = st,
                   marker = x, age = age,
                   grade = factor(sample(c("G1", "G2", "G3"), n, TRUE)),
                   N = factor(sample(c("N0", "N1"), n, TRUE)), stringsAsFactors = FALSE)
  env$.build_marker_df <- function(...) df
  env$.loc_clin <- function(acc) data.frame(ID = df$ID, age = df$age, grade = df$grade, N = df$N)
  rv <- shiny::reactiveValues(sel = list(acc = "TESTCOHORT", type = "Test", family = "OS", ep = "OS"))
  testServer(get("server_mod_cox", envir = env),
             args = list(id = "cox", rv = rv, dataset_info = CanPAS::dataset_info), {
    session$setInputs(kind = "gene", gene = "MARKER", ref = "__all__", rule = "max",
                      show_all = FALSE, clin = c("age", "grade", "N"), p_thr = 0.05)
    session$setInputs(go = 1)
    session$setInputs(go = 2)
    h <- paste(as.character(output$print_table), collapse = "")
    expect_gt(nchar(h), 1000L)
    # a three-line table: top rule, rule under the header, bottom rule
    expect_match(h, "border-top: 1.5pt solid", fixed = TRUE)
    expect_match(h, "border-bottom: 0.75pt solid", fixed = TRUE)
    expect_match(h, "border-bottom: 1.5pt solid", fixed = TRUE)
    # both model groups and the caption are present
    expect_match(h, "Univariate Cox", fixed = TRUE)
    expect_match(h, "Multivariate Cox", fixed = TRUE)
    expect_match(h, "TESTCOHORT", fixed = TRUE)
  })
})

# --- row ordering differs by page -----------------------------------------
# The COX page must keep the model order: sorting by HR splits a clinical
# covariate into rows scattered across the plot, so its levels no longer sit
# under their own name. The per-gene and multi-dataset pages have one row per
# gene/cohort, where sorting by HR is what is wanted.

test_that("the COX page plots in model order, the per-gene pages sort by HR", {
  skip_if_not(cox_app_ready(), "shiny app files unavailable")
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  cox <- readLines(file.path(ad, "mod_cox.R"), warn = FALSE)
  # both forest plots of this page ask for the unsorted model order
  expect_equal(sum(grepl('HR_order = "none"', cox, fixed = TRUE)), 2L)
  # and the other two result pages leave the default (sorted by HR) in place
  for (f in c("mod_cox_by_genes.R", "mod_cox_by_datasets.R")) {
    src <- readLines(file.path(ad, f), warn = FALSE)
    expect_false(any(grepl("HR_order", src)), info = f)
  }
})
