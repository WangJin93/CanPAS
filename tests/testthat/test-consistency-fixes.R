# The three consistency fixes of spec B10 ------------------------------------

test_that("B10a: cpas_meta_panel() defaults to the cpas_meta method, REML", {
  mk <- function(seed, beta) {
    set.seed(seed)
    n <- 140
    x <- stats::rnorm(n)
    data.frame(ID = sprintf("S%04d", seq_len(n)),
               OS_time = round(stats::rexp(n, 0.2) * exp(-beta * x), 3),
               OS_status = stats::rbinom(n, 1, 0.6),
               G1 = x, G2 = stats::rnorm(n), stringsAsFactors = FALSE)
  }
  # deliberately heterogeneous per-cohort effects, so tau^2 > 0 and the
  # estimator choice is visible
  m <- list(A = mk(1, -0.9), B = mk(2, -0.1), C = mk(3, 0.7))
  p <- cpas_meta_panel(c("A", "B", "C"), genes = c("G1", "G2"), type = "OS",
                       merged = m, progress = FALSE)
  expect_identical(p$input$method, "REML")
  expect_identical(cpas_manifest(p)$meta_method, "REML")
  # the two entry points now agree on their default
  d <- cpas_meta(c("A", "B", "C"), "G1", "OS", merged = m)
  expect_identical(d$input$method, "REML")
  expect_identical(p$input$method, d$input$method)
  # explicitly asking for DL still works and is not the default
  dl <- cpas_meta_panel(c("A", "B", "C"), genes = "G1", type = "OS",
                        method = "DL", merged = m, progress = FALSE)
  expect_identical(dl$input$method, "DL")
  expect_identical(cpas_manifest(dl)$meta_method, "DL")
  expect_false(isTRUE(all.equal(dl$table$tau2[p$table$gene == "G1"],
                                p$table$tau2[p$table$gene == "G1"])))
  expect_false(isTRUE(all.equal(dl$table$HR[p$table$gene == "G1"],
                                p$table$HR[p$table$gene == "G1"])))
  expect_error(cpas_meta_panel(c("A", "B"), genes = "G1", merged = m,
                               method = "nope", progress = FALSE),
               "should be one of")
})

test_that("B10a: the man page documents REML as the panel default", {
  db <- tryCatch(tools::Rd_db("CanPAS"), error = function(e) NULL)
  skip_if(is.null(db) || !"cpas_meta_panel.Rd" %in% names(db),
          "the installed Rd database is unavailable")
  txt <- paste(as.character(db[["cpas_meta_panel.Rd"]]), collapse = "\n")
  expect_match(txt, "REML")
  expect_false(grepl("DL.*default kept for continuity", txt))
  # and the pi_method argument is documented there too
  expect_match(txt, "pi_method")
})

test_that("B10b: the not-estimable result carries a top-level auto_repair", {
  set.seed(21)
  n <- 200
  x <- stats::rnorm(n)
  d <- data.frame(ID = sprintf("S%04d", seq_len(n)),
                  OS_time = stats::rexp(n, exp(0.4 * x) / 10),
                  OS_status = stats::rbinom(n, 1, 0.7),
                  marker = x, flat = 1, stringsAsFactors = FALSE)
  r <- suppressMessages(COX_analysis(d, type = "OS",
                                     cont_Variates = c("marker", "flat"),
                                     method = "multi"))
  expect_s3_class(r, "cpas_not_estimable")
  expect_identical(r$status, "not estimable")
  expect_false(is.null(r$auto_repair))
  expect_false(isTRUE(r$auto_repair))
  # consistent with the manifest, which previously carried it only in the notes
  expect_identical(r$auto_repair, r$manifest$auto_repair)
  expect_identical(r$auto_repair, r$input_params$auto_repair)
  expect_true(any(grepl("auto_repair = FALSE", r$manifest$notes)))
  df <- as.data.frame(cpas_manifest(r))
  expect_true("auto_repair" %in% df$field)
  expect_identical(df$value[df$field == "auto_repair"], "FALSE")
})

test_that("B10b: the manifest carries auto_repair for every analysis", {
  mk <- function(seed) {
    set.seed(seed)
    n <- 120
    x <- stats::rnorm(n)
    data.frame(ID = sprintf("S%04d", seq_len(n)),
               OS_time = round(stats::rexp(n, 0.2) * exp(-0.4 * x), 3),
               OS_status = stats::rbinom(n, 1, 0.6),
               G1 = x, stringsAsFactors = FALSE)
  }
  m <- list(A = mk(1), B = mk(2), C = mk(3))
  d <- cpas_meta(c("A", "B", "C"), "G1", "OS", merged = m)
  expect_false(isTRUE(cpas_manifest(d)$auto_repair))
  expect_identical(cpas_manifest(d)$auto_repair, d$input$auto_repair)
  expect_identical(cpas_manifest(d)$auto_repair, FALSE)
  p <- cpas_meta_panel(c("A", "B", "C"), genes = "G1", type = "OS",
                       merged = m, progress = FALSE)
  expect_identical(cpas_manifest(p)$auto_repair, FALSE)
  expect_error(cpas_meta(c("A", "B"), "G1", "OS", merged = m, auto_repair = NA),
               "'auto_repair' must be TRUE or FALSE")
})

test_that("B10c: the pipeline README states the verified 491-definition figure", {
  rd <- system.file("pipeline", "README.md", package = "CanPAS")
  expect_true(nzchar(rd))
  expect_true(file.exists(rd))
  txt <- paste(readLines(rd, warn = FALSE), collapse = "\n")
  expect_match(txt, "512 one-level function definitions, 0\\s+altered")
  expect_match(txt, "Counting rule")
  expect_false(grepl("270 definitions", txt))
    expect_match(txt, "534")             # the any-depth count, for contrast
})

test_that("B10c: the 491 figure is reproducible from the shipped scripts", {
  # Re-derive the count the README states, from the installed package alone:
  # one-level = a `name <- function` binding at the top level of a runner body
  # or one nesting level inside it.
  steps <- cpas_pipeline_steps()
  is_fun <- function(x) is.call(x) && is.symbol(x[[1]]) &&
    identical(x[[1]], as.name("function"))
  is_assign_fun <- function(x) is.call(x) && is.symbol(x[[1]]) &&
    (identical(x[[1]], as.name("<-")) || identical(x[[1]], as.name("="))) &&
    length(x) == 3 && is.symbol(x[[2]]) && is_fun(x[[3]])
  walk <- function(e, depth, acc) {
    if (is.call(e)) {
      if (depth %in% c(1L, 3L) && is_assign_fun(e)) acc <- acc + 1L
      for (i in seq_along(e)) acc <- walk(e[[i]], depth + 1L, acc)
    } else if (is.expression(e) || is.pairlist(e)) {
      for (i in seq_along(e)) acc <- walk(e[[i]], depth + 1L, acc)
    }
    acc
  }
  total <- 0L
  files <- system.file("pipeline", unique(steps$consolidated_file), package = "CanPAS")
  for (f in files) {
    if (!nzchar(f) || !file.exists(f)) next
    e <- parse(f, keep.source = FALSE)
    for (ex in e) {
      if (is.call(ex) && is.symbol(ex[[1]]) && identical(as.character(ex[[1]]), "<-") &&
          is.symbol(ex[[2]]) && grepl("^run_", as.character(ex[[2]])) && is_fun(ex[[3]]))
        total <- walk(as.expression(as.list(ex[[3]][[3]])), 0L, total)
    }
  }
  expect_equal(total, 512L)
})
