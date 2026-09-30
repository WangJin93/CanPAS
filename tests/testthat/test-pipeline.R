# Tests for the shipped curation pipeline (spec B6 / hardening item A10).
#
# The package ships the curation pipeline as nine consolidated scripts under
# inst/pipeline/.  This file checks the three properties the consolidation must
# keep:
#   (1) every consolidated file parses and defines its documented runner
#       inventory;
#   (2) source()-ing the nine files has no side effects - nothing is written;
#   (3) every function definition of the original pipeline/R scripts is present
#       in the consolidated copy with an identical deparsed body (checked only
#       when the (unshipped) source tree is reachable).

.pipeline_dir <- function() {
  p <- system.file("pipeline", package = "CanPAS")
  if (length(p) == 1L && nzchar(p) && dir.exists(p)) return(p)
  p2 <- file.path("..", "..", "inst", "pipeline")     # source tree fallback
  if (dir.exists(p2)) return(normalizePath(p2))
  ""
}

.pipeline_files <- function(dir) {
  sort(list.files(dir, pattern = "^[0-9]{2}_.*\\.R$"))
}

# Documented runner inventory (see inst/pipeline/README.md)
.pipeline_inventory <- c(
  "01_ingest_parse.R"          = 6L,
  "02_platform_maps.R"         = 5L,
  "03_survival_tables.R"       = 10L,
  "04_cohort_builds.R"         = 5L,
  "05_clinical_and_scale.R"    = 9L,
  "06_catalog_and_registry.R"  = 11L,
  "07_mirror_upload.R"         = 7L,
  "08_qc_screen_and_verify.R"  = 5L,
  "09_analysis_and_figures.R"  = 9L
)

# every function definition that is not nested inside another function body
.pipeline_defs <- function(e, out = list()) {
  if (!is.call(e)) return(out)
  if ((identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
      is.symbol(e[[2]]) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) {
    out[[as.character(e[[2]])]] <- paste(deparse(e[[3]]), collapse = "\n")
    return(out)
  }
  if (identical(e[[1]], as.name("{"))) {
    for (i in seq_along(e)[-1]) out <- .pipeline_defs(e[[i]], out)
    return(out)
  }
  if (identical(e[[1]], as.name("if"))) {
    for (i in 2:3) if (length(e) >= i) out <- .pipeline_defs(e[[i]], out)
    return(out)
  }
  if (identical(e[[1]], as.name("for")))   return(.pipeline_defs(e[[4]], out))
  if (identical(e[[1]], as.name("while"))) return(.pipeline_defs(e[[2]], out))
  out
}

.pipeline_defs_of <- function(text) {
  ee <- parse(text = text)
  out <- list()
  for (i in seq_along(ee)) {
    d <- .pipeline_defs(ee[[i]])
    for (nm in names(d)) out[[nm]] <- d[[nm]]
  }
  out
}

.pipeline_snapshot <- function(paths) {
  out <- character(0)
  for (r in paths) {
    if (!dir.exists(r)) next
    fs <- list.files(r, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    if (!length(fs)) next
    out <- c(out, sprintf("%s|%s|%s", fs, file.size(fs),
                          format(file.mtime(fs), "%Y-%m-%d %H:%M:%OS6")))
  }
  sort(out)
}

# the runner bodies of one consolidated file, keyed by runner name
.pipeline_bodies <- function(path) {
  e <- parse(path)
  out <- list()
  for (i in seq_along(e)) {
    ex <- e[[i]]
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.symbol(ex[[2]]) &&
        grepl("^run_", as.character(ex[[2]]))) {
      out[[sub("^run_", "", as.character(ex[[2]]))]] <- paste(deparse(ex[[3]][[3]]), collapse = "\n")
    }
  }
  out
}

test_that("the nine consolidated pipeline files are shipped and parse", {
  dir <- .pipeline_dir()
  skip_if(!nzchar(dir), "inst/pipeline/ is not available")
  files <- .pipeline_files(dir)
  expect_equal(files, names(.pipeline_inventory))
  for (f in files) {
    expect_error(parse(file.path(dir, f)), NA, info = f)
  }
})

test_that("every consolidated file defines its documented runner inventory", {
  dir <- .pipeline_dir()
  skip_if(!nzchar(dir), "inst/pipeline/ is not available")
  for (f in names(.pipeline_inventory)) {
    bodies <- .pipeline_bodies(file.path(dir, f))
    expect_equal(length(bodies), unname(.pipeline_inventory[[f]]), info = f)
    expect_true(all(grepl("^run_", paste0("run_", names(bodies)))), info = f)
    expect_true(all(nzchar(unlist(bodies))), info = f)
  }
})

test_that("sourcing the nine consolidated files has no side effects", {
  dir <- .pipeline_dir()
  skip_if(!nzchar(dir), "inst/pipeline/ is not available")
  # a fresh, empty directory used as the working directory, so a script that
  # wrote a relative path would leave a mark here
  sandbox <- file.path(tempdir(), "cpas-pipeline-source-smoke")
  dir.create(sandbox, showWarnings = FALSE, recursive = TRUE)
  leftover <- list.files(sandbox, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (length(leftover)) unlink(leftover, recursive = TRUE)
  watched <- c(sandbox, dir)
  root <- Sys.getenv("CPAS_DATA_ROOT", unset = "")
  if (nzchar(root) && dir.exists(file.path(root, "pipeline", "out")))
    watched <- c(watched, file.path(root, "pipeline", "out"))
  old <- setwd(sandbox)
  on.exit(setwd(old), add = TRUE)
  before <- .pipeline_snapshot(watched)
  env <- new.env(parent = globalenv())
  for (f in .pipeline_files(dir)) sys.source(file.path(dir, f), envir = env)
  after <- .pipeline_snapshot(watched)
  runners <- ls(env, pattern = "^run_")
  expect_equal(length(runners), 67L)
  expect_true(all(vapply(runners, function(n) is.function(get(n, envir = env)), logical(1))))
  expect_identical(after, before)
  expect_length(list.files(sandbox, recursive = TRUE, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("the original function bodies survive the consolidation unchanged", {
  dir <- .pipeline_dir()
  skip_if(!nzchar(dir), "inst/pipeline/ is not available")
  root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
  src <- file.path(root, "pipeline", "R")
  skip_if(!dir.exists(src), "the unshipped pipeline/R/ source tree is not available")

  checked <- 0L
  mismatches <- character(0)
  for (f in .pipeline_files(dir)) {
    bodies <- .pipeline_bodies(file.path(dir, f))
    for (nm in names(bodies)) {
      cand <- file.path(src, paste0(nm, ".R"))
      if (!file.exists(cand)) next
      a <- .pipeline_defs_of(paste(readLines(cand, warn = FALSE), collapse = "\n"))
      b <- .pipeline_defs_of(bodies[[nm]])
      expect_identical(names(a), names(b), info = nm)
      for (k in names(a)) {
        checked <- checked + 1L
        if (!identical(a[[k]], b[[k]])) mismatches <- c(mismatches, paste0(nm, "::", k))
      }
    }
  }
  expect_gt(checked, 0L)
  expect_identical(mismatches, character(0))
})
