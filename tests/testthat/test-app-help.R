# The app's Help page renders the HELP.md that ships with the package, and the
# examples printed there must call real package functions - not placeholders.

help_ready <- function() {
  for (p in c("shiny", "bs4Dash")) if (!requireNamespace(p, quietly = TRUE)) return(FALSE)
  f <- system.file("shiny", "CanPAS", "HELP.md", package = "CanPAS")
  nzchar(f) && file.exists(f)
}

test_that("HELP.md ships with the app and is renderable", {
  skip_if_not(help_ready(), "shiny app files unavailable")
  md <- system.file("shiny", "CanPAS", "HELP.md", package = "CanPAS")
  txt <- paste(readLines(md, warn = FALSE), collapse = "\n")
  expect_match(txt, "# CanPAS Shiny App", fixed = TRUE)
  expect_match(txt, "R examples (real functions, real cohorts)", fixed = TRUE)
  expect_gt(length(regmatches(txt, gregexpr("```r", txt, fixed = TRUE))[[1]]), 3L)
})

test_that("the help examples use real CanPAS functions and real cohorts", {
  skip_if_not(help_ready(), "shiny app files unavailable")
  md <- system.file("shiny", "CanPAS", "HELP.md", package = "CanPAS")
  txt <- paste(readLines(md, warn = FALSE), collapse = "\n")
  # the prerequisite chain and the analyses it feeds
  expect_match(txt, "dataset_info")            # the catalog itself, used as data
  for (fn in c("cohort_merged", "get_expr_data", "merge_surv_expr",
               "plot_km", "plot_roc", "COX_analysis", "COX_by_genes", "COX_by_datasets",
               "cpas_km_pooled", "cpas_meta", "cpas_meta_panel", "get_signature_value",
               "tcga_merged", "competing_risk_COX", "endpoint_resolve")) {
    expect_match(txt, paste0(fn, "("), fixed = TRUE, info = fn)
  }
  # real accessions from the catalog, and no simulated data left in the examples
  expect_match(txt, "GSE14814")
  expect_match(txt, "TCGA-LUAD")
  expect_false(grepl("rnorm(", txt, fixed = TRUE))
  expect_false(grepl("rexp(", txt, fixed = TRUE))
  # every accession quoted in the examples must exist in the catalog
  accs <- unique(regmatches(txt, gregexpr('"(GSE[0-9]+[A-Za-z0-9_]*|TCGA-[A-Z]+)"', txt))[[1]])
  accs <- gsub('"', "", accs)
  expect_gt(length(accs), 3L)
  expect_true(all(accs %in% CanPAS::dataset_info$Accession))
})

test_that("the app serves a Help page built from HELP.md", {
  skip_if_not(help_ready(), "shiny app files unavailable")
  suppressMessages({ library(shiny); library(bs4Dash) })
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  env <- new.env(parent = globalenv())
  sys.source(file.path(ad, "core.R"), envir = env)
  sys.source(file.path(ad, "mod_help.R"), envir = env)
  # the renderer turns the markdown into HTML
  p <- get(".help_md_path", envir = env)()
  expect_true(!is.null(p) && file.exists(p))
  html <- as.character(get(".md_to_html", envir = env)(p))
  expect_match(html, "cohort_merged", fixed = TRUE)
  expect_match(html, "CanPAS Shiny App", fixed = TRUE)
  # and the page itself is wired into the UI
  app_src <- system.file("shiny", "CanPAS", "app.R", package = "CanPAS")
  src <- paste(readLines(app_src, warn = FALSE), collapse = "\n")
  expect_match(src, 'tabName = "help"', fixed = TRUE)
  expect_match(src, "ui_mod_help(\"help\")", fixed = TRUE)
})
