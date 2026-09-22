# The Shiny app renders every figure through core.R's figure block, so that the
# size controls sit under the figure, resizing is live, and the downloaded file
# has exactly the displayed size. These tests pin that contract.

fig_block_ready <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) return(FALSE)
  if (!requireNamespace("bs4Dash", quietly = TRUE)) return(FALSE)
  if (!requireNamespace("shinycssloaders", quietly = TRUE)) return(FALSE)
  f <- system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS")
  nzchar(f) && file.exists(f)
}

png_size <- function(path) {
  b <- readBin(path, "raw", 24L)
  c(w = sum(as.integer(b[17:20]) * c(256^3, 256^2, 256, 1)),
    h = sum(as.integer(b[21:24]) * c(256^3, 256^2, 256, 1)))
}

test_that(".fig_dims returns the control values and falls back on bad input", {
  skip_if_not(fig_block_ready(), "shiny app files unavailable")
  env <- new.env(parent = globalenv())
  sys.source(system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS"), envir = env)
  fd <- get(".fig_dims", envir = env)
  expect_equal(unname(unlist(fd(list(k_width = 1150, k_height = 680), "k", 900, 520))), c(1150, 680))
  expect_equal(unname(unlist(fd(list(k_width = NA, k_height = 10), "k", 900, 520))), c(900, 520))
  expect_equal(unname(unlist(fd(list(), "k", 800, 500))), c(800, 500))
})

test_that("the downloaded figure keeps the displayed size", {
  skip_if_not(fig_block_ready(), "shiny app files unavailable")
  env <- new.env(parent = globalenv())
  sys.source(system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS"), envir = env)
  fd <- get(".fig_dims", envir = env)
  for (wh in list(c(900, 520), c(1150, 680), c(400, 300))) {
    d <- fd(list(k_width = wh[1], k_height = wh[2]), "k", 900, 520)
    png_f <- tempfile(fileext = ".png")
    grDevices::png(png_f, width = d$w, height = d$h, units = "px", res = 96)
    dev_px <- grDevices::dev.size("px"); graphics::plot(1:3); grDevices::dev.off()
    expect_equal(unname(png_size(png_f)), as.numeric(wh))      # PNG pixels = displayed pixels
    expect_equal(round(unname(dev_px)), as.numeric(wh))        # the device opened at that size

    pdf_f <- tempfile(fileext = ".pdf")
    grDevices::pdf(pdf_f, width = d$w / 96, height = d$h / 96)
    dev_in <- grDevices::dev.size("in"); graphics::plot(1:3); grDevices::dev.off()
    expect_equal(unname(dev_in), c(d$w / 96, d$h / 96), tolerance = 1e-3)  # same physical size
    expect_gt(file.size(pdf_f), 0)
  }
})

test_that("every analysis page carries the figure controls and drops the old PDF panel", {
  skip_if_not(fig_block_ready(), "shiny app files unavailable")
  # the app files assume app.R has attached these
  suppressMessages({ library(shiny); library(bs4Dash); library(DT) })
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  env <- new.env(parent = globalenv())
  for (f in c("core.R", "mod_welcome.R", "mod_datasets.R", "mod_km.R", "mod_cox.R",
              "mod_cox_by_genes.R", "mod_cox_by_datasets.R", "mod_pooled_km.R",
              "mod_meta.R", "mod_methods.R"))
    sys.source(file.path(ad, f), envir = env)
  pages <- list(km = "ui_mod_km", cox = "ui_mod_cox", genes = "ui_mod_cox_by_genes",
                bycohort = "ui_mod_cox_by_datasets", pooled = "ui_mod_pooled_km",
                meta = "ui_mod_meta")
  for (nm in names(pages)) {
    h <- paste(as.character(get(pages[[nm]], envir = env)(nm)), collapse = "")
    expect_true(grepl("Figure width (px)", h, fixed = TRUE))
    expect_true(grepl("Figure height (px)", h, fixed = TRUE))
    expect_true(grepl("Download PNG", h, fixed = TRUE))
    expect_true(grepl("Download PDF", h, fixed = TRUE))
    expect_false(grepl("PDF width (cm)", h, fixed = TRUE))
    expect_false(grepl("Download figure (PDF)", h, fixed = TRUE))
  }
})
