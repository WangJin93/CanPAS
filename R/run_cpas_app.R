#' @title Launch the CanPAS Shiny application
#' @description
#' Starts the Shiny app shipped inside the CanPAS package. The app implements
#' most of the package: GEO/TCGA cohort loading, single- and multi-variable
#' Cox analysis, Kaplan-Meier, time-dependent ROC, cross-cohort meta-analysis
#' with pooled Kaplan-Meier, weighted gene-signature scoring and TCGA clinical
#' exploration (see \code{HELP.md} inside the app).
#' @param host Host to bind, default \code{"127.0.0.1"}.
#' @param port Port to listen on (default: a free port chosen by Shiny).
#' @param launch.browser Launch the default browser (\code{TRUE}) or not.
#' @param ... Further arguments passed to \code{shiny::runApp}.
#' @return Starts a Shiny server (blocks until the app is stopped).
#' @details
#' Requires the suggested packages \code{shiny}, \code{DT}, \code{bs4Dash},
#' \code{shinyWidgets} and \code{shinycssloaders}. GEO queries need access to the
#' CanPAS public API;
#' TCGA expression queries additionally need \code{UCSCXenaShiny}. Cohort data are
#' fetched from the CanPAS public API at analysis time, so the analysis pages need
#' network access; only the Dashboard, Datasets and Methods pages work offline.
#' @examples
#' \dontrun{
#'    ## The app ships inside the package; the Help page documents every page.
#'    ## Launching it needs the suggested packages (shiny, bs4Dash, DT, ...).
#'    run_cpas_app()                 # free port, opens the browser
#'    run_cpas_app(port = 3838)      # fixed port
#' }
#' @export
run_cpas_app <- function(host = "127.0.0.1", port = NULL,
                         launch.browser = TRUE, ...) {
  needed <- c("shiny", "DT", "bs4Dash", "shinyWidgets", "shinycssloaders")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("The CanPAS Shiny app requires packages: ",
         paste(miss, collapse = ", "),
         ". Install with install.packages(c('shiny','DT','bs4Dash','shinyWidgets','shinycssloaders')).")
  app_dir <- system.file("shiny", "CanPAS", package = "CanPAS")
  if (!nzchar(app_dir) || !file.exists(file.path(app_dir, "app.R")))
    stop("The bundled Shiny app is missing; please re-install the CanPAS package.")
  shiny::runApp(appDir = app_dir, host = host, port = port,
                launch.browser = launch.browser, ...)
}
