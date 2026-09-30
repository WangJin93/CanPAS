# Launch the CanPAS Shiny application

Starts the Shiny app shipped inside the CanPAS package. The app
implements most of the package: GEO/TCGA cohort loading, single- and
multi-variable Cox analysis, Kaplan-Meier, time-dependent ROC,
cross-cohort meta-analysis with pooled Kaplan-Meier, weighted
gene-signature scoring and TCGA clinical exploration (see `HELP.md`
inside the app).

## Usage

``` r
run_cpas_app(host = "127.0.0.1", port = NULL, launch.browser = TRUE, ...)
```

## Arguments

- host:

  Host to bind, default `"127.0.0.1"`.

- port:

  Port to listen on (default: a free port chosen by Shiny).

- launch.browser:

  Launch the default browser (`TRUE`) or not.

- ...:

  Further arguments passed to
  [`shiny::runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

Starts a Shiny server (blocks until the app is stopped).

## Details

Requires the suggested packages `shiny`, `DT`, `bs4Dash`, `shinyWidgets`
and `shinycssloaders`. GEO queries need access to the CanPAS public API;
TCGA expression queries additionally need `UCSCXenaShiny`. Cohort data
are fetched from the CanPAS public API at analysis time, so the analysis
pages need network access; only the Dashboard, Datasets and Methods
pages work offline.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## The app ships inside the package; the Help page documents every page.
   ## Launching it needs the suggested packages (shiny, bs4Dash, DT, ...).
   run_cpas_app()                 # free port, opens the browser
   run_cpas_app(port = 3838)      # fixed port
} # }
```
