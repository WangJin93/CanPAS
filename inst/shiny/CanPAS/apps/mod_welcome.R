# Dashboard module --------------------------------------------------------
ui_mod_welcome <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(8,
      h2("Cancer Prognosis Analysis Suite (CanPAS)"),
      h5("Integrated survival analysis of GEO, CGGA and TCGA cancer cohorts"),
      hr(),
      h4("What this app does", class = "section-title"),
      p("CanPAS combines a curated GEO survival database (MySQL mirror), the CGGA
        glioma cohorts and TCGA cohorts (local clinical tables + on-demand
        UCSC Xena expression) behind one unified analysis workflow."),
      tags$ul(
        tags$li(tags$b("Datasets: "), "browse cohorts by cancer type and survival
                endpoint, inspect sample size and platform, and push a selected
                cohort to KM or COX analysis."),
        tags$li(tags$b("Single dataset analysis: "), "KM analysis (curves for a
                gene or signature), COX analysis (univariable + multivariable
                model restricted to the significant covariates, with clinical
                diagnostics) and COX by genes (one model per gene with FDR)."),
        tags$li(tags$b("Multi-datasets analysis: "), "COX by datasets (the same
                gene across datasets), Pooled KM (IPD or time-point pooling) and
                Meta-analysis (one gene, or a per-gene panel with FDR)."),
        tags$li(tags$b("Methods: "), "how the data were processed (parsing,
                cleaning, standardisation, mirroring) and how survival endpoints
                are grouped into families and mapped to each cohort.")
      ),
      h4("Common conventions", class = "section-title"),
      p("Survival times are expressed in years and status is coded 0/1;
        multiple probes per gene are collapsed with max/mean/median/min;
        signatures follow the form 0.5*GAPDH + 0.5*TNS1. Datasets and endpoints
        selected on any page are shared with the other pages.")
    ),
    column(4,
      bs4Dash::bs4Card(
        title = "Database at a glance",
        width = 12,
        uiOutput(ns("stats"))
      )
    )
  )
}

server_mod_welcome <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$stats <- renderUI({
      di <- .norm_catalog(dataset_info)
      n_ds   <- nrow(di)
      n_type <- length(unique(di$Type))
      fc     <- .family_counts(di)
      fams   <- .family_order[fc > 0]
      tagList(
        bs4Dash::bs4ValueBox(value = n_ds, subtitle = "cohorts in catalog",
                             color = "primary", width = 12),
        bs4Dash::bs4ValueBox(value = n_type, subtitle = "cancer types",
                             color = "warning", width = 12),
        bs4Dash::bs4ValueBox(value = length(fams), subtitle = "survival endpoint families",
                             color = "danger", width = 12),
        p(class = "note", "Endpoint families (cohorts): ",
          paste(sprintf("%s %d", fams, fc[fams]), collapse = " · "),
          br(), "Analyses are run per family; the token each cohort contributes
          (e.g. RFS for DFS, PFI for PFS) is resolved and reported automatically.")
      )
    })
  })
}
