# Datasets module: browse the catalog and push a cohort to KM / COX ----------
ui_mod_datasets <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(width = 3,
      h5("Filters"),
      selectInput(ns("category"), .lab("Endpoint", "Survival endpoint family (OS/DSS/DFS/PFS/MFS). Chosen first: the table below then lists only the cohorts that carry that endpoint, and a column shows the token each one contributes."),
                  choices = c("Any" = "__all__"), selected = "__all__"),
      selectInput(ns("type"), "Cancer type", choices = NULL),
      hr(),
      h5("How to use"),
      p("Select a row in the table to make that cohort the shared selection,",
        br(), "then use the buttons below:",
        br(), "· KM plot — jump to KM analysis",
        br(), "· COX analysis — jump to COX analysis",
        br(), "· Source page — NCBI GEO, CGGA or GDC Portal in a new tab.
        ", br(), br(),
        tags$b("Sources: "), "GEO (curated mirror), EMBL-EBI (ArrayExpress/BioStudies),
        CGGA (3 glioma cohorts), TCGA (on-demand UCSC Xena expression + local
        clinical tables) and two cBioPortal-hosted studies.",
        br(), br(),
        tags$b("Workflow: "), "pick the endpoint family first — the table then lists
        only the cohorts that carry it, and a column shows the token each one
        contributes at that endpoint (e.g. RFS for GSE31210 at DFS). Narrow by
        cancer type afterwards if you wish. 'PFS or MFS' combines the two families.",
        br(), br(),
        tags$b("Note: "), "flags cohorts with shared patients (same study on another
        platform, or the same series deposited twice) — pooling those counts patients
        twice. The multi-dataset pages show the same warning for the current selection.",
        br(), br(),
        tags$b("N and Events: "), "N is the analysable sample count — samples with
        expression data and a non-missing time and status for the cohort's primary
        endpoint — and Events is the number of those samples with an event. '\u2014'
        means the cohort has no expression table in the mirror, so gene-level
        analysis is not possible.")
    ),
    mainPanel(width = 9,
      h4("Cohort catalog", class = "section-title"),
      DT::DTOutput(ns("table")),
      br(),
      fluidRow(
        bs4Dash::column(4, shinyWidgets::actionBttn(ns("km_btn"), "KM plot",
          style = "gradient", color = "primary", icon = icon("chart-line"),
          block = TRUE)),
        bs4Dash::column(4, shinyWidgets::actionBttn(ns("cox_btn"), "COX analysis",
          style = "gradient", color = "success", icon = icon("cogs"),
          block = TRUE)),
        bs4Dash::column(4, shinyWidgets::actionBttn(ns("geo_btn"), "Source page",
          style = "gradient", color = "warning", icon = icon("external-link-alt"),
          block = TRUE))
      ),
      br(),
      uiOutput(ns("overlap_sel")),
      verbatimTextOutput(ns("selmsg"))
    )
  )
}

server_mod_datasets <- function(id, rv, on_jump, dataset_info) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    catalog <- .norm_catalog(dataset_info)
    fam_n <- reactive(.family_counts(catalog))

    observe({
      updateSelectInput(session, "type", choices = sort(unique(catalog$Type)))
      n <- fam_n()
      fams <- .family_order[n > 0]
      ch <- c(stats::setNames("__all__", "Any"))
      if (length(fams))
        ch <- c(ch, stats::setNames(fams, sprintf("%s (%d)", fams, n[fams])))
      if (all(c("PFS", "MFS") %in% fams))
        ch <- c(ch, stats::setNames("__pfsmfs__",
                 sprintf("PFS or MFS (%d)", n[["PFS"]] + n[["MFS"]])))
      updateSelectInput(session, "category", choices = ch, selected = "__all__")
    })

    # A cohort is shown for a family when the catalog annotation says so; when
    # the family annotation is missing the raw endpoint tokens are mapped to
    # families instead, so a family filter never hides an eligible cohort and
    # never shows a cohort that has no such endpoint at all.
    # "can this cohort answer this family" is defined once in core.R, the same
    # definition the multi-dataset pages now use, so filter and analysis agree
    family_hit <- function(d, fam) .family_hit(d, fam)

    df_sub <- reactive({
      d <- catalog
      if (!is.null(input$type) && nzchar(input$type))
        d <- d[d$Type == input$type, , drop = FALSE]
      if (!is.null(input$category) && nzchar(input$category) &&
          input$category != "__all__")
        d <- d[family_hit(d, input$category), , drop = FALSE]
      d
    })

    output$table <- DT::renderDT({
      d <- df_sub()
      show <- data.frame(
        Accession   = as.character(d$Accession),
        "Cancer type" = as.character(d$Type),
        "Endpoint families" = vapply(as.character(d$Accession),
                                     function(a) .family_display(catalog, a), character(1)),
        Endpoints   = as.character(d$SurvivalTypes),
        Platform    = as.character(d$GPL),
        "N (analysable)" = .dash(d$N),
        Events      = .dash(d$n_events),
        Note         = ifelse(is.na(d$Note), "", as.character(d$Note)),
        Source      = .cohort_source(as.character(d$Accession)),
        Link        = .cohort_link(as.character(d$Accession)),
        check.names = FALSE, stringsAsFactors = FALSE)
      ## when a single family is selected, the table answers "what token does
      ## this cohort contribute at that endpoint?" without the user having to
      ## read it out of the "Endpoint families" column
      fam <- input$category
      if (!is.null(fam) && nzchar(fam) && !identical(fam, "__all__")) {
        toks <- vapply(as.character(d$Accession),
                       function(a) .family_token_label(catalog, a, fam), character(1))
        col <- sprintf("Token for %s", if (identical(fam, "__pfsmfs__")) "PFS/MFS" else fam)
        show <- cbind(show[, 1:3, drop = FALSE],
                      stats::setNames(data.frame(ifelse(nzchar(toks), toks, "\u2014"),
                                                 stringsAsFactors = FALSE), col),
                      show[, 4:ncol(show), drop = FALSE])
      }
      DT::datatable(show, rownames = FALSE, escape = FALSE,
                    selection = "single", filter = "top",
                    options = list(pageLength = 12, scrollX = TRUE, dom = "ftp"))
    })

    sel_acc <- reactive({
      s <- input$table_rows_selected
      if (is.null(s) || !length(s)) return(NULL)
      df_sub()$Accession[s]
    })

    # selecting a row already updates the shared selection (modules inherit it)
    observeEvent(input$table_rows_selected, {
      acc <- sel_acc(); if (is.null(acc)) return()
      .set_shared_selection(rv, catalog, acc)
    }, ignoreNULL = TRUE)

    # note of the selected cohort (overlap, expression-only, ...)
    output$overlap_sel <- renderUI({
      row <- catalog[catalog$Accession == sel_acc(), , drop = FALSE]
      if (!nrow(row) || is.na(row$Note) || !nzchar(as.character(row$Note))) return(NULL)
      div(class = "alert alert-warning", style = "padding:8px; margin-top:8px;",
          tags$b("Note: "), as.character(row$Note))
    })

    output$selmsg <- renderPrint({
      acc <- sel_acc()
      if (is.null(acc)) return(cat("(no cohort selected — click a row above)"))
      row <- catalog[catalog$Accession == acc, , drop = FALSE]
      cat(sprintf("Selected: %s | cancer type: %s | families: %s | endpoints: %s | platform: %s",
                  acc, row$Type[1],
                  ifelse(is.na(as.character(row$EndpointFamilies[1])), "-",
                         as.character(row$EndpointFamilies[1])),
                  ifelse(is.na(as.character(row$SurvivalTypes[1])), "-",
                         as.character(row$SurvivalTypes[1])),
                  row$GPL[1]))
    })

    observeEvent(input$km_btn, {
      acc <- sel_acc()
      if (is.null(acc)) return(showNotification("Please select a cohort first", type = "error"))
      .set_shared_selection(rv, catalog, acc)
      on_jump("sd_km")
    })
    observeEvent(input$cox_btn, {
      acc <- sel_acc()
      if (is.null(acc)) return(showNotification("Please select a cohort first", type = "error"))
      .set_shared_selection(rv, catalog, acc)
      on_jump("sd_cox")
    })
    observeEvent(input$geo_btn, {
      acc <- sel_acc()
      if (is.null(acc)) return(showNotification("Please select a cohort first", type = "error"))
      url <- if (startsWith(acc, "TCGA-"))
        paste0("https://portal.gdc.cancer.gov/projects/", acc)
      else paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", acc)
      session$sendCustomMessage("open_url", url)
    })
  })
}
