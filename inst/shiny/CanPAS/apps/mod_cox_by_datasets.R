# Multi-datasets analysis · COX by datasets --------------------------------
# The same univariable Cox model (one gene) in several datasets: per-dataset
# estimates + FDR, forest plot above the table. Estimates are NOT pooled here;
# the Meta-analysis page does the pooling.

ui_mod_cox_by_datasets <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(width = 3,
      h5("Datasets", class = "section-title"),
      selectInput(ns("mtype"), "Cancer type", choices = NULL),
      selectInput(ns("mep"), .lab("Endpoint", "Survival endpoint family (OS/DSS/DFS/PFS/MFS). Chosen first: the dataset list below then offers only the cohorts that carry it; the concrete token is still resolved per cohort and reported with the results."), choices = NULL),
      selectizeInput(ns("mcohorts"), .lab("Datasets", "At least 1 cohort(s) that carry the endpoint family chosen above; each one contributes the token it actually has (e.g. RFS for the DFS family)."),
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Select datasets...", maxItems = 20)),
      hr(),
      h5("Marker", class = "section-title"),
      textInput(ns("mgene"), .lab("Gene", "HUGO gene symbol; a separate model is fitted in every cohort."), value = "GAPDH"),
      selectInput(ns("mrule_g"), .lab("Probe rule", "How several probes of one gene are combined into one value per sample: max (default), mean, median or min."),
                  choices = c("max", "mean", "median", "min"), selected = "max"),
      shinyWidgets::actionBttn(ns("btn_run"), "Run analysis",
                               style = "gradient", color = "primary",
                               icon = icon("play"), block = TRUE),
      br(),
      p(class = "note", "One univariable Cox model per dataset for the same gene.
        Hazard ratios are per raw expression unit and are NOT comparable across
        platforms; the FDR column adjusts over the dataset x gene combinations.
        Use the Meta-analysis page to pool the estimates.")
    ),
    mainPanel(width = 9,
      verbatimTextOutput(ns("status")),
      uiOutput(ns("overlap")),
      ui_fig_block(ns, "forest", w = 600, h = 520),
      br(), hr(),
      h5("Table (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("tbl")),
      br(),
      downloadButton(ns("dl_csv"), "Download table (CSV)", class = "btn-success")
    )
  )
}

server_mod_cox_by_datasets <- function(id, rv, dataset_info) {
  moduleServer(id, function(input, output, session) {
    output$overlap <- renderUI({
      .overlap_warning(dataset_info, input$mcohorts %||% character(0))
    })
    catalog <- .norm_catalog(dataset_info)
        # cancer type drives the dataset list; the shared cohort is preselected
    observeEvent(rv$sel$acc, {
      acc <- rv$sel$acc
      if (is.null(acc) || !acc %in% catalog$Accession) return()
      .init_type_family_cohorts(session, catalog, shared = acc)
    }, ignoreNULL = TRUE, ignoreInit = FALSE)
    ## cancer type changed: refresh the family choices for that type, then the
    ## cohort list for the family that is (still) selected
    observeEvent(input$mtype, {
      if (is.null(input$mtype) || !nzchar(input$mtype)) return()
      ch_ep <- .ep_choices_for_type(catalog, input$mtype)
      cur <- shiny::isolate(input$mep)
      sel_fam <- if (!is.null(cur) && cur %in% unname(ch_ep)) cur else unname(ch_ep)[1]
      updateSelectInput(session, "mep", choices = ch_ep, selected = sel_fam)
      r <- .refresh_cohorts(catalog, input$mtype, sel_fam,
                            selected = shiny::isolate(rv$sel$acc))
      updateSelectizeInput(session, "mcohorts", choices = r$choices, selected = r$selected)
    }, ignoreNULL = TRUE)
    ## family changed: the cohort list is restricted to cohorts that carry it
    observeEvent(input$mep, {
      if (is.null(input$mep) || !nzchar(input$mep)) return()
      r <- .refresh_cohorts(catalog, shiny::isolate(input$mtype), input$mep,
                            selected = shiny::isolate(input$mcohorts))
      if (length(r$dropped))
        showNotification(sprintf("%d selected cohort(s) do not carry %s and were removed: %s",
                                 length(r$dropped), input$mep,
                                 paste(r$dropped, collapse = ", ")),
                         type = "warning", duration = 8)
      updateSelectizeInput(session, "mcohorts", choices = r$choices, selected = r$selected)
    }, ignoreNULL = TRUE)
    # first render: fill the cancer type and dataset selectors
    session$onFlushed(function()
      .init_type_family_cohorts(session, catalog, shared = shiny::isolate(rv$sel$acc)),
      once = TRUE)

    ## A cohort change can no longer move the endpoint (the family is chosen
    ## first); it can only reveal that a cohort does not carry that family, which
    ## happens when the shared selection comes from another page.
    observeEvent(input$mcohorts, {
      accs <- input$mcohorts %||% character(0)
      if (!length(accs)) return()
      miss <- setdiff(accs, .cohorts_with_family(catalog, shiny::isolate(input$mtype),
                                                shiny::isolate(input$mep)))
      if (length(miss))
        showNotification(sprintf("%s does not carry %s and is skipped by the analysis",
                                 paste(miss, collapse = ", "), shiny::isolate(input$mep)),
                         type = "warning", duration = 8)
      cl <- tryCatch(.common_clinical(accs), error = function(e) character(0))
      updateSelectizeInput(session, "clin", choices = cl, selected = character(0))
    }, ignoreNULL = TRUE)

    run <- eventReactive(input$btn_run, {
      accs <- input$mcohorts %||% character(0)
      fam <- input$mep %||% "OS"
      gene <- input$mgene %||% "GAPDH"
      tryCatch({
        if (!length(accs)) stop("Please select at least one dataset.")
        if (!nzchar(trimws(gene))) stop("Please provide a gene symbol.")
        res <- COX_by_datasets(datasets = accs, gene = trimws(gene), type = fam,
                               process_duplicates = input$mrule_g %||% "max")
        list(res = res, accs = accs, family = fam, error = NULL, ts = Sys.time())
      }, error = function(e)
        list(res = NULL, accs = accs, family = fam, error = conditionMessage(e),
             ts = Sys.time()))
    })
    current <- reactive({ if (is.null(input$btn_run)) NULL else run() })

    table4 <- reactive({
      r <- current(); if (is.null(r) || is.null(r$res)) return(NULL)
      d <- r$res$combined_results
      data.frame(dataset = d$dataset, endpoint = d$endpoint, N = d$N,
                 HR = .fmt4(d$HR),
                 `95% CI` = paste0(.fmt4(d$HR95L), " - ", .fmt4(d$HR95H)),
                 p = .fmtp4(d$Pvalue), `FDR (BH)` = .fmtp4(d$P_adj),
                 check.names = FALSE, stringsAsFactors = FALSE)
    })

    output$status <- renderPrint({
      r <- current()
      if (is.null(r)) return(cat("Select the datasets, the endpoint family and the gene, then click 'Run analysis'."))
      if (!is.null(r$error)) return(cat("Analysis failed: ", r$error))
      md <- r$res$metadata; d <- r$res$combined_results
      cat(sprintf("Datasets requested: %d | analysed: %d | failed: %d | gene: %s | endpoint family: %s\n",
                  md$total_datasets, md$successful_datasets, md$failed_datasets,
                  md$gene_analyzed, r$family))
      cat("Endpoints used: ",
          paste(sprintf("%s=%s", d$dataset, d$endpoint), collapse = ", "), "\n", sep = "")
      if (length(r$res$errors))
        cat("Failed: ", paste(names(r$res$errors), unlist(r$res$errors), sep = ": ",
                              collapse = "; "), "\n", sep = "")
    })
    draw_forest <- function() {
      r <- current()
      if (is.null(r) || is.null(r$res)) {
        plot.new(); title(if (is.null(r)) "No analysis yet" else "Analysis failed")
        return(invisible(NULL))
      }
      d <- r$res$combined_results
      p <- forest_plot(COX_out = data.frame(
        Variates = sprintf("%s [%s]", d$dataset, d$endpoint),
        HR = d$HR, HR95L = d$HR95L, HR95H = d$HR95H, Pvalue = d$Pvalue),
        digits = 4) +
        ggplot2::labs(title = sprintf("Univariable Cox by dataset · %s · %s",
                                      r$res$metadata$gene_analyzed, r$family),
                      subtitle = "per raw expression unit (not comparable across platforms)")
      print(p)
    }
    server_fig_block(input, output, session, key = "forest", draw = draw_forest,
                     ready = function() { r <- current(); !is.null(r) && !is.null(r$res) },
                     w0 = 600, h0 = 520,
                     filename = function() paste0("COX_by_datasets_",
                       current()$res$metadata$gene_analyzed %||% "gene"))
    output$tbl <- DT::renderDT({
      d <- table4(); if (is.null(d)) return(NULL); .dtx4(d)
    })
    output$dl_csv <- downloadHandler(
      filename = function() paste0("COX_by_datasets_",
                                   current()$res$metadata$gene_analyzed %||% "gene", ".csv"),
      content = function(file) {
        d <- table4()
        if (is.null(d)) { showNotification("No result to download.", type = "warning"); return(NULL) }
        utils::write.csv(d, file, row.names = FALSE)
      })
  })
}
