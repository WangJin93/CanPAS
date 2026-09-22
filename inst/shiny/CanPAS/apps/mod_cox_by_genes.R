# Single dataset analysis · COX by genes ----------------------------------
# One univariable Cox model per gene in ONE dataset (COX_by_genes()), shown as a
# forest plot above the table, with Benjamini-Hochberg FDR across the genes.

ui_mod_cox_by_genes <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(width = 3,
      h5("Dataset", class = "section-title"),
      selectInput(ns("ctype"), "Cancer type", choices = NULL),
      selectizeInput(ns("ds"), "Dataset", choices = NULL,
                     options = list(placeholder = "Select a dataset...", maxOptions = 30)),
      selectInput(ns("ep"), .lab("Endpoint", "Survival endpoint family (OS/DSS/DFS/PFS/MFS); the concrete token is resolved per cohort and reported with the results."), choices = NULL),
      hr(),
      h5("Genes", class = "section-title"),
      textInput(ns("genes"), .lab("Genes", "Comma-separated gene list; every gene is collapsed on its own and the results are written into one table."),
                value = "GAPDH, TP53, TNS1, PTEN"),
      selectInput(ns("rule"), .lab("Probe rule", "How several probes of one gene are combined into one value per sample: max (default), mean, median or min."),
                  choices = c("max", "mean", "median", "min"), selected = "max"),
      shinyWidgets::actionBttn(ns("go"), "Run analysis",
                               style = "gradient", color = "primary",
                               icon = icon("play"), block = TRUE),
      br(),
      p(class = "note", "Each gene is fitted separately (univariable Cox) and the
        p-values are adjusted for multiple testing across the panel
        (Benjamini-Hochberg). For a joint adjusted model use the COX analysis page.")
    ),
    mainPanel(width = 9,
      verbatimTextOutput(ns("status")),
      ui_fig_block(ns, "forest", w = 600, h = 520),
      br(), hr(),
      h5("Table (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("tbl")),
      br(),
      downloadButton(ns("dl_csv"), "Download table (CSV)", class = "btn-success")
    )
  )
}

server_mod_cox_by_genes <- function(id, rv, dataset_info) {
  moduleServer(id, function(input, output, session) {
    catalog <- .norm_catalog(dataset_info)
    init_acc <- shiny::isolate(rv$sel$acc)
    init_fam <- shiny::isolate(rv$sel$family)
    # cancer type drives the dataset list (first render + later changes)
    session$onFlushed(function()
      .init_single_selectors(session, catalog, shared = shiny::isolate(rv$sel$acc)),
      once = TRUE)
    observeEvent(input$ctype, {
      if (is.null(input$ctype) || !nzchar(input$ctype)) return()
      ds <- .dataset_choices_of(catalog, input$ctype)
      accs <- unlist(unname(ds), use.names = FALSE)
      sh <- intersect(shiny::isolate(rv$sel$acc) %||% character(0), accs)
      updateSelectizeInput(session, "ds", choices = ds,
                           selected = if (length(sh)) sh else accs[1], server = TRUE)
    }, ignoreNULL = TRUE)
    if (!is.null(init_acc)) {
      eps <- .endpoint_choices(catalog, init_acc)
      updateSelectInput(session, "ep", choices = eps,
                        selected = if (!is.null(init_fam) && init_fam %in% eps)
                          init_fam else eps[1])
    }
    observeEvent(rv$sel$acc, {
      if (is.null(rv$sel$acc)) return()
      updateSelectizeInput(session, "ds", selected = rv$sel$acc)
      eps <- .endpoint_choices(catalog, rv$sel$acc)
      if (!is.null(rv$sel$family) && rv$sel$family %in% eps)
        updateSelectInput(session, "ep", choices = eps, selected = rv$sel$family)
    }, ignoreNULL = TRUE)
    observeEvent(input$ds, {
      if (is.null(input$ds) || !nzchar(input$ds)) return()
      if (!identical(rv$sel$acc, input$ds)) {
        .set_shared_selection(rv, catalog, input$ds)
        eps <- .endpoint_choices(catalog, input$ds)
        updateSelectInput(session, "ep", choices = eps,
                          selected = if (!is.null(rv$sel$ep) && rv$sel$ep %in% eps)
                            rv$sel$ep else eps[1])
      }
    }, ignoreInit = TRUE)
    observeEvent(input$ep, {
      if (is.null(input$ep) || !nzchar(input$ep)) return()
      if (!identical(rv$sel$family, input$ep)) {
        rv$sel$family <- input$ep
        rv$sel$ep <- .resolve_family(catalog, rv$sel$acc, input$ep)
      }
    }, ignoreInit = TRUE)

    run <- eventReactive(input$go, {
      acc <- shiny::isolate(rv$sel$acc)
      family <- shiny::isolate(rv$sel$family) %||% "OS"
      ep <- shiny::isolate(rv$sel$ep) %||% .resolve_family(catalog, acc, family)
      genes <- trimws(unlist(strsplit(input$genes %||% "", "[,;[:space:]]+")))
      genes <- genes[nzchar(genes)]
      tryCatch({
        if (is.null(acc)) stop("Please select a dataset first.")
        if (!length(genes)) stop("Please provide at least one gene symbol.")
        if (is.na(ep)) stop("Cohort ", acc, " has no ", family, " endpoint.")
        e <- get_expr_data(acc, genes, process_duplicates = input$rule %||% "max")
        df <- merge_surv_expr(acc, e)$merged_data
        gcols <- intersect(genes, colnames(df))
        if (!length(gcols))
          stop("None of the requested genes is measurable in this dataset.")
        res <- COX_by_genes(df, type = family, genes = gcols)
        list(res = res, acc = acc, family = family, ep = ep, rule = input$rule %||% "max",
             genes = gcols, error = NULL, ts = Sys.time())
      }, error = function(e)
        list(res = NULL, acc = acc, family = family, ep = ep, error = conditionMessage(e),
             ts = Sys.time()))
    })
    current <- reactive({ if (is.null(input$go)) NULL else run() })

    # presentation: 4 decimals everywhere (HR, CI, p)
    table4 <- reactive({
      r <- current(); if (is.null(r) || is.null(r$res)) return(NULL)
      tb <- r$res$results_table
      data.frame(gene = tb$gene,
                 HR = .fmt4(tb$HR), `95% CI` = paste0(.fmt4(tb$HR95L), " - ", .fmt4(tb$HR95H)),
                 p = .fmtp4(tb$Pvalue), `FDR (BH)` = .fmtp4(tb$P_adj),
                 N = r$res$metadata$sample_size, events = r$res$metadata$events,
                 check.names = FALSE, stringsAsFactors = FALSE)
    })

    output$status <- renderPrint({
      r <- current()
      if (is.null(r)) return(cat("Choose a dataset, the endpoint and the gene list, then click 'Run analysis'."))
      if (!is.null(r$error)) return(cat("Analysis failed: ", r$error))
      m <- r$res$metadata
      cat(sprintf("Dataset %s | endpoint %s | n = %d, events = %d\n",
                  r$acc, .endpoint_status_label(r$family, r$ep),
                  m$sample_size, m$events))
      cat(sprintf("Genes tested: %d (%s) | probe rule: %s\n",
                  m$gene_count, paste(m$genes_used %||% r$genes, collapse = ", "), r$rule))
      cat(sprintf("FDR < 0.05: %d gene(s)\n", m$n_significant_fdr005))
      if (length(m$failed_genes))
        cat("Not estimable: ", paste(m$failed_genes, collapse = ", "), "\n", sep = "")
    })

    draw_forest <- function() {
      r <- current()
      if (is.null(r) || is.null(r$res)) {
        plot.new(); title(if (is.null(r)) "No analysis yet" else "Analysis failed")
        return(invisible(NULL))
      }
      tb <- r$res$results_table
      d <- data.frame(Variates = tb$gene, HR = tb$HR, HR95L = tb$HR95L,
                      HR95H = tb$HR95H, Pvalue = tb$Pvalue, stringsAsFactors = FALSE)
      print(forest_plot(COX_out = d, digits = 4) +
              ggplot2::labs(title = sprintf("Univariable Cox by gene · %s · %s",
                                            r$acc, .endpoint_status_label(r$family, r$ep))))
    }
    server_fig_block(input, output, session, key = "forest", draw = draw_forest,
                     ready = function() { r <- current(); !is.null(r) && !is.null(r$res) },
                     w0 = 600, h0 = 520,
                     filename = function() paste0("COX_by_genes_",
                                                  current()$acc %||% "dataset"))
    output$tbl <- DT::renderDT({
      d <- table4(); if (is.null(d)) return(NULL); .dtx4(d)
    })
    output$dl_csv <- downloadHandler(
      filename = function() paste0("COX_by_genes_",
                                   current()$acc %||% "dataset", ".csv"),
      content = function(file) {
        d <- table4()
        if (is.null(d)) { showNotification("No result to download.", type = "warning"); return(NULL) }
        utils::write.csv(d, file, row.names = FALSE)
      })
  })
}
