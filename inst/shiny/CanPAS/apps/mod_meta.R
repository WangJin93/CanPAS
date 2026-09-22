# Multi-datasets analysis · Meta-analysis -------------------------------
# Two-stage meta-analysis (cpas_meta) for one gene/signature, plus a per-gene
# panel mode (cpas_meta_panel) that pools every gene of a list and adjusts the
# gene-level p-values. Forest plot above the tables throughout.

ui_mod_meta <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(width = 3,
      h5("Datasets", class = "section-title"),
      selectInput(ns("mtype"), "Cancer type", choices = NULL),
      selectInput(ns("mep"), .lab("Endpoint", "Survival endpoint family (OS/DSS/DFS/PFS/MFS). Chosen first: the dataset list below then offers only the cohorts that carry it; the concrete token is still resolved per cohort and reported with the results."), choices = NULL),
      selectizeInput(ns("mcohorts"), .lab("Datasets", "At least 2 cohort(s) that carry the endpoint family chosen above; each one contributes the token it actually has (e.g. RFS for the DFS family)."),
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Select datasets...", maxItems = 20)),
      hr(),
      h5("Mode", class = "section-title"),
      radioButtons(ns("mode"), NULL,
                   choices = c("Single gene / signature" = "single",
                               "Gene panel (per-gene meta)" = "panel"),
                   selected = "single"),
      conditionalPanel(condition = "input.mode == 'single'", ns = ns,
        radioButtons(ns("mkind"), .lab("Marker", "A single gene, or a weighted signature such as 0.5*GAPDH + 0.5*TNS1."),
                     choices = c("Gene" = "gene", "Signature" = "signature"),
                     selected = "gene"),
        conditionalPanel(condition = "input.mkind == 'gene'", ns = ns,
          textInput(ns("mgene"), .lab("Gene", "HUGO gene symbol; its probes are collapsed by the rule below."), value = "TP53"),
          selectInput(ns("mrule_g"), .lab("Probe rule", "How several probes of one gene are combined into one value per sample: max (default), mean, median or min."),
                      choices = c("max", "mean", "median", "min"), selected = "max")),
        conditionalPanel(condition = "input.mkind == 'signature'", ns = ns,
          textInput(ns("msig"), .lab("Signature", "Weighted signature formula, e.g. 0.5*GAPDH + 0.5*TNS1."), value = "0.5*GAPDH + 0.5*TNS1"),
          selectInput(ns("mrule"), .lab("Probe rule", "Each gene of the signature is first collapsed over its probes by this rule, then weighted and summed."),
                      choices = c("max", "mean", "median", "min"), selected = "max"),
          checkboxInput(ns("mallow"), .lab("Allow missing genes", "On: genes this platform cannot measure are skipped and the score is computed on the remaining ones."), value = FALSE))),
      conditionalPanel(condition = "input.mode == 'panel'", ns = ns,
        textInput(ns("mgenes"), .lab("Genes", "Comma-separated gene list; every gene is collapsed on its own and the results are written into one table."),
                  value = "TP53, GAPDH, TNS1, PTEN")),
      hr(),
      h5("Pooling", class = "section-title"),
      selectInput(ns("mmet"), .lab("Pooling", "RE = random effects (default), FE = fixed effect."), choices = c("RE", "FE"), selected = "RE"),
      numericInput(ns("mmin"), .lab("Min events", "Cohorts with fewer events than this are left out of the pooling, so a single tiny cohort cannot drive the result."), value = 5, min = 1),
      selectizeInput(ns("mclin"), .lab("Adjust for", "Optional: adjust for these covariates inside each cohort first, then pool the adjusted hazard ratios (the column must exist in that cohort)."),
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "none (unadjusted)")),
      shinyWidgets::actionBttn(ns("btn_meta"), "Run meta-analysis",
                               style = "gradient", color = "primary",
                               icon = icon("chart-bar"), block = TRUE),
      br(),
      p(class = "note", "Stage 1: one Cox model per dataset (marker standardised
        per SD within the dataset). Stage 2: inverse-variance pooling with
        heterogeneity statistics, a prediction interval and leave-one-out
        sensitivity analysis.")
    ),
    mainPanel(width = 9,
      verbatimTextOutput(ns("meta_status")),
      uiOutput(ns("overlap")),
      ui_fig_block(ns, "meta_forest", w = 600, h = 500),
      br(), hr(),
      h5("Pooled estimate (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("meta_pool")),
      br(),
      h5("Per-dataset estimates (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("meta_per")),
      br(),
      h5("Leave-one-out sensitivity (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("meta_loo")),
      br(),
      downloadButton(ns("dl_meta_csv"), "Download per-dataset table (CSV)",
                     class = "btn-success")
    )
  )
}

server_mod_meta <- function(id, rv, dataset_info) {
  moduleServer(id, function(input, output, session) {
    output$overlap <- renderUI({
      .overlap_warning(dataset_info, input$mcohorts %||% character(0))
    })
    catalog <- .norm_catalog(dataset_info)
    loc <- reactiveValues(meta = NULL, panel = NULL)

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
      if (length(accs) < 2) return()
      miss <- setdiff(accs, .cohorts_with_family(catalog, shiny::isolate(input$mtype),
                                                shiny::isolate(input$mep)))
      if (length(miss))
        showNotification(sprintf("%s does not carry %s and is skipped by the analysis",
                                 paste(miss, collapse = ", "), shiny::isolate(input$mep)),
                         type = "warning", duration = 8)
      cl <- tryCatch(.common_clinical(accs), error = function(e) character(0))
      updateSelectizeInput(session, "mclin", choices = cl, selected = character(0))
    }, ignoreNULL = TRUE)

    observeEvent(input$btn_meta, {
      accs <- input$mcohorts %||% character(0)
      fam <- input$mep %||% "OS"
      conf <- if (length(input$mclin)) input$mclin else NULL
      if (identical(input$mode %||% "single", "panel")) {
        genes <- trimws(unlist(strsplit(input$mgenes %||% "", "[,;[:space:]]+")))
        genes <- genes[nzchar(genes)]
        loc$panel <- tryCatch({
          if (length(accs) < 2) stop("Please select at least two datasets.")
          if (!length(genes)) stop("Please provide at least one gene symbol.")
          res <- cpas_meta_panel(datasets = accs, genes = genes, type = fam,
                                 method = input$mmet %||% "RE",
                                 confounders = conf, min_events = input$mmin %||% 5)
          list(res = res, accs = accs, family = fam, error = NULL)
        }, error = function(e) list(error = conditionMessage(e)))
        if (!is.null(loc$panel$error))
          showNotification(paste0("Per-gene meta panel failed: ", loc$panel$error),
                           type = "error", duration = 10)
        return(invisible(NULL))
      }
      spec <- .marker_spec(input)
      loc$meta <- tryCatch({
        b <- .multi_bundle(session, catalog, accs, fam, spec, min_n = 2)
        res <- cpas_meta(datasets = names(b$merged), marker = "marker", type = b$family,
                         method = input$mmet %||% "RE", confounders = conf,
                         min_events = input$mmin %||% 5, merged = b$merged)
        list(res = res, failed = b$failed, family = b$family, spec = spec)
      }, error = function(e) list(error = conditionMessage(e)))
      if (!is.null(loc$meta$error))
        showNotification(paste0("Meta-analysis failed: ", loc$meta$error),
                         type = "error", duration = 10)
    })

    panel_mode <- reactive(identical(input$mode %||% "single", "panel"))

    output$meta_status <- renderPrint({
      if (panel_mode()) {
        r <- loc$panel
        if (is.null(r)) return(cat("Provide the gene list, then click 'Run meta-analysis'."))
        if (!is.null(r$error)) return(cat("Failed: ", r$error))
        md <- r$res$metadata
        cat(sprintf("Per-gene meta panel | datasets %d/%d retrieved | genes %d | pooled %d | single-dataset %d | failed %d\n",
                    md$datasets_retrieved, md$datasets_requested, md$genes_requested,
                    md$genes_pooled, md$genes_single_dataset, md$genes_failed))
        cat(sprintf("FDR (BH) < 0.05: %d gene(s)\n", md$n_significant_fdr005))
        if (length(r$res$fetch_errors))
          cat("Retrieval failures: ",
              paste(names(r$res$fetch_errors), collapse = ", "), "\n", sep = "")
        if (length(r$res$errors))
          cat("Not pooled: ", substr(paste(names(r$res$errors), unlist(r$res$errors),
                                           sep = ": ", collapse = "; "), 1, 400), "\n", sep = "")
        return(invisible(NULL))
      }
      r <- loc$meta
      if (is.null(r)) return(cat("Select >= 2 datasets, the marker and the pooling method, then click 'Run meta-analysis'."))
      if (!is.null(r$error)) return(cat("Failed: ", r$error))
      po <- r$res$pooled
      cat(sprintf("Cohorts included: %d | endpoint %s | pooling %s | marker: %s\n",
                  nrow(r$res$per_dataset), r$res$input$type, r$res$input$method,
                  if (identical(r$spec$kind, "gene")) paste0("gene ", r$spec$gene)
                  else paste0("signature ", r$spec$sig)))
      cat("Endpoints used: ",
          paste(sprintf("%s=%s", r$res$per_dataset$dataset, r$res$per_dataset$endpoint),
                collapse = ", "), "\n", sep = "")
      cat(sprintf("Pooled HR = %s [%s, %s] | p = %s | I2 = %s | tau2 = %s\n",
                  .fmt4(po$HR), .fmt4(po$lower), .fmt4(po$upper), .fmtp4(po$p),
                  .fmt4(po$I2), .fmt4(po$tau2)))
      if (is.finite(po$pi_lower %||% NA_real_))
        cat(sprintf("95%% prediction interval for a new dataset: [%s, %s]\n",
                    .fmt4(po$pi_lower), .fmt4(po$pi_upper)))
      if (length(r$failed))
        cat("Excluded: ", paste(names(r$failed), r$failed, sep = ": ",
                                collapse = "; "), "\n", sep = "")
      if (length(r$res$errors))
        cat("Dropped by the meta step: ", substr(paste(names(r$res$errors),
            unlist(r$res$errors), sep = ": ", collapse = "; "), 1, 400), "\n", sep = "")
    })

    draw_forest <- function() {
      if (panel_mode()) {
        r <- loc$panel
        if (is.null(r) || !is.null(r$error)) {
          plot.new(); title(if (is.null(r)) "No analysis yet" else "Analysis failed")
          return(invisible(NULL))
        }
        print(plot_meta_panel(r$res, digits = 4)); return(invisible(NULL))
      }
      r <- loc$meta
      if (is.null(r) || !is.null(r$error)) {
        plot.new(); title(if (is.null(r)) "No analysis yet" else "Analysis failed")
        return(invisible(NULL))
      }
      print(plot_meta_forest(r$res, digits = 4))
    }
    server_fig_block(input, output, session, key = "meta_forest", draw = draw_forest,
                     ready = function() {
                       if (panel_mode()) { r <- loc$panel; return(!is.null(r) && is.null(r$error) && !is.null(r$res)) }
                       r <- loc$meta; !is.null(r) && is.null(r$error) && !is.null(r$res)
                     },
                     w0 = 600, h0 = 500,
                     filename = function() if (panel_mode()) "CPAS_meta_panel_forest" else "CPAS_meta_forest")

    output$meta_pool <- DT::renderDT({
      if (panel_mode()) {
        r <- loc$panel
        if (is.null(r) || !is.null(r$error)) return(NULL)
        d <- r$res$table
        return(.dtx4(data.frame(gene = d$gene, k = d$k, N = d$total_n,
                                events = d$total_events, HR = .fmt4(d$HR),
                                `95% CI` = paste0(.fmt4(d$lower), " - ", .fmt4(d$upper)),
                                p = .fmtp4(d$p), `FDR (BH)` = .fmtp4(d$P_adj),
                                I2 = .fmt4(d$I2),
                                `95% PI` = paste0(.fmt4(d$pi_lower), " - ", .fmt4(d$pi_upper)),
                                check.names = FALSE, stringsAsFactors = FALSE)))
      }
      r <- loc$meta
      if (is.null(r) || !is.null(r$error)) return(NULL)
      po <- r$res$pooled
      .dtx4(data.frame(method = po$method, k = po$k, N = po$total_n,
                       events = po$total_events, HR = .fmt4(po$HR),
                       lower = .fmt4(po$lower), upper = .fmt4(po$upper),
                       p = .fmtp4(po$p), Q = .fmt4(po$Q), df = po$df,
                       p_het = .fmtp4(po$p_heterogeneity), I2 = .fmt4(po$I2),
                       tau2 = .fmt4(po$tau2),
                       `95% PI lower` = .fmt4(po$pi_lower),
                       `95% PI upper` = .fmt4(po$pi_upper),
                       check.names = FALSE, stringsAsFactors = FALSE))
    })

    output$meta_per <- DT::renderDT({
      if (panel_mode()) {
        r <- loc$panel
        if (is.null(r) || !is.null(r$error) || is.null(r$res$per_dataset)) return(NULL)
        d <- r$res$per_dataset
        return(.dtx4(data.frame(gene = d$gene, dataset = d$dataset,
                                endpoint = d$endpoint, N = d$n, events = d$events,
                                HR = .fmt4(d$HR),
                                `95% CI` = paste0(.fmt4(d$lower), " - ", .fmt4(d$upper)),
                                p = .fmtp4(d$p), check.names = FALSE,
                                stringsAsFactors = FALSE)))
      }
      r <- loc$meta
      if (is.null(r) || !is.null(r$error)) return(NULL)
      d <- r$res$per_dataset
      .dtx4(data.frame(dataset = d$dataset, endpoint = d$endpoint, N = d$n,
                       events = d$events, HR = .fmt4(d$HR),
                       `95% CI` = paste0(.fmt4(d$lower), " - ", .fmt4(d$upper)),
                       p = .fmtp4(d$p), logHR = .fmt4(d$logHR), SE = .fmt4(d$se),
                       check.names = FALSE, stringsAsFactors = FALSE))
    })

    output$meta_loo <- DT::renderDT({
      if (panel_mode()) return(NULL)
      r <- loc$meta
      if (is.null(r) || !is.null(r$error)) return(NULL)
      lo <- tryCatch(loo_meta(r$res), error = function(e) NULL)
      if (is.null(lo) || !nrow(lo)) return(NULL)
      .dtx4(data.frame(left_out = lo$left_out, k = lo$k, HR = .fmt4(lo$HR),
                       lower = .fmt4(lo$lower), upper = .fmt4(lo$upper),
                       p = .fmtp4(lo$p), I2 = .fmt4(lo$I2),
                       check.names = FALSE, stringsAsFactors = FALSE))
    })

    output$dl_meta_csv <- downloadHandler(
      filename = function() "CPAS_meta_per_dataset.csv",
      content = function(file) {
        d <- if (panel_mode()) loc$panel$res$per_dataset else loc$meta$res$per_dataset
        if (is.null(d) || !nrow(d)) { showNotification("No result to download.", type = "warning"); return(NULL) }
        utils::write.csv(.round4(d), file, row.names = FALSE)
      })
  })
}
