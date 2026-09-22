# KM analysis module ------------------------------------------------------
ui_mod_km <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(width = 3,
      h5("Cohort", class = "section-title"),
      selectInput(ns("ctype"), "Cancer type", choices = NULL),
      selectizeInput(ns("ds"), "Dataset", choices = NULL,
                     options = list(placeholder = "Select a cohort...",
                                    maxOptions = 30)),
      selectInput(ns("ep"), .lab("Endpoint", "Survival endpoint family (OS/DSS/DFS/PFS/MFS); this cohort contributes the token it actually has, e.g. RFS for DFS."), choices = NULL),
      hr(),
      radioButtons(ns("kind"), .lab("Marker", "A single gene (probes are resolved for the platform) or a weighted signature such as 0.5*GAPDH + 0.5*TNS1."),
                   choices = c("Gene" = "gene", "Signature" = "signature"),
                   selected = "gene"),
      conditionalPanel(condition = "input.kind == 'gene'", ns = ns,
        textInput(ns("gene"), .lab("Gene", "HUGO gene symbol, e.g. GAPDH. Multiple probes of the chosen platform are resolved for it."), value = "GAPDH"),
        selectizeInput(ns("ref"), .lab("Probe", "Leave 'Auto: all probes' to use every probe of the gene (collapsed per sample by the rule below); pick one probe id to use that probe alone."),
                       choices = NULL,
                       options = list(placeholder = "Auto: all probes -> max",
                                      maxOptions = 500)),
        uiOutput(ns("ref_hint"))),
      conditionalPanel(condition = "input.kind == 'signature'", ns = ns,
        textInput(ns("sig"), .lab("Signature", "Weighted sum of genes, e.g. 0.5*GAPDH + 0.5*TNS1. Each gene is first collapsed over its probes."),
                  value = "0.5*GAPDH + 0.5*TNS1"),
        selectInput(ns("rule"), .lab("Probe rule", "How several probes of one gene are combined into one value per sample: max (default), mean, median or min."),
                    choices = c("max", "mean", "median", "min"),
                    selected = "max"),
        checkboxInput(ns("allow_missing"), .lab("Allow missing genes", "On: genes that this platform cannot measure are skipped and the signature is scored on the rest. Off: such a gene stops the analysis."),
                      value = FALSE)),
      hr(),
      radioButtons(ns("cutmode"), .lab("Cut point", "How the high/low groups are defined: median split (default), the cut point with the best log-rank statistic (search-adjusted p is shown when maxstat is installed), or your own value."),
                   choices = c("50 (top 50% high)" = "median",
                               "Auto (best log-rank)" = "auto",
                               "Custom" = "custom"),
                   selected = "median"),
      conditionalPanel(condition = "input.cutmode == 'custom'", ns = ns,
        radioButtons(ns("cutunit"), .lab("Unit", "Top percent: the highest x% of patients are 'high' (entering 25 gives the top quarter). Absolute threshold: an expression value on the log2 scale."),
                     choices = c("Top percent (%)" = "pct",
                                 "Absolute threshold" = "abs"),
                     selected = "pct"),
        numericInput(ns("cutval"), .lab("Value", "With 'Top percent' enter 1-99 (the highest x% becomes High); with 'Absolute threshold' enter the expression value."),
                     value = 50, min = 1, max = 99, step = 1)),
      shinyWidgets::actionBttn(ns("go"), "Analysis",
                               style = "gradient", color = "primary",
                               icon = icon("play"), block = TRUE),
      br(),
      p(class = "note", "Analysis downloads expression on demand, merges the
        survival table, draws the curve and lists the analysis data.")
    ),
    mainPanel(width = 9,
      verbatimTextOutput(ns("status")),
      ui_fig_block(ns, "plot", w = 600, h = 520),
      br(), hr(),
      h4("Analysis data (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("table")),
      br(),
      downloadButton(ns("csv"), "Download CSV")
    )
  )
}

server_mod_km <- function(id, rv, dataset_info) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    catalog <- .norm_catalog(dataset_info)
    ds_choices <- .dataset_choices(catalog)

    # ---------------- shared dataset / endpoint selection ----------------
    # initial values must be read with isolate(): reading a reactive value
    # outside a reactive consumer is an error in a real Shiny session
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
      if (!identical(input$ds, rv$sel$acc))
        updateSelectizeInput(session, "ds", selected = rv$sel$acc)
      eps <- .endpoint_choices(catalog, rv$sel$acc)
      sel <- if (!is.null(rv$sel$family) && rv$sel$family %in% eps)
        rv$sel$family else eps[1]
      updateSelectInput(session, "ep", choices = eps, selected = sel)
      if (is.na(rv$sel$ep %||% NA_character_))
        rv$sel$ep <- .resolve_family(catalog, rv$sel$acc, sel)
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

    # endpoint selector carries the FAMILY; the concrete token is resolved per cohort
    observeEvent(input$ep, {
      if (is.null(input$ep) || !nzchar(input$ep)) return()
      if (!identical(rv$sel$family, input$ep)) {
        rv$sel$family <- input$ep
        rv$sel$ep <- .resolve_family(catalog, rv$sel$acc, input$ep)
      }
    }, ignoreInit = TRUE)

    observeEvent(rv$sel$family, {
      if (!is.null(rv$sel$family) && !identical(input$ep, rv$sel$family))
        updateSelectInput(session, "ep", selected = rv$sel$family)
    }, ignoreNULL = TRUE)

    # ---------------- REF_ID choices ----------------
    # probe ids of the current gene, remembered for the hint below
    probes <- reactiveVal(character(0))
    observeEvent(list(rv$sel$acc, input$gene, input$kind), {
      if (is.null(input$kind) || input$kind != "gene") return()
      if (is.null(rv$sel$acc)) return()
      ch <- c("__all__" = "Auto: all probes -> max")
      tryCatch({
        pr <- .ref_probes(rv$sel$acc, input$gene)
        probes(pr)                                  # remember for the hint
        if (length(pr)) ch <- c(ch, stats::setNames(pr, pr))
        updateSelectizeInput(session, "ref", choices = ch, selected = "__all__",
                             server = FALSE)
      }, error = function(e) {
        updateSelectizeInput(session, "ref", choices = ch, selected = "__all__")
        showNotification(paste0("REF_ID lookup failed: ", conditionMessage(e)),
                         type = "warning")
      })
    })

    # what the probe box currently means (all probes -> max, one probe, or the
    # gene-level value for a TCGA cohort)
    output$ref_hint <- renderUI({
      if (!identical(input$kind %||% "gene", "gene")) return(NULL)
      txt <- .ref_hint(shiny::isolate(rv$sel$acc), input$gene %||% "GAPDH",
                       input$ref, length(probes()))
      if (is.null(txt)) return(NULL)
      div(class = "cpas-hint", txt)
    })

    # ---------------- analysis ----------------
    run_km <- eventReactive(input$go, {
      acc_run <- shiny::isolate(rv$sel$acc)
      tryCatch({
        if (is.null(acc_run)) stop("Please select a dataset first.")
        family <- rv$sel$family %||% "OS"
        ep <- rv$sel$ep %||% .resolve_family(catalog, rv$sel$acc, family)
        if (is.na(ep))
          stop("Cohort ", rv$sel$acc, " has no ", family, " endpoint.")
        kind <- input$kind %||% "gene"
        gene <- input$gene %||% "GAPDH"
        sig <- input$sig %||% "0.5*GAPDH + 0.5*TNS1"
        rule <- input$rule %||% "max"
        cutmode <- input$cutmode %||% "median"
        cutunit <- input$cutunit %||% "pct"
        cutval <- input$cutval %||% 50
        withProgress(message = "Fetching data and analysing...", value = 0.5, {
          df <- .build_marker_df(
            acc = acc_run, ep = ep, kind = kind,
            gene = gene, ref = if (kind == "gene") input$ref else NULL,
            formula_txt = if (kind == "signature") sig else NULL,
            rule = rule, allow_missing = isTRUE(input$allow_missing))
        })
        cut_desc <- switch(cutmode,
          median = "Top 50%", auto = "Auto",
          custom = if (cutunit == "pct") paste0("Top ", cutval, "%")
                   else paste0("abs ", cutval))
        cut_info <- NULL
        cut <- if (cutmode == "median") {
          stats::quantile(df$marker, 0.5, na.rm = TRUE)
        } else if (cutmode == "auto") {
          ac <- .cut_auto(df)
          cut_info <- ac
          if (is.finite(ac$cut)) ac$cut else stats::quantile(df$marker, 0.5, na.rm = TRUE)
        } else {
          if (cutunit == "pct") {
            if (cutval <= 0 || cutval >= 100)
              stop("Custom percentile must be between 1 and 99.")
            ## "Top x%" = the highest x% of patients are High, so the threshold is
            ## the (100 - x)th percentile: cutval = 25 -> quantile(0.75).
            stats::quantile(df$marker, 1 - cutval / 100, na.rm = TRUE)
          } else cutval
        }
        if (length(cut) == 0L || !is.finite(cut)) cut <- 0
        # analysis is called by FAMILY; the resolver picks the token present in df
        km <- plot_km(df, type = family, marker = "marker", cutpoint = cut,
                      pval = TRUE, risk.table = TRUE, xlab = "Time (years)")
        # the cut-point is shown in the figure itself, not only in the side panel
        leg <- switch(kind, gene = paste0("gene ", gene), signature = paste0("signature ", sig))
        km$plot <- km$plot + ggplot2::labs(
          subtitle = sprintf("%s | %s cut = %.4g (%s)", leg, toupper(cutmode), cut, cut_desc))
        list(df = df, km = km, cut = cut, cut_desc = cut_desc, ep = ep,
             family = family, cutmode = cutmode, cut_info = cut_info,
             acc = acc_run, kind = kind, gene = gene, sig = sig, rule = rule,
             ref = if (kind == "gene") input$ref else NULL,
             error = NULL, ts = Sys.time())
      }, error = function(e)
        list(df = NULL, km = NULL, cut = NULL, cut_desc = "", ep = NULL,
             family = NULL, cutmode = NULL, acc = acc_run,
             error = conditionMessage(e), ts = Sys.time()))
    })

    # NULL before the first click: outputs must not error out on a cold page
    current <- reactive({
      if (is.null(input$go)) return(NULL)
      run_km()
    })

    observeEvent(run_km(), {
      r <- run_km()
      if (!is.null(r$error))
        return(showNotification(paste0("Analysis failed: ", r$error), type = "error"))
      kind <- input$kind %||% "gene"
      rv$km_last <- r
    })

    output$status <- renderPrint({
      r <- current()
      if (is.null(r)) return(cat("Set the cohort, marker and cut point, then click 'Analysis'."))
      if (!is.null(r$error)) return(cat("Analysis failed: ", r$error))
      if (is.null(r$df)) return(cat("Set the marker and click Analysis."))
      n <- nrow(r$df); ev <- sum(r$df[[paste0(r$ep, "_status")]] == 1, na.rm = TRUE)
      cat(sprintf("Dataset %s | endpoint %s | n = %d, events = %d",
                  r$acc, .endpoint_status_label(r$family, r$ep), n, ev))
      grp <- ifelse(r$df$marker > r$cut, "High", "Low")
      cat(sprintf("\nGroups: High n = %d (%d events) | Low n = %d (%d events) | median follow-up %.4f y",
                  sum(grp == "High"), sum(r$df[[paste0(r$ep, "_status")]][grp == "High"] == 1, na.rm = TRUE),
                  sum(grp == "Low"), sum(r$df[[paste0(r$ep, "_status")]][grp == "Low"] == 1, na.rm = TRUE),
                  stats::median(r$df[[paste0(r$ep, "_time")]], na.rm = TRUE)))
      if (ev < 10L)
        cat(sprintf("\nNote: only %d event(s) in this cohort - the log-rank p-value and the HR are unstable.",
                    ev))
      if (identical(r$cutmode, "auto") && !is.null(r$cut_info)) {
        ci <- r$cut_info
        if (!is.finite(ci$p))
          cat("\nCut point: the best-log-rank search found no usable candidate cut point; the median was used.")
        else {
          cat(sprintf("\nCut point: best-log-rank split, searched over %d candidate cut points (naive log-rank p = %s)",
                      ci$n_cand, format.pval(ci$p, digits = 3)))
          if (!is.null(ci$adjusted_p) && is.finite(ci$adjusted_p))
            cat(sprintf(" | search-adjusted p = %s (maxstat %s)",
                        format.pval(ci$adjusted_p, digits = 3),
                        if (is.null(ci$adjusted_method) || is.na(ci$adjusted_method)) "" else ci$adjusted_method))
          else
            cat("\n  NOT adjusted for the cut-point search: treat this p-value as exploratory (install 'maxstat' for an adjusted p-value)")
        }
      }
      dr <- attr(r$df, "dropped_na")
      if (!is.null(dr) && sum(dr) > 0)
        cat(sprintf("\nRows excluded before analysis: %d without %s time/status, %d without marker",
                    dr[["surv"]], r$ep, dr[["marker"]]))
      cat(sprintf("\nMarker: %s | Cut: %s (threshold = %.4f)",
                  if (identical(r$kind, "gene"))
                    paste0("gene ", r$gene,
                           if (!is.null(r$ref) && nzchar(r$ref) && r$ref != "__all__")
                             paste0(" @probe ", r$ref, " [rule=", r$rule, "]")
                           else paste0(" (all probes, rule=", r$rule, ")"))
                  else paste0("signature ", r$sig, " [rule=", r$rule, "]"),
                  r$cut_desc, r$cut))
    })

    draw_km <- function() {
      r <- current()
      if (is.null(r)) { plot.new(); title("No analysis yet"); return(invisible(NULL)) }
      if (!is.null(r$error) || is.null(r$km)) { plot.new(); title("Analysis failed"); return(invisible(NULL)) }
      print(r$km)
    }
    server_fig_block(input, output, session, key = "plot", draw = draw_km,
                     ready = function() { r <- current(); !is.null(r) && is.null(r$error) && !is.null(r$km) },
                     w0 = 600, h0 = 520,
                     filename = function() paste0("KM_plot_",
                       if (is.null(current()) || is.null(current()$acc)) "cohort" else current()$acc))

    tbl_df <- reactive({
      r <- current()
      if (is.null(r) || !is.null(r$error) || is.null(r$df)) return(NULL)
      cols <- intersect(c("ID", paste0(r$ep, "_time"), paste0(r$ep, "_status"), "marker"),
                        colnames(r$df))
      r$df[, cols, drop = FALSE]
    })

    output$table <- DT::renderDT({
      req(tbl_df())
      .dtx4(tbl_df(), pageLength = 10)
    })

    output$csv <- downloadHandler(
      filename = function() {
        r <- current()
        paste0("KM_data_", if (is.null(r) || is.null(r$acc)) "cohort" else r$acc, ".csv")
      },
      content = function(file) {
        d <- tbl_df()
        if (is.null(d) || !nrow(d)) {
          showNotification("No completed analysis to download - click 'Analysis' first.",
                           type = "warning")
          return(NULL)
        }
        utils::write.csv(d, file, row.names = FALSE)
      })

  })
}

rounddf <- function(d) {
  for (cn in colnames(d)) if (is.numeric(d[[cn]])) d[[cn]] <- round(d[[cn]], 4)
  d
}
