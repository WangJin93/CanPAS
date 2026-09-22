# COX analysis module (univariable + multivariable; marker inherited from KM) --
ui_mod_cox <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(width = 3,
      h5("Cohort", class = "section-title"),
      selectInput(ns("ctype"), "Cancer type", choices = NULL),
      selectizeInput(ns("ds"), "Dataset", choices = NULL,
                     options = list(placeholder = "Select a cohort...",
                                    maxOptions = 30)),
      selectInput(ns("ep"), .lab("Endpoint", "Survival endpoint family (OS/DSS/DFS/PFS/MFS); this cohort contributes the token it has, e.g. RFS for DFS."), choices = NULL),
      hr(),
      h5("Marker", class = "section-title"),
      radioButtons(ns("kind"), .lab("Marker", "A single gene or a weighted signature such as 0.5*GAPDH + 0.5*TNS1."),
                   choices = c("Gene" = "gene", "Signature" = "signature"),
                   selected = "gene"),
      conditionalPanel(condition = "input.kind == 'gene'", ns = ns,
        textInput(ns("gene"), .lab("Gene", "HUGO gene symbol, e.g. GAPDH. The probes of the platform are resolved for it."), value = "GAPDH"),
        selectizeInput(ns("ref"), .lab("Probe", "Leave 'Auto: all probes' to use every probe of the gene (collapsed per sample by the maximum); pick one probe id to use that probe alone."),
                       choices = NULL,
                       options = list(placeholder = "Auto: all probes -> max", maxOptions = 500)),
        uiOutput(ns("ref_hint"))),
      conditionalPanel(condition = "input.kind == 'signature'", ns = ns,
        textInput(ns("sig"), .lab("Signature", "Weighted sum of genes, e.g. 0.5*GAPDH + 0.5*TNS1; each gene is collapsed over its probes first."), value = "0.5*GAPDH + 0.5*TNS1"),
        selectInput(ns("rule"), .lab("Probe rule", "How several probes of one gene are combined into one value per sample: max (default), mean, median or min."),
                    choices = c("max", "mean", "median", "min"), selected = "max"),
        checkboxInput(ns("allow_missing"), .lab("Allow missing genes", "On: genes this platform cannot measure are skipped and the signature is scored on the rest. Off: such a gene stops the analysis."),
                      value = FALSE)),
      hr(),
      h5("Clinical covariates (multi-select)", class = "section-title"),
      selectizeInput(ns("clin"), .lab("Covariates", "Clinical columns of this cohort. Numeric ones are fitted as continuous, the others as categorical (first level = reference); the marker is always included."),
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "e.g. age, sex, stage, grade",
                                    create = TRUE)),
      checkboxInput(ns("show_all"), .lab("Use all clinical covariates", "Fits every clinical column the cohort provides, instead of the ones selected above."),
                    value = FALSE),
      numericInput(ns("p_thr"), .lab("Multivariable p threshold", "Covariates whose univariable p is below this enter the multivariable model; the rest are shown in the univariable panel only."),
                   value = 0.05, min = 0, max = 1, step = 0.01),
      br(),
      shinyWidgets::actionBttn(ns("go"), "Run COX analysis",
                               style = "gradient", color = "success",
                               icon = icon("play"), block = TRUE)
    ),
    mainPanel(width = 9,
      verbatimTextOutput(ns("status")),
      tabsetPanel(
        id = ns("res_tabs"),
        tabPanel("Forest plots & tables",
          br(),
          h4("Univariable COX (all covariates; covariates kept in model order, levels grouped)", class = "section-title"),
          ui_fig_block(ns, "forest_uni", w = 600, h = 460),
          DT::DTOutput(ns("tbl_uni")),
          downloadButton(ns("csv_uni"), "Download CSV (univariable)"),
          br(), hr(),
          h4("Multivariable COX (significant covariates only; covariates kept in model order, levels grouped)", class = "section-title"),
          uiOutput(ns("multi_note")),
          ui_fig_block(ns, "forest_multi", w = 600, h = 460),
          DT::DTOutput(ns("tbl_multi")),
          downloadButton(ns("csv_multi"), "Download CSV (multivariable)"),
          br(), hr(),
          h4("Marker estimate: univariable vs multivariable", class = "section-title"),
          DT::DTOutput(ns("marker_cmp"))
        ),
        tabPanel("Print results",
          br(),
          h4("Univariable and multivariable Cox regression", class = "section-title"),
          uiOutput(ns("print_note")),
          uiOutput(ns("print_table")),
          br(),
          downloadButton(ns("dl_print_docx"), "Download table (Word)"),
          downloadButton(ns("dl_print_csv"), "Download table (CSV)")
        )
      )
    )
  )
}

server_mod_cox <- function(id, rv, dataset_info) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    catalog <- .norm_catalog(dataset_info)
    ds_choices <- .dataset_choices(catalog)

    # ---------------- shared dataset / endpoint selection ----------------
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

    # ---------------- marker inherited from KM page ----------------
    marker_info <- reactive({
      kind <- input$kind %||% "gene"
      list(kind = kind,
           gene = if (kind == "gene") (input$gene %||% "GAPDH") else NULL,
           ref  = if (kind == "gene") (input$ref %||% "__all__") else NULL,
           sig  = if (kind == "signature") (input$sig %||% "0.5*GAPDH + 0.5*TNS1") else NULL,
           rule = input$rule %||% "max",
           allow_missing = isTRUE(input$allow_missing))
    })
    # probe list for the selected gene (loaded on demand)
    # probe ids of the current gene, remembered for the hint below
    probes <- reactiveVal(character(0))
    observeEvent(list(input$ds, input$gene, input$kind), {
      if (!identical(input$kind %||% "gene", "gene")) return()
      acc <- shiny::isolate(rv$sel$acc); g <- input$gene %||% "GAPDH"
      if (is.null(acc) || !nzchar(g)) return()
      pr <- tryCatch(.ref_probes(acc, g), error = function(e) character(0))
      probes(pr)
      updateSelectizeInput(session, "ref",
                           choices = c("Auto: all probes -> max" = "__all__", pr),
                           selected = "__all__")
    }, ignoreNULL = TRUE)

    # what the probe box currently means for the marker
    output$ref_hint <- renderUI({
      if (!identical(input$kind %||% "gene", "gene")) return(NULL)
      txt <- .ref_hint(shiny::isolate(rv$sel$acc), input$gene %||% "GAPDH",
                       input$ref, length(probes()))
      if (is.null(txt)) return(NULL)
      div(class = "cpas-hint", txt)
    })

    # ---------------- clinical covariates ----------------
    observeEvent(rv$sel$acc, {
      acc <- rv$sel$acc
      ch <- character(0)
      if (!is.null(acc)) {
        cl <- tryCatch(.loc_clin(acc), error = function(e) NULL)
        if (!is.null(cl)) ch <- sort(setdiff(colnames(cl), "ID"))
      }
      updateSelectizeInput(session, "clin", choices = ch, selected = character(0))
    }, ignoreNULL = TRUE)

    run_cox <- eventReactive(input$go, {
      acc_run <- shiny::isolate(rv$sel$acc)
      # multivariable models whose repair failed are reported separately, so a
      # single bad covariate set cannot discard the univariable panel
      multi_error <- NULL
      tryCatch({
        if (is.null(acc_run)) stop("Please select a dataset first.")
        family <- shiny::isolate(rv$sel$family) %||% "OS"
        ep <- shiny::isolate(rv$sel$ep) %||% .resolve_family(catalog, acc_run, family)
        if (is.na(ep)) stop("Cohort ", acc_run, " has no ", family, " endpoint.")
        m <- marker_info()
        withProgress(message = "Fetching data and fitting models...", value = 0.5, {
          df <- .build_marker_df(
            acc = acc_run, ep = ep, kind = m$kind,
            gene = m$gene, ref = if (m$kind == "gene") m$ref else NULL,
            formula_txt = if (m$kind == "signature") m$sig else NULL,
            rule = m$rule %||% "max", allow_missing = m$allow_missing)
        })
        clin_chosen <- if (isTRUE(input$show_all)) {
          cl <- .loc_clin(acc_run)
          if (is.null(cl)) character(0) else setdiff(colnames(cl), "ID")
        } else input$clin %||% character(0)
        # a typed covariate that is not a column must not vanish silently
        clin_missing <- setdiff(clin_chosen, colnames(df))
        clin_chosen <- intersect(clin_chosen, colnames(df))
        if (length(clin_missing))
          showNotification(paste0("Covariate(s) not available in this cohort and skipped: ",
                                  paste(clin_missing, collapse = ", ")), type = "warning",
                           duration = 8)
        conts <- "marker"; cates <- character(0)
        for (v in clin_chosen) {
          if (is.numeric(df[[v]])) conts <- c(conts, v) else cates <- c(cates, v)
        }
        uni <- COX_analysis(df, type = family, cont_Variates = conts,
                            cate_Variates = cates, method = "uni")
        # the multivariable model keeps the covariates that reached the
        # univariable p threshold (COX_analysis(method = "uni") rows per variable)
        thr <- suppressWarnings(as.numeric(input$p_thr %||% 0.05))
        if (!is.finite(thr)) thr <- 0.05
        ut <- uni$results_table
        keep_vars <- unique(ut$Var1[!is.na(ut$Pvalue) & ut$Pvalue < thr])
        if (!length(keep_vars)) {
          multi <- NULL
          sel_note <- sprintf("No covariate reached p < %.4f; the multivariable model was not fitted.", thr)
        } else if (length(keep_vars) < 2L) {
          multi <- NULL
          sel_note <- sprintf("Only one covariate (%s) reached p < %.4f; a multivariable model needs at least two.",
                              keep_vars, thr)
        } else {
          # the multivariable step is isolated: a model that cannot be repaired
          # must not hide the univariable results and figures
          multi <- tryCatch(COX_analysis(df, type = family,
                                cont_Variates = intersect(keep_vars, conts),
                                cate_Variates = intersect(keep_vars, cates), method = "multi"),
                            error = function(e) structure(list(message = conditionMessage(e)),
                                                          class = "cpas_multi_error"))
          multi_error <- if (inherits(multi, "cpas_multi_error")) multi$message else NULL
          if (inherits(multi, "cpas_multi_error")) multi <- NULL
          sel_note <- sprintf("Covariates with univariable p < %.4f: %s", thr,
                              paste(keep_vars, collapse = ", "))
          if (!is.null(multi_error))
            sel_note <- paste0(sel_note,
              "\nThe multivariable model could not be estimated; the univariable results above are unaffected.")
        }
        # the print-ready summary: the same univariable and multivariable models
        # combined by COX_screen_adjust() (verified to give numerically identical
        # rows to the two panels above) and formatted as a three-line table; its
        # own failure must not take the page down
        screen <- tryCatch(suppressMessages(suppressWarnings(
                    COX_screen_adjust(df, type = family,
                                      cont_Variates = conts, cate_Variates = cates,
                                      p.threshold = thr))),
                  error = function(e) structure(list(message = conditionMessage(e)),
                                                class = "cpas_screen_error"))
        screen_error <- if (inherits(screen, "cpas_screen_error")) screen$message else NULL
        if (inherits(screen, "cpas_screen_error")) screen <- NULL
        list(df = df, uni = uni, multi = multi, multi_error = multi_error,
             screen = screen, screen_error = screen_error, ep = ep, family = family,
             multi_drops = tryCatch(multi$metadata$dropped_covariates, error = function(e) NULL),
             multi_reduced = isTRUE(tryCatch(multi$metadata$reduced, error = function(e) FALSE)),
             multi_note_txt = tryCatch(as.character(multi$metadata$reduced_note), error = function(e) character(0)),
             clinical = clin_chosen, missing_clinical = clin_missing,
             selected = keep_vars, selection_note = sel_note, p_thr = thr,
             acc = acc_run, kind = m$kind, gene = m$gene, sig = m$sig,
             error = NULL, ts = Sys.time())
      }, error = function(e)
        list(df = NULL, uni = NULL, multi = NULL, ep = NULL, family = NULL,
             clinical = character(0), acc = acc_run,
             error = conditionMessage(e), ts = Sys.time()))
    })

    current <- reactive({
      if (is.null(input$go)) return(NULL)
      run_cox()
    })

    # When the multivariable model needed covariates to be dropped, merged or
    # excluded so that it stayed estimable, say so where the model is shown.
    output$multi_note <- renderUI({
      r <- current()
      if (is.null(r) || is.null(r$multi_error) && is.null(tryCatch(r$multi, error = function(e) NULL)))
        return(NULL)
      dc <- tryCatch(r$multi_drops, error = function(e) NULL)
      blocks <- list()
      if (!is.null(r$multi_error)) {
        blocks <- c(blocks, list(div(class = "alert alert-danger", style = "padding:8px; margin-top:6px;",
          tags$b("The multivariable model could not be estimated. "),
          "The univariable results above are unaffected. Reasons reported by the ",
          "estimator and by the automatic screening:", tags$br(),
          tags$code(r$multi_error), tags$br(),
          "Try fewer covariates, or use the univariable panel / the 'COX by genes' page.")))
      }
      if (isTRUE(r$multi_reduced)) {
        blocks <- c(blocks, list(div(class = "alert alert-danger", style = "padding:8px; margin-top:6px;",
          tags$b("The multivariable model was reduced and is NOT adjusted for confounding. "),
          paste(r$multi_note_txt, collapse = " "))))
      }
      if (!is.null(dc) && nrow(dc)) {
        blocks <- c(blocks, list(div(class = "alert alert-warning", style = "padding:8px; margin-top:6px;",
          tags$b("Multivariable model adjusted automatically: "),
          "the covariates below could not be estimated as supplied and were ",
          "handled as stated, so that the remaining model could be fitted. ",
          "The univariable panel above still shows every covariate.",
          tags$ul(lapply(seq_len(nrow(dc)), function(i)
            tags$li(tags$b(as.character(dc$variable[i])), " — ",
                    as.character(dc$reason[i])))))))
      }
      if (!length(blocks)) return(NULL)
      tagList(blocks)
    })

    ph_text <- function(r) {
      ph <- tryCatch(r$multi$metadata$ph_test, error = function(e) NULL)
      if (is.null(ph)) return("")
      p <- ph$table["GLOBAL", "p"]
      sprintf("\nProportional-hazards check (cox.zph, GLOBAL): p = %s%s",
              format.pval(p, digits = 3),
              if (is.finite(p) && p < 0.05) " - the PH assumption is questionable; interpret HRs with care" else "")
    }

    output$status <- renderPrint({
      r <- current()
      if (is.null(r)) return(cat("Choose the marker (KM page), covariates and click 'Run COX analysis'."))
      if (!is.null(r$error)) return(cat("Analysis failed: ", r$error))
      if (is.null(r$df)) return(cat("Click 'Run COX analysis' to start."))
      ev <- sum(r$df[[paste0(r$ep, "_status")]] == 1, na.rm = TRUE)
      cat(sprintf("Dataset %s | endpoint %s | n = %d, events = %d | covariates: %s",
                  r$acc, .endpoint_status_label(r$family, r$ep), nrow(r$df), ev,
                  ifelse(length(r$clinical) == 0, "(none, marker only)",
                         paste(r$clinical, collapse = ", "))))
      cat(sprintf("\nMarker: %s", if (identical(r$kind, "gene"))
        paste0("gene ", r$gene) else paste0("signature ", r$sig)))
      cat("\n", r$selection_note %||% "", "\n", sep = "")
      mi <- tryCatch(r$multi$metadata, error = function(e) NULL)
      if (!is.null(r$multi_error)) {
        cat("\nMultivariable model: NOT ESTIMATED. Reasons: ",
            r$multi_error, "\n", sep = "")
        cat("The univariable results above are unaffected.\n")
      }
      if (!is.null(mi)) {
        cat(sprintf("Multivariable model: complete cases = %s, events in model = %s (%s events per coefficient)",
                    ifelse(is.na(mi$complete_cases), "NA", mi$complete_cases),
                    ifelse(is.na(mi$events_in_model), "NA", mi$events_in_model),
                    if (!is.null(mi$complete_cases) && !is.na(mi$complete_cases) &&
                        length(mi$cate_vars) + length(mi$cont_vars) > 0)
                      sprintf("%.4f", mi$events_in_model /
                                max(1, length(mi$cont_vars) + length(mi$cate_vars))) else "NA"))
        if (length(mi$skipped_variables))
          cat("\nVariable(s) skipped: ", paste(mi$skipped_variables, collapse = "; "), sep = "")
        dc <- mi$dropped_covariates
        if (!is.null(dc) && nrow(dc)) {
          cat("\nCovariates handled automatically to keep the model estimable:\n")
          for (i in seq_len(nrow(dc)))
            cat(sprintf("  - %s [%s]: %s\n", dc$variable[i], dc$detail[i], dc$reason[i]))
        }
        if (isTRUE(mi$reduced))
          cat("\nWARNING: the multivariable model was reduced to ",
              length(mi$final_covariates), " covariate(s) (",
              paste(mi$final_covariates, collapse = ", "),
              ") and is NOT adjusted for confounding.\n", sep = "")
        cat(ph_text(r))
      }
      dr <- attr(r$df, "dropped_na")
      if (!is.null(dr) && sum(dr) > 0)
        cat(sprintf("\nRows excluded before analysis: %d without %s time/status, %d without marker",
                    dr[["surv"]], r$ep, dr[["marker"]]))
    })

    multi_ok <- function() { r <- current(); !is.null(r) && !is.null(r$multi) }
    draw_uni <- function() {
      r <- current()
      if (is.null(r) || !is.null(r$error) || is.null(r$uni)) {
        plot.new(); title(if (is.null(r)) "No analysis yet" else
          if (!is.null(r$error)) "Analysis failed" else "No univariable result"); return(invisible(NULL))
      }
      # model order, NOT sorted by HR: the levels of a clinical covariate must
      # stay together under their own name (see ?forest_plot, group_levels)
      print(forest_plot(COX_out = r$uni$results_table, digits = 4,
                        HR_order = "none"))
    }
    draw_multi <- function() {
      r <- current()
      if (multi_ok()) {
        p <- forest_plot(COX_out = r$multi$results_table, digits = 4,
                         HR_order = "none")
        # a model that the repair had to reduce must be labelled on the figure
        # itself, so a reduced fit cannot be mistaken for an adjusted model
        if (isTRUE(r$multi_reduced)) {
          fc <- tryCatch(as.character(r$multi$metadata$final_covariates), error = function(e) character(0))
          p <- p + ggplot2::ggtitle(sprintf("REDUCED model: %d covariate(s) (%s) - not adjusted for confounding",
                                            length(fc), paste(fc, collapse = ", "))) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 11, colour = "#B00020",
                                                             face = "bold", hjust = 0.5))
        }
        print(p); return(invisible(NULL))
      }
      # the model could not be fitted: draw the reason instead of an empty frame
      plot.new()
      msg <- if (is.null(r)) "No analysis yet" else
        if (!is.null(r$multi_error)) strwrap(paste0("Multivariable model not estimable. ",
             "Reasons: ", r$multi_error), width = 90) else "No multivariable result"
      title(main = "Multivariable COX model could not be estimated",
            cex.main = 1, col.main = "#B00020")
      text(0.5, 0.55, paste(msg, collapse = "\n"), cex = 0.8, col = "grey20")
    }
    server_fig_block(input, output, session, key = "forest_uni", draw = draw_uni,
                     ready = function() { r <- current(); !is.null(r) && is.null(r$error) && !is.null(r$uni) },
                     w0 = 600, h0 = 460,
                     filename = function() paste0("COX_univariable_",
                       if (is.null(current()) || is.null(current()$acc)) "dataset" else current()$acc))
    server_fig_block(input, output, session, key = "forest_multi", draw = draw_multi,
                     ready = function() { r <- current(); !is.null(r) && (multi_ok() || !is.null(r$multi_error)) },
                     w0 = 600, h0 = 460,
                     filename = function() paste0("COX_multivariable_",
                       if (is.null(current()) || is.null(current()$acc)) "dataset" else current()$acc))
    output$tbl_uni <- DT::renderDT({
      r <- current(); if (is.null(r) || !is.null(r$error)) return(NULL)
      if (!is.null(r$uni)) .dtx4(r$uni$results_table)
    })
    output$tbl_multi <- DT::renderDT({
      r <- current(); if (is.null(r) || !is.null(r$error)) return(NULL)
      if (!is.null(r$multi)) .dtx4(r$multi$results_table)
    })
    output$marker_cmp <- DT::renderDT({
      r <- current()
      if (is.null(r) || !is.null(r$error)) return(NULL)
      if (is.null(r$uni) || is.null(r$multi)) return(NULL)
      take <- function(tb) tb[tb$Var1 == "marker" & !is.na(tb$HR),
                              c("Level", "HR", "HR95L", "HR95H", "Pvalue"), drop = FALSE]
      a <- take(r$uni$results_table); b <- take(r$multi$results_table)
      if (!nrow(a) || !nrow(b)) return(NULL)
      cm <- data.frame(model = c("univariable", "multivariable"), rbind(a[1, ], b[1, ]))
      .dtx4(cm)
    })

    # ---- Print results tab -------------------------------------------------
    # COX_screen_adjust() returns the same two models as one summary table,
    # formatted as a three-line table (top rule, rule under the header, bottom
    # rule) so it can go straight into a manuscript.
    print_ft <- function() {
      r <- current()
      if (is.null(r) || is.null(r$screen)) return(NULL)
      ft <- r$screen$print_result
      if (is.null(ft)) return(NULL)
      cap <- sprintf("%s (%s, %s), n = %d, events = %d. Univariable models for all covariates; multivariable model for those with p < %.4f (%s).",
                     r$acc, .endpoint_status_label(r$family, r$ep),
                     if (identical(r$kind, "gene")) paste0("marker: ", r$gene) else paste0("signature: ", r$sig),
                     r$screen$result$N[1] %||% nrow(r$df),
                     sum(r$df[[paste0(r$ep, "_status")]] == 1, na.rm = TRUE), r$p_thr,
                     if (is.null(r$multi)) "none reached the threshold" else
                       paste(intersect(r$selected, c(r$screen$cont_Variates, r$screen$cate_Variates)), collapse = ", "))
      flextable::set_caption(ft, caption = cap)
    }

    output$print_note <- renderUI({
      r <- current()
      if (is.null(r)) return(div(class = "note", "Run the analysis to build the summary table."))
      if (!is.null(r$screen_error))
        return(div(class = "alert alert-danger", style = "padding:8px;",
          tags$b("The summary table could not be built. "), tags$code(r$screen_error)))
      if (is.null(r$screen)) return(NULL)
      msgs <- list()
      if (isTRUE(r$multi_reduced))
        msgs <- c(msgs, list(div(class = "alert alert-danger", style = "padding:8px;",
          tags$b("The multivariable model was reduced and is NOT adjusted for confounding. "),
          paste(r$multi_note_txt, collapse = " "))))
      dc <- tryCatch(r$multi_drops, error = function(e) NULL)
      if (!is.null(dc) && nrow(dc))
        msgs <- c(msgs, list(div(class = "alert alert-warning", style = "padding:8px;",
          tags$b("Covariates handled automatically (see the Forest tab for details): "),
          paste(sprintf("%s - %s", dc$variable, dc$reason), collapse = "; "))))
      if (!length(msgs)) return(NULL)
      tagList(msgs)
    })

    output$print_table <- renderUI({
      r <- current()
      if (is.null(r)) return(NULL)
      ft <- print_ft()
      if (is.null(ft)) return(div(class = "note", "No table yet."))
      div(style = "overflow-x:auto;", flextable::htmltools_value(ft))
    })

    output$dl_print_docx <- downloadHandler(
      filename = function() {
        r <- current(); paste0("COX_table_", if (is.null(r) || is.null(r$acc)) "cohort" else r$acc, ".docx")
      },
      content = function(file) {
        ft <- print_ft()
        if (is.null(ft)) {
          showNotification("No table to download - run the analysis first.", type = "warning"); return(NULL)
        }
        flextable::save_as_docx(ft, path = file)
      })
    output$dl_print_csv <- downloadHandler(
      filename = function() {
        r <- current(); paste0("COX_table_", if (is.null(r) || is.null(r$acc)) "cohort" else r$acc, ".csv")
      },
      content = function(file) {
        r <- current()
        if (is.null(r) || is.null(r$screen)) {
          showNotification("No table to download - run the analysis first.", type = "warning"); return(NULL)
        }
        utils::write.csv(r$screen$result, file, row.names = FALSE)
      })

    output$csv_uni <- downloadHandler(
      filename = function() {
        r <- current(); paste0("COX_uni_", if (is.null(r) || is.null(r$acc)) "cohort" else r$acc, ".csv")
      },
      content = function(file) {
        r <- current()
        if (is.null(r) || !is.null(r$error) || is.null(r$uni)) {
          showNotification("No completed analysis to download - run the analysis first.",
                           type = "warning"); return(NULL)
        }
        utils::write.csv(r$uni$results_table, file, row.names = FALSE)
      })
    output$csv_multi <- downloadHandler(
      filename = function() {
        r <- current(); paste0("COX_multi_", if (is.null(r) || is.null(r$acc)) "cohort" else r$acc, ".csv")
      },
      content = function(file) {
        r <- current()
        if (is.null(r) || !is.null(r$error) || is.null(r$multi)) {
          showNotification("No completed analysis to download - run the analysis first.",
                           type = "warning"); return(NULL)
        }
        utils::write.csv(r$multi$results_table, file, row.names = FALSE)
      })
  })
}

dtx <- function(d) {
  if ("Pvalue" %in% colnames(d))
    d$Pvalue <- ifelse(is.na(d$Pvalue), "",
                       ifelse(d$Pvalue < 0.001, "<0.001", sprintf("%.4f", d$Pvalue)))
  DT::datatable(d, rownames = FALSE, filter = "top",
                options = list(pageLength = 12, scrollX = TRUE, dom = "ftp"))
}
