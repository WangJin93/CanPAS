# Multi-datasets analysis · Pooled KM -------------------------------------
# Two ways of integrating several datasets into one KM analysis:
#   "ipd"  = pool the patients (median split inside each dataset) and fit one
#            survival curve, reported with the unstratified and the
#            dataset-stratified log-rank p;
#   "meta" = pool the survival probabilities S(t) at landmark times through
#            log(-log S) inverse variance;
#   "both" = both routes (default).

ui_mod_pooled_km <- function(id) {
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
      h5("Marker", class = "section-title"),
      radioButtons(ns("mkind"), .lab("Marker", "A single gene, or a weighted signature such as 0.5*GAPDH + 0.5*TNS1."),
                   choices = c("Gene" = "gene", "Signature" = "signature"),
                   selected = "gene"),
      conditionalPanel(condition = "input.mkind == 'gene'", ns = ns,
        textInput(ns("mgene"), .lab("Gene", "HUGO gene symbol; its probes are collapsed by the rule below."), value = "GAPDH"),
        selectInput(ns("mrule_g"), .lab("Probe rule", "How several probes of one gene are combined into one value per sample: max (default), mean, median or min."),
                    choices = c("max", "mean", "median", "min"), selected = "max")),
      conditionalPanel(condition = "input.mkind == 'signature'", ns = ns,
        textInput(ns("msig"), .lab("Signature", "Weighted signature formula, e.g. 0.5*GAPDH + 0.5*TNS1."), value = "0.5*GAPDH + 0.5*TNS1"),
        selectInput(ns("mrule"), .lab("Probe rule", "Each gene of the signature is first collapsed over its probes by this rule, then weighted and summed."),
                    choices = c("max", "mean", "median", "min"), selected = "max"),
        checkboxInput(ns("mallow"), .lab("Allow missing genes", "On: genes this platform cannot measure are skipped and the score is computed on the remaining ones."), value = FALSE)),
      hr(),
      h5("Integration", class = "section-title"),
      selectInput(ns("km_method"), .lab("Pooling", "ipd = pool the patients and estimate one curve; meta = estimate KM inside each cohort and pool S(t) at each time point; both = output the two."),
                  choices = c("Both (IPD + time-point meta)" = "both",
                              "IPD (pool patients)" = "ipd",
                              "Time-point meta (pool S(t))" = "meta"),
                  selected = "both"),
      ## the meta-pooling estimator and the landmark years only matter when the
      ## time-point route is part of the run, so they appear with it
      conditionalPanel(condition = "input.km_method != 'ipd'", ns = ns,
        selectInput(ns("mmet"), .lab("Meta pooling", "Estimator for the time-point pooling of S(t): RE = random effects DerSimonian-Laird (default), FE = fixed effect. Only used when the time-point route runs."), choices = c("RE", "FE"), selected = "RE"),
        fluidRow(
          column(4, numericInput(ns("km_l1"), .lab("1y", "Landmark time(s) in years at which S(t) is pooled across cohorts; edit any of the three boxes. Times beyond a cohort's follow-up are carried forward and reported as such."), 1, min = 0)),
          column(4, numericInput(ns("km_l3"), "3y", 3, min = 0)),
          column(4, numericInput(ns("km_l5"), "5y", 5, min = 0))
        )
      ),
      hr(),
      h5("Cut point", class = "section-title"),
      radioButtons(ns("km_cut"), .lab("Cut point",
        paste("How the high/low groups are formed, always inside each cohort separately so the",
              "groups stay comparable across platforms. '50 (top 50% high)' splits at that",
              "cohort's own median. 'Top percent (%)' takes the highest x% of that cohort by",
              "expression as High and the rest as Low: entering 25 means the top quarter of",
              "patients are High. A per-cohort search for the best cut point is deliberately",
              "not offered: that search inflates the p-value (see the cut-point simulation in",
              "the paper) and pooling it over cohorts compounds the inflation.")),
                   choices = c("50 (top 50% high)" = "median",
                               "Top percent (%)" = "top_pct"),
                   selected = "median"),
      conditionalPanel(condition = "input.km_cut == 'top_pct'", ns = ns,
        numericInput(ns("km_top_pct"), .lab("Top percent (%)",
          "The highest x% of patients in each cohort by expression are High, the rest Low. Enter 25 for the top quarter. The threshold is that cohort's (100 - x)th percentile; where several patients tie at the threshold slightly more than x% can fall into High, and the actual count is reported with the result."), 25, min = 1, max = 99, step = 1)
      ),
      shinyWidgets::actionBttn(ns("btn_km"), "Run pooled KM",
                               style = "gradient", color = "primary",
                               icon = icon("chart-line"), block = TRUE),
      br(),
      p(class = "note", "Each dataset is grouped inside itself — at its own median
        marker level, or at the top-percent threshold you set — so the comparison
        is 'high vs low within each dataset'. The dataset-stratified log-rank p is
        the one to report.")
    ),
    mainPanel(width = 9,
      verbatimTextOutput(ns("km_status")),
      uiOutput(ns("overlap")),
      ui_fig_block(ns, "km_plot", w = 680, h = 780),
      br(), hr(),
      h5("Landmark survival (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("km_land")),
      br(),
      h5("Per-dataset KM curves", class = "section-title"),
      ui_fig_block(ns, "km_grid", w = 680, h = 620),
      br(), hr(),
      h5("Per-dataset estimates (4 decimals)", class = "section-title"),
      DT::DTOutput(ns("km_per")),
      br(),
      downloadButton(ns("dl_km_csv"), "Download landmark table (CSV)", class = "btn-success")
    )
  )
}

server_mod_pooled_km <- function(id, rv, dataset_info) {
  moduleServer(id, function(input, output, session) {
    output$overlap <- renderUI({
      .overlap_warning(dataset_info, input$mcohorts %||% character(0))
    })
    catalog <- .norm_catalog(dataset_info)
    loc <- reactiveValues(km = NULL)

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

    observeEvent(input$btn_km, {
      spec <- .marker_spec(input)
      accs <- input$mcohorts %||% character(0)
      loc$km <- tryCatch({
        b <- .multi_bundle(session, catalog, accs, input$mep %||% "OS", spec, min_n = 2)
        lm <- unique(c(input$km_l1 %||% 1, input$km_l3 %||% 3, input$km_l5 %||% 5))
        lm <- lm[is.finite(lm) & lm > 0]
        if (!length(lm)) lm <- c(1, 3, 5)
        res <- cpas_km_pooled(b$merged, marker = "marker", type = b$family,
                              method = input$km_method %||% "both",
                              landmarks = lm, meta_method = input$mmet %||% "RE",
                              cut = input$km_cut %||% "median",
                              top_pct = if (identical(input$km_cut, "top_pct"))
                                input$km_top_pct else NULL,
                              cut_value = if (identical(input$km_cut, "custom"))
                                input$km_cut_value else NULL)
        list(res = res, failed = b$failed, family = b$family, spec = spec)
      }, error = function(e) list(error = conditionMessage(e)))
      if (!is.null(loc$km$error))
        showNotification(paste0("Pooled KM failed: ", loc$km$error),
                         type = "error", duration = 10)
    })

    output$km_status <- renderPrint({
      r <- loc$km
      if (is.null(r)) return(cat("Select >= 2 datasets, the marker and the method, then click 'Run pooled KM'."))
      if (!is.null(r$error)) return(cat("Failed: ", r$error))
      k <- r$res
      cat(sprintf("Datasets pooled: %d | samples %d (High %d / Low %d) | method %s | family %s\n",
                  length(k$datasets), nrow(k$df), k$n_high, k$n_low, k$method, r$family))
      cat("Marker: ", if (identical(r$spec$kind, "gene"))
        paste0("gene ", r$spec$gene) else paste0("signature ", r$spec$sig), "\n", sep = "")
      cat("Cut: ", k$cutpoint, "\n", sep = "")
      if (length(k$empty_cohorts))
        cat("Left out by the threshold: ", paste(k$empty_cohorts, collapse = ", "),
            " (", paste(unlist(k$empty_reasons), collapse = "; "), ")\n", sep = "")
      if (length(k$skipped_cohorts))
        cat("Skipped (insufficient or constant data): ",
            paste(k$skipped_cohorts, collapse = ", "),
            " (", paste(unlist(k$skipped_reasons), collapse = "; "), ")\n", sep = "")
      if (!is.null(k$logrank_p))
        cat(sprintf("Pooled log-rank p = %s | dataset-stratified p = %s\n",
                    .fmtp4(k$logrank_p), .fmtp4(k$logrank_p_stratified)))
      if (!is.null(k$dataset_endpoints))
        cat("Endpoints used: ",
            paste(sprintf("%s=%s", names(k$dataset_endpoints), k$dataset_endpoints),
                  collapse = ", "), "\n", sep = "")
      if (length(r$failed))
        cat("Excluded: ", paste(names(r$failed), r$failed, sep = ": ",
                                collapse = "; "), "\n", sep = "")
    })
    draw_km <- function() {
      r <- loc$km
      if (is.null(r) || !is.null(r$error)) {
        plot.new(); title(if (is.null(r)) "No analysis yet" else "Analysis failed")
        return(invisible(NULL))
      }
      k <- r$res
      has_ipd <- !is.null(k$fit); has_meta <- !is.null(k$meta_curve)
      if (has_ipd && has_meta) {
        ## method = "both" holds BOTH routes. A single plot_cpas_km() call returns
        ## only the IPD curve when both are present, so the two are composed here.
        ## The composition uses patchwork::wrap_plots() rather than the / and +
        ## patchwork operators: the app only loads patchwork's namespace (it is a
        ## Suggests dependency), so operator dispatch is not guaranteed.
        if (!requireNamespace("patchwork", quietly = TRUE)) {
          print(plot_cpas_km(k, which = "ipd"))
          showNotification("Install the 'patchwork' package to see both pooling routes in one figure.",
                           type = "warning", duration = 8)
          return(invisible(NULL))
        }
        pi <- plot_cpas_km(k, which = "ipd")
        pm <- plot_cpas_km(k, which = "meta")
        curve_a <- pi$plot +
          ggplot2::ggtitle("(A) IPD: patients of all cohorts pooled into one Kaplan-Meier") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
        panels <- list(curve_a)
        hs <- 3
        if (!is.null(pi$table)) {
          panels <- c(panels, list(pi$table +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 9))))
          hs <- c(hs, 1)
        }
        panels <- c(panels, list(pm +
          ggplot2::ggtitle("(B) Time-point meta: S(t) pooled across cohorts (points = the landmark years)") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
                         legend.position = "bottom")))
        hs <- c(hs, 2)
        print(patchwork::wrap_plots(panels, ncol = 1, heights = hs))
      } else {
        print(plot_cpas_km(k))
      }
    }
    server_fig_block(input, output, session, key = "km_plot", draw = draw_km,
                     ready = function() { r <- loc$km; !is.null(r) && is.null(r$error) && !is.null(r$res) },
                     ## taller than a single curve: with method = "both" this
                     ## canvas holds the IPD curve, its risk table and the
                     ## time-point meta curve
                     w0 = 680, h0 = 780, filename = function() "CPAS_pooled_KM")
    output$km_land <- DT::renderDT({
      r <- loc$km
      if (is.null(r) || !is.null(r$error) || is.null(r$res$meta_landmarks)) return(NULL)
      d <- r$res$meta_landmarks
      if (is.null(d) || !nrow(d)) return(NULL)
      .dtx4(data.frame(group = d$group, time = d$time, k = d$k,
                       S = .fmt4(d$S),
                       `95% CI` = paste0(.fmt4(d$lower), " - ", .fmt4(d$upper)),
                       I2 = .fmt4(d$I2), p_het = .fmtp4(d$p_het),
                       check.names = FALSE, stringsAsFactors = FALSE))
    })
    draw_grid <- function() {
      r <- loc$km
      if (is.null(r) || !is.null(r$error)) { plot.new(); title(""); return(invisible(NULL)) }
      p <- tryCatch(plot_cpas_km_perdataset(r$res, ncol = 3), error = function(e) NULL)
      if (!is.null(p)) grid::grid.draw(p) else { plot.new(); title("Per-dataset curves unavailable") }
    }
    server_fig_block(input, output, session, key = "km_grid", draw = draw_grid,
                     ready = function() { r <- loc$km; !is.null(r) && is.null(r$error) && !is.null(r$res) },
                     w0 = 600, h0 = 620, filename = function() "CPAS_pooled_KM_per_dataset")
    output$km_per <- DT::renderDT({
      r <- loc$km
      if (is.null(r) || !is.null(r$error)) return(NULL)
      d <- r$res$df
      if (is.null(d) || !nrow(d)) return(NULL)
      agg <- do.call(rbind, lapply(sort(unique(d$dataset)), function(cf) {
        x <- d[d$dataset == cf, ]
        fit <- tryCatch(survival::survfit(survival::Surv(time, status) ~ group, data = x),
                        error = function(e) NULL)
        hr <- tryCatch({
          m <- survival::coxph(survival::Surv(time, status) ~ (group == "High"), data = x)
          s <- summary(m)$coefficients
          c(HR = exp(s[1, 1]), lower = exp(s[1, 1] - 1.96 * s[1, 3]),
            upper = exp(s[1, 1] + 1.96 * s[1, 3]), p = s[1, 5])
        }, error = function(e) c(HR = NA, lower = NA, upper = NA, p = NA))
        data.frame(dataset = cf, n = nrow(x), events = sum(x$status == 1),
                   HR = hr[["HR"]], lower = hr[["lower"]], upper = hr[["upper"]],
                   p = hr[["p"]], stringsAsFactors = FALSE)
      }))
      .dtx4(data.frame(dataset = agg$dataset, N = agg$n, events = agg$events,
                       HR = .fmt4(agg$HR),
                       `95% CI` = paste0(.fmt4(agg$lower), " - ", .fmt4(agg$upper)),
                       p = .fmtp4(agg$p), check.names = FALSE, stringsAsFactors = FALSE))
    })
    output$dl_km_csv <- downloadHandler(
      filename = function() "CPAS_pooled_KM_landmarks.csv",
      content = function(file) {
        r <- loc$km
        if (is.null(r) || !is.null(r$error) || is.null(r$res$meta_landmarks)) {
          showNotification("No result to download.", type = "warning"); return(NULL)
        }
        utils::write.csv(.round4(r$res$meta_landmarks), file, row.names = FALSE)
      })
  })
}
