#' @title Per-gene meta-analysis panel across datasets
#' @description
#' Runs the two-stage meta-analysis of \code{\link{cpas_meta}} for every gene in
#' \code{genes} over the same set of \code{datasets} and collects the pooled
#' estimates, heterogeneity statistics, prediction intervals and coverage into
#' one table, with Benjamini-Hochberg FDR across the genes. Expression and
#' survival tables are fetched **once per dataset for all genes**, so a panel of
#' G genes over D datasets costs D retrievals rather than G x D.
#'
#' This is the pooled counterpart of the per-dataset screening table: for each
#' gene the per-dataset estimates are pooled (inverse variance, random or fixed
#' effects) instead of merely being listed, and the gene-level p-values are
#' adjusted for multiple testing.
#' @param datasets Character vector of dataset accessions or TCGA projects.
#' @param genes Character vector of gene symbols (the panel).
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or a concrete token; the
#'   family is resolved per dataset and the token used is reported in
#'   \code{per_dataset$endpoint} and \code{endpoints}.
#' @param method Pooling method for \code{\link{cpas_meta}}: \code{"RE"}
#'   (random effects, DerSimonian-Laird, default) or \code{"FE"} (fixed effect).
#' @param confounders Optional covariates adjusted for inside every dataset.
#' @param min_events Minimum number of events a dataset must contribute
#'   (default 5); datasets below it are excluded with a recorded reason.
#' @param process_duplicates Probe collapsing rule used when retrieving
#'   expression (see \code{\link{get_expr_data}}).
#' @param fdr_method Method passed to \code{\link[stats]{p.adjust}} for the
#'   gene-level FDR (default \code{"BH"}).
#' @param merged Optional named list of already merged data frames
#'   (dataset -> data.frame), which skips all retrieval (used by the test suite
#'   and by callers who cache their own data).
#' @param max_try Number of retrieval attempts per dataset.
#' @param progress Print a progress message per gene.
#' @return List of class \code{cpas_meta_panel}:
#'   \item{\code{table}:}{one row per gene: \code{gene}, \code{k} (datasets
#'     pooled), \code{datasets_used}, \code{total_n}, \code{total_events},
#'     \code{HR}, \code{lower}, \code{upper}, \code{p}, \code{P_adj},
#'     \code{P_adj_text}, \code{I2}, \code{tau2}, \code{pi_lower}, \code{pi_upper},
#'     \code{p_heterogeneity}}
#'   \item{\code{per_dataset}:}{long table: \code{gene}, \code{dataset},
#'     \code{endpoint}, \code{n}, \code{events}, \code{HR}, \code{lower},
#'     \code{upper}, \code{p}}
#'   \item{\code{endpoints}:}{named list (gene -> named vector of tokens used)}
#'   \item{\code{errors}:}{reasons why a gene or a dataset was not pooled}
#'   \item{\code{fetch_errors}:}{datasets that could not be retrieved}
#'   \item{\code{input}:}{echoed inputs}
#' @details
#' Each gene is standardised within each dataset (per SD) before pooling, so the
#' hazard ratios are comparable across platforms, and the endpoint token each
#' dataset contributes is reported per gene (a dataset reporting RFS and another
#' reporting DFI are both usable under the DFS family; the mixed tokens are
#' listed). A gene with fewer than two usable datasets is reported with
#' \code{k = 1} and no pooled p-value, and the reason is recorded in
#' \code{errors}.
#' @export
#' @examples
#' \dontrun{
#'    ## Two real lung cohorts with a DFS-family endpoint (RFS where available).
#'    p <- cpas_meta_panel(c("GSE31210", "GSE37745"), genes = c("TP53", "GAPDH"),
#'                          type = "DFS")
#'    p$table
#'    plot_meta_panel(p)
#' }
cpas_meta_panel <- function(datasets, genes, type = "OS",
                            method = c("RE", "FE"),
                            confounders = NULL, min_events = 5,
                            process_duplicates = "max", fdr_method = "BH",
                            merged = NULL, max_try = 3, progress = TRUE) {
  method <- match.arg(method)
  if (missing(datasets) || !length(datasets))
    stop("'datasets' must contain at least one accession.")
  datasets <- unique(as.character(datasets))
  if (missing(genes) || !length(genes))
    stop("'genes' must contain at least one gene symbol.")
  genes <- unique(as.character(genes))

  # ---- retrieve once per dataset (all genes together) ----------------------
  fetch_errors <- list()
  if (is.null(merged)) {
    merged <- list()
    for (d in datasets) {
      got <- NULL
      for (i in seq_len(max_try)) {
        got <- tryCatch({
          e <- get_expr_data(d, genes, process_duplicates = process_duplicates)
          merge_surv_expr(d, e)$merged_data
        }, error = function(e) e)
        if (!inherits(got, "error")) break
        Sys.sleep(1 + i)
      }
      if (inherits(got, "error")) {
        fetch_errors[[d]] <- conditionMessage(got)
        next
      }
      merged[[d]] <- got
    }
  } else {
    if (is.null(names(merged)) || any(!nzchar(names(merged))))
      stop("'merged' must be a named list (dataset -> data.frame).")
    missing_ds <- setdiff(datasets, names(merged))
    if (length(missing_ds)) {
      datasets <- intersect(datasets, names(merged))
      for (d in missing_ds) fetch_errors[[d]] <- "not supplied in 'merged'"
    }
  }
  if (!length(merged)) {
    stop("No dataset could be retrieved: ",
         paste(sprintf("%s (%s)", names(fetch_errors), unlist(fetch_errors)),
               collapse = "; "))
  }

  # ---- one meta-analysis per gene ------------------------------------------
  rows <- list(); long <- list(); eps <- list(); errs <- list(); notes <- list()
  for (g in genes) {
    if (progress) message("cpas_meta_panel: ", g, " (",
                          match(g, genes), "/", length(genes), ")")
    res <- NULL
    withCallingHandlers(
      res <- tryCatch(
        cpas_meta(datasets = names(merged), marker = g, type = type,
                  method = method, confounders = confounders,
                  min_events = min_events, merged = merged,
                  max_try = max_try),
        error = function(e) e),
      warning = function(w) {
        notes[[g]] <<- unique(c(notes[[g]], conditionMessage(w)))
        tryCatch(invokeRestart("muffleWarning"), error = function(e) NULL)
      })
    if (inherits(res, "error")) {
      errs[[g]] <- conditionMessage(res)
      rows[[length(rows) + 1L]] <- data.frame(
        gene = g, k = 0L, datasets_used = "", total_n = NA_integer_,
        total_events = NA_integer_, HR = NA_real_, lower = NA_real_,
        upper = NA_real_, p = NA_real_, I2 = NA_real_, tau2 = NA_real_,
        pi_lower = NA_real_, pi_upper = NA_real_, p_heterogeneity = NA_real_,
        stringsAsFactors = FALSE)
      next
    }
    po <- res$pooled
    rows[[length(rows) + 1L]] <- data.frame(
      gene = g, k = po$k,
      datasets_used = paste(res$per_dataset$dataset, collapse = ","),
      total_n = po$total_n, total_events = po$total_events,
      HR = po$HR, lower = po$lower, upper = po$upper, p = po$p,
      I2 = po$I2, tau2 = po$tau2,
      pi_lower = po$pi_lower, pi_upper = po$pi_upper,
      p_heterogeneity = po$p_heterogeneity, stringsAsFactors = FALSE)
    long[[length(long) + 1L]] <- cbind(
      gene = g,
      res$per_dataset[, c("dataset", "endpoint", "n", "events", "HR",
                          "lower", "upper", "p")])
    eps[[g]] <- stats::setNames(res$per_dataset$endpoint, res$per_dataset$dataset)
    if (length(res$errors))
      errs[[g]] <- paste(names(res$errors), unlist(res$errors), sep = ": ",
                         collapse = "; ")
  }

  tbl <- do.call(rbind, rows)
  tbl$P_adj <- NA_real_
  # the FDR family contains only genes that were actually pooled (k >= 2): a
  # single-dataset row is that dataset's estimate, not a meta-analysis, and must
  # not enter the multiplicity correction
  ok <- !is.na(tbl$p) & tbl$k >= 2L
  if (any(ok)) tbl$P_adj[ok] <- stats::p.adjust(tbl$p[ok], method = fdr_method)
  tbl$P_adj_text <- ifelse(is.na(tbl$P_adj), "",
    ifelse(tbl$P_adj < 1e-4, formatC(tbl$P_adj, format = "g", digits = 3),
           formatC(tbl$P_adj, format = "f", digits = 4)))
  tbl <- tbl[order(tbl$p, na.last = TRUE), , drop = FALSE]
  rownames(tbl) <- NULL

  out <- list(input = list(datasets = datasets, genes = genes, type = type,
                           method = method, confounders = confounders,
                           min_events = min_events,
                           process_duplicates = process_duplicates,
                           fdr_method = fdr_method, time = Sys.time()),
              table = tbl,
              per_dataset = if (length(long)) do.call(rbind, long) else NULL,
              endpoints = eps,
              errors = errs,
              notes = notes,
              fetch_errors = fetch_errors,
              metadata = list(datasets_requested = length(datasets),
                              datasets_retrieved = length(merged),
                              genes_requested = length(genes),
                              genes_pooled = sum(tbl$k >= 2L, na.rm = TRUE),
                              genes_single_dataset = sum(tbl$k == 1L, na.rm = TRUE),
                              genes_failed = sum(tbl$k == 0L, na.rm = TRUE),
                              n_significant_fdr005 =
                                sum(tbl$P_adj < 0.05, na.rm = TRUE)))
  class(out) <- "cpas_meta_panel"
  out
}

#' @title Print method for cpas_meta_panel objects
#' @param x An object of class \code{cpas_meta_panel}.
#' @param n Number of genes to print (the table is ordered by p-value).
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_meta_panel <- function(x, n = 20, ...) {
  cat("CanPAS per-gene meta-analysis panel\n")
  cat(sprintf("datasets: %d/%d retrieved | genes: %d | pooled (k>=2): %d | single dataset: %d | failed: %d | FDR<0.05: %d\n",
              x$metadata$datasets_retrieved, x$metadata$datasets_requested,
              x$metadata$genes_requested, x$metadata$genes_pooled,
              x$metadata$genes_single_dataset, x$metadata$genes_failed,
              x$metadata$n_significant_fdr005))
  if (length(x$fetch_errors))
    cat("retrieval failures: ",
        paste(names(x$fetch_errors), collapse = ", "), "\n", sep = "")
  d <- utils::head(x$table, n)
  print(d[, c("gene", "k", "total_n", "total_events", "HR", "lower", "upper",
              "p", "P_adj", "I2")], row.names = FALSE)
  if (nrow(x$table) > n) cat("... ", nrow(x$table) - n, " more gene(s)\n", sep = "")
  invisible(x)
}

#' @title Plot a per-gene meta-analysis panel
#' @description Forest-style plot of the pooled hazard ratios of a
#' \code{\link{cpas_meta_panel}} result: one row per gene, point size scaled by
#' the number of pooled datasets, colour by FDR significance, with the
#' heterogeneity (I2) shown next to each estimate.
#' @param x An object returned by \code{\link{cpas_meta_panel}}.
#' @param digits Decimal places for the printed estimates (default 4).
#' @param annotate Add the per-gene I2 (and the FDR-adjusted p for the
#'   significant genes) to the axis labels.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#'    p <- cpas_meta_panel(c("GSE31210", "GSE37745"), genes = c("TP53", "GAPDH"),
#'                          type = "DFS")
#'    plot_meta_panel(p, annotate = TRUE, digits = 4)
#' }
plot_meta_panel <- function(x, annotate = TRUE, digits = 4) {
  if (!inherits(x, "cpas_meta_panel"))
    stop("'x' must be a 'cpas_meta_panel' object (see cpas_meta_panel()).")
  d <- x$table[!is.na(x$table$HR), , drop = FALSE]
  if (!nrow(d)) stop("No gene with an estimable pooled hazard ratio to plot.")
  d$sig <- ifelse(!is.na(d$P_adj) & d$P_adj < 0.05, "FDR < 0.05", "not significant")
  d$label <- if (isTRUE(annotate))
    sprintf("%s%s", d$gene,
            ifelse(is.na(d$I2), "", sprintf(" (I2=%.0f%%, k=%d)", 100 * d$I2, d$k)))
  else d$gene
  ggplot2::ggplot(d, ggplot2::aes(x = HR, y = stats::reorder(label, HR))) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                            height = 0.2, colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(size = k, colour = sig)) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "datasets pooled") +
    ggplot2::labs(title = sprintf("Per-gene meta-analysis (%s, %s)",
                                  x$input$type, x$input$method),
                  subtitle = sprintf("%d dataset(s) | per-SD within dataset | BH-FDR across genes",
                                     x$metadata$datasets_retrieved),
                  x = "Pooled hazard ratio (95% CI, log scale)", y = NULL,
                  colour = NULL) +
    ggplot2::theme_bw()
}
