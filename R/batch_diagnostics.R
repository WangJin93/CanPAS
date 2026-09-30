# batch_diagnostics.R --------------------------------------------------------
# Per-cohort expression diagnostics for a multi-cohort (batch) analysis (spec B5).
#
# Why this exists.  A pooled analysis over several GEO/TCGA cohorts mixes
# platforms, so the absolute expression of a gene is not comparable between
# cohorts.  Two facts about the pooled machinery make that harmless, and both
# are visible here:
#
#   * the two-stage meta-analysis of cpas_meta() pools WITHIN-COHORT effect
#     sizes: every cohort is standardised inside itself (a per-SD log-HR) before
#     the inverse-variance pool, so a constant shift of the whole cohort enters
#     neither the point estimate nor its variance;
#   * the pooled Kaplan-Meier of cpas_km_pooled() splits each cohort at its OWN
#     percentile (median or top-percent), never at one absolute value.
#
# batch_diagnostics() shows, per cohort and gene, the distribution that the
# standardisation and the cut-points are computed from, and plot_batch_diagnostics()
# puts those distributions side by side on one axis so the cross-platform shift
# is visible next to the cut-points that absorb it.
#
# Reads cohorts through the existing accessors/caches (get_expr_data() /
# merge_surv_expr() via cohort_merged(), all cached), so a previously fetched
# cohort is analysed offline; a caller may also pass 'merged' and touch nothing.

#' @title Per-cohort expression diagnostics for a multi-cohort analysis
#' @description
#' Summarises the expression distribution of one gene (or a gene vector) in each
#' cohort, together with the cut-point values the pooled analyses use, so a
#' batch effect can be inspected before anything is pooled.
#'
#' \strong{Multi-cohort pooling uses within-cohort standardised effect sizes.}
#' \code{\link{cpas_meta}} is a two-stage meta-analysis: stage 1 fits the Cox
#' model inside each cohort, standardises the marker within that cohort and
#' returns a log-HR and its SE; stage 2 pools those per-cohort estimates by
#' inverse variance.  A cohort-wide additive shift of the expression scale
#' (the usual platform batch effect) cancels inside the cohort and therefore
#' \emph{does not enter the pooled estimate}.  The pooled \strong{Kaplan-Meier}
#' of \code{\link{cpas_km_pooled}} splits each cohort at a \strong{within-cohort
#' percentile} (the cohort's own median, or its own \eqn{(100-p)}th percentile),
#' so again no absolute cross-platform threshold is used.  The summaries this
#' function returns make both facts checkable: \code{median} and the two
#' cut-point columns are the values used, and the standardised columns are what
#' actually enters the pool.
#' @param datasets Character vector of dataset accessions or TCGA projects.
#' @param genes One or more gene symbols.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or a concrete token. It only
#'   selects which survival columns are joined and which token is reported; the
#'   expression summaries are independent of it.
#' @param top_pct Percentile used by the \code{cut = "top_pct"} rule of
#'   \code{\link{cpas_km_pooled}}: the cut-point reported is the cohort's
#'   \eqn{(100 - top\_pct)}th percentile, i.e. \code{top_pct = 25} means the top
#'   25\% of the cohort is "High" (default 25).
#' @param merged Optional named list (cohort -> merged data.frame) that skips all
#'   retrieval; useful offline and in tests.
#' @param process_duplicates Probe collapsing rule passed to
#'   \code{\link{get_expr_data}} (default \code{"max"}).
#' @param max_try Retrieval attempts per cohort.
#' @return A list of class \code{cpas_batch_diagnostics}:
#'   \item{\code{table}:}{one row per cohort and gene: \code{dataset},
#'     \code{gene}, \code{endpoint} (the token resolved for that cohort),
#'     \code{n} (non-missing expression values), \code{n_missing},
#'     \code{median}, \code{q1}, \code{q3}, \code{iqr}, \code{min}, \code{max},
#'     \code{mean}, \code{sd}, \code{cut_median} (the median-split cut-point
#'     used), \code{cut_top_pct} (the percentile cut-point used),
#'     \code{n_high_median}, \code{n_low_median}, \code{n_high_top_pct},
#'     \code{n_low_top_pct}, and \code{median_z} (the cohort median on the
#'     within-cohort standardised scale, i.e. what the pool compares)}
#'   \item{\code{cohorts}:}{one row per cohort: how many genes were found, the
#'     token resolved and the sample count}
#'   \item{\code{values}:}{the long per-sample values the plot and the summaries
#'     are built from: \code{dataset}, \code{gene}, \code{value}}
#'   \item{\code{batch}:}{the cross-cohort spread: per gene, \code{min_median},
#'     \code{max_median}, \code{median_shift} (max minus min, the absolute
#'     cross-platform difference in expression units), and
#'     \code{median_z_shift} (the same shift after within-cohort
#'     standardisation)}
#'   \item{\code{errors}, \code{input}, \code{manifest}:}{cohorts that could not
#'     be read, the echoed inputs, and the analysis manifest}
#' @details
#' The cut-point columns are computed exactly as \code{\link{cpas_km_pooled}}
#' computes them, on the same rows: the median-split value is
#' \code{stats::median(value)} and the percentile value is
#' \code{stats::quantile(value, 1 - top_pct/100)} with the default type-7
#' interpolation.  They are reported per cohort, never pooled, because pooling an
#' absolute expression threshold across platforms is exactly the mistake the
#' within-cohort rules avoid.
#'
#' Cohorts are read through \code{\link{cohort_merged}}, which goes through the
#' same cached accessors as the analyses (\code{\link{get_expr_data}},
#' \code{\link{merge_surv_expr}}, \code{\link{get_data}},
#' \code{\link{tcga_surv_table}}).  A cohort whose tables are already in the
#' local cache is therefore analysed without a network call; passing
#' \code{merged} skips retrieval entirely.
#' @seealso \code{\link{plot_batch_diagnostics}}, \code{\link{cpas_meta}},
#'   \code{\link{cpas_km_pooled}}, \code{\link{cohort_merged}}
#' @export
#' @examples
#' \dontrun{
#'    ## Two real lung cohorts; the second run reads the cache offline.
#'    b <- batch_diagnostics(c("GSE31210", "GSE37745"), genes = "GAPDH",
#'                           type = "DFS", top_pct = 25)
#'    b$table[, c("dataset", "gene", "n", "median", "cut_median", "cut_top_pct")]
#'    b$batch                      # how far apart the absolute scales are
#'    plot_batch_diagnostics(b)    # the distributions side by side
#' }
batch_diagnostics <- function(datasets, genes, type = "OS", top_pct = 25,
                              merged = NULL, process_duplicates = "max",
                              max_try = 3) {
  if (missing(datasets) || !length(datasets))
    stop("'datasets' must contain at least one accession.", call. = FALSE)
  datasets <- unique(as.character(datasets))
  if (missing(genes) || !length(genes))
    stop("'genes' must contain at least one gene symbol.", call. = FALSE)
  genes <- unique(as.character(genes))
  if (!is.numeric(top_pct) || length(top_pct) != 1L || !is.finite(top_pct) ||
      top_pct <= 0 || top_pct >= 100)
    stop("'top_pct' must be a single finite number between 1 and 99.", call. = FALSE)
  prob <- 1 - top_pct / 100

  errors <- list()
  if (is.null(merged)) {
    merged <- list()
    for (d in datasets) {
      got <- NULL
      for (i in seq_len(max_try)) {
        got <- tryCatch(cohort_merged(d, genes, type = type,
                                      process_duplicates = process_duplicates),
                        error = function(e) e)
        if (!inherits(got, "error")) break
        Sys.sleep(min(1 + i, 3))
      }
      if (inherits(got, "error")) {
        errors[[d]] <- conditionMessage(got)
        next
      }
      merged[[d]] <- got
    }
  } else {
    if (is.null(names(merged)) || any(!nzchar(names(merged))))
      stop("'merged' must be a named list (dataset -> data.frame).", call. = FALSE)
    for (d in setdiff(datasets, names(merged)))
      errors[[d]] <- "not supplied in 'merged'"
    merged <- merged[intersect(datasets, names(merged))]
  }
  if (!length(merged))
    stop("No cohort could be read: ",
         paste(sprintf("%s (%s)", names(errors), unlist(errors)), collapse = "; "),
         call. = FALSE)

  rows <- list(); vals <- list(); cohort_rows <- list()
  for (d in names(merged)) {
    dd <- merged[[d]]
    if (!is.data.frame(dd)) {
      errors[[d]] <- "not a data.frame"
      next
    }
    tok <- tryCatch(endpoint_resolve(d, type), error = function(e) NA_character_)
    for (g in genes) {
      if (!g %in% colnames(dd)) {
        errors[[paste(d, g, sep = "/")]] <- sprintf("gene column '%s' absent", g)
        next
      }
      v <- suppressWarnings(as.numeric(dd[[g]]))
      n <- sum(!is.na(v))
      if (n < 1L) {
        errors[[paste(d, g, sep = "/")]] <- "no non-missing expression value"
        next
      }
      vv <- v[!is.na(v)]
      q <- stats::quantile(vv, probs = c(0.25, 0.5, 0.75, prob), names = FALSE)
      med <- q[2]
      thr <- q[4]
      mu <- mean(vv); s <- stats::sd(vv)
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = d, gene = g,
        endpoint = if (is.null(tok) || is.na(tok)) NA_character_ else as.character(tok),
        n = as.integer(n), n_missing = as.integer(length(v) - n),
        median = med, q1 = q[1], q3 = q[3], iqr = q[3] - q[1],
        min = min(vv), max = max(vv), mean = mu, sd = s,
        cut_median = med, cut_top_pct = thr,
        n_high_median = as.integer(sum(vv > med)),
        n_low_median = as.integer(sum(vv <= med)),
        n_high_top_pct = as.integer(sum(vv > thr)),
        n_low_top_pct = as.integer(sum(vv <= thr)),
        median_z = if (is.finite(s) && s > 0) (med - mu) / s else NA_real_,
        stringsAsFactors = FALSE)
      vals[[length(vals) + 1L]] <- data.frame(dataset = d, gene = g, value = vv,
                                              stringsAsFactors = FALSE)
    }
    cohort_rows[[length(cohort_rows) + 1L]] <- data.frame(
      dataset = d, n_samples = nrow(dd),
      genes_found = sum(genes %in% colnames(dd)),
      endpoint = if (is.null(tok) || is.na(tok)) NA_character_ else as.character(tok),
      stringsAsFactors = FALSE)
  }
  if (!length(rows))
    stop("No gene could be summarised in any cohort; see the errors.", call. = FALSE)
  tbl <- do.call(rbind, rows); rownames(tbl) <- NULL
  vdf <- do.call(rbind, vals); rownames(vdf) <- NULL
  cdf <- do.call(rbind, cohort_rows); rownames(cdf) <- NULL

  batch <- do.call(rbind, lapply(split(tbl, tbl$gene), function(x) {
    data.frame(gene = x$gene[1],
               n_cohorts = nrow(x),
               min_median = min(x$median), max_median = max(x$median),
               median_shift = max(x$median) - min(x$median),
               median_z_shift = max(x$median_z, na.rm = TRUE) -
                 min(x$median_z, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(batch) <- NULL

  out <- list(input = list(datasets = datasets, genes = genes, type = type,
                           top_pct = top_pct, process_duplicates = process_duplicates,
                           time = Sys.time()),
              table = tbl, cohorts = cdf, values = vdf, batch = batch,
              errors = errors,
              manifest = .cpas_manifest_new(
                analysis = "batch_diagnostics",
                cohorts = names(merged),
                family = tryCatch(endpoint_family(type), error = function(e) NA_character_),
                token = as.character(type)[1],
                selection_rule = sprintf(paste0("expression was summarised per cohort for %d gene(s) ",
                                                "over %d cohort(s); the cut-point columns are the ",
                                                "within-cohort values cpas_km_pooled() would use"),
                                         length(genes), length(merged)),
                dropped_rows = if (length(errors)) data.frame(
                  cohort = names(errors), n_dropped = NA_integer_,
                  reason = unlist(errors, use.names = FALSE), stringsAsFactors = FALSE) else NULL,
                cut_rule = sprintf(paste0("within-cohort percentile cut-points only: median split ",
                                          "(High = value > cohort median) and top-%g%% split ",
                                          "(High = value > cohort's %gth percentile); no absolute ",
                                          "cross-cohort threshold is used"),
                                   top_pct, 100 - top_pct),
                cut_points_searched = 0L,
                auto_repair = FALSE,
                notes = c(paste0("pooling is over within-cohort standardised effect sizes ",
                                 "(two-stage meta-analysis), so a constant cohort-wide shift ",
                                 "of the expression scale does not enter the pooled estimate"),
                          paste0("pooled Kaplan-Meier splits each cohort at its own percentile, ",
                                 "so absolute cross-platform expression differences do not enter ",
                                 "the pooled curve"),
                          if (length(errors)) sprintf("%d cohort/gene combination(s) could not be summarised",
                                                      length(errors)) else NULL)))
  class(out) <- "cpas_batch_diagnostics"
  out
}

#' @title Plot per-cohort expression distributions
#' @description
#' Draws the per-cohort expression distribution of every gene in a
#' \code{\link{batch_diagnostics}} result as side-by-side box plots on one
#' shared axis, with the within-cohort median-split and percentile cut-points
#' marked, so the cross-platform shift and the cut-points that absorb it are
#' visible in the same panel.
#' @param x An object returned by \code{\link{batch_diagnostics}}.
#' @param genes Optional subset of genes to draw (default: all).
#' @param show_cut_points Draw the per-cohort median (solid segment) and
#'   percentile (dashed segment) cut-points (default \code{TRUE}).
#' @param free_y Let each gene panel keep its own y scale (default \code{FALSE}:
#'   one shared axis, which is what shows the cross-cohort differences).
#' @param ... Reserved, currently unused.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{batch_diagnostics}}
#' @export
#' @examples
#' \dontrun{
#'    b <- batch_diagnostics(c("GSE31210", "GSE37745"), genes = "GAPDH", type = "DFS")
#'    plot_batch_diagnostics(b)
#' }
plot_batch_diagnostics <- function(x, genes = NULL, show_cut_points = TRUE,
                                   free_y = FALSE, ...) {
  if (!inherits(x, "cpas_batch_diagnostics"))
    stop("'x' must be a 'cpas_batch_diagnostics' object (see batch_diagnostics()).",
         call. = FALSE)
  d <- x$values
  if (!is.null(genes)) d <- d[d$gene %in% as.character(genes), , drop = FALSE]
  if (!nrow(d)) stop("No value to plot for the requested gene(s).")
  d$dataset <- factor(d$dataset, levels = unique(d$dataset))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = dataset, y = value)) +
    ggplot2::geom_boxplot(outlier.size = 0.6, width = 0.7,
                          fill = "grey92", colour = "grey25")
  if (isTRUE(show_cut_points)) {
    cut <- x$table[x$table$gene %in% unique(d$gene), , drop = FALSE]
    cut$dataset <- factor(cut$dataset, levels = levels(d$dataset))
    p <- p +
      ggplot2::geom_segment(data = cut,
                            ggplot2::aes(x = as.numeric(dataset) - 0.35,
                                         xend = as.numeric(dataset) + 0.35,
                                         y = cut_median, yend = cut_median),
                            inherit.aes = FALSE, colour = "firebrick", linewidth = 0.6) +
      ggplot2::geom_segment(data = cut,
                            ggplot2::aes(x = as.numeric(dataset) - 0.35,
                                         xend = as.numeric(dataset) + 0.35,
                                         y = cut_top_pct, yend = cut_top_pct),
                            inherit.aes = FALSE, colour = "steelblue",
                            linetype = 2, linewidth = 0.6)
  }
  if (length(unique(d$gene)) > 1L)
    p <- p + ggplot2::facet_wrap(~gene, scales = if (isTRUE(free_y)) "free_y" else "fixed")
  else
    p <- p + ggplot2::labs(subtitle = unique(d$gene))
  p +
    ggplot2::labs(
      title = "Per-cohort expression distributions",
      x = NULL, y = "Expression (log2 scale)",
      caption = paste0("segments: within-cohort median split (solid) and top-percentile ",
                       "cut-point (dashed); pooling is over within-cohort standardised ",
                       "effect sizes, so the absolute shift between cohorts does not ",
                       "enter the pooled estimate")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' @title Print a batch-diagnostics result
#' @param x An object of class \code{cpas_batch_diagnostics}.
#' @param n Number of cohort/gene rows to print (default 20).
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_batch_diagnostics <- function(x, n = 20, ...) {
  cat("CanPAS batch diagnostics\n")
  cat(sprintf("cohorts: %d | genes: %d | cohort/gene rows: %d\n",
              nrow(x$cohorts), length(x$input$genes), nrow(x$table)))
  d <- utils::head(x$table, n)
  print(d[, c("dataset", "gene", "n", "median", "iqr", "min", "max",
              "cut_median", "cut_top_pct")], row.names = FALSE)
  if (nrow(x$table) > n) cat("... ", nrow(x$table) - n, " more row(s)\n", sep = "")
  print(x$batch, row.names = FALSE)
  if (length(x$errors))
    cat("errors: ", paste(names(x$errors), collapse = ", "), "\n", sep = "")
  invisible(x)
}
