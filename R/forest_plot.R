#' @title Forest Plot of Cox Regression Results
#' @description
#' Draws a forest plot (point = HR, error bar = 95\% CI) from a tidy Cox
#' results table as produced by \code{\link{COX_analysis}} or
#' \code{\link{COX_screen_adjust}}. Rows are coloured by significance and
#' risk direction; categorical variable names are shown as group headers.
#' @param COX_out A data.frame with (at least) columns \code{Variates} (or
#'   \code{gene}, as returned by \code{\link{COX_by_genes}}),
#' \code{HR}, \code{HR95L}, \code{HR95H}, \code{Pvalue}; optionally columns
#' \code{Type} (\code{cont} / \code{cate_header} / \code{cate_level}) and
#' \code{Level}. If the columns \code{HR_uni} ... \code{Pvalue_uni} from
#' \code{\link{COX_screen_adjust}} are present instead, use
#' \code{which = "uni"} or \code{"multi"}.
#' @param which Which set of columns to plot when \code{COX_out} contains both
#' univariate and multivariate estimates (\code{"uni"} or \code{"multi"}).
#' @param HR_threshold Reference value for risk direction (dashed vertical
#' line). Default 1.
#' @param p_threshold Significance threshold for colouring. Default 0.05.
#' @param colors Vector of three colours: [low-risk, non-significant,
#' high-risk]. Defaults \code{c("#00AFBB", "black", "#FC4E07")}.
#' @param log_x Whether to plot the x-axis on a log2 scale. Default TRUE.
#' @param HR_order Ordering of estimator rows: \code{"none"} (input order),
#' \code{"decrease"} (default) or \code{"increase"} by HR.
#' @param group_levels How the rows of a categorical covariate stay
#'   identifiable once the rows are reordered (only used when \code{HR_order}
#'   is not \code{"none"}). \code{"label"} (default) sorts strictly by HR: a
#'   categorical level is then no longer adjacent to its covariate name, so
#'   every level is labelled \code{"covariate: level"} and the separate header
#'   rows are dropped. \code{"block"} keeps the levels of each covariate
#'   together, sorts them inside their block and orders the blocks by HR (by the
#'   most extreme level in the requested direction), so the bold covariate
#'   header still sits directly above its own levels.
#' @param point_size Point size.
#' @param digits Decimal places for the printed estimates (default 4).
#' @param bar_height Height of the error bars.
#' @return A ggplot object. The plotted rows and their order are in the
#'   returned object's \code{data}, so a caller can verify that each label and
#'   its estimate belong together.
#' @details
#' A categorical covariate contributes one row per level, and a level name alone
#' ("G3") is only meaningful next to its covariate name. Sorting rows by HR
#' therefore interleaves covariates; the previous implementation re-attached the
#' covariate headers by moving whole position ranges, which could place a level
#' under another covariate's header. Now either the level labels carry the
#' covariate name (\code{group_levels = "label"}, default, strictly sorted) or
#' the covariate blocks stay contiguous (\code{group_levels = "block"}).
#' \code{HR_order = "none"} leaves the input order untouched, header rows
#' included.
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#'    ## Fit a real model first (expression + survival + clinical covariates
#'    ## come from the CanPAS mirror through cohort_merged()).
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS", clin = TRUE)
#'    r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
#'                       cate_Variates = c("sex", "stage"), method = "uni")
#' 
#'    forest_plot(r$results_table)
#' 
#'    ## The same figure for the multivariable model, on a linear scale.
#'    r2 <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
#'                        cate_Variates = c("sex", "stage"), method = "multi")
#'    forest_plot(r2$results_table, log_x = FALSE, p_threshold = 0.01)
#' }
forest_plot <- function(COX_out,
                        which = c("uni", "multi"),
                        HR_threshold = 1,
                        p_threshold = 0.05,
                        colors = c("#00AFBB", "black", "#FC4E07"),
                        log_x = TRUE,
                        HR_order = c("decrease", "increase", "none"),
                        group_levels = c("label", "block"),
                        point_size = 2.2,
                        bar_height = 0.3,
                        digits = 4) {
  which <- match.arg(which); HR_order <- match.arg(HR_order)
  group_levels <- match.arg(group_levels)
  if (!is.data.frame(COX_out)) stop("'COX_out' must be a data.frame.")

  d <- COX_out
  if (!"Variates" %in% colnames(d) && "gene" %in% colnames(d))
    d$Variates <- d$gene                     # COX_by_genes() output
  if (!"HR" %in% colnames(d)) {              # COX_screen_adjust long output
    suf <- if (which == "multi") "_multi" else "_uni"
    need <- c("HR", "HR95L", "HR95H", "Pvalue")
    cols <- paste0(need, suf)
    if (!all(cols %in% colnames(d)))
      stop("'COX_out' must contain HR/HR95L/HR95H/Pvalue (or *_uni/*_multi variants).")
    d$HR <- d[[cols[1]]]; d$HR95L <- d[[cols[2]]]
    d$HR95H <- d[[cols[3]]]; d$Pvalue <- d[[cols[4]]]
    if (!"Variates" %in% colnames(d)) {
      d$Variates <- if ("gene" %in% colnames(d)) d$gene else d$Var1
    }
    if (!"Type" %in% colnames(d)) d$Type <- ifelse(is.na(d$Level), "cont", "cate_level")
  } else {
    if (!all(c("Variates", "HR", "HR95L", "HR95H", "Pvalue") %in% colnames(d)))
      stop("'COX_out' must contain columns Variates, HR, HR95L, HR95H, Pvalue.")
    if (!"Type" %in% colnames(d)) d$Type <- "cont"
    if (!"Level" %in% colnames(d)) d$Level <- NA_character_
  }

  for (cn in c("HR", "HR95L", "HR95H", "Pvalue"))
    d[[cn]] <- suppressWarnings(as.numeric(d[[cn]]))

  est_rows <- !is.na(d$HR) & is.finite(d$HR)
  if (!any(est_rows)) stop("No rows with a numeric HR were found.")
  if (any(d$HR <= 0, na.rm = TRUE))
    stop("HR must be strictly positive (log scale).")

  keep_rows <- est_rows | (!is.na(d$Type) & d$Type == "cate_header")
  d <- d[keep_rows, , drop = FALSE]

  est <- d[is.na(d$Type) | d$Type != "cate_header", , drop = FALSE]
  headers <- d[!is.na(d$Type) & d$Type == "cate_header", , drop = FALSE]

  sort_rows <- function(x, decreasing) {
    if (HR_order == "none") return(x)
    x[order(x$HR, decreasing = decreasing), , drop = FALSE]
  }

  if (HR_order == "none") {
    # rows exactly as supplied: the caller's ordering is authoritative and the
    # header rows stay where they are
    d <- d
  } else if (identical(group_levels, "label")) {
    # strict HR order. The rows of a categorical covariate are no longer
    # adjacent, so relying on a header row above them would mislabel them
    # (the old code moved a positional RANGE and could push a level under
    # another variable's header). Every label therefore names its variable,
    # and the separate header rows are dropped.
    d <- sort_rows(est, HR_order == "decrease")
  } else {
    # "block": the levels of a categorical covariate stay together, the blocks
    # are ordered by HR (by the most extreme level in the requested direction)
    # and the levels are sorted inside their block, so each header still sits
    # directly above its own levels
    dec <- HR_order == "decrease"
    key <- vapply(split(est$HR, est$Var1), function(v)
      if (dec) max(v, na.rm = TRUE) else min(v, na.rm = TRUE), numeric(1))
    var_order <- names(key)[order(key, decreasing = dec)]
    blocks <- list()
    for (v in var_order) {
      h <- headers[headers$Var1 == v, , drop = FALSE]
      if (nrow(h)) blocks[[length(blocks) + 1L]] <- h[1, , drop = FALSE]
      rows_v <- est[est$Var1 == v, , drop = FALSE]
      if (nrow(rows_v)) blocks[[length(blocks) + 1L]] <- sort_rows(rows_v, dec)
    }
    left <- headers[!headers$Var1 %in% var_order, , drop = FALSE]
    if (nrow(left)) blocks[[length(blocks) + 1L]] <- left
    d <- do.call(rbind, blocks)
  }
  d$row <- seq_len(nrow(d))

  # level labels name their variable whenever the rows are reordered
  d$label <- ifelse(is.na(d$Level) | d$Level == "" | d$Type == "cate_header",
                    d$Variates, d$Level)
  if (HR_order != "none" && identical(group_levels, "label"))
    d$label <- ifelse(d$Type == "cate_level",
                      paste0(d$Variates, ": ", d$label), d$label)
  d$label <- ifelse(d$Type == "cate_header",
                    paste0("**", d$label, "**"), d$label)

  d$col_group <- ifelse(d$Type == "cate_header", "header",
    ifelse(!is.na(d$Pvalue) & d$Pvalue < p_threshold & d$HR > HR_threshold,
           "high",
    ifelse(!is.na(d$Pvalue) & d$Pvalue < p_threshold & d$HR < HR_threshold,
           "low", "ns")))

  p <- ggplot(d, aes(x = HR, y = row)) +
    geom_vline(xintercept = HR_threshold, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmax = HR95H, xmin = HR95L, colour = col_group),
                   height = bar_height, show.legend = FALSE) +
    geom_point(aes(colour = col_group), size = point_size, show.legend = FALSE) +
    scale_colour_manual(values = c(header = "transparent", high = colors[3],
                                   low = colors[1], ns = colors[2])) +
    scale_y_continuous(breaks = d$row, labels = d$label, limits = c(0.5, nrow(d) + 0.5),
                       expand = expansion(mult = c(0.02, 0.08))) +
    theme_bw() +
    labs(x = paste0("Hazard Ratio (", if (log_x) "log2" else "linear", " scale)"),
         y = NULL) +
    theme(axis.text.y = ggtext::element_markdown(size = 12),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 11))
  if (log_x) p <- p + coord_trans(x = "log2")
  p
}
