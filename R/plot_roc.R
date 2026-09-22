#' @title Time-dependent ROC Curve for a Marker
#' @description
#' Computes a time-dependent ROC curve at a given prediction horizon with
#' \code{survivalROC::survivalROC} and draws it with ggplot2.
#' @param df Data.frame with columns \code{<type>_time}, \code{<type>_status}
#' and \code{marker}.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token (default \code{"OS"}).
#' @param marker Name of the numeric marker column.
#' @param predict.time Prediction horizon (same unit as \code{<type>_time},
#' typically years). Default 1.
#' @param method ROC method passed to \code{survivalROC}: \code{"KM"} (default)
#' or \code{"NNE"}.
#' @param span Span for \code{method = "NNE"} (ignored for KM).
#' @return A ggplot object; the estimated AUC is attached as attribute
#' \code{"AUC"} of the returned object.
#' @details Rows with missing time/status/marker are removed first. The
#' \code{<type>_time} values are interpreted as the same unit as
#' \code{predict.time}.
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
#' 
#'    g <- plot_roc(d, type = "OS", marker = "GAPDH", predict.time = 3)
#'    g$auc
#' }
plot_roc <- function(df, type = "OS", marker,
                     predict.time = 1, method = c("KM", "NNE"), span = NULL) {
  method <- match.arg(method)
  type <- .resolve_type_in_df(df, type)          # family or token
  tc <- paste0(type, "_time"); sc <- paste0(type, "_status")
  if (!all(c(tc, sc, marker) %in% colnames(df)))
    stop("'df' must contain columns ", tc, ", ", sc, " and ", marker, ".")
  d <- df[c(tc, sc, marker)]
  d[[tc]] <- suppressWarnings(as.numeric(d[[tc]]))
  d[[sc]] <- suppressWarnings(as.numeric(d[[sc]]))
  d[[marker]] <- suppressWarnings(as.numeric(d[[marker]]))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (any(!is.na(d[[sc]]) & !d[[sc]] %in% c(0, 1)))
    stop("The status column must be coded 0 (censored) / 1 (event).")
  if (any(!is.na(d[[tc]]) & d[[tc]] < 0))
    stop("The time column contains negative values.")
  if (nrow(d) < 10L) stop("Too few complete rows for a ROC analysis.")
  if (sum(d[[sc]] == 1) < 1L) stop("No events available for ROC analysis.")
  if (predict.time <= 0) stop("'predict.time' must be positive.")
  # a cumulative/dynamic ROC needs events before, and controls after, the horizon
  tmax <- max(d[[tc]], na.rm = TRUE)
  if (predict.time > tmax)
    stop(sprintf(paste0("'predict.time' = %.3g years exceeds the longest follow-up in this ",
                        "cohort (%.3g years): the time-dependent AUC is not estimable."),
                 predict.time, tmax))
  n_case <- sum(d[[sc]] == 1 & d[[tc]] <= predict.time)
  n_ctrl <- sum(d[[tc]] > predict.time)
  if (n_case < 1L)
    stop(sprintf(paste0("No event occurs before predict.time = %.3g years (n = %d events ",
                        "in total); the time-dependent AUC is not estimable."),
                 predict.time, sum(d[[sc]] == 1)))
  if (n_ctrl < 1L)
    stop(sprintf(paste0("No patient is still at risk beyond predict.time = %.3g years ",
                        "(no controls); the time-dependent AUC is not estimable."), predict.time))

  args <- list(Stime = d[[tc]], status = d[[sc]], marker = d[[marker]],
               predict.time = predict.time, method = method)
  if (method == "NNE") args$span <- if (is.null(span)) 0.25 * nrow(d)^(-0.20) else span
  roc_obj <- tryCatch(do.call(survivalROC::survivalROC, args),
                      error = function(e) stop("survivalROC failed: ",
                                               conditionMessage(e), call. = FALSE))
  if (!is.finite(roc_obj$AUC))
    stop("The time-dependent AUC could not be estimated at predict.time = ",
         predict.time, " (no usable cases/controls at this horizon).")
  msg <- sprintf("t = %.3g y: %d cases (event <= t), %d controls (follow-up > t), AUC = %.3f",
                 predict.time, n_case, n_ctrl, roc_obj$AUC)

  roc_data <- data.frame(FPR = roc_obj$FP, TPR = roc_obj$TP)
  roc_data <- roc_data[order(roc_data$FPR, roc_data$TPR), , drop = FALSE]

  p <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)",
         title = paste0(marker, " (", type, ", t = ", predict.time, ")")) +
    annotate("text", x = 0.65, y = 0.15,
             label = sprintf("AUC = %.3f", roc_obj$AUC), size = 4.5) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  attr(p, "AUC") <- roc_obj$AUC
  p
}
