## ---------------------------------------------------------------------------
## Group colours for the high/low Kaplan-Meier panels.
##
## Convention: HIGH expression is red, LOW expression is green. survminer's
## default hue palette assigns its first colour to the first factor level, and
## the levels are c("Low", "High"), so without an explicit palette the LOW group
## is drawn in the warm colour and the HIGH group in the cool one — the reverse
## of the convention above. Users can override the two colours with
## options(CanPAS.km_palette = c(low = "...", high = "...")) or per call by
## passing palette = c(...) through ... to plot_km().
## ---------------------------------------------------------------------------
.cpas_group_colours <- function() {
  p <- getOption("CanPAS.km_palette", NULL)
  if (is.null(p)) p <- c("#2CA02C", "#D62728")            # green, red
  if (length(p) != 2L) stop("'CanPAS.km_palette' must hold exactly two colours ",
                            "(low, high).", call. = FALSE)
  p <- as.character(p)
  names(p) <- c("Low", "High")
  p
}

#' @title Kaplan-Meier Survival Plot by Marker
#' @description
#' Splits samples into \code{High} / \code{Low} groups around a marker
#' cut-point (default: median) and draws a Kaplan-Meier curve with
#' \code{survminer::ggsurvplot}.
#' @param df Data.frame with columns \code{ID}, \code{<type>_time},
#' \code{<type>_status} and \code{marker}.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token (default \code{"OS"}).
#' @param marker Name of the numeric marker column used for stratification.
#' @param cutpoint Cut-point rule: a number (fixed threshold), \code{"median"}
#' (default) or \code{"mean"}.
#' @param ... Further arguments passed to \code{survminer::ggsurvplot}
#' (e.g. \code{pval = TRUE}, \code{risk.table = TRUE}, \code{legend = "bottom"}).
#' @return A \code{ggsurvplot} object (access the ggplot via \code{$plot}).
#' @details Rows with missing time/status/marker are removed; groups are
#' \code{High = marker > cutpoint} and \code{Low = marker <= cutpoint}.
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#'    ## The prerequisite chain returns one analysis-ready table per cohort.
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
#' 
#'    plot_km(d, type = "OS", marker = "GAPDH", pval = TRUE)
#' 
#'    ## A signature or a specific cutpoint works on the same table.
#'    plot_km(d, type = "OS", marker = "GAPDH", cutpoint = "median")
#' }
plot_km <- function(df, type = "OS", marker, cutpoint = "median", ...) {
  type <- .resolve_type_in_df(df, type)          # family or token
  tc <- paste0(type, "_time"); sc <- paste0(type, "_status")
  if (!all(c("ID", tc, sc, marker) %in% colnames(df)))
    stop("'df' must contain columns ID, ", tc, ", ", sc, " and ", marker, ".")
  dat <- df[c("ID", tc, sc, marker)]
  dat[[tc]] <- suppressWarnings(as.numeric(dat[[tc]]))
  dat[[sc]] <- suppressWarnings(as.numeric(dat[[sc]]))
  dat[[marker]] <- suppressWarnings(as.numeric(dat[[marker]]))
  colnames(dat) <- c("ID", "time", "status", "marker")
  if (any(!is.na(dat$status) & !dat$status %in% c(0, 1)))
    stop("The status column must be coded 0 (censored) / 1 (event).")
  if (any(!is.na(dat$time) & dat$time < 0))
    stop("The time column contains negative values.")
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < 5L) stop("Too few complete rows for a KM plot.")
  if (length(unique(dat$status)) < 2L || sum(dat$status == 1) < 1L)
    stop("Need both censored and event observations to draw a KM curve.")

  cp <- if (is.numeric(cutpoint)) {
    if (length(cutpoint) != 1L || !is.finite(cutpoint))
      stop("'cutpoint' must be a single finite number, 'median' or 'mean'.")
    as.numeric(cutpoint)
  } else if (identical(cutpoint, "median")) stats::median(dat$marker)
    else if (identical(cutpoint, "mean")) mean(dat$marker)
    else stop("'cutpoint' must be a number, 'median' or 'mean'.")
  dat$group <- factor(ifelse(dat$marker > cp, "High", "Low"),
                      levels = c("Low", "High"))
  if (length(unique(dat$group)) < 2L)
    stop("The chosen cutpoint yields only one group.")
  grp_n <- table(dat$group)
  grp_ev <- tapply(dat$status, dat$group, sum)
  if (min(grp_n) < 5L)
    warning(sprintf(paste0("The smallest group has only %d patient(s) (%d event(s)); ",
                           "the log-rank p-value and the curve are unstable."),
                    min(grp_n), min(grp_ev)), call. = FALSE)

  fit <- survfit(Surv(time, status) ~ group, data = dat)
  ## high = red, low = green unless the caller supplies palette = ... explicitly
  dots <- list(...)
  if (!"palette" %in% names(dots)) dots$palette <- unname(.cpas_group_colours())
  p <- do.call(ggsurvplot, c(list(fit = fit, data = dat), dots))
  p$plot <- p$plot +
    ggtitle(marker) +
    theme(plot.title = element_text(hjust = 0.5))
  p
}
