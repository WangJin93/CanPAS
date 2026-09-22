#' @title Univariate or Multivariate Cox Regression
#' @description
#' Fits Cox proportional-hazards models on a data.frame of survival plus
#' continuous and/or categorical covariates, and formats a tidy results table
#' (hazard ratios, 95\% CIs, Wald p-values, sample sizes per variable).
#' @param df Data.frame. If \code{type} is given it must contain \code{ID},
#' \code{<type>_time} and \code{<type>_status}; if \code{type = NULL} the
#' columns must already be named \code{time} and \code{status}.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token; the family
#' is resolved to the token available in \code{df}. \code{NULL} when \code{df}
#' already contains \code{time}/\code{status}. Default \code{NULL}.
#' @param cont_Variates Character vector of continuous covariate columns.
#' @param cate_Variates Character vector of categorical covariate columns
#' (the first factor level is the reference level).
#' @param method Either \code{"uni"} (one model per covariate) or
#' \code{"multi"} (single multivariable model with all covariates).
#' @param precision Decimal places used to format the display columns
#' (\code{HR_text}, \code{P_text}); numeric columns keep full precision.
#' @param min_level_n Integer. Minimum number of patients required for a
#' categorical level to keep its own coefficient in the multivariable model.
#' Rarer levels are merged into an \code{"Other"} level when that leaves the
#' covariate usable; if merging is not possible the patients belonging to those
#' levels are excluded from the multivariable model and the covariate is kept.
#' Default \code{5}.
#' @param min_covariates Integer. Number of covariates below which the
#' multivariable model is no longer regarded as an adjusted model. When the
#' automatic repair can only keep fewer than \code{min_covariates} covariates,
#' the model is still fitted and returned (with a warning) so that the
#' analysis can be inspected, but \code{metadata$reduced} is \code{TRUE} and
#' \code{metadata$reduced_note} states which covariate(s) survived and why the
#' others could not be co-estimated. Default \code{2}. Only a model with no
#' estimable covariate at all is refused.
#' @param drop_nonestimable Logical. If \code{TRUE} (default) a non-estimable
#' multivariable model is repaired instead of being refused: constant or nearly
#' constant covariates and perfectly separated or very rare levels are screened
#' out first, covariates that still keep the joint model from converging are
#' then removed one at a time (or, when the offender cannot be isolated, by
#' trying each removal in turn, down to a single covariate). Every decision is
#' reported in \code{metadata$dropped_covariates} and
#' \code{metadata$drop_notes}. If \code{FALSE} the model is fitted on the raw
#' complete cases and an error is raised when it is not estimable.
#' @return List of class \code{cpas_COX}:
#'   \item{\code{input_params}:}{echoed inputs}
#'   \item{\code{processed_data}:}{analysis data after numeric/factor coercion}
#'   \item{\code{models}:}{list of fitted \code{coxph} models}
#'   \item{\code{summaries}:}{corresponding \code{summary.coxph} objects}
#'   \item{\code{results_table}:}{tidy results with columns \code{Var1},
#'     \code{Variates}, \code{Level}, \code{Type} (\code{cont},
#'     \code{cate_header}, \code{cate_level}), \code{N}, \code{HR},
#'     \code{HR95L}, \code{HR95H}, \code{Pvalue}, \code{HR_text},
#'     \code{P_text}}
#'   \item{\code{metadata}:}{sample sizes, events, and, for multivariable
#'     models, \code{dropped_covariates} and \code{drop_notes} describing every
#'     covariate that the automatic repair removed and why, plus
#'     \code{final_covariates}, \code{reduced} and \code{reduced_note} which
#'     flag a model that had to be reduced below \code{min_covariates}}
#' @details
#' Status is treated as \code{1 = event, 0 = censored}; rows with missing
#' time/status are dropped. For categorical covariates the first factor level
#' is the reference (no row is printed for it); header rows carry
#' \code{HR = NA} and are intended for plot annotations. Univariate models are
#' estimated on the complete cases of \code{(time, status, variable)};
#' multivariate models on the complete cases of all covariates together.
#'
#' A multivariable model can fail not because the data are unusable but because
#' one covariate is: a constant or nearly constant continuous covariate, a
#' categorical level with very few patients, a level in which all (or no)
#' patients have an event, or a covariate collinear with the others. With
#' \code{drop_nonestimable = TRUE} such covariates are located and removed, and
#' the remaining model is fitted and reported together with the reason for each
#' removal; coefficients are never reported for a model that did not converge,
#' because those estimates are meaningless rather than merely imprecise. The
#' search continues down to a single covariate: if no pair of covariates can be
#' co-estimated, the largest estimable model is returned with
#' \code{metadata$reduced = TRUE}, a warning, and \code{metadata$reduced_note}
#' explaining that the result is not adjusted for confounding. Only when nothing
#' at all can be estimated does the function stop, and it then reports both the
#' covariates it had already removed and the reason the final fit failed.
#' @export
#' @examples
#' \dontrun{
#'    ## A real cohort from the CanPAS mirror: Lung Cancer, OS + DSS, n = 133.
#'    di <- dataset_info[dataset_info$Accession == "GSE14814", ]
#'    di[, c("Accession", "Type", "N", "EndpointFamilies")]
#' 
#'    ## Prerequisite chain: expression (probes collapsed) then survival, and
#'    ## with clin = TRUE the clinical covariates this cohort really carries.
#'    d <- cohort_merged("GSE13507", c("GAPDH", "ACTB"), type = "OS", clin = TRUE)
#'    colnames(d)
#' 
#'    ## One model per covariate.
#'    r <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
#'                       cate_Variates = c("grade", "N"), method = "uni")
#'    head(r$results_table)
#' 
#'    ## Joint model: covariates that cannot be co-estimated are located and
#'    ## removed, and every decision is recorded together with its reason.
#'    r2 <- COX_analysis(d, type = "OS", cont_Variates = c("GAPDH", "age"),
#'                        cate_Variates = c("grade", "N"), method = "multi")
#'    r2$metadata$dropped_covariates
#'    r2$metadata$reduced
#' }
COX_analysis <- function(df,
                         type = NULL,
                         cont_Variates = NULL,
                         cate_Variates = NULL,
                         method = c("uni", "multi"),
                         precision = 3,
                         min_level_n = 5L,
                         min_covariates = 2L,
                         drop_nonestimable = TRUE) {
  method <- match.arg(method)
  if (!is.data.frame(df)) stop("'df' must be a data.frame.")
  cont_Variates <- if (is.null(cont_Variates)) character(0) else as.character(cont_Variates)
  cate_Variates <- if (is.null(cate_Variates)) character(0) else as.character(cate_Variates)
  all_vars <- unique(c(cont_Variates, cate_Variates))
  if (length(all_vars) == 0L) stop("No covariates supplied.")
  miss <- setdiff(all_vars, colnames(df))
  if (length(miss)) stop("Covariate columns not found in 'df': ",
                         paste(miss, collapse = ", "), ".")
  if (length(intersect(cont_Variates, cate_Variates)))
    stop("A variable cannot be both continuous and categorical.")

  if (!is.null(type)) {
    type <- .resolve_type_in_df(df, type)      # family or token
    tc <- paste0(type, "_time"); sc <- paste0(type, "_status")
    if (!all(c("ID", tc, sc) %in% colnames(df)))
      stop("'df' must contain columns ID, ", tc, ", ", sc, ".")
    dat <- df[c("ID", tc, sc, all_vars)]
    dat[[tc]] <- suppressWarnings(as.numeric(dat[[tc]]))
    dat[[sc]] <- suppressWarnings(as.numeric(dat[[sc]]))
    colnames(dat)[2:3] <- c("time", "status")
  } else {
    if (!all(c("time", "status") %in% colnames(df)))
      stop("With type = NULL, 'df' must contain columns 'time' and 'status'.")
    dat <- df[c("time", "status", all_vars)]
  }
  dat$time <- suppressWarnings(as.numeric(dat$time))
  dat$status <- suppressWarnings(as.numeric(dat$status))
  if (any(!is.na(dat$status) & !dat$status %in% c(0, 1)))
    stop("'status' must be coded 0 (censored) / 1 (event).")
  if (any(!is.na(dat$time) & dat$time < 0))
    stop("'time' contains negative values.")
  if (sum(dat$status == 1, na.rm = TRUE) < 2L)
    stop("Fewer than 2 events available.")

  processed_data <- dat
  for (v in cont_Variates) processed_data[[v]] <- suppressWarnings(as.numeric(processed_data[[v]]))
  for (v in cate_Variates) processed_data[[v]] <- factor(processed_data[[v]])

  fmt <- function(x) sprintf(paste0("%.", precision, "f"), x)

  # ---- extract one tidy row per model term (exact coef/vcov, no summary rounding) ----
  plain <- function(x) gsub("^`(.*)`$", "\\1", x)   # non-syntactic names come back-quoted
  coef_rows <- function(fit, dd, cont_vars, cate_vars) {
    b <- stats::coef(fit); V <- stats::vcov(fit)
    terms <- attr(fit$terms, "term.labels")
    xl <- fit$xlevels                       # named list, only categorical vars
    rows <- list()
    for (tlab in terms) {
      tl <- plain(tlab)
      if (tl %in% cont_vars) {
        nm <- if (tlab %in% names(b)) tlab else if (tl %in% names(b)) tl else names(b)[1]
        if (!nm %in% names(b)) next
        est <- unname(b[nm]); se <- sqrt(unname(V[nm, nm]))
        if (!is.finite(est) || !is.finite(se) || se <= 0) next
        rows[[length(rows) + 1L]] <- data.frame(
          Var1 = tl, Variates = tl, Level = NA_character_, Type = "cont",
          HR = exp(est), HR95L = exp(est - 1.96 * se),
          HR95H = exp(est + 1.96 * se),
          Pvalue = 2 * stats::pnorm(-abs(est / se)), stringsAsFactors = FALSE)
      } else if (tl %in% cate_vars) {
        lv <- xl[[tlab]]; if (is.null(lv)) lv <- xl[[tl]]
        if (is.null(lv)) next
        for (lvl in lv[-1L]) {
          nm <- if (paste0(tlab, lvl) %in% names(b)) paste0(tlab, lvl) else paste0(tl, lvl)
          if (!nm %in% names(b)) next
          est <- unname(b[nm]); se <- sqrt(unname(V[nm, nm]))
          if (!is.finite(est) || !is.finite(se) || se <= 0) next
          rows[[length(rows) + 1L]] <- data.frame(
            Var1 = tl, Variates = tl, Level = lvl, Type = "cate_level",
            HR = exp(est), HR95L = exp(est - 1.96 * se),
            HR95H = exp(est + 1.96 * se),
            Pvalue = 2 * stats::pnorm(-abs(est / se)), stringsAsFactors = FALSE)
        }
      }
    }
    if (!length(rows)) stop("No estimable coefficients were returned by coxph.")
    do.call(rbind, rows)
  }

  header_row <- function(var, n) {
    data.frame(Var1 = var, Variates = var, Level = NA_character_,
               Type = "cate_header", HR = NA_real_, HR95L = NA_real_,
               HR95H = NA_real_, Pvalue = NA_real_, N = as.integer(n),
               stringsAsFactors = FALSE)
  }

  models <- list(); summaries <- list(); parts <- list()
  failed <- character(0)
  # covariates handled automatically in the multivariable branch (empty for the
  # univariable branch, but always present so the metadata can report it)
  drops <- data.frame(variable = character(0), detail = character(0),
                      reason = character(0), stringsAsFactors = FALSE)
  # TRUE when the multivariable repair had to keep fewer than 'min_covariates'
  # covariates, so the fitted model is not an adjusted model
  reduced <- FALSE

  if (method == "uni") {
    for (v in cont_Variates) {
      dd <- processed_data[stats::complete.cases(processed_data[c("time", "status", v)]), , drop = FALSE]
      fml <- stats::as.formula(paste0("Surv(time, status) ~ `", v, "`"))
      if (nrow(dd) == 0L) { failed <- c(failed, paste0(v, ": no complete cases")); next }
      diag <- .cpas_COX_diag(fml, dd)
      if (!diag$ok) { failed <- c(failed, paste0(v, ": ", diag$reason)); next }
      fit <- diag$fit
      models[[v]] <- fit; summaries[[v]] <- summary(fit)
      r <- tryCatch(coef_rows(fit, dd, v, character(0)), error = function(e) e)
      if (inherits(r, "error")) { failed <- c(failed, paste0(v, ": ", conditionMessage(r))); next }
      r$N <- nrow(dd)
      parts[[v]] <- r
    }
    for (v in cate_Variates) {
      dd <- processed_data[stats::complete.cases(processed_data[c("time", "status", v)]), , drop = FALSE]
      fml <- stats::as.formula(paste0("Surv(time, status) ~ `", v, "`"))
      if (nrow(dd) == 0L) { failed <- c(failed, paste0(v, ": no complete cases")); next }
      diag <- .cpas_COX_diag(fml, dd)
      if (!diag$ok) { failed <- c(failed, paste0(v, ": ", diag$reason)); next }
      fit <- diag$fit
      models[[v]] <- fit; summaries[[v]] <- summary(fit)
      hdr <- header_row(v, nrow(dd))
      r <- tryCatch(coef_rows(fit, dd, character(0), v), error = function(e) e)
      if (inherits(r, "error")) { failed <- c(failed, paste0(v, ": ", conditionMessage(r))); next }
      if (is.null(r)) r <- data.frame()
      r$N <- nrow(dd)
      parts[[v]] <- rbind(hdr, r)
    }
  } else {
    dd <- processed_data[stats::complete.cases(processed_data[c("time", "status", all_vars)]), ,
                         drop = FALSE]
    if (nrow(dd) == 0L)
      stop("No complete cases for the requested covariates: ",
           paste(all_vars, collapse = ", "), ".")
    if (sum(dd$status == 1) < 2L) stop("Fewer than 2 events in complete cases.")

    # ---- screen the covariates, then drop what still breaks the joint model ----
    # Rather than refusing the whole model, locate the offending covariates
    # (constant / nearly constant, rare or perfectly separated levels, then
    # whatever keeps the joint fit non-estimable), drop them and report why.
    screen <- .cpas_COX_screen(dd, status = "status",
                               cont = intersect(cont_Variates, all_vars),
                               cate = intersect(cate_Variates, all_vars),
                               min_level_n = min_level_n)
    drops <- screen$drops
    vars <- screen$keep
    if (!length(vars))
      stop("Multivariate Cox model is not estimable: ",
           paste(unique(drops$reason), collapse = "; "), ".", call. = FALSE)
    if (drop_nonestimable) {
      repeat {
        # take the screened data: categorical levels removed above must stay
        # removed, otherwise the rare/separated levels return and the model
        # fails for a reason we already handled
        dat <- screen$data[stats::complete.cases(screen$data[c("time", "status", vars)]), ,
                           drop = FALSE]
        if (!nrow(dat) || sum(dat$status == 1) < 2L) {
          drops <- rbind(drops, data.frame(
            variable = paste(vars, collapse = "+"), detail = "complete cases",
            reason = "fewer than 2 events after dropping the other covariates",
            stringsAsFactors = FALSE))
          vars <- character(0); break
        }
        fml <- stats::as.formula(paste0("Surv(time, status) ~ ",
                                        paste0("`", vars, "`", collapse = " + ")))
        diag <- .cpas_COX_diag(fml, dat[, c("time", "status", vars), drop = FALSE])
        if (diag$ok) break
        off <- .cpas_COX_drop_worst(diag, vars)
        if (length(off) != 1L || is.na(off) || !(off %in% vars)) {
          # the offending variable cannot be isolated (typically a collinear
          # pair): try each variable in turn and keep the first subset that is
          # estimable, dropping the one that costs the least information
          trial <- NULL
          for (v in rev(vars)) {
            sub_vars <- setdiff(vars, v)
            # subsets may shrink to a single covariate: a reduced model that is
            # reported and flagged beats refusing the analysis altogether
            if (!length(sub_vars)) next
            f <- stats::as.formula(paste0("Surv(time, status) ~ ",
                                          paste0("`", sub_vars, "`", collapse = " + ")))
            d2 <- .cpas_COX_diag(f, dat[, c("time", "status", sub_vars), drop = FALSE])
            if (d2$ok) { trial <- list(vars = sub_vars, diag = d2, dropped = v); break }
          }
          if (is.null(trial)) break
          vars <- trial$vars
          diag <- trial$diag
          drops <- rbind(drops, data.frame(
            variable = trial$dropped, detail = "collinear / joint model",
            reason = "the joint model was not estimable and the offending covariate could not be isolated; this covariate was removed by trying each removal in turn",
            stringsAsFactors = FALSE))
          break
        }
        vars <- setdiff(vars, off)
        drops <- rbind(drops, data.frame(
          variable = off, detail = "joint model",
          reason = paste0("non-estimable in the multivariable model: ", diag$reason),
          stringsAsFactors = FALSE))
        if (!length(vars)) break
      }
    } else {
      dat <- dd
    }
    if (is.null(diag) || !diag$ok || !length(vars)) {
      # report the whole sequence, not just the last fit that failed: the
      # covariates already removed and their reasons are what the user needs
      why <- unique(c(if (!is.null(diag) && !diag$ok) diag$reason else NULL,
                      if (!is.null(drops) && nrow(drops))
                        paste0(drops$variable, " (", drops$detail, ": ", drops$reason, ")") else NULL))
      stop("Multivariate Cox model is not estimable: ",
           if (length(why)) paste(why, collapse = "; ") else "no covariate could be estimated",
           ". Complete separation, collinear or constant covariates produce ",
           "meaningless coefficients; remove or combine the offending variable.",
           call. = FALSE)
    }
    # the repair may have had to remove so much that the "multivariable" model
    # is no longer multivariable; it is still returned (the caller asked for a
    # result), but it is flagged so it cannot be reported as an adjusted model
    reduced <- length(vars) < min_covariates
    if (reduced)
      warning(sprintf(paste0("Multivariable model reduced to %d covariate(s) (%s): ",
                             "%s. This is not an adjusted model; the covariate(s) ",
                             "listed above could not be co-estimated."),
                      length(vars), paste(vars, collapse = ", "),
                      paste(unique(drops$reason), collapse = "; ")), call. = FALSE)
    dd <- dat
    cont_Variates <- intersect(cont_Variates, vars)
    cate_Variates <- intersect(cate_Variates, vars)
    fit <- diag$fit
    models$multivariate <- fit; summaries$multivariate <- summary(fit)
    n_par <- length(stats::coef(fit)); events_multi <- sum(dd$status == 1)
    if (n_par > 0 && events_multi / n_par < 10)
      warning(sprintf(paste0("Multivariable model has %.1f events per variable ",
                             "(%d events, %d coefficients); estimates may be unstable ",
                             "(a common rule of thumb is >= 10 events per variable)."),
                      events_multi / n_par, events_multi, n_par), call. = FALSE)
    ph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
    if (!is.null(ph)) summaries$ph_test <- ph
    tab <- coef_rows(fit, dd, cont_Variates, cate_Variates)
    if (is.null(tab)) tab <- data.frame()
    # categorical header rows, one per categorical variable, placed above its levels
    ord <- unique(c(intersect(cont_Variates, tab$Var1),
                    intersect(cate_Variates, tab$Var1)))
    out <- list()
    for (v in ord) {
      if (v %in% cate_Variates) out[[length(out) + 1L]] <- header_row(v, nrow(dd))
      lv_rows <- tab[tab$Var1 == v & tab$Type == "cate_level", , drop = FALSE]
      if (v %in% cont_Variates) lv_rows <- tab[tab$Var1 == v, , drop = FALSE]
      if (nrow(lv_rows)) { lv_rows$N <- nrow(dd); out[[length(out) + 1L]] <- lv_rows }
    }
    if (!length(out)) out[[1L]] <- tab
    parts[["multivariate"]] <- do.call(rbind, out)
  }

  new_failed <- setdiff(failed, names(models))
  if (length(new_failed))
    warning("Variable(s) skipped (not estimable): ", paste(new_failed, collapse = "; "),
            call. = FALSE)
  if (!length(parts))
    stop("No estimable covariate remained.",
         if (length(failed)) paste0(" Skipped: ", paste(failed, collapse = "; ")) else "",
         call. = FALSE)
  results_table <- do.call(rbind, parts)
  if (is.null(results_table) || !nrow(results_table))
    stop("No results could be produced.")
  results_table$N <- as.integer(results_table$N)
  results_table$HR_text <- ifelse(is.na(results_table$HR), "",
    paste0(fmt(results_table$HR), " [", fmt(results_table$HR95L), ", ",
           fmt(results_table$HR95H), "]"))
  p_thr <- 10^(-precision)
  results_table$P_text <- ifelse(is.na(results_table$Pvalue), "",
    ifelse(results_table$Pvalue < p_thr,
           paste0("<", format(p_thr, scientific = FALSE)),
           sprintf(paste0("%.", precision, "f"), results_table$Pvalue)))
  rownames(results_table) <- NULL

  result_obj <- list(
    input_params = list(df_name = deparse(substitute(df)), type = type,
                        cont_Variates = cont_Variates,
                        cate_Variates = cate_Variates,
                        method = method, precision = precision,
                        analysis_time = Sys.time()),
    processed_data = processed_data,
    models = models,
    summaries = summaries,
    results_table = results_table,
    metadata = list(sample_size = nrow(processed_data),
                    events = sum(processed_data$status == 1, na.rm = TRUE),
                    complete_cases = if (method == "multi")
                      nrow(processed_data[stats::complete.cases(processed_data[c("time", "status",
                        intersect(all_vars, unique(c(cont_Variates, cate_Variates))))]), ,
                                          drop = FALSE]) else NA_integer_,
                    events_in_model = if (method == "multi" && !is.null(models$multivariate))
                      tryCatch(sum(models$multivariate$y[, 2] == 1), error = function(e) NA_integer_) else NA_integer_,
                    skipped_variables = failed,
                    dropped_covariates = drops,
                    drop_notes = if (nrow(drops))
                      paste0(unique(drops$variable), ": ", unique(drops$reason)) else character(0),
                    final_covariates = if (method == "multi") vars else all_vars,
                    reduced = reduced,
                    reduced_note = if (reduced)
                      paste0("the multivariable model could only keep ",
                             length(vars), " covariate(s) (",
                             paste(vars, collapse = ", "), "), because the others ",
                             "could not be co-estimated: ",
                             paste(unique(drops$reason), collapse = "; "),
                             ". It is therefore not adjusted for confounding.")
                      else character(0),
                    ph_test = if (!is.null(summaries$ph_test)) summaries$ph_test else NULL,
                    cont_vars = cont_Variates, cate_vars = cate_Variates,
                    analysis_type = ifelse(method == "uni", "Univariate COX",
                                           "Multivariate COX"),
                    survival_type = type))
  class(result_obj) <- "cpas_COX"
  result_obj
}
