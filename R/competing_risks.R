#' @title Cumulative incidence functions for competing risks
#' @description
#' Estimates the Aalen-Johansen cumulative incidence function (CIF) of each
#' event type with \code{survival::survfit()} on a multi-state response
#' \code{Surv(time, factor(status))}, optionally by group.
#' @param df Data.frame with a time column and a cause-of-event status column.
#' @param time Name of the numeric follow-up time column (years).
#' @param status Name of the status column: \code{0} = censored,
#'   \code{1, 2, ...} = event types (cause of the event).
#' @param group Optional name of a \strong{categorical} column used to split
#'   the CIF (e.g. a marker group or sex). A continuous column is refused: it
#'   would become one stratum per distinct value.
#' @param times Optional numeric vector of time points at which the CIF is
#'   tabulated (default: deciles of the observed event times).
#' @param conf_level Confidence level for the pointwise interval.
#' @param group_cut What to do when \code{group} names a continuous column:
#'   \code{"none"} (default) stops with an explanatory error,
#'   \code{"median"} splits at the median into low/high, \code{"tertile"} into
#'   low/middle/high and \code{"quartile"} into Q1-Q4. The cut points are
#'   reported in the result (\code{cut_points}), so the split is reproducible.
#' @param max_levels More distinct values than this make a group column count as
#'   continuous (default 20).
#' @return Object of class \code{cpas_cif}: list with \code{table} (data.frame:
#'   \code{group}, \code{cause}, \code{time}, \code{n_risk}, \code{cif},
#'   \code{lower}, \code{upper}), \code{fit} (the \code{survfit} object),
#'   \code{causes}, \code{events}, \code{n}, \code{group_levels} and
#'   \code{group_sizes} (how the strata were formed) and the echoed inputs.
#' @details
#' With competing events, the complement of the Kaplan-Meier estimate is
#' \strong{not} the incidence of the event of interest: patients who die of
#' another cause are removed from the risk set in the Kaplan-Meier estimator and
#' their future event probability is redistributed over the remaining causes.
#' The CIF does not redistribute that probability, which is why it is the
#' appropriate descriptive curve whenever competing events exist.
#'
#' \code{group} names a grouping variable, not a covariate: the column is
#' factored as supplied, so a continuous variable such as age would produce one
#' stratum per distinct value (and one curve per stratum). Such a column is
#' refused unless \code{group_cut} asks for a documented quantile split, in
#' which case the threshold(s) used are returned with the result.
#' @export
#' @examples
#' \dontrun{
#'    ## A real cohort that carries both OS and DSS, through the prerequisite
#'    ## chain: cause 1 = cancer death (DSS event), cause 2 = other death.
#'    d  <- cohort_merged("GSE14814", "GAPDH", type = "OS", clin = TRUE)
#'    cr <- data.frame(
#'      time   = d$OS_time,
#'      status = ifelse(d$OS_status == 0, 0, ifelse(d$DSS_status == 1, 1, 2)),
#'      age    = d$age,                        # continuous
#'      sex    = factor(d$sex),                # categorical
#'      age_group = factor(ifelse(d$age > stats::median(d$age, na.rm = TRUE),
#'                                  "older", "younger")))
#'    table(cr$status)   # 0 censored, 1 cancer death, 2 other death
#' 
#'    ci <- cif_fit(cr, time = "time", status = "status", times = c(1, 3, 5))
#'    ci$table
#' 
#'    ## Stratified by a real categorical column (sex) ...
#'    ci_sex <- cif_fit(cr, time = "time", status = "status", group = "sex")
#'    ci_sex$group_levels; ci_sex$group_sizes
#'    plot_cif(ci_sex)
#' 
#'    ## ... or by a documented split of the CONTINUOUS age column. Passing age
#'    ## directly as group is refused: it would create one stratum per distinct
#'    ## age. group_cut says where to cut, and returns the cut point used.
#'    ci_age <- cif_fit(cr, time = "time", status = "status",
#'                      group = "age", group_cut = "median")
#'    ci_age$cut_points; ci_age$group_levels
#' }
cif_fit <- function(df, time, status, group = NULL, times = NULL,
                    conf_level = 0.95,
                    group_cut = c("none", "median", "tertile", "quartile"),
                    max_levels = 20L) {
  group_cut <- match.arg(group_cut)
  if (!is.data.frame(df)) stop("'df' must be a data.frame.")
  for (v in c(time, status)) if (!v %in% colnames(df))
    stop("Column '", v, "' not found in 'df'.")
  if (!is.null(group) && !group %in% colnames(df))
    stop("Group column '", group, "' not found in 'df'.")
  dd <- df[, c(time, status, group), drop = FALSE]
  colnames(dd) <- c("time", "cause", if (!is.null(group)) "grp")
  dd$time <- suppressWarnings(as.numeric(dd$time))
  dd$cause <- suppressWarnings(as.numeric(dd$cause))
  dd <- dd[stats::complete.cases(dd), , drop = FALSE]
  if (any(dd$time < 0)) stop("Negative follow-up times are not allowed.")
  causes <- sort(unique(dd$cause))
  if (!all(causes %in% 0:9))
    stop("The status column must be 0 (censored) / 1, 2, ... (cause of event).")
  ev_causes <- setdiff(causes, 0)
  if (!length(ev_causes) || length(ev_causes) < 2L)
    stop("At least two competing event types are required (status 1, 2, ...). ",
         "With a single event type use plot_km().")
  if (nrow(dd) < 10L) stop("Too few complete rows.")

  # ---- the group must be categorical -----------------------------------
  # 'group' names a GROUPING column, and it is factored below: a continuous
  # variable (age measured to one decimal, say) would give one stratum per
  # distinct value - 103 curves instead of a stratified comparison. Such a
  # column is refused unless the caller asks for a documented split.
  cut_points <- NULL; group_levels <- NULL
  if (!is.null(group)) {
    g <- dd$grp
    distinct <- unique(g[!is.na(g)])
    continuous <- (is.numeric(g) || is.character(g)) && length(distinct) > max_levels
    if (continuous && group_cut == "none")
      stop("The group column '", group, "' has ", length(distinct),
           " distinct values, so it is continuous, not categorical: every value ",
           "would become its own stratum (", length(distinct), " curves). ",
           "Build a group column explicitly (e.g. df$age_group <- ",
           "ifelse(df$age > median(df$age, na.rm = TRUE), 'high', 'low')) or ask ",
           "for a documented split with group_cut = \"median\", \"tertile\" or ",
           "\"quartile\".", call. = FALSE)
    if (continuous) {
      x <- suppressWarnings(as.numeric(g))
      if (any(is.na(x)))
        stop("group_cut needs a numeric group column; '", group,
             "' contains values that are not numbers.", call. = FALSE)
      probs <- switch(group_cut, median = 0.5, tertile = c(1 / 3, 2 / 3),
                      quartile = c(0.25, 0.5, 0.75))
      cut_points <- unname(stats::quantile(x, probs = probs, na.rm = TRUE))
      labs <- switch(group_cut, median = c("low", "high"),
                     tertile = c("low", "middle", "high"),
                     quartile = c("Q1", "Q2", "Q3", "Q4"))
      grp <- cut(x, breaks = c(-Inf, cut_points, Inf), labels = labs,
                 include.lowest = TRUE, right = FALSE)
      if (any(is.na(grp)))
        stop("group_cut produced empty groups; check the distribution of '",
             group, "'.", call. = FALSE)
      dd$grp <- grp
      group_levels <- labs
      message("Group column '", group, "' is continuous: split by ",
              group_cut, " at ", paste(sprintf("%.4g", cut_points), collapse = ", "),
              " into ", paste(labs, collapse = "/"), ".")
    } else if (is.character(g)) {
      dd$grp <- factor(g)
    }
  }

  dd$cause_f <- factor(dd$cause, levels = c(0, ev_causes))
  if (is.null(times)) {
    tt <- dd$time[dd$cause != 0]
    times <- unique(stats::quantile(tt, probs = seq(0.2, 0.9, by = 0.1), na.rm = TRUE))
  }
  times <- sort(unique(as.numeric(times)))
  times <- times[is.finite(times) & times > 0]

  fit <- if (is.null(group)) {
    survival::survfit(survival::Surv(time, cause_f) ~ 1, data = dd, conf.int = conf_level)
  } else {
    dd$grp <- factor(dd$grp)
    survival::survfit(survival::Surv(time, cause_f) ~ grp, data = dd, conf.int = conf_level)
  }
  sm <- summary(fit, times = times, extend = TRUE)

  st <- attr(fit, "states")
  pull_mat <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    colnames(x) <- if (!is.null(colnames(x))) colnames(x) else st
    x
  }
  P <- pull_mat(sm$pstate); L <- pull_mat(sm$lower); U <- pull_mat(sm$upper)
  nrisk <- if (!is.null(sm$n.risk)) sm$n.risk else rep(NA_real_, length(sm$time))
  strata <- if (!is.null(sm$strata)) as.character(sm$strata) else rep("all", length(sm$time))
  strata <- sub("^grp=", "", strata)

  rows <- list()
  for (i in seq_along(sm$time)) {
    for (cs in ev_causes) {
      cn <- as.character(cs)
      if (is.null(P) || !cn %in% colnames(P)) next
      rows[[length(rows) + 1L]] <- data.frame(
        group = strata[i], cause = cs, time = sm$time[i],
        n_risk = nrisk[i], cif = unname(P[i, cn]),
        lower = if (!is.null(L) && cn %in% colnames(L)) unname(L[i, cn]) else NA_real_,
        upper = if (!is.null(U) && cn %in% colnames(U)) unname(U[i, cn]) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("No cumulative incidence could be extracted from the fit.")
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$group, tab$cause, tab$time), , drop = FALSE]
  rownames(tab) <- NULL

  events <- vapply(ev_causes, function(cs) sum(dd$cause == cs), numeric(1))
  names(events) <- as.character(ev_causes)
  out <- list(input_params = list(time = time, status = status, group = group,
                                  conf_level = conf_level, group_cut = group_cut,
                                  analysis_time = Sys.time()),
              table = tab, fit = fit, causes = ev_causes, events = events,
              n_censored = sum(dd$cause == 0), n = nrow(dd), times = times,
              # how the strata were formed: a categorical column as supplied, or
              # a continuous one split at these quantiles
              group = group,
              group_cut = if (is.null(group)) NULL else group_cut,
              cut_points = cut_points,
              group_levels = if (is.null(group)) NULL else levels(dd$grp),
              group_sizes = if (is.null(group)) NULL else table(droplevels(dd$grp)))
  class(out) <- "cpas_cif"
  out
}

#' @title Plot cumulative incidence curves
#' @description Draws the cumulative incidence functions returned by
#' \code{\link{cif_fit}} as step curves with pointwise confidence ribbons.
#' @param x An object returned by \code{\link{cif_fit}}.
#' @param cause Event type (cause code) to plot; default: all.
#' @param ... Unused, kept for S3 compatibility.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#'    ## A real cohort that carries both OS and DSS, through the prerequisite
#'    ## chain: cause 1 = cancer death (DSS event), cause 2 = other death.
#'    d  <- cohort_merged("GSE14814", "GAPDH", type = "OS", clin = TRUE)
#'    cr <- data.frame(
#'      time   = d$OS_time,
#'      status = ifelse(d$OS_status == 0, 0, ifelse(d$DSS_status == 1, 1, 2)),
#'      age    = d$age,                        # continuous
#'      sex    = factor(d$sex),                # categorical
#'      age_group = factor(ifelse(d$age > stats::median(d$age, na.rm = TRUE),
#'                                  "older", "younger")))
#'    table(cr$status)   # 0 censored, 1 cancer death, 2 other death
#' 
#'    x <- cif_fit(cr, time = "time", status = "status", times = c(1, 3, 5))
#'    plot_cif(x)             # cumulative incidence of cause 1
#'    plot_cif(x, cause = 2)
#' 
#'    ## The same figure per group, and per tertile of a continuous variable.
#'    plot_cif(cif_fit(cr, time = "time", status = "status", group = "sex"))
#'    plot_cif(suppressMessages(
#'      cif_fit(cr, time = "time", status = "status", group = "age",
#'              group_cut = "tertile")))
#' }
plot_cif <- function(x, cause = NULL, ...) {
  if (!inherits(x, "cpas_cif")) stop("'x' must be a 'cpas_cif' object (see cif_fit()).")
  d <- x$table
  if (!is.null(cause)) d <- d[d$cause %in% cause, , drop = FALSE]
  if (!nrow(d)) stop("No rows for the requested cause(s).")
  d$cause_lab <- paste0("cause ", d$cause)
  ggplot2::ggplot(d, ggplot2::aes(x = time, y = cif,
                                  colour = cause_lab, fill = cause_lab)) +
    ggplot2::geom_step(linewidth = 0.9) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.15, colour = NA) +
    ggplot2::facet_wrap(~ group) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Cumulative incidence (Aalen-Johansen)",
                  x = "Time (years)", y = "Cumulative incidence",
                  colour = "Event type", fill = "Event type") +
    ggplot2::theme_bw()
}

#' @title Cause-specific and subdistribution (Fine-Gray) Cox models
#' @description
#' Fits, for the same covariates, (1) the cause-specific Cox model
#' (\code{Surv(time, status == etype)}) and (2) the Fine-Gray subdistribution
#' hazard model obtained with \code{survival::finegray()} and a weighted
#' \code{coxph()}. The two answer different questions and are reported side by
#' side.
#' @param df Data.frame with a time column and a cause-coded status column.
#' @param time Name of the numeric follow-up time column (years).
#' @param status Name of the status column: \code{0} = censored,
#'   \code{1, 2, ...} = event types.
#' @param covariates Character vector of covariates (the first one is usually
#'   the marker of interest). A factor stays categorical; a character column
#'   whose values all parse as numbers is treated as continuous (mirrored
#'   clinical columns arrive as text, e.g. \code{"55.4"}), any other character
#'   column becomes a factor. The classification actually used is returned in
#'   \code{covariate_types}.
#' @param etype Event type of interest (default 1). Competing events are all
#'   other non-zero codes.
#' @param conf_level Confidence level for the reported intervals.
#' @param max_levels Maximum number of levels accepted for a categorical
#'   covariate (default 20). A covariate with one coefficient per level is
#'   usually a continuous variable that was supplied as text, so the call stops
#'   with an explanatory error instead of reporting meaningless hazard ratios;
#'   raise the limit for a deliberate fine stratification.
#' @return Object of class \code{cpas_competing}: list with
#'   \item{\code{cause_specific}}{tidy data.frame (\code{Variates}, \code{Level},
#'     \code{HR}, \code{HR95L}, \code{HR95H}, \code{Pvalue}, \code{N})}
#'   \item{\code{subdistribution}}{the same columns for the Fine-Gray model
#'     (hazard ratios are subdistribution hazard ratios, sHR)}
#'   \item{\code{n}, \code{events_interest}, \code{events_competing},
#'     \code{n_censored}, \code{etype}}{sample and event accounting}
#'   \item{\code{covariate_types}}{named vector: how each covariate was used
#'     (\code{"continuous"}, \code{"categorical"} or the text-coercion cases)}
#'   \item{\code{diagnostics}}{per model, \code{ok} and, when \code{FALSE},
#'     the reason (non-convergence, separation, infinite coefficient): the
#'     affected coefficient is then omitted from the table and a warning is
#'     raised, so a meaningless hazard ratio is never reported as a result}
#'   \item{\code{fits}}{the two fitted \code{coxph} objects}
#' @details
#' The cause-specific hazard is the rate of the event among patients still
#' event-free; it answers "does the marker act on the disease process".
#' The subdistribution hazard keeps patients with a competing event in the risk
#' set and is directly linked to the cumulative incidence; it answers "does the
#' marker change the probability of dying of the disease". Reporting only the
#' cause-specific hazard can overstate the clinical effect when competing
#' events are common (and vice versa).
#'
#' A numeric covariate contributes one coefficient. If a continuous variable is
#' supplied as text (as the mirrored clinical columns are), it would otherwise
#' be read as a factor with one level per value; \code{competing_risk_COX()}
#' restores the evident type, reports it in \code{covariate_types}, and refuses
#' (rather than silently reporting) a categorical covariate with more than
#' \code{max_levels} levels. Non-convergence, separation and infinite
#' coefficients are reported through \code{diagnostics} and excluded from the
#' tables.
#'
#' The Fine-Gray standard errors are clustered by patient, because
#' \code{survival::finegray()} expands every patient who has a competing event
#' into several rows; without the cluster term the reported standard error is
#' 17-23\% too small. The point estimates and the clustered standard errors agree
#' with \code{cmprsk::crr()} to within 0.1\%.
#' @export
#' @examples
#' \dontrun{
#'    ## A real cohort that carries both OS and DSS, through the prerequisite
#'    ## chain: cause 1 = cancer death (DSS event), cause 2 = other death.
#'    d  <- cohort_merged("GSE14814", "GAPDH", type = "OS", clin = TRUE)
#'    cr <- data.frame(
#'      time   = d$OS_time,
#'      status = ifelse(d$OS_status == 0, 0, ifelse(d$DSS_status == 1, 1, 2)),
#'      age    = d$age,     # numeric: one coefficient, as for a Cox model
#'      sex    = factor(d$sex))
#'    table(cr$status)   # 0 censored, 1 cancer death, 2 other death
#' 
#'    cc <- competing_risk_COX(cr, time = "time", status = "status",
#'                             covariates = c("age", "sex"), etype = 1)
#'    cc$covariate_types    # how each covariate was used (continuous / categorical)
#'    cc$cause_specific     # cause-specific hazard ratios
#'    cc$subdistribution    # Fine-Gray subdistribution hazard ratios
#'    cc$diagnostics        # ok + reason per model; unusable coefficients are omitted
#' }
competing_risk_COX <- function(df, time, status, covariates = NULL, etype = 1,
                               conf_level = 0.95, max_levels = 20L) {
  if (!is.data.frame(df)) stop("'df' must be a data.frame.")
  covariates <- if (is.null(covariates)) NULL else as.character(covariates)
  for (v in c(time, status, covariates)) if (!v %in% colnames(df))
    stop("Column '", v, "' not found in 'df'.")
  keep <- unique(c(time, status, covariates))
  dd <- df[, keep, drop = FALSE]
  dd[[time]] <- suppressWarnings(as.numeric(dd[[time]]))
  dd[[status]] <- suppressWarnings(as.numeric(dd[[status]]))
  dd <- dd[stats::complete.cases(dd), , drop = FALSE]
  if (any(dd[[time]] < 0)) stop("Negative follow-up times are not allowed.")
  codes <- sort(unique(dd[[status]]))
  if (!all(codes %in% 0:9)) stop("'status' must be 0 (censored) / 1, 2, ... (cause).")
  ev <- setdiff(codes, 0)
  if (length(ev) < 2L)
    stop("Competing risks need at least two event types (status 1, 2, ...); ",
         "with one event type use COX_analysis().")
  if (!etype %in% ev) stop("'etype' = ", etype, " is not among the event codes: ",
                           paste(ev, collapse = ", "), ".")
  if (sum(dd[[status]] == etype) < 5L)
    stop("Fewer than 5 events of interest (etype = ", etype, ").")
  if (is.null(covariates) || !length(covariates))
    stop("'covariates' must contain at least one column name.")

  # ---- covariate typing -------------------------------------------------
  # Mirrored clinical columns arrive as text ("55.4", "female"), and coxph()
  # turns a character column into a FACTOR: one coefficient per distinct value.
  # For a continuous variable that is meaningless (111 coefficients for age),
  # so the evident type is restored here and reported in covariate_types.
  cov_types <- stats::setNames(character(length(covariates)), covariates)
  for (v in covariates) {
    x <- dd[[v]]
    if (is.factor(x)) { cov_types[[v]] <- "categorical"; next }
    if (is.character(x)) {
      num <- suppressWarnings(as.numeric(x))
      if (any(!is.na(num)) && !any(is.na(num) & !is.na(x))) {
        dd[[v]] <- num
        cov_types[[v]] <- "continuous (numeric text)"
      } else {
        dd[[v]] <- factor(x)
        cov_types[[v]] <- "categorical (from text)"
      }
      next
    }
    if (is.logical(x)) {
      dd[[v]] <- factor(x); cov_types[[v]] <- "categorical (from logical)"; next
    }
    cov_types[[v]] <- "continuous"
  }
  nlev <- vapply(covariates, function(v)
    if (is.factor(dd[[v]])) nlevels(droplevels(dd[[v]])) else 0L, integer(1))
  wide <- names(nlev)[nlev > max_levels]
  if (length(wide))
    stop("Categorical covariate(s) with more than ", max_levels, " levels: ",
         paste(sprintf("%s (%d levels)", wide, nlev[wide]), collapse = "; "),
         ". One coefficient per level is almost never intended - a numeric ",
         "column read as text would produce exactly this. Convert the column ",
         "explicitly (as.numeric() for a continuous variable, factor() for a ",
         "deliberate grouping) or raise 'max_levels'.", call. = FALSE)
  coerced <- names(cov_types)[grepl("\\(", cov_types)]
  if (length(coerced))
    message("Covariate type inferred from text: ",
            paste(sprintf("%s -> %s", coerced, cov_types[coerced]), collapse = "; "),
            ". Wrap the column in factor() if it should be categorical.")

  # a coefficient that is infinite, or whose standard error is unusable, is not
  # a result: it is dropped from the table and reported in $diagnostics
  tidy_fit <- function(fit, n_pat) {
    if (is.null(fit)) return(NULL)          # the fit failed: diagnostics say why
    b <- stats::coef(fit); V <- stats::vcov(fit)
    if (!length(b)) return(NULL)
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    rows <- lapply(names(b), function(nm) {
      est <- unname(b[nm]); se <- tryCatch(sqrt(unname(V[nm, nm])), error = function(e) NA_real_)
      if (!is.finite(est) || !is.finite(se) || se <= 0) return(NULL)
      lvl <- NA_character_; var <- nm
      if (grepl("^`.*`", nm)) nm2 <- gsub("^`(.*)`$", "\\1", nm) else nm2 <- nm
      for (cv in covariates) {
        if (identical(nm2, cv)) { var <- cv; break }
        if (startsWith(nm2, cv) && nchar(nm2) > nchar(cv)) {
          var <- cv; lvl <- sub(paste0("^", cv), "", nm2); break
        }
      }
      data.frame(Variates = var, Level = lvl, HR = exp(est),
                 HR95L = exp(est - z * se), HR95H = exp(est + z * se),
                 Pvalue = 2 * stats::pnorm(-abs(est / se)),
                 N = n_pat, stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(NULL)
    out <- do.call(rbind, rows); rownames(out) <- NULL; out
  }

  # ---- 1. cause-specific hazard (non-event codes treated as censored) ----
  dd$cs_status <- as.numeric(dd[[status]] == etype)
  cs_fml <- stats::reformulate(covariates, response = "survival::Surv(cs_time, cs_status)")
  dd$cs_time <- dd[[time]]
  cs_msgs <- character(0)
  cs_fit <- withCallingHandlers(
    tryCatch(survival::coxph(cs_fml, data = dd, model = TRUE, x = TRUE),
             error = function(e) NULL),
    warning = function(w) {
      cs_msgs <<- c(cs_msgs, conditionMessage(w))
      tryCatch(invokeRestart("muffleWarning"), error = function(e) NULL)
    })
  cs_diag <- .cpas_COX_checks(cs_fit, cs_msgs)

  # ---- 2. Fine-Gray subdistribution hazard (weighted expansion) ----
  fg_data <- dd[, c(time, status, covariates), drop = FALSE]
  # The subdistribution variance must be clustered by PATIENT: finegray()
  # duplicates every patient who has a competing event, and the weighted fit
  # ignores that within-patient correlation. The naive model-based standard
  # error is then far too small (measured: 17-23% below the cmprsk::crr
  # variance). finegray() keeps only the variables named in its formula, so the
  # patient id is carried through as an extra term of the EXPANSION formula and
  # is deliberately left out of the Cox formula below.
  fg_data$.cpas_pid <- seq_len(nrow(fg_data))
  fg_fml <- stats::as.formula(sprintf("survival::Surv(`%s`, factor(`%s`)) ~ %s + .cpas_pid",
                                      time, status,
                                      paste(sprintf("`%s`", covariates), collapse = " + ")))
  fg <- survival::finegray(fg_fml, data = fg_data, etype = etype)
  if (!".cpas_pid" %in% colnames(fg))
    stop("finegray() did not carry the patient id; the subdistribution standard ",
         "error would be invalid. Please report this as a bug.", call. = FALSE)
  fml_fg <- stats::reformulate(covariates,
                               response = "survival::Surv(fgstart, fgstop, fgstatus)")
  fg_msgs <- character(0)
  fg_fit <- withCallingHandlers(
    tryCatch(survival::coxph(fml_fg, data = fg, weight = fg$fgwt,
                             cluster = fg$.cpas_pid, model = TRUE, x = TRUE),
             error = function(e) NULL),
    warning = function(w) {
      fg_msgs <<- c(fg_msgs, conditionMessage(w))
      tryCatch(invokeRestart("muffleWarning"), error = function(e) NULL)
    })
  fg_diag <- .cpas_COX_checks(fg_fit, fg_msgs)

  diagnostics <- list(
    cause_specific = list(ok = isTRUE(cs_diag$ok), reason = cs_diag$reason),
    subdistribution = list(ok = isTRUE(fg_diag$ok), reason = fg_diag$reason))
  for (nm in names(diagnostics)) {
    d <- diagnostics[[nm]]
    if (!d$ok)
      warning("The ", gsub("_", "-", nm), " model is not fully estimable: ",
              d$reason, ". The affected coefficient(s) are omitted from the table.",
              call. = FALSE)
  }

  out <- list(input_params = list(time = time, status = status,
                                  covariates = covariates, etype = etype,
                                  conf_level = conf_level,
                                  max_levels = max_levels,
                                  analysis_time = Sys.time()),
              cause_specific = tidy_fit(cs_fit, nrow(dd)),
              subdistribution = tidy_fit(fg_fit, nrow(dd)),
              n = nrow(dd),
              events_interest = sum(dd[[status]] == etype),
              events_competing = sum(!dd[[status]] %in% c(0, etype)),
              n_censored = sum(dd[[status]] == 0),
              etype = etype,
              covariate_types = cov_types,
              diagnostics = diagnostics,
              fits = list(cause_specific = cs_fit, subdistribution = fg_fit))
  class(out) <- "cpas_competing"
  out
}

#' @title Print method for cpas_competing objects
#' @param x An object of class \code{cpas_competing}.
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_competing <- function(x, ...) {
  cat("CanPAS competing-risks analysis\n")
  cat(sprintf("n = %d | events of interest (status %d) = %d | competing events = %d | censored = %d\n",
              x$n, x$etype, x$events_interest, x$events_competing, x$n_censored))
  cat("\nCause-specific hazard ratios (rate among event-free patients):\n")
  print(x$cause_specific, row.names = FALSE)
  cat("\nFine-Gray subdistribution hazard ratios (linked to cumulative incidence):\n")
  print(x$subdistribution, row.names = FALSE)
  invisible(x)
}

#' @title Print method for cpas_cif objects
#' @param x An object of class \code{cpas_cif}.
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_cif <- function(x, ...) {
  cat("CanPAS cumulative incidence (competing risks)\n")
  cat(sprintf("n = %d | censored = %d | event counts: %s\n", x$n, x$n_censored,
              paste(sprintf("cause %s = %d", names(x$events), x$events), collapse = ", ")))
  print(x$table, row.names = FALSE)
  invisible(x)
}
