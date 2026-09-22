# Internal utilities shared across the package (not exported).

# Fetch expression + survival of a dataset in one merged data.frame
.cpas_merged <- function(dataset, genes, process_duplicates = "max") {
  e <- get_expr_data(dataset, genes, process_duplicates = process_duplicates)
  merge_surv_expr(dataset, e)$merged_data
}

# Pick the first endpoint whose *_time/_status columns exist in `df`
.cpas_endpoint <- function(df, prefer = c("OS", "DSS", "DFS", "RFS", "PFS",
                                          "MFS", "DRFS", "EFS", "DFI", "PFI")) {
  have <- sub("_(time|status)$", "",
              colnames(df)[grepl("_(time|status)$", colnames(df))])
  have <- unique(have)
  hit <- intersect(prefer, have)
  if (length(hit)) hit[1] else NA_character_
}

# Marker (gene/probe) columns of a merged CanPAS data.frame
.cpas_markers <- function(df) {
  setdiff(colnames(df), c("ID", grep("_(status|time)$", colnames(df), value = TRUE)))
}

# Resolve requested type: use it if present, otherwise the first endpoint
.cpas_resolve_type <- function(df, type) {
  if (!is.null(type) && length(type) == 1L &&
      all(c(paste0(type, "_time"), paste0(type, "_status")) %in% colnames(df)))
    return(type)
  ep <- .cpas_endpoint(df)
  if (is.na(ep))
    stop("No survival endpoint columns (*_time/*_status) found in the data.")
  if (!is.null(type) && length(type) == 1L && !is.na(type))
    message("Requested endpoint '", type, "' is unavailable; using '", ep,
            "' instead.")
  ep
}

# Fit a Cox model and decide whether the estimate is usable ----------------
# Complete separation, non-convergence and near-singular designs do NOT make
# coxph fail: it returns a warning ("Loglik converged before variable 1;
# coefficient may be infinite") together with a finite but meaningless
# coefficient (e.g. HR 1.8e+09 with a 95% CI of [0, Inf]). The helper captures
# those warnings and the numeric red flags so callers can drop the estimate and
# report why.
# Inspect an already-fitted coxph model (and the warnings it produced) and say
# whether its coefficients are usable. Shared by the Cox, competing-risks and
# Fine-Gray paths so that all of them refuse the same meaningless results.
.cpas_COX_checks <- function(fit, msgs = character(0)) {
  if (is.null(fit))
    return(list(fit = NULL, ok = FALSE, reason = "the model could not be fitted",
                warnings = msgs))
  b <- tryCatch(stats::coef(fit), error = function(e) numeric(0))
  se <- tryCatch(sqrt(diag(stats::vcov(fit))), error = function(e) numeric(0))
  reason <- character(0)
  bad_msg <- grepl(paste0("converged before variable|may be infinite|did not converge|",
                          "Ran out of iterations|essentially perfect fit|singular|",
                          "coefficient may be infinite"),
                   msgs, ignore.case = TRUE)
  if (any(bad_msg)) reason <- c(reason, unique(msgs[bad_msg]))
  if (!length(b) || any(!is.finite(b)) || any(!is.finite(se)))
    reason <- c(reason, "non-finite coefficient or standard error")
  if (length(se) && any(se <= 0, na.rm = TRUE))
    reason <- c(reason, "zero or negative standard error")
  if (!is.null(fit$fail) && any(fit$fail > 0, na.rm = TRUE))
    reason <- c(reason, "iterations did not converge")
  if (length(b) && any(abs(b) > 20, na.rm = TRUE))
    reason <- c(reason, "extreme coefficient (|log HR| > 20): separation or unidentifiable effect")
  if (length(se) && any(se > 10, na.rm = TRUE))
    reason <- c(reason, "extremely large standard error (> 10 on the log-HR scale)")
  list(fit = fit, ok = !length(reason),
       reason = paste(unique(reason), collapse = "; "), warnings = msgs)
}

.cpas_COX_diag <- function(fml, data, ...) {
  msgs <- character(0)
  # model/x are kept so that post-fit routines which reconstruct the data from
  # the fit (e.g. survival::cox.zph) still work outside this helper frame
  fit <- withCallingHandlers(
    tryCatch(survival::coxph(fml, data = data, model = TRUE, x = TRUE, ...),
             error = function(e) e),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      tryCatch(invokeRestart("muffleWarning"), error = function(e) NULL)
    })
  if (inherits(fit, "error"))
    return(list(fit = NULL, ok = FALSE, reason = conditionMessage(fit), warnings = msgs))
  .cpas_COX_checks(fit, msgs)
}

# Screen covariates before fitting a multivariable Cox model.
#
# The package refuses to publish a non-estimable multivariable model (complete
# separation, collinearity or a constant covariate produce meaningless hazard
# ratios). Refusing the whole model is not helpful on its own: this helper
# locates the offending covariates so the caller can drop them, fit the
# remaining model and tell the user exactly what was removed and why.
#
#   * a continuous covariate that is constant or nearly constant (>= 95% of
#     patients share one value) is dropped;
#   * categorical levels with fewer than `min_level_n` patients, with no events,
#     or with events for every patient, are dropped as levels (they are the
#     usual cause of separation); a covariate that keeps fewer than two levels
#     is dropped as a whole;
#   * whatever still makes the joint model non-estimable is dropped one
#     variable at a time by .cpas_COX_drop_worst() in the caller.
#
# Returns the screened data, the covariates that remain, and a table of drops
# (variable / detail / reason) for reporting.
.cpas_COX_screen <- function(data, status = "status",
                            cont = character(0), cate = character(0),
                            min_level_n = 5L) {
  drops <- data.frame(variable = character(0), detail = character(0),
                      reason = character(0), stringsAsFactors = FALSE)
  add <- function(v, detail, reason) {
    drops[nrow(drops) + 1L, ] <<- list(as.character(v), as.character(detail),
                                       as.character(reason))
  }
  dat <- data
  keep <- character(0)

  for (v in cont) {
    x <- suppressWarnings(as.numeric(dat[[v]]))
    if (length(unique(x[!is.na(x)])) <= 1L) {
      add(v, "continuous", "constant across the complete cases"); next
    }
    if (mean(x == stats::median(x, na.rm = TRUE), na.rm = TRUE) >= 0.95) {
      add(v, "continuous", "nearly constant (\u2265 95% of patients share one value)"); next
    }
    keep <- c(keep, v)
  }

  for (v in cate) {
    x <- dat[[v]]
    if (!is.factor(x)) x <- factor(as.character(x))
    dat[[v]] <- x
    tb <- table(x)
    ev <- tapply(dat[[status]], x, function(s) sum(s == 1, na.rm = TRUE))
    ev[is.na(ev)] <- 0L
    n_bad <- names(tb)[tb < min_level_n]
    s_bad <- names(tb)[ev == 0L | ev == as.integer(tb)]
    bad <- unique(c(n_bad, s_bad))
    if (!length(bad)) { keep <- c(keep, v); next }
    rest <- setdiff(names(tb), bad)
    why <- character(0)
    if (length(n_bad))
      why <- c(why, sprintf("level(s) %s have fewer than %d patients (%s)",
                            paste(n_bad, collapse = "/"), min_level_n,
                            paste(as.integer(tb[n_bad]), collapse = "/")))
    sep_only <- setdiff(s_bad, n_bad)
    if (length(sep_only))
      why <- c(why, sprintf("level(s) %s have no events or events for every patient (%s of %s)",
                            paste(sep_only, collapse = "/"),
                            paste(as.integer(ev[sep_only]), collapse = "/"),
                            paste(as.integer(tb[sep_only]), collapse = "/")))
    if (length(rest) < 2L) {
      add(v, paste0("levels: ", paste(bad, collapse = ", ")),
          paste0(paste(why, collapse = "; "),
                 "; fewer than two usable levels remain, so the covariate was dropped"))
      next
    }
    # collapse the problematic levels instead of dropping their patients: a
    # merged "Other" level usually stays estimable, whereas a 1-patient level
    # or an all-events level pushes the coefficient to infinity
    lv <- as.character(x)
    lv[lv %in% bad] <- "Other"
    dat[[v]] <- factor(lv, levels = c(rest, "Other"))
    o_n  <- sum(dat[[v]] == "Other")
    o_ev <- sum(dat[[status]][dat[[v]] == "Other"] == 1, na.rm = TRUE)
    if (o_ev == 0L || o_ev == o_n) {
      # merging did not help (the merged level is still perfectly separated);
      # fall back to excluding the patients in those levels, which keeps the
      # covariate usable whenever two or more levels remain
      keep_n <- sum(tb[rest])
      if (length(rest) >= 2L && keep_n > 0) {
        lv2 <- as.character(x)
        lv2[lv2 %in% bad] <- NA_character_
        dat[[v]] <- factor(lv2, levels = rest)
        add(v, paste0("levels excluded from the model: ", paste(bad, collapse = ", ")),
            paste0(paste(why, collapse = "; "),
                   sprintf("; merging them into 'Other' stays perfectly separated (Other: %d patients, %d events), so those %d patients were excluded from the multivariable model and the covariate was kept with its remaining levels (%s)",
                           o_n, o_ev, length(bad), paste(rest, collapse = "/"))))
        keep <- c(keep, v)
        next
      }
      add(v, paste0("levels: ", paste(bad, collapse = ", ")),
          paste0(paste(why, collapse = "; "),
                 sprintf("; merged into 'Other' (%d patients, %d events) which stays perfectly separated and fewer than two levels would remain, so the covariate was dropped",
                         o_n, o_ev)))
      next
    }
    add(v, paste0("levels collapsed into 'Other': ", paste(bad, collapse = ", ")),
        paste0(paste(why, collapse = "; "),
               sprintf("; the patients were kept and merged into one 'Other' level (%d patients, %d events)",
                       o_n, o_ev)))
    keep <- c(keep, v)
  }

  list(data = dat, keep = keep, drops = drops)
}

# Identify the covariate to drop when a fitted multivariable model is not
# estimable: the one behind a non-finite coefficient, else the one with an
# extreme coefficient or standard error. Returns NA when no single covariate can
# be blamed (typically a collinear pair), leaving the caller to try removals in
# turn.
.cpas_COX_drop_worst <- function(diag, vars) {
  b  <- tryCatch(stats::coef(diag$fit), error = function(e) numeric(0))
  se <- tryCatch(sqrt(diag(stats::vcov(diag$fit))), error = function(e) numeric(0))
  if (!length(b) || is.null(names(b))) return(NA_character_)
  nm <- names(b)
  owner <- function(x) {
    if (length(x) != 1L || is.na(x)) return(NA_character_)
    hit <- vars[vapply(vars, function(v) identical(x, v) || startsWith(x, v), logical(1))]
    if (!length(hit)) return(NA_character_)
    as.character(hit[[which.max(nchar(hit))]])
  }
  bad <- which(!is.finite(b))
  if (!length(bad) && length(se) == length(b)) bad <- which(!is.finite(se) | se > 10)
  if (!length(bad)) bad <- which(abs(b) > 20)
  if (!length(bad)) return(NA_character_)
  for (i in bad[order(abs(b[bad]), decreasing = TRUE)]) {
    o <- owner(nm[i])
    if (length(o) == 1L && !is.na(o)) return(o)
  }
  NA_character_
}
