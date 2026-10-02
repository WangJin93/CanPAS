#' @title Two-stage cross-cohort integrative (meta) analysis
#' @description
#' Runs a two-stage integrative analysis over several cohorts of one cancer type:
#'   stage 1 fits \code{coxph(Surv(time, status) ~ marker + confounders)} in every
#'   cohort, giving log-HR +/- SE per cohort; stage 2 pools them by inverse
#'   variance and reports the pooled HR (95\% CI), Z/p, the heterogeneity
#'   statistics Q/df/I\eqn{^2}/tau\eqn{^2} and two prediction intervals.
#'
#' The pooling method is chosen by \code{method}: \code{"REML"} (default,
#' iterative restricted maximum likelihood for tau\eqn{^2}), \code{"DL"}
#' (DerSimonian-Laird, the previous default, kept for continuity), \code{"HK"}
#' (Hartung-Knapp-Sidik-Jonkman adjusted variance of the pooled estimate, with
#' its characteristic wider interval) or \code{"FE"} (fixed effect).
#' \code{"RE"} is still accepted as a synonym of \code{"DL"}.
#' @param datasets Character vector of dataset accessions (matching
#'   \code{dataset_info$Accession}).
#' @param marker One gene symbol (e.g. \code{"TP53"}) or a signature formula
#'   string (e.g. \code{"0.5*TP53+0.3*GAPDH"}).
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or a raw endpoint token
#'   (OS/DSS/DFS/RFS/PFS/MFS/DFI/PFI/DRFS/EFS, ...). A family is resolved per
#'   cohort to the concrete token that cohort provides (see
#'   \code{\link{endpoint_resolve}}) and recorded in \code{per_dataset$endpoint}
#'   and in \code{per_dataset$pooling_class}.
#' @param method Pooling method: \code{"REML"} (default), \code{"DL"},
#'   \code{"HK"} or \code{"FE"}; \code{"RE"} is accepted as \code{"DL"}.
#' @param pi_method Construction rule of the 95\% prediction interval for a new
#'   cohort (spec B6): \code{"t"} (default) is the documented
#'   \eqn{\hat\mu \pm t_{0.975,k-2}\sqrt{se_{pooled}^2+\tau^2}} interval
#'   (needs at least 3 cohorts); \code{"normal"} is the normal approximation
#'   \eqn{\hat\mu \pm 1.96\sqrt{se_{pooled}^2+\tau^2}} (needs at least 2);
#'   \code{"HK"} uses the Hartung-Knapp-Sidik-Jonkman adjusted standard error of
#'   the pooled estimate with \eqn{t_{0.975,k-1}} (needs at least 2). The rule
#'   actually used is recorded in \code{$pooled$pi_method} and
#'   \code{$pooled$pi_rule}, and the other construction is returned beside it,
#'   clearly labelled, in \code{$pooled$pi_alt_lower}/\code{pi_alt_upper} with
#'   \code{$pooled$pi_alt_rule}. Only the prediction interval is affected:
#'   \code{pi_method} does not change the pooled HR, its confidence interval, or
#'   any heterogeneity statistic.
#' @param confounders Optional character vector of covariate names to adjust for
#'   (they must exist in the merged data). With the default
#'   \code{auto_repair = FALSE}, a cohort in which a requested covariate is
#'   missing or constant is reported as
#'   \code{not estimable} instead of being adjusted for a reduced covariate set.
#' @param min_events Minimum number of events for a cohort to enter (default 5).
#' @param max_try Retries per cohort when fetching fails (API jitter).
#' @param merged Optional named list (cohort name -> already merged data.frame).
#'   When supplied nothing is fetched from the network.
#' @param pooling \code{"family"} (default) pools the Exact-equivalent and the
#'   Clinically-related rows of the endpoint family, i.e. the cohorts that
#'   contribute RFS/EFS/DFI to DFS, CSS/BCSS to DSS, PFI to PFS and DRFS to MFS
#'   as well as those whose token is the family itself. \code{"exact"} pools only
#'   rows whose pooling class is \code{Exact-equivalent}. Both modes are
#'   available because \code{"family"} is what makes the cross-token claim
#'   (the DFS pool would fall from 98 cohorts to the ~20 whose token is
#'   literally DFS); the mode actually used is recorded in
#'   \code{$pooled$pooling}, in \code{$per_dataset$pooling_class} and in
#'   \code{$manifest}.
#' @param overlap What to do about cohort pairs known to share patients:
#'   \code{"warn"} (default) proceeds but warns and records the pairs in
#'   \code{$pooled$overlap_pairs}, \code{"refuse"} stops with an actionable
#'   error, \code{"dedupe"} keeps one member of each overlapping group (the
#'   larger cohort, ties by accession sort) and records what it dropped in
#'   \code{$excluded} and \code{$manifest$overlap_dropped}.
#' @param auto_repair \code{FALSE} (default) is the fail-safe: a cohort whose
#'   model would need covariates dropped, levels merged or rows excluded is
#'   reported as \code{not estimable} with its reasons and offending terms in
#'   \code{$excluded}, and no model is returned for it. \code{TRUE} restores the
#'   previous behaviour (repair the model) and records every modification in
#'   \code{$modifications}.
#' @param class_table Optional endpoint-semantics table (the frozen schema of
#'   \code{\link{endpoint_semantics}}) used instead of the shipped copy.
#' @param overlap_table Optional shared-patient register (the frozen schema of
#'   \code{\link{cohort_overlap}}) used instead of the shipped copy.
#' @return A list of class \code{cpas_meta}:
#'   \item{\code{per_dataset}:}{one row per pooled cohort with \code{dataset},
#'     \code{endpoint} (the token used), \code{pooling_class}, \code{n},
#'     \code{events}, \code{HR}, \code{lower}, \code{upper}, \code{logHR},
#'     \code{se}, \code{p}}
#'   \item{\code{pooled}:}{one row: the method used, \code{pi_method} (the
#'     prediction-interval rule used), \code{k}, \code{HR},
#'     \code{lower}, \code{upper}, \code{p}, \code{Q}, \code{df},
#'     \code{p_heterogeneity}, \code{I2}, \code{tau2}, the primary prediction
#'     interval \code{pi_lower}/\code{pi_upper} (by default t with k-2 df) and
#'     the alternative interval \code{pi_alt_lower}/\code{pi_alt_upper} (by
#'     default the normal approximation), plus \code{pooling},
#'     \code{pooling_classes}, \code{overlap_pairs} and \code{overlap_mode}}
#'   \item{\code{excluded}:}{cohorts that did not enter and why
#'     (\code{cohort}, \code{stage}, \code{reason}); \code{stage = "model"} rows
#'     are the fail-safe refusals}
#'   \item{\code{modifications}:}{every modification \code{auto_repair = TRUE}
#'     made (\code{cohort}, \code{variable}, \code{detail}, \code{reason})}
#'   \item{\code{errors}, \code{input}, \code{manifest}:}{fetch/endpoint errors
#'     per cohort, the echoed inputs, and the analysis manifest
#'     (\code{\link{cpas_manifest}})}
#' @examples
#' \dontrun{
#'    ## Real cohorts, real endpoint: each cohort contributes the DFS-family
#'    ## token it actually has, which the function resolves and reports.
#'    m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
#'    m$per_dataset[, c("dataset", "endpoint", "pooling_class", "n", "events", "HR")]
#'    m$pooled
#'    plot_meta_forest(m)
#'    loo_meta(m)               # leave-one-out sensitivity
#'    cpas_manifest(m)          # what was pooled, and under which defaults
#'    print(m)
#'
#'    ## Strict mode: only cohorts whose token IS the family definition.
#'    m2 <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS",
#'                    pooling = "exact")
#'    m2$per_dataset$pooling_class
#'
#'    ## Prediction-interval rule (spec B6): the default t(k-2) interval, the
#'    ## normal approximation, or the HK-based interval. The interval used is
#'    ## recorded in $pooled$pi_method.
#'    m_normal <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53",
#'                          type = "DFS", pi_method = "normal")
#'    m_normal$pooled[, c("pi_method", "pi_lower", "pi_upper", "pi_alt_lower")]
#' }
#' @seealso \code{\link{endpoint_semantics}}, \code{\link{cohort_overlap}},
#'   \code{\link{cpas_manifest}}
#' @export
cpas_meta <- function(datasets, marker, type = "OS",
                      method = c("REML", "DL", "HK", "FE"),
                      pi_method = c("t", "normal", "HK"),
                      confounders = NULL,
                      min_events = 5,
                      max_try = 3,
                      merged = NULL,
                      pooling = c("family", "exact"),
                      overlap = c("warn", "refuse", "dedupe"),
                      auto_repair = FALSE,
                      class_table = NULL,
                      overlap_table = NULL) {
  method_req <- as.character(method)[1]
  method <- .cpas_meta_method(method_req)
  pi_method <- .cpas_pi_method(pi_method)
  pooling <- match.arg(pooling)
  overlap <- match.arg(overlap)
  if (!is.logical(auto_repair) || length(auto_repair) != 1L || is.na(auto_repair))
    stop("'auto_repair' must be TRUE or FALSE.", call. = FALSE)
  if (missing(datasets) || is.null(datasets))
    stop("'datasets' must contain at least one accession.", call. = FALSE)
  datasets <- as.character(datasets)
  if (!length(datasets)) stop("'datasets' must contain at least one accession.", call. = FALSE)
  datasets_req <- datasets
  genes <- if (grepl("[+*:]", marker)) {
    unique(unlist(regmatches(marker, gregexpr("[A-Za-z][A-Za-z0-9._]*", marker))))
  } else marker

  if (anyDuplicated(datasets))
    stop("'datasets' contains duplicated accessions: ",
         paste(unique(datasets[duplicated(datasets)]), collapse = ", "),
         ". Each dataset can be used once.", call. = FALSE)
  if (!is.null(merged) && anyDuplicated(names(merged)))
    stop("'merged' contains duplicated names: ",
         paste(unique(names(merged)[duplicated(names(merged))]), collapse = ", "),
         ". Use unique dataset names.", call. = FALSE)

  family <- endpoint_family(type)
  if (is.na(family)) family <- as.character(type)[1]

  # ---- shared-patient overlap --------------------------------------------
  ov_tab <- overlap_table
  if (is.null(ov_tab))
    ov_tab <- tryCatch(cohort_overlap(), error = function(e) NULL)
  ov_source <- if (is.null(ov_tab)) "unavailable (companion table not shipped)" else
    "inst/extdata/cohort_overlap.csv"
  ov_in <- .cpas_overlap_pairs_in(datasets, table = ov_tab)
  ov_pairs <- ov_in$table
  ov_drop <- .cpas_empty_df(c("cohort", "kept", "shared_patients", "group", "reason"))
  excluded <- .cpas_empty_df(c("cohort", "stage", "reason"))
  add_excluded <- function(cohort, stage, reason) {
    excluded[nrow(excluded) + 1L, ] <<- list(as.character(cohort), as.character(stage),
                                             as.character(reason))
  }
  if (!is.null(ov_pairs) && nrow(ov_pairs)) {
    if (overlap == "refuse")
      stop("Refusing to pool cohorts that share patients (overlap = \"refuse\"): ",
           paste(sprintf("%s / %s share %s patients", ov_pairs$AccessionA,
                         ov_pairs$AccessionB, ov_pairs$SharedPatients),
                 collapse = "; "),
           ". Drop one member of each pair, or use overlap = \"dedupe\" to let ",
           "the function keep the larger cohort of every overlapping group.",
           call. = FALSE)
    if (overlap == "warn")
      warning("Cohorts sharing patients are being pooled (overlap = \"warn\"): ",
              paste(sprintf("%s / %s (%s shared)", ov_pairs$AccessionA,
                            ov_pairs$AccessionB, ov_pairs$SharedPatients),
                    collapse = "; "),
              ". Their patients are counted twice, which narrows the pooled ",
              "interval without adding information; use overlap = \"dedupe\" to ",
              "keep one member of each group, or overlap = \"refuse\" to stop.",
              call. = FALSE)
  } else if (is.null(ov_tab) && overlap != "warn") {
    stop("overlap = \"", overlap, "\" needs the shared-patient register, but ",
         "it is not shipped with this build (inst/extdata/cohort_overlap.csv ",
         "is missing; it is produced by the curation pipeline). Pass the pairs ",
         "through overlap_table=, or use overlap = \"warn\" to proceed without ",
         "the check.", call. = FALSE)
  }
  if (overlap == "dedupe" && !is.null(ov_pairs) && nrow(ov_pairs)) {
    n_of <- list()
    di <- .cpas_dataset_info()
    for (t in datasets) {
      n <- NA_real_
      if (t %in% di$Accession) {
        v <- suppressWarnings(as.numeric(di$N[match(t, di$Accession)]))
        if (length(v) && is.finite(v)) n <- v
      }
      if (!is.finite(n) && !is.null(merged[[t]])) n <- nrow(merged[[t]])
      n_of[[t]] <- n
    }
    dd <- .cpas_overlap_dedupe(datasets, ov_pairs, n_of)
    ov_drop <- dd$dropped
    for (i in seq_len(nrow(ov_drop)))
      add_excluded(ov_drop$cohort[i], "overlap", ov_drop$reason[i])
    datasets <- dd$keep
    ov_in <- .cpas_overlap_pairs_in(datasets, table = ov_tab)
    ov_pairs <- ov_in$table
  }

  per <- list(); errors <- list(); mods <- list(); rows_dropped <- list()
  ph_list <- list(); lookup <- list()
  received_counts <- integer(0); analyzed_inputs <- list()
  for (t in datasets) {
    r <- NULL
    if (!is.null(merged) && !is.null(merged[[t]])) {
      df <- merged[[t]]
      ok <- TRUE
    } else {
      ok <- FALSE; err <- NULL
      if (startsWith(t, "TCGA-")) {
        r <- tryCatch(list(df = tcga_merged(sub("^TCGA-", "", t), genes, type)),
                      error = function(e) e)
        if (!inherits(r, "error")) ok <- TRUE else err <- r
      } else {
        for (i in seq_len(max_try)) {
          r <- tryCatch({
            ex <- get_expr_data(t, genes, process_duplicates = "max")
            es <- merge_surv_expr(t, ex)
            list(df = es$merged_data, geneCols = setdiff(colnames(ex$expr_data), "ID"))
          }, error = function(e) e)
          if (!inherits(r, "error")) { ok <- TRUE; break }
          err <- r; Sys.sleep(2 + 1.5 * i)
        }
      }
      if (!ok) {
        errors[[t]] <- paste0("fetch/merge: ", substr(err$message, 1, 100))
        add_excluded(t, "fetch", errors[[t]])
        next
      }
      df <- r$df
    }
    received_counts[t] <- nrow(df)
    if ("ID" %in% names(df) && (anyNA(df$ID) || anyDuplicated(df$ID)))
      stop("Dataset ", t, ": sample ID values must be non-missing and unique.", call. = FALSE)
    # family (OS/DSS/DFS/PFS/MFS) -> the concrete endpoint this cohort has; a
    # raw token is also accepted
    tok <- endpoint_resolve(t, type)
    if (is.na(tok)) {
      # fallback: a cohort outside the catalog (user-supplied data) whose data
      # already carry the endpoint columns is handled under the raw token
      if (all(c(paste0(type, "_time"), paste0(type, "_status")) %in% colnames(df)))
        tok <- as.character(type)[1]
      else {
        errors[[t]] <- paste0("no ", type, " endpoint")
        add_excluded(t, "endpoint", errors[[t]])
        next
      }
    }
    tc <- paste0(tok, "_time"); sc <- paste0(tok, "_status")
    if (!all(c(tc, sc) %in% colnames(df))) {
      errors[[t]] <- paste0("no ", type, " endpoint (resolved ", tok,
                            " not in merged data)")
      add_excluded(t, "endpoint", errors[[t]])
      next
    }
    cls <- .cpas_pooling_lookup(t, family, table = class_table, tokens = tok)
    if (!is.na(cls$token) && !identical(as.character(cls$token), as.character(tok))) {
      add_excluded(t, "pooling", sprintf("endpoint token conflict: analysis uses %s but classification table records %s",
                                        tok, cls$token))
      next
    }
    cls$family <- family
    lookup[[t]] <- cls
    allowed_classes <- if (identical(pooling, "exact")) "Exact-equivalent" else
      c("Exact-equivalent", "Clinically-related")
    if (!cls$pooling_class %in% allowed_classes) {
      add_excluded(t, "pooling",
                   sprintf("pooling = \"%s\" keeps only %s endpoint classes; this cohort's token is %s (%s in family %s)",
                           pooling, paste(allowed_classes, collapse = " or "), tok,
                           cls$pooling_class, family))
      next
    }
    df[[tc]] <- suppressWarnings(as.numeric(df[[tc]]))
    df[[sc]] <- suppressWarnings(as.numeric(df[[sc]]))
    if (any(!is.na(df[[sc]]) & (!is.finite(df[[sc]]) | !df[[sc]] %in% c(0, 1))))
      stop("Dataset ", t, ": status must be coded 0 (censored) / 1 (event).", call. = FALSE)
    if (any(!is.na(df[[tc]]) & df[[tc]] < 0))
      stop("Dataset ", t, ": negative survival time is not allowed.", call. = FALSE)
    for (gn in genes) if (gn %in% colnames(df)) df[[gn]] <- suppressWarnings(as.numeric(df[[gn]]))

    if (grepl("[+*:]", marker)) {
      gv <- intersect(genes, colnames(df))
      missing_sig <- setdiff(genes, colnames(df))
      if (!length(gv)) {
        errors[[t]] <- paste0("signature genes absent: ",
                              paste(missing_sig, collapse = ","))
        add_excluded(t, "marker", errors[[t]])
        next
      }
      score <- tryCatch(eval(parse(text = marker), envir = as.list(df[gv])),
                        error = function(e) NULL)
      if (is.null(score) || all(is.na(score))) {
        errors[[t]] <- paste0("signature eval failed (missing on platform: ",
                              paste(missing_sig, collapse = ","), ")")
        add_excluded(t, "marker", errors[[t]])
        next
      }
      df$marker <- score
    } else {
      if (!marker %in% colnames(df)) {
        errors[[t]] <- paste0("gene missing: ", marker)
        add_excluded(t, "marker", errors[[t]])
        next
      }
      df$marker <- suppressWarnings(as.numeric(df[[marker]]))
    }

    # analysis sample = complete cases of (time, status, marker) and the
    # covariates used, so that standardisation and the model use one sample
    missing_cov <- setdiff(confounders, colnames(df))
    if (length(missing_cov)) {
      if (!isTRUE(auto_repair)) {
        add_excluded(t, "model", paste0("not estimable (auto_repair = FALSE): requested covariate(s) absent: ",
                                         paste(missing_cov, collapse = ", ")))
        next
      }
      for (v in missing_cov)
        mods[[length(mods) + 1L]] <- data.frame(cohort = t, variable = v,
          detail = "absent from supplied data", reason = "requested covariate dropped with explicit auto_repair = TRUE",
          stringsAsFactors = FALSE)
    }
    conf_all <- intersect(confounders, colnames(df))
    conf_avail <- conf_all
    if (length(conf_all)) {
      okv <- vapply(conf_all, function(cn)
        length(unique(df[[cn]][!is.na(df[[cn]])])) > 1L, logical(1))
      conf_avail <- conf_all[okv]
    }
    const_cov <- setdiff(conf_all, conf_avail)
    if (length(const_cov)) {
      detail <- "constant in this cohort"
      reason <- sprintf(paste0("covariate(s) %s take a single value in this cohort, ",
                               "so they cannot be co-estimated"),
                        paste(const_cov, collapse = ", "))
      for (v in const_cov)
        mods[[length(mods) + 1L]] <- data.frame(cohort = t, variable = v,
                                                detail = detail, reason = reason,
                                                stringsAsFactors = FALSE)
      if (!isTRUE(auto_repair)) {
        add_excluded(t, "model", paste0("not estimable (auto_repair = FALSE): ", reason,
                                        ". The covariate(s) ", paste(const_cov, collapse = ", "),
                                        " would have to be dropped; set auto_repair = TRUE ",
                                        "to repair the model and record the modification."))
        next
      }
    }
    for (cn in conf_avail)
      if (is.character(df[[cn]])) {
        num <- suppressWarnings(as.numeric(df[[cn]]))
        df[[cn]] <- if (all(is.na(num) == is.na(df[[cn]]))) num else factor(df[[cn]])
      } else if (is.logical(df[[cn]])) df[[cn]] <- factor(df[[cn]])
    req <- c(tc, sc, "marker", conf_avail)
    n_before <- nrow(df)
    keep <- stats::complete.cases(df[req]) & is.finite(df[[tc]]) &
      is.finite(df$marker) & df[[tc]] >= 0
    df <- df[keep, , drop = FALSE]
    if (n_before - nrow(df) > 0)
      rows_dropped[[length(rows_dropped) + 1L]] <- data.frame(
        cohort = t, n_dropped = n_before - nrow(df),
        reason = "incomplete or non-finite (time, status, marker, covariates)",
        stringsAsFactors = FALSE)
    if (nrow(df) < 10) {
      errors[[t]] <- sprintf("insufficient after cleaning (n=%d)", nrow(df))
      add_excluded(t, "sample", errors[[t]])
      next
    }
    sdx <- stats::sd(df$marker)
    if (is.na(sdx) || sdx <= 0) {
      errors[[t]] <- "marker constant"
      add_excluded(t, "marker", errors[[t]])
      next
    }
    df$marker <- (df$marker - mean(df$marker)) / sdx   # per-SD scaling
    events <- sum(df[[sc]] == 1, na.rm = TRUE)
    if (events < min_events || length(unique(df[[sc]])) < 2) {
      errors[[t]] <- sprintf("insufficient (n=%d, events=%d)", nrow(df), events)
      add_excluded(t, "sample", errors[[t]])
      next
    }
    vars <- c("marker", conf_avail)
    diag <- .cpas_COX_diag(stats::as.formula("survival::Surv(time, status) ~ ."),
                           data.frame(time = df[[tc]], status = df[[sc]],
                                      df[, vars, drop = FALSE]))
    if (!diag$ok) {
      errors[[t]] <- paste0("cox: ", substr(diag$reason, 1, 120))
      add_excluded(t, "model", errors[[t]])
      next
    }
    fit <- diag$fit
    analyzed_inputs[[t]] <- df[, unique(c("ID", tc, sc, vars))[unique(c("ID", tc, sc, vars)) %in% names(df)], drop = FALSE]
    ph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
    if (!is.null(ph)) {
      tb <- as.data.frame(ph$table)
      ph_list[[length(ph_list) + 1L]] <- data.frame(
        cohort = t, term = rownames(tb), p = tb$p, stringsAsFactors = FALSE)
    }
    b <- stats::coef(fit)[["marker"]]; se <- sqrt(diag(stats::vcov(fit)))[["marker"]]
    per[[t]] <- data.frame(dataset = t, endpoint = tok,
                           pooling_class = cls$pooling_class,
                           n = nrow(df), events = events,
                           covariates_used = paste(conf_avail, collapse = ", "),
                           n_parameters = length(stats::coef(fit)),
                           events_per_parameter = events / length(stats::coef(fit)),
                           HR = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
                           logHR = b, se = se,
                           p = 2 * stats::pnorm(-abs(b / se)), stringsAsFactors = FALSE)
  }

  if (!length(per)) {
    # keep the per-cohort reasons: the caller needs to know whether a gene was
    # missing on the platform, the marker was constant, the events were few, the
    # token was not exact under pooling = "exact", or the model needed repair
    why <- if (nrow(excluded))
      paste(sprintf("%s: %s", excluded$cohort, excluded$reason), collapse = "; ")
    else if (length(errors))
      paste(sprintf("%s: %s", names(errors), unlist(errors)), collapse = "; ")
    else "no reason recorded"
    stop("no dataset produced estimable results. Reasons: ", why, call. = FALSE)
  }
  if (length(per) == 1L)
    message("Only one dataset produced estimable results; the 'pooled' row is that dataset's estimate.")
  pc <- do.call(rbind, per)
  pc <- pc[order(match(pc$dataset, datasets)), , drop = FALSE]
  rownames(pc) <- NULL
  tok_used <- unique(pc$endpoint)
  if (length(tok_used) > 1L)
    warning("Cohorts contribute different endpoint tokens within family '", type, "': ",
            paste(sprintf("%s=%s", pc$dataset, pc$endpoint), collapse = ", "),
            ". Tokens of one family are pooled by design; the per-cohort token is kept in ",
            "per_dataset$endpoint and should be reported.", call. = FALSE)
  pooled <- meta_pool(pc$logHR, pc$se, p = pc$p, method = method,
                      pi_method = pi_method)
  pooled$total_n <- sum(pc$n)
  pooled$total_events <- sum(pc$events)
  pooled$pooling <- pooling
  pooled$pooling_classes <- paste(sort(unique(pc$pooling_class)), collapse = ", ")
  pooled$overlap_mode <- overlap
  # The pairs stay a data.frame, so they are stored as a one-element list column
  # (the standard way to keep a table inside a rectangular data.frame) and the
  # plain data.frame is also returned as $overlap_pairs on the result.
  pooled[["overlap_pairs"]] <- list(ov_pairs)
  pooled$overlap_pairs_text <- if (nrow(ov_pairs))
    paste(sprintf("%s/%s (%s shared)", ov_pairs$AccessionA, ov_pairs$AccessionB,
                  ov_pairs$SharedPatients), collapse = "; ") else ""
  pooled$overlap_source <- ov_source
  mods_tb <- if (length(mods)) do.call(rbind, mods) else
    .cpas_empty_df(c("cohort", "variable", "detail", "reason"))
  rd_tb <- if (length(rows_dropped)) do.call(rbind, rows_dropped) else
    .cpas_empty_df(c("cohort", "n_dropped", "reason"))
  rd_tb$n_dropped <- as.integer(rd_tb$n_dropped)
  lookup_tb <- if (length(lookup)) {
    tb <- do.call(rbind, lookup)
    tb <- tb[, c("accession", "family", "token", "token_role",
                 "pooling_class", "pooling_class_source"), drop = FALSE]
    rownames(tb) <- NULL
    tb[tb$accession %in% pc$dataset, , drop = FALSE]
  } else .cpas_empty_df(c("accession", "family", "token", "token_role",
                          "pooling_class", "pooling_class_source"))
  ph_tb <- if (length(ph_list)) do.call(rbind, ph_list) else NULL
  notes <- character(0)
  if (!isTRUE(auto_repair))
    notes <- c(notes, "auto_repair = FALSE (fail-safe): a model that would need covariates dropped, levels merged or rows excluded is reported as not estimable instead of being repaired")
  if (!identical(pi_method, "t"))
    notes <- c(notes, sprintf("pi_method = \"%s\": the primary prediction interval uses %s instead of the default t(k-2) construction (the default interval is returned as the clearly labelled alternative)",
                              pi_method, pooled$pi_rule))
  if (identical(pooling, "exact"))
    notes <- c(notes, sprintf("pooling = \"exact\": only Exact-equivalent rows were pooled (%d of %d requested cohorts entered)",
                              nrow(pc), length(datasets)))
  if (nrow(ov_drop))
    notes <- c(notes, sprintf("overlap = \"dedupe\" removed %d cohort(s): %s",
                              nrow(ov_drop), paste(ov_drop$cohort, collapse = ", ")))
  if (is.na(ov_source) || grepl("^unavailable", ov_source))
    notes <- c(notes, "the shared-patient register is not shipped in this build, so overlap could not be checked")
  manifest <- .cpas_manifest_new(
    analysis = "cpas_meta",
    cohorts = pc$dataset, accession = pc$dataset,
    raw_token = pc$endpoint, pooling_class = pc$pooling_class,
    n_input = sum(received_counts), n_analyzed = sum(pc$n),
    n_excluded = sum(received_counts) - sum(pc$n),
    hash_scope = "serialized analyzed per-cohort model input list; not upstream raw dataset",
    events = sum(pc$events), marker_requested = marker, marker_definition = marker,
    covariates_requested = confounders,
    fitted_covariates = stats::setNames(pc$covariates_used, pc$dataset),
    estimator = paste0("survival::coxph; random-effects ", pooled$method),
    inference = paste0("Wald cohort estimates; ", pooled$se_method, "; PI: ", pooled$pi_rule),
    analyzed_data = analyzed_inputs,
    dataset_hash = .cpas_hash_analyzed_inputs(analyzed_inputs),
    coverage_status = "measured received-input and analyzed-row audit; unfetched inputs have unknown counts",
    family = family,
    token = paste(sort(unique(pc$endpoint)), collapse = ", "),
    tokens = stats::setNames(pc$endpoint, pc$dataset),
    token_role = paste(sort(unique(lookup_tb$token_role)), collapse = ", "),
    pooling = pooling,
    pooling_classes = strsplit(pooled$pooling_classes, ", ", fixed = TRUE)[[1]],
    pooling_table = lookup_tb,
    overlap_mode = overlap, overlap_pairs = ov_pairs, overlap_dropped = ov_drop,
    overlap_source = ov_source,
    selection_rule = sprintf(paste0("the %d cohort(s) supplied in 'datasets' were requested for family %s ",
                                    "and the concrete token each one provides was resolved from the catalog; ",
                                    "%d cohort(s) entered the pool and %d were excluded (reasons in $excluded)"),
                             length(datasets), family, nrow(pc), nrow(excluded)),
    dropped_rows = rd_tb, dropped_covariates = mods_tb,
    cut_rule = "not applicable (meta-analysis of a continuous marker; no cut-point is searched)",
    meta_method = pooled$method,
    auto_repair = isTRUE(auto_repair),
    pi_method = pi_method,
    tau2 = pooled$tau2, I2 = pooled$I2,
    pi_primary = c(lower = pooled$pi_lower, upper = pooled$pi_upper),
    pi_primary_rule = pooled$pi_rule,
    pi_alt = c(lower = pooled$pi_alt_lower, upper = pooled$pi_alt_upper),
    pi_alt_rule = pooled$pi_alt_rule,
    ph_test = if (is.null(ph_tb)) NULL else list(table = ph_tb),
    notes = notes)
  out <- list(input = list(datasets = datasets, datasets_requested = datasets_req,
                           marker = marker, type = type, family = family,
                           method = method, method_requested = method_req,
                           pi_method = pi_method,
                           pooling = pooling, overlap = overlap,
                           auto_repair = auto_repair, confounders = confounders,
                           time = Sys.time()),
              settings = list(method = method, method_requested = method_req,
                              pi_method = pi_method,
                              pooling = pooling, overlap = overlap,
                              auto_repair = auto_repair),
              per_dataset = pc, pooled = pooled, errors = errors,
              excluded = excluded,
              not_estimable = excluded[excluded$stage == "model", , drop = FALSE],
              modifications = mods_tb,
              overlap_pairs = ov_pairs,
              manifest = manifest)
  class(out) <- "cpas_meta"
  out
}

# The four pooling methods; "RE" is the historical spelling of "DL".
.cpas_meta_methods <- c("REML", "DL", "HK", "FE")

.cpas_meta_method <- function(method) {
  m <- toupper(as.character(method)[1])
  if (is.na(m) || !nzchar(m))
    stop("'method' must be one of ", paste(.cpas_meta_methods, collapse = ", "), ".",
         call. = FALSE)
  if (identical(m, "RE")) m <- "DL"          # today's estimator, legacy spelling
  if (!m %in% .cpas_meta_methods)
    stop("'method' must be one of ", paste(.cpas_meta_methods, collapse = ", "),
         " (\"RE\" is accepted as \"DL\"); got \"", method[1], "\".", call. = FALSE)
  m
}

# DerSimonian-Laird tau^2 (the estimator this package used before REML)
.cpas_tau2_dl <- function(b, se, w = 1 / se^2) {
  k <- length(b)
  bfe <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bfe)^2)
  max(0, (Q - k + 1) / (sum(w) - sum(w^2) / sum(w)))
}

# REML tau^2: the root of the REML estimating equation of the standard
# random-effects model,
#   g(tau2) = sum(w^2 * ((b - mu)^2 - v)) / sum(w^2) + 1 / sum(w) - tau2 = 0,
#   w = 1/(v + tau2), mu = sum(w b)/sum(w)
# solved by bisection (only a sign change is needed, so it cannot stall the way
# a fixed-point iteration does when the map's derivative approaches 1). A root
# at the boundary means the estimate is truncated at tau^2 = 0. Falls back to DL
# (with a recorded note) when no bracket can be found.
.cpas_tau2_reml <- function(b, se, tol = 1e-12, maxit = 300L) {
  v <- se^2
  start <- .cpas_tau2_dl(b, se)
  g <- function(t2) {
    w <- 1 / (v + t2)
    mu <- sum(w * b) / sum(w)
    sum(w^2 * ((b - mu)^2 - v)) / sum(w^2) + 1 / sum(w) - t2
  }
  g0 <- g(0)
  if (!is.finite(g0))
    return(list(tau2 = start, converged = FALSE, truncated = FALSE,
                iterations = 0L, start = start))
  if (g0 <= 0)
    return(list(tau2 = 0, converged = TRUE, truncated = TRUE,
                iterations = 0L, start = start))
  hi <- max(start, 1e-8)
  it <- 0L
  while (g(hi) > 0 && it < maxit) { hi <- hi * 2 + 1e-6; it <- it + 1L }
  if (g(hi) > 0)
    return(list(tau2 = start, converged = FALSE, truncated = FALSE,
                iterations = it, start = start))
  lo <- 0
  for (i in seq_len(200L)) {
    mid <- (lo + hi) / 2
    if (g(mid) > 0) lo <- mid else hi <- mid
    if (hi - lo <= tol * (1 + hi)) break
  }
  list(tau2 = (lo + hi) / 2, converged = TRUE, truncated = FALSE,
       iterations = it + 200L, start = start)
}

# Hartung-Knapp-Sidik-Jonkman adjusted variance of the pooled estimate:
#   se = sqrt( (1/(k-1)) * sum(w (b - mu)^2) / sum(w) ),  95% CI = mu +/- t(k-1) * se
.cpas_hk <- function(b, w, mu, se_re) {
  k <- length(b)
  if (k < 2L) return(list(se = se_re, note = "HK needs at least 2 cohorts; the random-effects SE was used"))
  denom <- (k - 1) * sum(w)
  num <- sum(w * (b - mu)^2)
  se <- if (denom > 0 && num > 0) sqrt(num / denom) else se_re
  list(se = se, note = NA_character_)
}

meta_pool <- function(b, se, p = NULL, method = c("REML", "DL", "HK", "FE"),
                      pi_method = c("t", "normal", "HK")) {
  method_req <- as.character(method)[1]
  method <- .cpas_meta_method(method_req)
  pi_method <- .cpas_pi_method(pi_method)
  k <- length(b)
  if (k == 0L) stop("meta_pool(): no dataset estimates were supplied.")
  if (length(se) != k) stop("meta_pool(): 'b' and 'se' must have the same length.")
  w <- 1 / se^2
  bfe <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bfe)^2)
  df <- k - 1
  zi <- if (!is.null(p) && length(p) == k) sign(b) * stats::qnorm(1 - p / 2) else NULL
  pi_lower <- NA_real_; pi_upper <- NA_real_
  pi_alt_lower <- NA_real_; pi_alt_upper <- NA_real_
  # The prediction interval of a NEW cohort.  Three construction rules, selected
  # by pi_method; the rule actually used is recorded in pi_rule / pi_alt_rule:
  #   "t"      bm +/- t(0.975, k-2) * sqrt(se_pooled^2 + tau2)   (the documented
  #            default; needs k >= 3)
  #   "normal" bm +/- 1.96         * sqrt(se_pooled^2 + tau2)     (needs k >= 2)
  #   "HK"     bm +/- t(0.975, k-1) * sqrt(se_HK^2 + tau2), where se_HK is the
  #            Hartung-Knapp-Sidik-Jonkman adjusted SE of the pooled estimate
  #            (needs k >= 2)
  # Whichever is the primary rule, the other construction is returned beside it
  # as a clearly labelled alternative (pi_alt_*) so a paper can report both.
  pi_rules <- .cpas_pi_rule_text(pi_method)
  pi_rule <- pi_rules$primary
  pi_alt_rule <- pi_rules$alt
  if (k == 1L) {
    # A single cohort is not a meta-analysis: report that cohort's own estimate
    # and leave the heterogeneity statistics undefined (NA, never NaN).
    return(data.frame(
      method = method, method_requested = method_req, tau2_method = method,
      pi_method = pi_method,
      note = if (identical(method, "HK")) "HK needs at least 2 cohorts; a single cohort's own estimate is reported" else NA_character_,
      k = 1L,
      total_n = NA_integer_, total_events = NA_integer_,
      HR = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
      logHR = b, se = se, p = 2 * stats::pnorm(-abs(b / se)),
      Q = 0, df = 0L, p_heterogeneity = NA_real_, I2 = NA_real_, tau2 = NA_real_,
      pi_lower = NA_real_, pi_upper = NA_real_,
      pi_alt_lower = NA_real_, pi_alt_upper = NA_real_,
      pi_rule = pi_rule, pi_alt_rule = pi_alt_rule,
      z_stouffer = if (is.null(zi)) NA_real_ else sum(zi),
      p_stouffer = if (is.null(zi)) NA_real_ else 2 * stats::pnorm(-abs(sum(zi)))))
  }
  note <- NA_character_
  tau2_method <- method
  if (identical(method, "FE")) {
    tau2 <- 0
    tau2_method <- "FE (tau2 fixed at 0)"
  } else if (identical(method, "DL")) {
    tau2 <- .cpas_tau2_dl(b, se, w)
  } else {
    # REML and HK share the REML tau^2 (HK is a variance adjustment on top of
    # the random-effects model, not a separate tau^2 estimator)
    r <- .cpas_tau2_reml(b, se)
    tau2 <- r$tau2
    tau2_method <- "REML"
    if (!isTRUE(r$converged)) {
      tau2 <- r$start
      tau2_method <- "DL (REML did not converge)"
      note <- sprintf("REML tau^2 did not converge within %d iterations; the DerSimonian-Laird estimate was used instead", 200L)
    } else if (isTRUE(r$truncated)) {
      note <- "the REML estimate was truncated at tau^2 = 0"
    }
  }
  w2 <- 1 / (se^2 + tau2)
  bm <- sum(w2 * b) / sum(w2)
  seM <- sqrt(1 / sum(w2))
  se_re <- seM
  se_method <- "inverse-variance random effects"
  crit <- 1.96
  pv <- 2 * stats::pnorm(-abs(bm / seM))
  if (identical(method, "HK")) {
    hk <- .cpas_hk(b, w2, bm, seM)
    seM <- hk$se
    se_method <- "Hartung-Knapp-Sidik-Jonkman adjusted"
    if (!is.na(hk$note)) note <- if (is.na(note)) hk$note else paste(note, hk$note, sep = "; ")
    crit <- stats::qt(0.975, df = k - 1)
    pv <- 2 * stats::pt(-abs(bm / seM), df = k - 1)
  }
  I2 <- if (Q > 0) max(0, (Q - df) / Q) else 0
  # The two ingredients every prediction-interval rule uses: the RE variance of
  # the pooled estimate (se_re, before HK adjustment) plus tau^2, and
  # the same sum with the HK-adjusted SE, used by pi_method = "HK".
  se_hk <- .cpas_hk(b, w2, bm, se_re)$se
  se_pi <- sqrt(se_re^2 + tau2)
  se_pi_hk <- sqrt(se_hk^2 + tau2)
  pi_use <- .cpas_pi_bounds(bm, pi_method, se_pi, se_pi_hk, k)
  pi_lower <- pi_use$lower; pi_upper <- pi_use$upper
  alt_use <- .cpas_pi_bounds(bm, pi_use$alt, se_pi, se_pi_hk, k)
  pi_alt_lower <- alt_use$lower; pi_alt_upper <- alt_use$upper
  zs <- NA_real_; ps <- NA_real_
  if (!is.null(zi)) {
    zs <- sum(zi) / sqrt(k)
    ps <- 2 * stats::pnorm(-abs(zs))
  }
  data.frame(method = method, method_requested = method_req, tau2_method = tau2_method,
             pi_method = pi_method,
             note = note, k = k, total_n = NA_integer_, total_events = NA_integer_,
             HR = exp(bm), lower = exp(bm - crit * seM), upper = exp(bm + crit * seM),
             logHR = bm, se = seM, p = pv,
             Q = Q, df = df, p_heterogeneity = stats::pchisq(Q, df, lower.tail = FALSE),
             I2 = I2, tau2 = tau2,
             pi_lower = pi_lower, pi_upper = pi_upper,
             pi_alt_lower = pi_alt_lower, pi_alt_upper = pi_alt_upper,
             pi_rule = pi_rule, pi_alt_rule = pi_alt_rule,
             se_method = se_method,
             z_stouffer = zs, p_stouffer = ps)
}

# The prediction-interval construction rules (spec B6). "t" is the documented
# default inherited from the frozen 1.0.0 numbers; "normal" is the normal
# approximation; "HK" uses the Hartung-Knapp-Sidik-Jonkman adjusted SE of the
# pooled estimate with t(k-1). Case-insensitive.
.cpas_pi_methods <- c("t", "normal", "HK")

.cpas_pi_method <- function(x) {
  v <- as.character(x)[1]
  if (is.na(v)) v <- "t"
  hit <- .cpas_pi_methods[tolower(.cpas_pi_methods) == tolower(v)]
  if (!length(hit))
    stop("'pi_method' must be one of ", paste(.cpas_pi_methods, collapse = ", "),
         " (\"t\" is the default); got \"", x[1], "\".", call. = FALSE)
  hit[1]
}

# Human-readable primary / alternative rule labels for a pi_method.
.cpas_pi_rule_text <- function(pi_method) {
  t_lab <- "t distribution with k-2 df on sqrt(se_pooled^2 + tau2)"
  n_lab <- "normal approximation (1.96) on sqrt(se_pooled^2 + tau2)"
  h_lab <- paste0("Hartung-Knapp-Sidik-Jonkman adjusted SE of the pooled estimate ",
                  "with t(k-1) on sqrt(se_HK^2 + tau2)")
  switch(pi_method,
         t = list(primary = t_lab, alt = n_lab),
         normal = list(primary = n_lab, alt = t_lab),
         HK = list(primary = h_lab,
                   alt = paste0("t distribution with k-2 df on sqrt(se_pooled^2 + tau2) ",
                                "(unadjusted random-effects SE)")))
}

# One interval under a named rule.  Returns exp() bounds, NA when the rule needs
# more cohorts than were supplied (t needs k >= 3, the other two k >= 2), so a
# two-cohort pool yields NA for "t" and a finite normal / HK interval rather
# than a NaN.
.cpas_pi_bounds <- function(bm, rule, se_pi, se_pi_hk, k) {
  rule <- .cpas_pi_method(rule)
  if (identical(rule, "t")) {
    if (k < 3L) return(list(lower = NA_real_, upper = NA_real_, alt = "normal",
                            rule = "t"))
    crit <- stats::qt(0.975, df = k - 2)
    return(list(lower = exp(bm - crit * se_pi), upper = exp(bm + crit * se_pi),
                alt = "normal", rule = rule))
  }
  if (k < 2L) return(list(lower = NA_real_, upper = NA_real_,
                          alt = rule, rule = rule))
  if (identical(rule, "normal"))
    return(list(lower = exp(bm - 1.96 * se_pi), upper = exp(bm + 1.96 * se_pi),
                alt = "t", rule = rule))
  crit <- stats::qt(0.975, df = k - 1)
  list(lower = exp(bm - crit * se_pi_hk), upper = exp(bm + crit * se_pi_hk),
       alt = "t", rule = rule)
}


#' @title Print method for cpas_meta objects
#' @description Prints a concise summary of a \code{cpas_meta} analysis result.
#' @param x An object of class \code{cpas_meta}.
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_meta <- function(x, ...) {
  cat("CanPAS integrative (meta) analysis\n")
  cat("marker:", x$input$marker, "| type:", x$input$type,
      "| method:", x$input$method,
      "| pi_method:", if (is.null(x$input$pi_method)) "t" else x$input$pi_method,
      "| datasets:", nrow(x$per_dataset), "(failed:", length(x$errors), ")\n")
  if (!is.null(x$input$pooling))
    cat("pooling:", x$input$pooling, "| pooling classes pooled:",
        if (is.null(x$pooled$pooling_classes) || is.na(x$pooled$pooling_classes))
          "none" else x$pooled$pooling_classes,
        "| auto_repair:", isTRUE(x$input$auto_repair), "\n")
  ovp <- x$overlap_pairs
  if (is.null(ovp) && is.list(x$pooled$overlap_pairs)) ovp <- x$pooled$overlap_pairs[[1]]
  if (is.data.frame(ovp) && nrow(ovp))
    cat("shared-patient pairs in this pool (overlap = \"", x$pooled$overlap_mode,
        "\"): ",
        paste(sprintf("%s/%s (%s)", ovp$AccessionA, ovp$AccessionB,
                      ovp$SharedPatients), collapse = "; "),
        "\n", sep = "")
  if (!is.null(x$not_estimable) && nrow(x$not_estimable))
    cat("not estimable (fail-safe, auto_repair = FALSE): ",
        paste(sprintf("%s (%s)", x$not_estimable$cohort,
                      x$not_estimable$reason), collapse = "; "), "\n", sep = "")
  # the list column (the overlapping pairs themselves) is not printed as a
  # column; it is summarised in the lines above and available as $overlap_pairs
  po_show <- x$pooled[, !vapply(x$pooled, is.list, logical(1)), drop = FALSE]
  print(po_show, row.names = FALSE)
  if (is.data.frame(ovp) && nrow(ovp))
    cat("$pooled$overlap_pairs[[1]] and $overlap_pairs hold ", nrow(ovp),
        " shared-patient pair(s); $pooled$overlap_pairs_text holds the same as text.\n",
        sep = "")
  pm <- if (is.null(x$pooled$pi_method)) "t" else x$pooled$pi_method
  pi_lab <- switch(pm, t = "t k-2", normal = "normal", HK = "HK t(k-1)", pm)
  alt_lab <- switch(pm, t = "normal", normal = "t k-2", HK = "t k-2, unadjusted SE", "alternative")
  if (!is.null(x$pooled$pi_lower) && is.finite(x$pooled$pi_lower))
    cat(sprintf("95%% prediction interval for a new dataset (primary, %s): [%.3f, %.3f]\n",
                pi_lab, x$pooled$pi_lower, x$pooled$pi_upper))
  if (!is.null(x$pooled$pi_alt_lower) && is.finite(x$pooled$pi_alt_lower))
    cat(sprintf("95%% prediction interval (alternative, %s): [%.3f, %.3f]\n",
                alt_lab, x$pooled$pi_alt_lower, x$pooled$pi_alt_upper))
  cat("\nPer-dataset:\n"); print(x$per_dataset[, c("dataset","n","events","HR","lower","upper","p")], row.names = FALSE)
  invisible(x)
}

#' @title Forest plot of an integrative (meta) analysis
#' @description Draws the per-cohort HR (95\% CI) together with the pooled HR of
#' a \code{\link{cpas_meta}} result as a forest plot.
#' @param x An object of class \code{cpas_meta}.
#' @param digits Decimal places used for the annotations in the figure
#'   (default 4).
#' @param show_stars Whether to draw the per-cohort significance stars at the
#'   left edge of the panel (\code{*} p<0.05, \code{**} p<0.01,
#'   \code{***} p<0.001); default \code{TRUE}.
#' @param label_size Font size of the pooled HR / 95\% PI annotation (default 4.2).
#' @param star_size Font size of the stars and of the \code{Overall} label;
#'   defaults to \code{0.85 * label_size}.
#' @param y_headroom Head-room at the top of the panel so the \code{Overall}
#'   diamond and its confidence interval are not clipped (default 3.6).
#' @param x_frac Horizontal position, as a fraction from the left edge of the
#'   panel, used to anchor the stars and the pooled annotation (default 0.02).
#' @param label_lift Vertical offset of the pooled annotation relative to the
#'   \code{Overall} row (default 1).
#' @param label_where Where to place the pooled annotation: \code{"inside"}
#'   puts it in the top-left corner of the panel and \code{"subtitle"} above the
#'   panel (default \code{"inside"}).
#' @param ... Reserved, currently unused.
#' @return A \code{ggplot} object.
#' @details
#' The stars and the pooled annotation are anchored to a \strong{finite} data
#' value derived from the fixed axis expansion, not to \code{x = -Inf}: an
#' infinite value becomes \code{NaN} through the log10 scale, and the text layer
#' is then silently dropped with a warning. The per-cohort significance stars sit
#' at the left edge inside the panel, and the pooled HR / 95\% PI annotation is
#' placed in the top-left corner of the panel by default.
#' @examples
#' \dontrun{
#'    m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
#'    plot_meta_forest(m, digits = 4)
#' }
#' @export
plot_meta_forest <- function(x, digits = 4,
                             show_stars = TRUE,
                             label_size = 4.2,
                             star_size  = NULL,
                             y_headroom = 3.6,
                             x_frac = 0.02,
                             label_lift = 1.00,
                             label_where = c("inside", "subtitle"), ...) {
  label_where <- match.arg(label_where)
  pc <- x$per_dataset
  po <- x$pooled
  pc$se <- ifelse(pc$se <= 0, NA, pc$se)
  w <- ifelse(is.na(pc$se), 0, 1 / pc$se^2)
  pc$w <- w
  order_t <- pc$dataset
  pc$dataset <- factor(pc$dataset, levels = rev(order_t))
  k <- nrow(pc)
  y_overall <- k + 1                      # Overall 行所在的离散 y
  if (is.null(star_size)) star_size <- 0.85 * label_size
  PT <- 2.845                             # ggplot size -> pt
  ## 面板 x 边界（对数空间 5% 比例留白）：lo/hi 为面板左右缘对应的数据值
  lo_d <- min(pc$lower, na.rm = TRUE); hi_d <- max(pc$upper, na.rm = TRUE)
  Rlog <- log10(hi_d) - log10(lo_d)
  lo <- 10^(log10(lo_d) - 0.05 * Rlog); hi <- 10^(log10(hi_d) + 0.05 * Rlog)
  frac_x <- function(f) 10^(log10(lo) + f * (log10(hi) - log10(lo)))   # 面板内固定比例 -> 有限数据值
  x_anchor <- frac_x(x_frac)

  ## 合并标注文本（两行；内容与包内一致）
  pooled_lab <- if (is.finite(po$pi_lower))
    sprintf(paste0("Pooled HR %.", digits, "f [%.", digits, "f, %.", digits, "f] | I2=%.0f%%\n",
                   "95%% PI [%.", digits, "f, %.", digits, "f]"),
            po$HR, po$lower, po$upper, 100 * po$I2, po$pi_lower, po$pi_upper)
  else sprintf(paste0("Pooled HR %.", digits, "f [%.", digits, "f, %.", digits, "f], I2=%.0f%%"),
               po$HR, po$lower, po$upper, 100 * po$I2)

  ## ③ 星号文本（按行；x 由面板比例决定，与 HR 数值无关）
  st <- rep("", k)
  if (isTRUE(show_stars) && "p" %in% names(pc)) {
    pv <- suppressWarnings(as.numeric(pc$p))
    st <- ifelse(is.na(pv), "",
          ifelse(pv < 0.001, "***",
          ifelse(pv < 0.01,  "**",
          ifelse(pv < 0.05,  "*", ""))))
  }
  ## 行对齐：直接复用 pc$dataset 这个**同一个因子**作 y 美学，星号便与它自己的队列行
  ## 一同落在离散轴的同一位置上（v7 原用 y = seq_len(k) 的数值 1..k，而 levels = rev(order_t)
  ## 使数据第 i 行位于面板第 k-i+1 位，于是星号整列被上下镜像到错误队列上）。
  star_df <- data.frame(dataset = pc$dataset, st = st, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(pc, ggplot2::aes(x = HR, y = dataset)) +
    ggplot2::geom_point(size = 3 * sqrt(pc$w) / max(sqrt(pc$w)) + 1, color = "steelblue") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
                           orientation = "y", width = 0.25, color = "grey30") +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
    ## Overall 行（红菱形 + CI），靠顶部留白保证在框内
    ggplot2::geom_point(data = data.frame(HR = po$HR), ggplot2::aes(x = HR, y = y_overall),
                        shape = 18, size = 6, color = "firebrick", inherit.aes = FALSE) +
    ggplot2::geom_errorbar(data = data.frame(HR = po$HR, lower = po$lower, upper = po$upper),
                           ggplot2::aes(xmin = lower, xmax = upper, y = y_overall),
                           orientation = "y", width = 0.15, color = "firebrick", inherit.aes = FALSE) +
    ## x 轴只留右侧少量空白（左侧不再需要为文本让位）
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.6, y_headroom))) +
    ggplot2::labs(title = paste0("Meta forest: ", x$input$marker, " (", x$input$type, ", ",
                                 x$input$method, ")"),
                  x = "Hazard ratio (95% CI, log scale)", y = NULL) +
    ggplot2::theme_bw() +
    ## ③ 星号（x = 面板左缘起固定比例，y = 各队列行）
    ggplot2::geom_text(data = star_df, ggplot2::aes(x = x_anchor, y = dataset, label = st),
                       inherit.aes = FALSE, hjust = 0, size = star_size,
                       fontface = "bold", colour = "grey20") +
    ## Overall 文字（同一列）
    ggplot2::annotate("text", x = x_anchor, y = y_overall, label = "Overall",
                      hjust = 0, size = star_size, fontface = "bold", colour = "firebrick") +
    ## ② 合并标注：inside = 面板内（同一列、Overall 行之上）；subtitle = 面板上方
    { if (label_where == "inside")
        ggplot2::annotate("text", x = x_anchor, y = y_overall + label_lift, label = pooled_lab,
                          hjust = 0, vjust = 0, size = label_size, colour = "grey15") } +
    { if (label_where == "subtitle")
        ggplot2::labs(subtitle = pooled_lab) } +
    { if (label_where == "subtitle")
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = label_size,
                                                             hjust = 0, colour = "grey15",
                                                             margin = ggplot2::margin(b = 4))) }
  p
}

#' @title Sensitivity analysis: leave-one-out
#' @description Refits the pooled estimate once per cohort, each time leaving one
#' cohort out, so the influence of a single cohort on the overall result can be
#' assessed.
#' @param x An object of class \code{cpas_meta}.
#' @return data.frame: \code{left_out} names the cohort that was removed and the
#'   remaining columns are the re-pooled summary for that omission, computed with
#'   the same pooling method (\code{x$input$method}) as the original analysis.
#' @examples
#' \dontrun{
#'    m <- cpas_meta(c("GSE31210", "GSE37745", "GSE42127"), marker = "TP53", type = "DFS")
#'    loo_meta(m)
#' }
#' @export
loo_meta <- function(x) {
  pc <- x$per_dataset
  if (nrow(pc) < 2L)
    stop("loo_meta(): leave-one-out sensitivity analysis needs at least 2 datasets.")
  out <- lapply(seq_len(nrow(pc)), function(i) {
    d <- pc[-i, , drop = FALSE]
    # with 2 cohorts the remaining estimate is that cohort's own result (k = 1)
    po <- meta_pool(d$logHR, d$se, p = d$p, method = x$input$method,
                    pi_method = if (is.null(x$input$pi_method)) "t" else x$input$pi_method)
    data.frame(left_out = pc$dataset[i], po)
  })
  do.call(rbind, out)
}

#' @title Pooled Kaplan-Meier (marker split inside every cohort)
#' @description Splits each cohort into High/Low at the marker and then pools by
#' \code{method}:
#'   \code{"ipd"}  pools the patients directly into one \code{survfit} (the
#'     log-rank test is reported both unstratified and stratified by cohort);
#'   \code{"meta"} two-stage: at every time point the cohorts' KM survival
#'     probabilities S(t) are transformed with log(-log S) and pooled by inverse
#'     variance (RE/FE), giving a pooled survival curve and the 1/3/5-year table;
#'   \code{"both"} returns both routes.
#' @param merged Named list of the cohorts' \code{merge_surv_expr} merged
#'   data.frames. The names must be catalog accessions so the concrete endpoint
#'   token of each cohort can be resolved. The result carries
#'   \code{dataset_endpoints}, the token every dataset actually used.
#' @param marker Gene column name.
#' @param type Endpoint family or token (OS/RFS/...).
#' @param method \code{"ipd"}, \code{"meta"} or \code{"both"}.
#' @param landmarks Time points, in years, at which the pooled survival
#'   probabilities are tabulated and drawn (default 1, 3 and 5 years).
#' @param meta_method Pooling method for the time-point survival probabilities:
#'   \code{"RE"}/\code{"DL"} (DerSimonian-Laird random effects, the default kept
#'   for continuity with the published curves), \code{"REML"}, \code{"HK"} or
#'   \code{"FE"} (fixed effect).
#' @param cut High/low split rule, one of three: \code{"median"} (default) takes
#'   the top 50\% inside every cohort, i.e. High = marker above \strong{that
#'   cohort's own} median; \code{"top_pct"} sorts each cohort by expression and
#'   takes the highest \code{top_pct}\% as High (threshold = that cohort's
#'   \eqn{100 - top_pct} percentile, so \code{top_pct = 25} means the top 25\%
#'   are high expressors); \code{"custom"} uses the absolute threshold
#'   \code{cut_value}, with High above it and Low at or below it. An absolute
#'   threshold requires comparable measurement and normalization across cohorts;
#'   a log2 transformation alone does not establish that comparability.
#'   Percentile splits are nevertheless the
#'   conventional reading of "high versus low expression", and the application
#'   offers only median and top_pct. A per-cohort search for the best cut point
#'   is deliberately \strong{not} offered: it would repeat, cohort by cohort,
#'   the cut-point search quantified in Section 3.4 and compound its inflation
#'   when pooled.
#' @param top_pct Percentage for \code{cut = "top_pct"}: a single finite value
#'   between 1 and 99; the patients with the highest \code{top_pct}\% of
#'   expression are High. The threshold is computed with the default
#'   \code{stats::quantile()} interpolation (type 7); when several patients sit
#'   exactly on the threshold the High group can be slightly larger than
#'   \code{top_pct}\%, and the actual threshold and counts are returned
#'   (\code{cohort_thresholds}, \code{cutpoint}). A cohort left one-sided by
#'   the rule is recorded in \code{empty_cohorts} with its reason and does not
#'   enter the pool.
#' @param cut_value Threshold for \code{cut = "custom"}: a single finite value on
#'   the log2 expression scale (or on the signature score). A cohort with no
#'   patient on one side of the threshold is recorded in \code{empty_cohorts}
#'   with its reason and does not enter the pool.
#' @return An ordinary \code{list} (no class attribute; not an S3 object, so
#'   there is no method dispatch):
#'   \item{\code{df}:}{the pooled analysis data (time/status/marker/dataset/group)}
#'   \item{\code{datasets}, \code{n_high}, \code{n_low}, \code{method}, \code{cutpoint}:}{the cohorts included and the size of each group}
#'   \item{\code{dataset_endpoints}:}{the endpoint token every dataset actually used (named vector)}
#'   \item{\code{cut}, \code{top_pct}, \code{cut_value}, \code{cutpoint}:}{the split
#'     rule actually used, its parameter or threshold, and a printable description
#'     of the rule; \code{cohort_thresholds} gives the per-cohort threshold when
#'     \code{cut = "top_pct"}}
#'   \item{\code{n_dropped}, \code{empty_cohorts}, \code{empty_reasons}:}{the split
#'     \code{n_dropped} counts incomplete rows removed from included cohorts;
#'     when any split rule leaves one side empty in a cohort,
#'     that cohort is recorded in \code{empty_cohorts} with its reason and does not
#'     enter the pool}
#'   \item{\code{skipped_cohorts}, \code{skipped_reasons}:}{cohorts excluded for
#'     insufficient data and why - the endpoint could not be resolved, a column is
#'     missing, fewer than 10 complete (time, status, marker) rows, or the marker
#'     takes a single value in that cohort. These cohorts do not enter the pool but
#'     are \strong{never silently dropped}}
#'   \item{\code{fit}, \code{logrank_p}, \code{logrank_p_stratified}:}{the
#'     \code{survfit} and the log-rank p (pooled / stratified by cohort) when
#'     \code{method} is ipd or both}
#'   \item{\code{meta_landmarks}, \code{meta_curve}:}{the time-point pooled
#'     survival table and the fine-grid curve when \code{method} is meta or both}
#'   \item{\code{manifest}:}{the analysis manifest
#'     (\code{\link{cpas_manifest}}): the cut-point rule, the number of cut-points
#'     searched (0: no search), the cohorts skipped and why, the versions and the
#'     timestamp}
#' @details
#' Each cohort is split on its own marker median (High = marker > that cohort's
#' median), so what is compared is the risk above versus below each cohort's own
#' median and not an absolute threshold across cohorts; when an absolute
#' threshold is wanted, use \code{plot_km()} cohort by cohort. The IPD route
#' reports the log-rank p both unstratified and stratified by cohort; reporting
#' the stratified one is recommended. In the time-point route, a time point
#' beyond a cohort's longest follow-up excludes that cohort at that time
#' (\code{extend = FALSE}), and the number of cohorts actually contributing at
#' each time point is recorded in \code{meta_landmarks$k}.
#' @examples
#' \dontrun{
#'    ## Pooled Kaplan-Meier from an explicitly built list of merged cohorts
#'    ## (names must be catalog accessions so the endpoint is resolved per cohort).
#'    merged <- list(GSE31210 = cohort_merged("GSE31210", "GAPDH", type = "RFS"),
#'                    GSE37745 = cohort_merged("GSE37745", "GAPDH", type = "RFS"))
#'    km <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS", method = "both")
#'    km$datasets; km$dataset_endpoints; km$logrank_p
#'
#'    plot_cpas_km(km)                       # pooled curve + landmark table
#'    plot_cpas_km_perdataset(km, ncol = 2)  # one panel per cohort
#' }
#' @export
cpas_km_pooled <- function(merged, marker, type = "OS",
                           method = c("both", "ipd", "meta"),
                           landmarks = c(1, 3, 5),
                           meta_method = "RE",
                           cut = c("median", "top_pct", "custom"),
                           top_pct = NULL,
                           cut_value = NULL) {
  method <- match.arg(method)
  cut <- match.arg(cut)
  empty <- list()                      # cohorts a threshold left one-sided
  skipped <- list()                    # cohorts dropped for insufficient/constant data
  med <- list()                        # per-cohort marker median (reported on failure)
  thr <- list()                        # per-cohort top-percent threshold
  dropped_n <- stats::setNames(integer(length(merged)), names(merged))
  if (identical(cut, "top_pct") &&
      (is.null(top_pct) || length(top_pct) != 1L || !is.finite(top_pct) ||
       top_pct <= 0 || top_pct >= 100))
    stop("cut = \"top_pct\" needs a single finite 'top_pct' between 1 and 99 ",
         "(the highest top_pct% of each cohort is High; 25 means the top 25%).",
         call. = FALSE)
  if (identical(cut, "custom") &&
      (is.null(cut_value) || length(cut_value) != 1L || !is.finite(cut_value)))
    stop("cut = \"custom\" needs a single finite 'cut_value' (a threshold on ",
         "the log2 expression scale, or on the signature score).", call. = FALSE)
  ep_used <- character(0)
  parts <- lapply(names(merged), function(t) {
    d <- merged[[t]]
    tok <- endpoint_resolve(t, type)
    if (is.na(tok)) {
      if (all(c(paste0(type, "_time"), paste0(type, "_status")) %in% colnames(d)))
        tok <- as.character(type)[1]
      else {
        skipped[[t]] <<- sprintf("no %s endpoint for this cohort and no %s_time/%s_status columns",
                                 type, type, type)
        return(NULL)
      }
    }
    tc <- paste0(tok, "_time"); sc <- paste0(tok, "_status")
    if (!all(c(tc, sc, marker) %in% colnames(d))) {
      skipped[[t]] <<- sprintf("column(s) absent: %s",
                               paste(setdiff(c(tc, sc, marker), colnames(d)), collapse = ", "))
      return(NULL)
    }
    ep_used[[t]] <<- tok
    st <- suppressWarnings(as.numeric(d[[sc]]))
    tm <- suppressWarnings(as.numeric(d[[tc]]))
    if (any(!is.na(st) & !st %in% c(0, 1)))
      stop("Cohort ", t, ": status column ", sc,
           " must be coded 0 (censored) / 1 (event).", call. = FALSE)
    if (any(!is.na(tm) & tm < 0))
      stop("Cohort ", t, ": negative survival time in ", tc, ".", call. = FALSE)
    v <- suppressWarnings(as.numeric(d[[marker]]))
    keep <- !is.na(v) & !is.na(suppressWarnings(as.numeric(d[[tc]]))) &
      !is.na(suppressWarnings(as.numeric(d[[sc]])))
    if (sum(keep) < 10) {
      skipped[[t]] <<- sprintf("only %d rows with a complete (time, status, marker) triple (minimum 10)",
                               sum(keep))
      return(NULL)
    }
    if (length(unique(v[keep])) < 2) {
      skipped[[t]] <<- sprintf("the marker takes a single value (%g) in this cohort, so no split is possible",
                               unique(v[keep])[1])
      return(NULL)
    }
    dropped_n[[t]] <<- nrow(d) - sum(keep)
    dd <- data.frame(time = suppressWarnings(as.numeric(d[[tc]][keep])),
                     status = suppressWarnings(as.numeric(d[[sc]][keep])),
                     marker = v[keep], dataset = t)
    ## Three rules:
    ##   "median"   the top 50% of THIS cohort (High = marker > cohort median);
    ##   "top_pct"  the top top_pct% of THIS cohort by expression, i.e. High =
    ##              marker > this cohort's (100 - top_pct)th percentile;
    ##   "custom"   an absolute threshold: High = marker > cut_value, Low otherwise.
    ## Absolute thresholds require a harmonized measurement scale; log2 alone
    ## does not make thresholds comparable across platforms. A per-cohort search
    ## for the best cut point is deliberately not offered: that search inflates
    ## the p-value (quantified in Section 3.4) and pooling it over k cohorts
    ## compounds the inflation.
    med[[t]] <<- stats::median(dd$marker)
    if (identical(cut, "top_pct")) {
      thr <- stats::quantile(dd$marker, probs = 1 - top_pct / 100, names = FALSE)
      thr[[t]] <<- thr
      dd$group <- ifelse(dd$marker > thr, "High", "Low")
      if (length(unique(dd$group)) < 2L) {
        ## the top-percent rule can still leave one side empty when every value
        ## in the cohort is identical (or the ties span the whole cohort)
        empty[[t]] <<- sprintf("a top-%g%% split left no patient %s the cut (%g) (n = %d)",
                               top_pct,
                               if (all(dd$marker > thr)) "at or below" else "above",
                               thr, nrow(dd))
        return(NULL)
      }
    } else if (identical(cut, "custom")) {
      dd$group <- ifelse(dd$marker > cut_value, "High", "Low")
      ## one side missing entirely (min(table()) < 1 never fires: it means "no
      ## rows at all", not "no rows on that side")
      if (length(unique(dd$group)) < 2L) {
        ## the threshold leaves one side empty in THIS cohort: record it with the
        ## reason and leave the cohort out of the pool instead of failing later
        empty[[t]] <<- sprintf("no patient %s the threshold %g (n = %d)",
                               if (all(dd$marker > cut_value)) "at or below" else "above",
                               cut_value, nrow(dd))
        return(NULL)
      }
    } else {
      med_cut <- stats::median(dd$marker)
      dd$group <- ifelse(dd$marker > med_cut, "High", "Low")
      if (length(unique(dd$group)) < 2L) {
        empty[[t]] <<- sprintf("the median split at %g left one group empty (n = %d)",
                               med_cut, nrow(dd))
        return(NULL)
      }
    }
    dd
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    if (length(empty)) {
      rule_desc <- switch(cut,
        top_pct = sprintf("A top-%g%% split", top_pct),
        custom  = sprintf("A custom cut of %g", cut_value),
        sprintf("Cut rule \"%s\"", cut))
      stop(rule_desc, " left every cohort one-sided (",
           paste(sprintf("%s: %s", names(empty), unlist(empty)), collapse = "; "),
           "). Pick a threshold inside the marker's observed range; the cohort ",
           "medians were ", paste(sprintf("%s = %.3g", names(med), unlist(med)),
                                  collapse = ", "), ".",
           if (length(skipped))
             paste0(" Cohorts skipped for insufficient or constant data: ",
                    paste(sprintf("%s (%s)", names(skipped), unlist(skipped)),
                          collapse = "; "), ".")
           else "", call. = FALSE)
    }
    stop("No cohort provided usable (time, status, marker) data for endpoint ",
         type, ".", call. = FALSE)
  }
  df <- do.call(rbind, parts)
  n_dropped <- sum(dropped_n[names(merged) %in% unique(df$dataset)])
  df$group <- factor(df$group, levels = c("Low", "High"))
  df <- df[stats::complete.cases(df[c("time", "status", "marker", "group")]), , drop = FALSE]
  if (!sum(df$status == 1))
    stop("No event was observed in any cohort for endpoint ", type,
         "; a log-rank test / pooled KM cannot be estimated.", call. = FALSE)
  n_hi <- sum(df$group == "High"); n_lo <- sum(df$group == "Low")
  if (min(n_hi, n_lo) < 1L)
    stop("Median split produced an empty group.")
  cut_label <- switch(cut,
    median = "50% split within each cohort: High = marker > that cohort's median",
    top_pct = sprintf(paste0("top %g%% of each cohort by expression: ",
                             "High = marker > that cohort's %gth percentile"),
                      top_pct, 100 - top_pct),
    custom = sprintf(paste0("custom absolute threshold: High = marker > %g, ",
                            "Low = marker <= %g (log2 scale, the same value in every cohort)"),
                     cut_value, cut_value))
  out <- list(df = df, datasets = unique(df$dataset), method = method,
              n_high = n_hi, n_low = n_lo,
              dataset_endpoints = ep_used, cut = cut, top_pct = top_pct,
              cut_value = cut_value,
              n_dropped = n_dropped, cohort_medians = unlist(med),
              cohort_thresholds = unlist(thr),
              empty_cohorts = names(empty),
              empty_reasons = unlist(empty),
              skipped_cohorts = names(skipped),
              skipped_reasons = unlist(skipped), cutpoint = cut_label)
  if (method %in% c("ipd", "both")) {
    fit <- survival::survfit(survival::Surv(time, status) ~ group, data = df)
    lr  <- survival::survdiff(survival::Surv(time, status) ~ group, data = df)
    p   <- 1 - stats::pchisq(lr$chisq, length(lr$n) - 1)
    lrs <- tryCatch(survival::survdiff(survival::Surv(time, status) ~ group + strata(dataset), data = df),
                    error = function(e) NULL)
    p_strat <- if (!is.null(lrs)) 1 - stats::pchisq(lrs$chisq, length(lrs$n) - 1) else NA_real_
    out$fit <- fit; out$logrank_p <- p; out$logrank_p_stratified <- p_strat
  }
  if (method %in% c("meta", "both")) {
    ms <- km_surv_meta(df, times = landmarks, meta_method = meta_method)
    out$meta_landmarks <- ms$landmarks
    out$meta_curve <- ms$curve
  }
  out$manifest <- .cpas_manifest_new(
    analysis = "cpas_km_pooled",
    cohorts = unique(df$dataset), accession = unique(df$dataset),
    family = endpoint_family(type), raw_token = unname(ep_used),
    n_input = sum(vapply(merged, nrow, integer(1))), n_analyzed = nrow(df),
    events = sum(df$status == 1), marker_requested = marker, marker_definition = marker,
    fitted_covariates = "group", estimator = paste0("survival::survfit; ", method),
    inference = paste0("Greenwood KM SE; log-rank; landmark ", meta_method),
    analyzed_data = df, coverage_status = "measured pooled input and analysis rows",
    tokens = ep_used,
    token_role = "resolved per cohort",
    selection_rule = sprintf(paste0("the %d cohort(s) in the supplied 'merged' list that provided a usable ",
                                    "(time, status, marker) triple for endpoint family '%s'; skipped: %s"),
                             length(unique(df$dataset)), type,
                             if (length(skipped)) paste(sprintf("%s (%s)", names(skipped),
                                                               unlist(skipped)), collapse = "; ") else "none"),
    dropped_rows = data.frame(
      cohort = names(dropped_n)[names(dropped_n) %in% unique(df$dataset)],
      n_dropped = as.integer(dropped_n[names(dropped_n) %in% unique(df$dataset)]),
      reason = "incomplete time/status/marker rows within included cohorts",
      stringsAsFactors = FALSE),
    cut_rule = cut_label,
    cut_points_searched = 0L,
    search_adjusted_p = NA_real_,
    meta_method = meta_method,
    tau2 = if (!is.null(out$meta_landmarks) && "tau2" %in% colnames(out$meta_landmarks))
      stats::median(out$meta_landmarks$tau2, na.rm = TRUE) else NA_real_,
    I2 = if (!is.null(out$meta_landmarks)) stats::median(out$meta_landmarks$I2, na.rm = TRUE) else NA_real_,
    notes = c("no cut-point is searched: the split is the requested rule applied inside every cohort, so the reported p-value is not search-adjusted",
               paste0("pooled KM is an IPD descriptive mixture or independent landmark meta-estimates; it is not a standardized clinical survival curve"),
              if (identical(cut, "top_pct")) "top-percent split" else NULL,
              if (length(empty)) paste0("cohort(s) left one-sided by the rule and not pooled: ",
                                        paste(names(empty), collapse = ", ")) else NULL))
  invisible(out)
}

# Time-point survival meta-analysis on log(-log) survival probabilities. The
# returned curve is a set of independent landmark estimates, not a guaranteed
# monotone Kaplan-Meier curve for a target population.
.cpas_km_meta_critical <- function(meta_method, k) {
  m <- toupper(as.character(meta_method)[1])
  if (identical(m, "HK") && k >= 2L) stats::qt(.975, df = k - 1L) else 1.96
}
km_surv_meta <- function(df, times, meta_method = "RE") {
  grp <- sort(unique(df$group))
  res_l <- lapply(grp, function(g) {
    dg <- df[df$group == g, ]
    cg <- unique(dg$dataset)
    fits <- lapply(cg, function(cf) {
      dd <- dg[dg$dataset == cf, ]
      s <- survival::survfit(survival::Surv(time, status) ~ 1, data = dd)
      sm <- summary(s, times = times, extend = FALSE)
      data.frame(dataset = rep(cf, length(sm$time)), time = as.numeric(sm$time),
                 S = as.numeric(sm$surv), SE = as.numeric(sm$std.err))
    })
    m <- do.call(rbind, fits); m <- m[is.finite(m$SE) & m$SE > 0 & m$S > 0 & m$S < 1, ]
    theta <- log(-log(m$S)); se_t <- m$SE / (m$S * abs(log(m$S)))
    tab <- lapply(unique(m$time), function(tm) {
      mm <- m[m$time == tm, ]
      if (nrow(mm) < 1) return(NULL)
      po <- meta_pool(theta[m$time == tm], se_t[m$time == tm], method = meta_method)
      data.frame(group = g, time = tm, k = nrow(mm),
                 S = exp(-exp(po$logHR)),
                 lower = exp(-exp(po$logHR + .cpas_km_meta_critical(meta_method, nrow(mm)) * po$se)),
                 upper = exp(-exp(po$logHR - .cpas_km_meta_critical(meta_method, nrow(mm)) * po$se)),
                 I2 = po$I2, p_het = po$p_heterogeneity)
    })
    do.call(rbind, tab)
  })
  land <- do.call(rbind, res_l)
  # 细网格曲线（每 0.5 年一步）
  grid_t <- seq(0.5, max(times) + 0.5, 0.5)
  res_c <- lapply(grp, function(g) {
    dg <- df[df$group == g, ]
    fits <- lapply(unique(dg$dataset), function(cf) {
      dd <- dg[dg$dataset == cf, ]
      s <- survival::survfit(survival::Surv(time, status) ~ 1, data = dd)
      sm <- summary(s, times = grid_t, extend = FALSE)
      data.frame(dataset = rep(cf, length(sm$time)), time = as.numeric(sm$time),
                 S = as.numeric(sm$surv), SE = as.numeric(sm$std.err))
    })
    m <- do.call(rbind, fits)
    m <- m[is.finite(m$SE) & m$SE > 0 & m$S > 0 & m$S < 1, ]
    tab <- lapply(unique(m$time), function(tm) {
      mm <- m[m$time == tm, ]; if (nrow(mm) < 1) return(NULL)
      po <- meta_pool(log(-log(mm$S)), mm$SE / (mm$S * abs(log(mm$S))), method = meta_method)
      data.frame(group = g, time = tm, S = exp(-exp(po$logHR)),
                 lower = exp(-exp(po$logHR + .cpas_km_meta_critical(meta_method, nrow(mm)) * po$se)),
                 upper = exp(-exp(po$logHR - .cpas_km_meta_critical(meta_method, nrow(mm)) * po$se)))
    })
    do.call(rbind, tab)
  })
  list(landmarks = land, curve = do.call(rbind, res_c))
}

#' @title Grid of per-cohort Kaplan-Meier panels
#' @description Splits a \code{\link{cpas_km_pooled}} result back into one
#' Kaplan-Meier panel per dataset (each cohort split at its own median) and lays
#' them out as a grid.
#' @param km An object returned by \code{\link{cpas_km_pooled}}.
#' @param ncol Number of columns in the grid.
#' @param draw \code{TRUE} (default) draws the grid on the current device with
#'   \code{gridExtra::grid.arrange}; \code{FALSE} returns a composable
#'   \code{patchwork} object instead, so the caller can combine it with other
#'   figures (passing the return value of \code{grid.arrange} to patchwork does
#'   not error but is silently dropped, which is why the composed case needs
#'   \code{draw = FALSE}).
#' @return Invisibly \code{NULL} when \code{draw = TRUE} (the figure has been
#'   drawn); a \code{patchwork} object when \code{draw = FALSE}.
#' @export
plot_cpas_km_perdataset <- function(km, ncol = 4, draw = TRUE) {
  df <- km$df
  plist <- lapply(sort(unique(df$dataset)), function(cf) {
    d <- df[df$dataset == cf, ]
    f <- tryCatch(survival::survfit(survival::Surv(time, status) ~ group, data = d),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    sp <- survminer::ggsurvplot(f, data = d, pval = TRUE, legend = "none",
                                title = cf, xlab = "Time (years)", conf.int = FALSE,
                                palette = unname(.cpas_group_colours()))
    sp$plot
  })
  plist <- plist[!vapply(plist, is.null, logical(1))]
  if (isTRUE(draw)) return(do.call(gridExtra::grid.arrange, c(plist, list(ncol = ncol))))
  patchwork::wrap_plots(plist, ncol = ncol)
}

#' @title Plot pooled / meta Kaplan-Meier result
#' @description Visualizes a \code{\link{cpas_km_pooled}} result: when
#' method \code{"ipd"}/\code{"both"} was used a classical stratified
#' Kaplan-Meier plot is drawn; otherwise (method \code{"meta"}) the merged
#' meta survival curve with confidence band is plotted.
#' @param km An object returned by \code{\link{cpas_km_pooled}}.
#' @param which Which route to draw: \code{"auto"} (default) draws the
#'   individual-patient-data curve when that route was run and the time-point meta
#'   curve otherwise, exactly as before; \code{"ipd"} or \code{"meta"} selects a
#'   route explicitly. A \code{method = "both"} result holds both, so a caller that
#'   wants to show them together (as the Shiny application does) can draw each in
#'   turn instead of only ever getting the IPD curve.
#' @return A \code{ggsurvplot} object for the IPD route (with its risk table) or a
#' ggplot object for the meta route, carrying the requested landmark estimates and
#' their confidence intervals as points on the curve.
#' @export
plot_cpas_km <- function(km, which = c("auto", "ipd", "meta")) {
  which <- match.arg(which)
  has_ipd <- !is.null(km$fit)
  has_meta <- !is.null(km$meta_curve)
  if (which == "auto") which <- if (has_ipd) "ipd" else "meta"
  if (which == "ipd") {
    if (!has_ipd)
      stop("This result has no individual-patient-data route (method = \"meta\"). ",
           "Run with method = \"ipd\" or \"both\".", call. = FALSE)
    return(survminer::ggsurvplot(km$fit, data = km$df, pval = TRUE, risk.table = TRUE,
                                 legend.title = "marker",
                                 legend.labs = levels(km$df$group),
                                 palette = unname(.cpas_group_colours()),
                                 xlab = "Time (years)", conf.int = TRUE))
  }
  if (!has_meta)
    stop("This result has no time-point meta route (method = \"ipd\"). ",
         "Run with method = \"meta\" or \"both\".", call. = FALSE)
  ## the meta curve, with the requested landmark estimates drawn on it so the
  ## 1/3/5-year values that the landmark table reports are visible on the figure
  d <- km$meta_curve
  p <- ggplot(d, aes(x = time, y = S, color = group, fill = group)) +
    geom_step(linewidth = 0.7) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    scale_colour_manual(values = .cpas_group_colours()) +
    scale_fill_manual(values = .cpas_group_colours()) +
    scale_y_continuous(limits = c(0, 1)) + theme_bw() +
    labs(x = "Time (years)", y = "Survival probability",
         title = "Time-point meta curve (S(t) pooled across cohorts)")
  lm <- km$meta_landmarks
  if (!is.null(lm) && nrow(lm)) {
    lm <- lm[stats::complete.cases(lm[, c("time", "S", "lower", "upper")]), , drop = FALSE]
    if (nrow(lm)) {
      p <- p +
        geom_point(data = lm, aes(x = time, y = S, colour = group),
                   size = 2.1, inherit.aes = FALSE) +
        geom_errorbar(data = lm, aes(x = time, ymin = lower, ymax = upper, colour = group),
                      width = 0.12, linewidth = 0.45, inherit.aes = FALSE) +
        geom_text(data = lm, aes(x = time, y = lower, colour = group,
                                label = sprintf("%s\n%.3f", group, S)),
                  vjust = 1.5, size = 2.6, show.legend = FALSE, inherit.aes = FALSE)
    }
  }
  p
}
