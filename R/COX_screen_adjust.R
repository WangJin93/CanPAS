#' @title Univariate -> Multivariate Cox Analysis Workflow
#' @description
#' Runs univariate Cox models for all supplied covariates, selects covariates
#' significant at \code{p.threshold} (main effect of a categorical variable =
#' any of its levels), then - if at least two covariates qualify - runs one
#' multivariable Cox model and returns both univariate and multivariate
#' estimates side by side, optionally as a formatted \code{flextable}.
#' @param df Data.frame (see \code{\link{COX_analysis}} for the required
#' layout, including \code{type}).
#' @param type Survival endpoint prefix, default \code{"OS"}.
#' @param cont_Variates Continuous covariate column names.
#' @param cate_Variates Categorical covariate column names.
#' @param precision Decimal places for formatted values.
#' @param p.threshold Significance threshold for entering the multivariate step.
#' @return List with:
#'   \item{\code{result}:}{long data.frame with univariate and (if available)
#'     multivariate HR/CI/p per covariate and level}
#'   \item{\code{print_result}:}{a \code{flextable} of \code{result}, formatted
#'     as a three-line table (top rule, rule under the header, bottom rule) with
#'     the univariable and multivariable columns grouped under their own header
#'     row, ready for \code{flextable::save_as_docx()} or printing}
#'   \item{\code{sig_variates}:}{covariates that entered the multivariate step}
#'   \item{\code{uni_table}, \code{multi_table}:}{raw tables of each step}
#'   \item{\code{multi_metadata}:}{metadata of the multivariate step, including
#'     \code{reduced} and \code{reduced_note} when the model had to be reduced
#'     below \code{min_covariates} (see \code{\link{COX_analysis}}); a message
#'     states the same}
#'   \item{\code{cont_Variates}, \code{cate_Variates}:}{echoed inputs}
#' @export
#' @examples
#' \dontrun{
#'    ## Real cohort with real covariates: expression + survival + clinical
#'    ## columns are fetched and merged by the prerequisite functions.
#'    d <- cohort_merged("GSE13507", "GAPDH", type = "OS", clin = TRUE)
#' 
#'    q <- COX_screen_adjust(d, type = "OS", cont_Variates = c("GAPDH", "age"),
#'                            cate_Variates = c("grade", "N"), p.threshold = 0.05)
#'    q$sig_variates
#'    head(q$result)
#'    ## multi_metadata$reduced is TRUE when the joint model had to be reduced
#'    q$multi_metadata$reduced_note
#' }
COX_screen_adjust <- function(df, type = "OS",
                                cont_Variates = NULL,
                                cate_Variates = NULL,
                                precision = 3,
                                p.threshold = 0.05) {
  if (!is.numeric(p.threshold) || p.threshold <= 0 || p.threshold >= 1)
    stop("'p.threshold' must lie in (0, 1).")
  cont_Variates <- if (is.null(cont_Variates)) character(0) else as.character(cont_Variates)
  cate_Variates <- if (is.null(cate_Variates)) character(0) else as.character(cate_Variates)

  uni <- COX_analysis(df, type = type, cont_Variates = cont_Variates,
                      cate_Variates = cate_Variates, method = "uni",
                      precision = precision)
  u <- uni$results_table
  u$seq <- seq_len(nrow(u))

  # a covariate is "significant" if the minimal p over its rows is < threshold
  sig_rows <- !is.na(u$Pvalue) & u$Pvalue < p.threshold
  sig_vars <- unique(u$Var1[sig_rows])
  if (!length(sig_vars)) {
    message("No covariate reached p < ", p.threshold,
            "; only univariate results are returned.")
    multi <- NULL
    out <- u
  } else if (length(sig_vars) < 2L) {
    message("Only one covariate is significant (", sig_vars,
            "); a multivariate model needs at least two covariates. ",
            "Only univariate results are returned.")
    multi <- NULL
    out <- u
  } else {
    multi <- COX_analysis(df, type = type,
                          cont_Variates = intersect(sig_vars, cont_Variates),
                          cate_Variates = intersect(sig_vars, cate_Variates),
                          method = "multi", precision = precision)
    m <- multi$results_table
    key_cols <- c("Var1", "Level")
    mm <- merge(m[, c(key_cols, "HR", "HR95L", "HR95H", "Pvalue")],
                u[, c(key_cols, "HR", "HR95L", "HR95H", "Pvalue")],
                by = key_cols, suffixes = c("_multi", "_uni"), all = FALSE)
    out <- merge(u, mm[, c(key_cols, grep("_multi$", colnames(mm), value = TRUE))],
                 by = key_cols, all.x = FALSE)
    out <- out[order(out$seq), , drop = FALSE]
    # a model the repair had to reduce is no longer an adjusted model: say so
    # here as well, since this function is exported on its own
    if (isTRUE(multi$metadata$reduced))
      message("The multivariate model could only keep ",
              length(multi$metadata$final_covariates), " covariate(s) (",
              paste(multi$metadata$final_covariates, collapse = ", "),
              ") and is not adjusted for confounding: ",
              paste(multi$metadata$reduced_note, collapse = " "))
  }
  if ("seq" %in% colnames(out)) out$seq <- NULL

  # ---- formatted flextable ----
  fmt3 <- function(x) ifelse(is.na(x), "",
                             sprintf(paste0("%.", precision, "f"), x))
  display <- data.frame(
    Variates = ifelse(is.na(out$Level) | out$Level == "" | out$Type == "cate_header",
                      out$Variates, paste0("  ", out$Level)),
    N = out$N,
    stringsAsFactors = FALSE)
  uni_cols <- c("Univariate.HR", "Univariate.95CI", "Univariate.P")
  display$Univariate.HR  <- fmt3(out$HR)
  display$Univariate.95CI <- paste0(fmt3(out$HR95L), "-", fmt3(out$HR95H))
  display$Univariate.P   <- ifelse(is.na(out$Pvalue), "",
                                   ifelse(out$Pvalue < 10^(-precision),
                                          paste0("<", 10^(-precision)),
                                          fmt3(out$Pvalue)))
  if (!is.null(multi)) {
    display$Multivariate.HR  <- fmt3(out$HR_multi)
    display$Multivariate.95CI <- paste0(fmt3(out$HR95L_multi), "-",
                                        fmt3(out$HR95H_multi))
    display$Multivariate.P   <- ifelse(is.na(out$Pvalue_multi), "",
                                       ifelse(out$Pvalue_multi < 10^(-precision),
                                              paste0("<", 10^(-precision)),
                                              fmt3(out$Pvalue_multi)))
  }
  key <- colnames(display)
  what <- c("Variates", "N",
            rep("Univariate Cox", 3), if (!is.null(multi)) rep("Multivariate Cox", 3) else NULL)
  measure <- c("Variates", "N", "HR", "95% CI", "P-value",
               if (!is.null(multi)) c("HR", "95% CI", "P-value") else NULL)
  typology <- data.frame(key = key, what = what, measure = measure,
                         stringsAsFactors = FALSE)

  ft <- flextable::flextable(display)
  ft <- flextable::set_header_df(ft, mapping = typology, key = "key")
  ft <- flextable::merge_h(ft, part = "header")
  ft <- flextable::merge_v(ft, j = "Variates", part = "header")
  ft <- flextable::align(ft, align = "center", part = "all")
  ft <- flextable::autofit(ft)
  # three-line table (三线表): a top rule, a rule under the header block and a
  # bottom rule, and nothing else. The flextable default left every border at
  # '0 solid', i.e. no rules at all, which is not a printable table.
  ft <- flextable::border_remove(ft)
  ft <- flextable::hline_top(ft, part = "header",
                             border = flextable::fp_border_default(width = 1.5, color = "black"))
  ft <- flextable::hline_bottom(ft, part = "header",
                                border = flextable::fp_border_default(width = 0.75, color = "black"))
  ft <- flextable::hline_bottom(ft, part = "body",
                                border = flextable::fp_border_default(width = 1.5, color = "black"))
  ft <- flextable::align(ft, j = "Variates", align = "left", part = "all")
  ft <- flextable::padding(ft, padding = 4, part = "all")

  list(result = out,
       print_result = ft,
       sig_variates = sig_vars,
       uni_table = u,
       multi_table = if (!is.null(multi)) multi$results_table else NULL,
       multi_metadata = if (!is.null(multi)) multi$metadata else NULL,
       cont_Variates = cont_Variates,
       cate_Variates = cate_Variates)
}
