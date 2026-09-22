#' @title Cox Regression for Many Genes (Panel Screening)
#' @description
#' Fits one univariable Cox proportional-hazards model per gene
#' (\code{Surv(time, status) ~ gene}) and returns hazard ratios, 95\%
#' confidence intervals, Wald p-values and Benjamini-Hochberg FDRs for every
#' gene, together with the fitted model objects. This is a panel screen: each
#' gene is estimated separately and the genes are never adjusted for one
#' another. Use \code{\link{COX_analysis}(method = "multi")} for a joint model
#' with adjusted hazard ratios, and \code{\link{COX_by_datasets}} for the same
#' analysis across several datasets.
#' @param df Data.frame whose first column is \code{ID} and that contains the
#' columns \code{<type>_time} (years), \code{<type>_status} (0/1) and one
#' numeric column per marker.
#' @param type Endpoint family (\code{"OS"}, \code{"DSS"}, \code{"DFS"},
#' \code{"PFS"}, \code{"MFS"}) or a concrete token (\code{"RFS"}, \code{"DFI"},
#' \code{"PFI"} ...). Families are resolved to the token available in \code{df}.
#' @param genes Character vector of gene (expression) columns to test.
#' @return List of class \code{cpas_COX_by_genes}:
#'   \item{\code{input_params}:}{echoed inputs and analysis time}
#'   \item{\code{processed_data}:}{analysis data frame after cleaning}
#'   \item{\code{individual_models}:}{named list of \code{coxph} fits}
#'   \item{\code{individual_summaries}:}{named list of model summaries}
#'   \item{\code{results_table}:}{data.frame with columns \code{gene},
#'     \code{HR}, \code{HR95L}, \code{HR95H}, \code{Pvalue} and \code{P_adj}
#'     (Benjamini-Hochberg FDR across all genes tested) plus \code{P_adj_text}}
#'   \item{\code{metadata}:}{sample size, events, gene count, and
#'     \code{failed_genes} / \code{failure_reasons} for genes that could not be
#'     estimated (separation, non-convergence, constant gene)}
#' @details
#' HRs are per one raw unit of the marker. For cross-platform comparability
#' consider per-SD standardization (see \code{\link{cpas_meta}}).
#' Every gene is tested separately, so the table also carries the
#' Benjamini-Hochberg FDR (\code{P_adj}): report it when the genes are a
#' screening panel rather than a pre-specified hypothesis. Rows with
#' missing time/status or missing gene values are dropped; genes that
#' produce a degenerate model are skipped with a warning and reported as
#' \code{NA} rows in \code{results_table} plus a note in the metadata.
#' @importFrom survival coxph Surv
#' @export
#' @examples
#' \dontrun{
#'    ## One real cohort, several genes: the prerequisite chain gives the table.
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB", "TP53"), type = "OS")
#' 
#'    r <- COX_by_genes(d, type = "OS", genes = c("GAPDH", "ACTB", "TP53"))
#'    head(r$results_table)
#' }
COX_by_genes <- function(df, type = "OS", genes) {
  if (!is.data.frame(df)) stop("'df' must be a data.frame.")
  if (missing(genes) || length(genes) == 0L)
    stop("'genes' must contain at least one column name.")
  genes <- unique(as.character(genes))

  # family (OS/DSS/DFS/PFS/MFS) or concrete token -> token present in the data
  type <- .resolve_type_in_df(df, type)
  tc <- paste0(type, "_time"); sc <- paste0(type, "_status")
  need <- c("ID", tc, sc)
  missing_cols <- setdiff(need, colnames(df))
  if (length(missing_cols))
    stop("Columns missing from 'df': ", paste(missing_cols, collapse = ", "), ".")
  miss_gene <- setdiff(genes, colnames(df))
  if (length(miss_gene))
    stop("Gene columns not found in 'df': ", paste(miss_gene, collapse = ", "), ".")

  dat <- df[c("ID", tc, sc, genes)]
  for (cn in setdiff(colnames(dat), "ID"))
    dat[[cn]] <- suppressWarnings(as.numeric(dat[[cn]]))
  if (any(!is.na(dat[[sc]]) & !dat[[sc]] %in% c(0, 1)))
    stop("Status column '", sc, "' must be coded 0 (censored) / 1 (event).")
  if (any(!is.na(dat[[tc]]) & dat[[tc]] < 0))
    stop("Time column '", tc, "' contains negative values.")

  keep <- stats::complete.cases(dat[c(tc, sc, genes)])
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) == 0L) stop("No complete rows available for analysis.")
  events <- sum(dat[[sc]] == 1, na.rm = TRUE)
  if (events < 2L) stop("Fewer than 2 events; Cox model cannot be fitted.")
  if (length(unique(dat[[sc]])) < 2L)
    stop("The status column contains only one level; Cox model cannot be fitted.")

  colnames(dat)[2:3] <- c("time", "status")

  individual_models <- list(); individual_summaries <- list()
  rows <- vector("list", length(genes)); notes <- character(0)
  reasons <- list()
  for (i in seq_along(genes)) {
    n <- genes[i]
    diag <- .cpas_COX_diag(stats::as.formula("Surv(time, status) ~ dat[[n]]"), dat)
    if (!diag$ok) {
      # separation / non-convergence / non-estimable effect: drop the marker and
      # record the reason instead of publishing a meaningless hazard ratio
      notes <- c(notes, n)
      reasons[[n]] <- diag$reason
      rows[[i]] <- data.frame(gene = n, HR = NA_real_, HR95L = NA_real_,
                              HR95H = NA_real_, Pvalue = NA_real_)
      next
    }
    f <- diag$fit
    individual_models[[n]] <- f
    s <- summary(f)
    individual_summaries[[n]] <- s
    cf <- s$coefficients[1, ]
    est <- cf["coef"]; se <- cf["se(coef)"]
    rows[[i]] <- data.frame(gene = n, HR = exp(est),
                            HR95L = exp(est - 1.96 * se),
                            HR95H = exp(est + 1.96 * se),
                            Pvalue = cf[["Pr(>|z|)"]])
  }
  res <- do.call(rbind, rows)
  # multiple testing: a panel of many genes needs an adjusted p-value
  res$P_adj <- stats::p.adjust(res$Pvalue, method = "BH")
  res$P_adj_text <- ifelse(is.na(res$P_adj), "",
    ifelse(res$P_adj < 1e-16, "<1e-16", formatC(res$P_adj, format = "g", digits = 4)))
  res <- res[order(res$Pvalue, na.last = TRUE), ]
  rownames(res) <- NULL
  if (length(notes))
    warning("The following genes could not be modelled (NA rows added): ",
            paste(notes, collapse = ", "), call. = FALSE)

  metadata <- list(sample_size = nrow(dat), events = events,
                   fdr_method = "BH",
                   n_significant_fdr005 = sum(res$P_adj < 0.05, na.rm = TRUE),
                   var_number = length(genes), gene_count = length(genes),
                   failed_genes = notes,
                   failure_reasons = reasons,
                   analysis_type = "Univariate Cox Regression",
                   survival_type = type)

  result_obj <- list(input_params = list(df_name = deparse(substitute(df)),
                                         type = type, genes = genes,
                                         analysis_time = Sys.time()),
                     processed_data = dat,
                     individual_models = individual_models,
                     individual_summaries = individual_summaries,
                     results_table = res,
                     metadata = metadata)
  class(result_obj) <- "cpas_COX_by_genes"
  result_obj
}
