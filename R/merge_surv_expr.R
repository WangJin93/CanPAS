#' @title Merge Survival and Expression Data of a Dataset
#' @description
#' Pulls the survival table of a dataset from the CanPAS API and joins it with
#' gene expression data produced by \code{\link{get_expr_data}} on the common
#' sample identifier (\code{ID}).
#'
#' For datasets whose accession starts with \code{"TCGA-"} the survival table is
#' not read from the mirror: the local TCGA clinical table is used instead
#' (\code{\link{tcga_surv_table}}, the same source the App and
#' \code{\link{tcga_merged}} use).
#' @param dataset Dataset accession in \code{dataset_info} (e.g. \code{"GSE14814"}
#'   or \code{"TCGA-LUAD"}).
#' @param expr_data Result object of class \code{cpas_get_expr} returned by
#' \code{\link{get_expr_data}} (an object with an \code{expr_data} element whose
#' first column is \code{ID}) or a plain data.frame in that format.
#' @return List of class \code{cpas_merge}:
#'   \item{\code{input_params}:}{echoed inputs and merge time}
#'   \item{\code{raw_surv_data}:}{survival table exactly as returned by the API
#'     (or the local TCGA clinical table)}
#'   \item{\code{raw_expr_data}:}{expression input}
#'   \item{\code{merged_data}:}{samples x (ID, \code{<TYPE>_time},
#'     \code{<TYPE>_status}, gene columns). All survival times are years and all
#'     status columns are binary (0 censored / 1 event); rows with missing time
#'     or status are kept (see \code{NA}) so that callers may decide how to treat
#'     them, while gene columns are numeric}
#'   \item{\code{metadata}:}{sample counts and matched columns}
#' @details
#' Time columns are interpreted as years; status columns must be binary
#' 0/1 (the format used by the CanPAS database). Columns of the survival table
#' other than \code{ID}, \code{*_status} and \code{*_time} (e.g. clinical
#' covariates) are intentionally not carried into \code{merged_data}: clinical
#' data live in the local survival RDS mirrors and can be joined by the caller.
#' @examples
#' \dontrun{
#'    ## The two prerequisite steps and the merge: expression first, survival
#'    ## second, giving the analysis-ready table used by every other function.
#'    e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
#'    m <- merge_surv_expr("GSE14814", e)
#'    colnames(m$merged_data)
#' 
#'    ## raw_surv_data keeps the clinical columns of that cohort (age, sex, ...);
#'    ## cohort_merged(clin = TRUE) adds them to the merged table for you.
#'    colnames(m$raw_surv_data)
#'    head(m$merged_data)
#' }
#' @export
merge_surv_expr <- function(dataset, expr_data) {
  if (inherits(expr_data, "cpas_get_expr")) {
    raw_expr <- expr_data$expr_data
  } else if (is.data.frame(expr_data)) {
    raw_expr <- expr_data
  } else {
    stop("'expr_data' must be a cpas_get_expr object or a data.frame ",
         "(first column named 'ID').")
  }
  if (!"ID" %in% colnames(raw_expr))
    stop("Expression data must contain an 'ID' column (samples).")

  # ---- survival table: local TCGA clinical table, or the mirror for GEO/CGGA ----
  if (startsWith(as.character(dataset)[1], "TCGA-")) {
    surv_df <- tcga_surv_table(sub("^TCGA-", "", as.character(dataset)[1]))
  } else {
    raw_surv <- get_data(dataset = dataset, action = "surv_data")
    surv_df <- raw_surv$response
  }
  if (is.null(surv_df) || !is.data.frame(surv_df) || nrow(surv_df) == 0L)
    stop("Survival query returned no rows for dataset ", dataset, ".")
  colnames(surv_df)[1] <- "ID"
  surv_df$ID <- as.character(surv_df$ID)

  keep <- c("ID", colnames(surv_df)[grepl("_(status|time)$", colnames(surv_df))])
  surv <- surv_df[, unique(keep), drop = FALSE]

  # ---- sanity checks / coercion ----
  for (tc in grep("_time$", colnames(surv), value = TRUE)) {
    v <- suppressWarnings(as.numeric(surv[[tc]]))
    if (any(!is.na(v) & v < 0))
      stop("Dataset ", dataset, ": negative survival time found in column ", tc, ".")
    surv[[tc]] <- v
  }
  for (sc in grep("_status$", colnames(surv), value = TRUE)) {
    v <- suppressWarnings(as.numeric(surv[[sc]]))
    bad <- !is.na(v) & !v %in% c(0, 1)
    if (any(bad))
      stop("Dataset ", dataset, ": status column ", sc,
           " contains values other than 0/1 (please standardize first).")
    surv[[sc]] <- v
  }
  if (anyDuplicated(as.character(raw_expr$ID)))
    stop("Expression data contain duplicated sample IDs (",
         sum(duplicated(as.character(raw_expr$ID))), " duplicated); merging them into the ",
         "survival table would duplicate patients and pseudo-replicate the analysis. ",
         "Please make sample IDs unique first.", call. = FALSE)

  # expression: numeric gene columns
  gcols <- setdiff(colnames(raw_expr), "ID")
  expr <- raw_expr
  expr[gcols] <- lapply(expr[gcols], function(x) suppressWarnings(as.numeric(x)))

  # ---- join ----
  merged <- merge(surv, expr, by = "ID", all.x = FALSE, all.y = FALSE, sort = TRUE)
  if (nrow(merged) == 0L)
    stop("No sample ID in common between the survival table (", nrow(surv),
         " samples) and the expression table (", nrow(expr),
         " samples) for dataset ", dataset, ". Check that both use the same sample IDs.",
         call. = FALSE)
  merged <- merged[, colSums(is.na(merged)) < nrow(merged), drop = FALSE]

  metadata <- list(dataset_accession = dataset,
                   samples_in_survival = nrow(surv),
                   samples_in_expression = nrow(expr),
                   samples_merged = nrow(merged),
                   survival_columns = setdiff(colnames(surv), "ID"),
                   genes_in_expression = length(gcols),
                   merged_columns = colnames(merged))

  result_obj <- list(input_params = list(dataset = dataset,
                                         expr_data_source = deparse(substitute(expr_data)),
                                         merge_time = Sys.time()),
                     raw_surv_data = surv_df,
                     raw_expr_data = expr_data,
                     merged_data = merged,
                     metadata = metadata)
  class(result_obj) <- "cpas_merge"
  result_obj
}
