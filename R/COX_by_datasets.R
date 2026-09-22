#' @title Cox Regression for Many Datasets (Single Gene)
#' @description
#' Runs the same univariable Cox analysis (one gene) on several datasets
#' (GEO accessions or \code{TCGA-<PROJECT>} projects) and returns one combined
#' tidy table plus the per-dataset objects. The estimates are reported per
#' dataset and are NOT pooled: use \code{\link{cpas_meta}} for an
#' inverse-variance meta-analysis with heterogeneity statistics. Expected
#' hazard ratios are per raw expression unit, so they are not comparable
#' across platforms; \code{cpas_meta} standardises the gene within each dataset
#' (per SD) before pooling.
#' @param datasets Character vector of dataset accessions or TCGA projects
#'   (e.g. \code{c("GSE14814", "GSE31210")} or \code{"TCGA-LUAD"}).
#' @param gene Single gene symbol to test.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token; families
#'   are resolved per cohort and the token used is reported in \code{endpoint}.
#' @param precision Decimal places for formatted display columns.
#' @param process_duplicates Probe collapsing rule
#' (see \code{\link{get_expr_data}}). Default \code{"max"}.
#' @return List of class \code{cpas_multi_cox}:
#'   \item{\code{input_params}:}{echoed inputs}
#'   \item{\code{individual_results}:}{per successful dataset, a list with
#'     \code{expr_data}, \code{merged_data} and \code{cox_analysis}}
#'   \item{\code{combined_results}:}{data.frame with columns \code{dataset},
#'     \code{endpoint} (the token actually used), \code{Variates} (gene),
#'     \code{Level}, \code{N}, \code{HR}, \code{HR95L}, \code{HR95H},
#'     \code{Pvalue} and a Benjamini-Hochberg \code{P_adj} across all tested
#'     dataset x gene combinations}
#'   \item{\code{metadata}, \code{errors}:}{bookkeeping}
#' @export
#' @examples
#' \dontrun{
#'    ## The same gene in real cohorts of one cancer type (lung, OS).
#'    di <- dataset_info[dataset_info$Type == "Lung Cancer", "Accession"]
#'    head(di)
#' 
#'    r <- COX_by_datasets(c("GSE14814", "GSE31210"), gene = "GAPDH", type = "OS")
#'    head(r$results_table)
#' }
COX_by_datasets <- function(datasets, gene, type = "OS", precision = 3,
                            process_duplicates = "max") {
  if (missing(datasets) || length(datasets) == 0L)
    stop("'datasets' must contain at least one accession.")
  datasets <- as.character(datasets)
  if (length(gene) != 1L) stop("'gene' must be a single symbol.")
  individual_results <- list(); errors <- list(); rows <- list()

  for (dataset in datasets) {
    one <- tryCatch({
      tok <- endpoint_resolve(dataset, type)
      if (is.na(tok)) tok <- as.character(type)[1]   # cohort-specific token
      e  <- get_expr_data(dataset, gene, process_duplicates = process_duplicates)
      m  <- merge_surv_expr(dataset, e)
      cx <- COX_analysis(m$merged_data, type = tok, cont_Variates = gene,
                         cate_Variates = NULL, method = "uni",
                         precision = precision)
      r <- cx$results_table
      list(expr_data = e, merged_data = m, cox_analysis = cx,
           row = data.frame(dataset = dataset, endpoint = tok,
                            Variates = gene, Level = NA_character_,
                            N = r$N[1], HR = r$HR[1], HR95L = r$HR95L[1],
                            HR95H = r$HR95H[1], Pvalue = r$Pvalue[1],
                            stringsAsFactors = FALSE))
    }, error = function(e) structure(list(error = conditionMessage(e)),
                                     class = "try-error"))
    if (inherits(one, "try-error")) {
      warning("Dataset ", dataset, " failed: ", one$error, call. = FALSE)
      errors[[dataset]] <- one$error
    } else {
      individual_results[[dataset]] <- one[c("expr_data", "merged_data", "cox_analysis")]
      rows[[dataset]] <- one$row
    }
  }

  combined_results <- if (length(rows)) dplyr::bind_rows(rows) else
    data.frame(dataset = character(0), endpoint = character(0),
               Variates = character(0),
               Level = character(0), N = integer(0), HR = numeric(0),
               HR95L = numeric(0), HR95H = numeric(0), Pvalue = numeric(0))
  # every dataset x gene combination is tested separately: report the FDR
  if ("Pvalue" %in% colnames(combined_results))
    combined_results$P_adj <- stats::p.adjust(combined_results$Pvalue, method = "BH")

  result_obj <- list(
    input_params = list(datasets = datasets, gene = gene, type = type,
                        precision = precision, analysis_time = Sys.time()),
    individual_results = individual_results,
    combined_results = combined_results,
    metadata = list(total_datasets = length(datasets),
                    successful_datasets = length(individual_results),
                    failed_datasets = length(errors),
                    gene_analyzed = gene,
                    survival_type = type,
                    precision = precision),
    errors = errors)
  class(result_obj) <- "cpas_COX_by_datasets"
  result_obj
}
