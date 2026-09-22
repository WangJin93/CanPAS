#' @title Quick Time-dependent ROC Plots for a Dataset
#' @description
#' One-stop helper: retrieves expression of the requested genes for a dataset,
#' merges it with the survival table, and draws a grid of time-dependent ROC
#' curves for every gene that maps on the platform.
#' @param dataset Dataset accession in \code{dataset_info} (GEO), or
#' \code{"TCGA-<PROJECT>"} (e.g. \code{"TCGA-LUAD"}).
#' @param genes Character vector of gene symbols.
#' @param process_duplicates Collapsing rule for multiple probes per gene
#' (see \code{\link{get_expr_data}}). Default \code{"max"}.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token. A family
#' is resolved to the token available in the cohort (e.g. DFS -> RFS/DFI/EFS).
#' @param predict.time Prediction horizon (same time unit as the survival
#' column, typically years).
#' @param method ROC method passed to \code{\link{plot_roc}}.
#' @param ncol Number of columns in the plot grid.
#' @param ... Further arguments passed to \code{\link{plot_roc}}.
#' @return Invisibly, the list of ggplot objects; the combined grid is drawn
#' to the active graphics device.
#' @export
#' @examples
#' \dontrun{
#'    quick_roc("GSE14814", c("GAPDH", "ACTB"), predict.time = 3)
#' 
#'    ## Equivalent, using the prerequisite chain for a single cohort.
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
#'    plot_roc(d, type = "OS", marker = "GAPDH", predict.time = 3)$auc
#' }
quick_roc <- function(dataset, genes,
                      process_duplicates = "max",
                      type = "OS",
                      predict.time = 1,
                      method = c("KM", "NNE"),
                      ncol = 3,
                      ...) {
  method <- match.arg(method)
  merged <- .cpas_merged(dataset, genes, process_duplicates)
  type   <- .resolve_type_in_df(merged, type)     # family or token
  markers <- intersect(unique(genes), .cpas_markers(merged))
  if (!length(markers))
    stop("None of the requested genes were mapped in dataset ", dataset, ".")
  markers <- markers[vapply(markers, function(m) {
    v <- suppressWarnings(as.numeric(merged[[m]])); !all(is.na(v))
  }, logical(1))]
  if (!length(markers)) stop("No numeric expression values available for the requested genes.")

  plots <- lapply(markers, function(marker)
    plot_roc(merged, type = type, marker = marker,
             predict.time = predict.time, method = method, ...))

  gridExtra::grid.arrange(grobs = plots, ncol = max(1, min(ncol, length(markers))))
  invisible(plots)
}
