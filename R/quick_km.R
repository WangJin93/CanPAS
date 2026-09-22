#' @title Quick Kaplan-Meier Plots for a Dataset
#' @description
#' One-stop helper: retrieves expression of the requested genes for a dataset,
#' merges it with the survival table, and draws a grid of Kaplan-Meier plots
#' (median split) for every gene that maps on the platform.
#' @param dataset Dataset accession in \code{dataset_info} (GEO), or
#' \code{"TCGA-<PROJECT>"} (e.g. \code{"TCGA-LUAD"}).
#' @param genes Character vector of gene symbols.
#' @param process_duplicates Collapsing rule for multiple probes per gene
#' (see \code{\link{get_expr_data}}). Default \code{"max"}.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token. A family
#' is resolved to the token available in the cohort (e.g. DFS -> RFS/DFI/EFS).
#' @param ncol Number of columns in the plot grid.
#' @param cutpoint Cut-point rule passed to \code{\link{plot_km}}.
#' @param ... Further arguments passed to \code{\link{plot_km}} /
#' \code{survminer::ggsurvplot}.
#' @return Invisibly, the list of \code{ggsurvplot} objects; the combined grid
#' is drawn to the active graphics device.
#' @export
#' @examples
#' \dontrun{
#'    ## quick_km() runs the same prerequisite chain internally, for one cohort
#'    ## and one or more genes, and draws the curves.
#'    quick_km("GSE14814", c("GAPDH", "ACTB"), pval = TRUE, legend = "bottom")
#' 
#'    ## The same data through the explicit chain, for customisation.
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
#'    plot_km(d, type = "OS", marker = "GAPDH", cutpoint = 12)   # absolute cut-off
#' }
quick_km <- function(dataset, genes,
                     process_duplicates = "max",
                     type = "OS",
                     ncol = 3,
                     cutpoint = "median",
                     ...) {
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
    plot_km(merged, type = type, marker = marker, cutpoint = cutpoint, ...)$plot)

  gridExtra::grid.arrange(grobs = plots, ncol = max(1, min(ncol, length(markers))))
  invisible(plots)
}
