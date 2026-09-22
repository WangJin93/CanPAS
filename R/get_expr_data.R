#' @title Retrieve Expression Data for Genes of a Dataset
#' @description
#' Resolves gene symbols to probes on the dataset's platform (via the ID_map
#' table and the platform's probe-to-gene mapping) and downloads the matching
#' expression values from the CanPAS API. Multiple probes mapping to the same
#' gene can be collapsed with a user-selected rule.
#'
#' Datasets whose accession starts with \code{"TCGA-"} are handled differently:
#' the mirror holds no TCGA expression table, so values are fetched per gene
#' from UCSC Xena on demand (\code{\link{tcga_gene_expr_df}}, the same source the
#' App uses). Xena returns gene-level log2 values, so there is no probe
#' resolution or collapsing and \code{process_duplicates} is ignored.
#' @param dataset Dataset accession in \code{dataset_info}, e.g. \code{"GSE14814"}
#'   or \code{"TCGA-LUAD"}.
#' @param gene Character vector of gene symbols.
#' @param genes One or more gene symbols (character); equivalent to \code{gene}
#'   (supply only one of the two).
#' @param process_duplicates How to collapse several probes that map to the same
#'   gene: \code{"max"} (default), \code{"mean"}, \code{"median"},
#'   \code{"min"} - or \code{"no"} to keep probe-level rows (columns then
#'   carry probe ids instead of symbols). Ignored for TCGA datasets.
#' @return List of class \code{cpas_get_expr}:
#'   \item{\code{input_params}:}{echoed inputs and analysis time}
#'   \item{\code{raw_ids}:}{rows of ID_map matching the requested symbols}
#'   \item{\code{platform_info}:}{dataset_info row of the dataset}
#'   \item{\code{ref_ids}:}{probe (row_names) to gene mapping actually used;
#'     for TCGA one row per gene, since there are no probes}
#'   \item{\code{expr_data}:}{data.frame, first column \code{ID} (samples),
#'     remaining columns are genes (or probes) with numeric expression}
#'   \item{\code{metadata}:}{genes requested / found / matched / sample count,
#'     plus \code{source} (\code{"mirror"} or \code{"xena"})}
#' @details
#' The internal tables \code{ID_map} and \code{dataset_info} ship with the
#' package and are loaded lazily; they never need to be attached manually.
#' Genes that exist in \code{ID_map} but have no probe on the dataset's
#' platform are skipped with a message. For TCGA datasets the symbol is passed
#' to Xena directly, so a gene missing from \code{ID_map} is still attempted,
#' and genes absent from Xena are skipped with a message.
#' @examples
#' \dontrun{
#'    ## GEO: probe ids are resolved through the platform map and collapsed to
#'    ## one value per gene.
#'    e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
#'    e$metadata
#'    head(e$expr_data)
#' 
#'    ## TCGA: the mirror holds no expression table, so values are fetched per
#'    ## gene from UCSC Xena (gene-level log2, no probe collapsing).
#'    tcga <- get_expr_data("TCGA-LUAD", "TP53")
#'    head(tcga$expr_data)
#' }
#' @export

get_expr_data <- function(dataset, gene = NULL,
                          process_duplicates = c("max", "mean", "median", "min", "no"),
                          genes = NULL) {
  process_duplicates <- match.arg(process_duplicates)
  if (is.null(gene)) {
    if (is.null(genes))
      stop("Supply the gene symbols as 'genes' (preferred) or 'gene'.")
    gene <- genes
  } else if (!is.null(genes) &&
             !identical(as.character(gene), as.character(genes))) {
    stop("'gene' and 'genes' were both supplied with different values; use only 'genes'.")
  }
  if (length(gene) == 0L)
    stop("'genes' must contain at least one gene symbol.")
  gene <- unique(as.character(gene))
  if (anyNA(gene) || any(!grepl("^[A-Za-z0-9._\\-]+$", gene)))
    stop("'gene' contains unsupported characters.")

  # ---- internal tables (lazy data) ----
  pkg_data <- function(name) {
    if (exists(name, envir = parent.frame(), inherits = TRUE)) {
      return(get(name, envir = parent.frame(), inherits = TRUE))
    }
    e <- new.env(parent = emptyenv())
    utils::data(list = name, package = "CanPAS", envir = e)
    if (exists(name, envir = e, inherits = FALSE)) get(name, envir = e) else NULL
  }
  ID_map       <- pkg_data("ID_map")
  dataset_info <- pkg_data("dataset_info")
  if (is.null(dataset_info) || !all(c("Accession", "GPL") %in% colnames(dataset_info)))
    stop("Internal table dataset_info (Accession/GPL) is unavailable.")

  # ---- genes in ID_map ----
  raw_ids <- ID_map[ID_map$Symbol %in% gene, , drop = FALSE]
  raw_ids <- raw_ids[!is.na(raw_ids$gene_id), , drop = FALSE]
  if (nrow(raw_ids) == 0L)
    stop("None of the requested genes were found in ID_map. Check symbols: ",
         paste(gene, collapse = ", "), ".")
  missing_sym <- setdiff(gene, raw_ids$Symbol)
  if (length(missing_sym))
    message("Gene(s) not present in ID_map, skipped: ",
            paste(missing_sym, collapse = ", "), ".")

  # ---- platform of the dataset ----
  platform_info <- dataset_info[dataset_info$Accession == dataset, , drop = FALSE]
  if (nrow(platform_info) == 0L)
    stop("Dataset ", sQuote(dataset), " not found in dataset_info.")
  gpl <- platform_info$GPL[1]

  # ============ TCGA: no mirror expression table, fetch from Xena =============
  if (startsWith(as.character(dataset)[1], "TCGA-")) {
    project <- sub("^TCGA-", "", as.character(dataset)[1])
    got <- list(); failed <- character(0)
    for (g in gene) {
      d <- tryCatch(tcga_gene_expr_df(project, g), error = function(e) NULL)
      # Xena answers with a column of NaN (517 rows, no error) for genes it does
      # not carry, so an all-non-finite result counts as "not available".
      if (is.null(d) || !nrow(d) || all(!is.finite(suppressWarnings(as.numeric(d[[2]]))))) {
        failed <- c(failed, g); next
      }
      d$sample <- as.character(d$sample)
      colnames(d)[2] <- g
      got[[g]] <- d
    }
    if (!length(got))
      stop("No expression returned for ", dataset, " (UCSC Xena). Requested genes: ",
           paste(gene, collapse = ", "), ".")
    if (length(failed))
      message("Gene(s) not available in ", dataset, " (Xena), skipped: ",
              paste(failed, collapse = ", "), ".")
    m <- Reduce(function(a, b) merge(a, b, by = "sample", all = TRUE), got)
    m <- m[order(m$sample), , drop = FALSE]
    gcols <- setdiff(colnames(m), "sample")
    m[gcols] <- lapply(m[gcols], function(x) suppressWarnings(as.numeric(x)))
    expr_data <- data.frame(ID = as.character(m$sample), m[, gcols, drop = FALSE],
                            stringsAsFactors = FALSE, check.names = FALSE)
    rownames(expr_data) <- NULL
    # one row per gene keeps the object shape identical to the GEO path
    ref_ids <- data.frame(row_names = gcols, gene_id = NA_character_,
                          Symbol = gcols, Ensembl_gene = NA_character_,
                          stringsAsFactors = FALSE)
    id_hit <- if (!is.null(ID_map) && all(c("gene_id", "Symbol") %in% colnames(ID_map)))
      ID_map[ID_map$Symbol %in% gcols, , drop = FALSE] else NULL
    if (!is.null(id_hit) && nrow(id_hit))
      ref_ids$gene_id <- id_hit$gene_id[match(ref_ids$Symbol, id_hit$Symbol)]
    metadata <- list(dataset_accession = dataset, platform = gpl, source = "xena",
                     total_genes_requested = length(gene),
                     genes_found = length(gcols),
                     genes_matched_platform = length(gcols),
                     duplicate_handling = NA_character_,
                     sample_count = nrow(expr_data),
                     gene_columns = gcols)
    out <- list(input_params = list(dataset = dataset, gene = gene,
                                    process_duplicates = process_duplicates,
                                    analysis_time = Sys.time()),
                raw_ids = id_hit, platform_info = platform_info,
                ref_ids = ref_ids, expr_data = expr_data, metadata = metadata)
    class(out) <- "cpas_get_expr"
    return(out)
  }

  # ============ GEO / CGGA: mirror API =======================================
  if (is.null(ID_map) || !all(c("gene_id", "Symbol") %in% colnames(ID_map)))
    stop("Internal table ID_map (gene_id/Symbol) is unavailable.")


  # ---- probe map on this platform ----
  ref <- get_data(dataset = gpl, action = "gpl",
                  ids = unique(as.character(raw_ids$gene_id)))
  if (is.null(ref$response) || !is.data.frame(ref$response) || nrow(ref$response) == 0L)
    stop("Platform ", gpl, " returned no probe-to-gene mapping for the requested genes.")
  ref_ids <- merge(ref$response, raw_ids, by = "gene_id")
  ref_ids <- ref_ids[!duplicated(ref_ids$row_names), , drop = FALSE]
  if (nrow(ref_ids) == 0L)
    stop("No probes for the requested genes were found on platform ", gpl, ".")
  if (nrow(ref_ids) < nrow(raw_ids))
    message("Some genes have no probe on platform ", gpl, " and are skipped.")

  # ---- expression ----
  expr <- get_data(dataset = dataset, action = "expression",
                   ids = unique(as.character(ref_ids$row_names)))
  if (is.null(expr$response) || !is.data.frame(expr$response) || nrow(expr$response) == 0L)
    stop("Expression query returned no rows for dataset ", dataset, ".")

  first_col <- colnames(expr$response)[1]
  expr_raw <- merge(ref_ids, expr$response, by.x = "row_names", by.y = first_col,
                    all.x = FALSE, all.y = FALSE, sort = FALSE)
  expr_raw <- expr_raw[!duplicated(expr_raw$row_names), , drop = FALSE]

  sample_cols <- setdiff(colnames(expr_raw),
                         c("row_names", "gene_id", "Symbol", "Ensembl_gene"))
  expr_raw[sample_cols] <- lapply(expr_raw[sample_cols],
                                  function(x) suppressWarnings(as.numeric(x)))

  # ---- collapse probes -> one value per gene ----
  if (process_duplicates != "no") {
    aggf <- switch(process_duplicates,
                   max    = function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE),
                   mean   = function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE),
                   median = function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE),
                   min = function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
    agg <- dplyr::summarise(dplyr::group_by(expr_raw, Symbol),
                           dplyr::across(dplyr::all_of(sample_cols), aggf),
                           .groups = "drop")
    expr_out <- as.data.frame(agg, stringsAsFactors = FALSE)
    colnames(expr_out)[1] <- "ID"
    # order columns by the order genes were requested
    ord <- match(gene, expr_out$ID)
    ord <- ord[!is.na(ord)]
    expr_out <- expr_out[order(ord), , drop = FALSE]
  } else {
    expr_out <- expr_raw[, c("row_names", sample_cols), drop = FALSE]
    colnames(expr_out)[1] <- "ID"       # ID column keeps probe ids when no collapse
  }

  # transpose: rows = samples, columns = genes (or probes); base-R only
  genes_vec <- as.character(expr_out[[1L]])
  vals <- as.matrix(expr_out[, -1L, drop = FALSE])
  tvals <- t(vals)
  expr_data <- data.frame(ID = rownames(tvals), tvals,
                          stringsAsFactors = FALSE, check.names = FALSE)
  colnames(expr_data)[-1L] <- genes_vec
  rownames(expr_data) <- NULL

  metadata <- list(dataset_accession = dataset, source = "mirror",
                   platform = gpl,
                   total_genes_requested = length(gene),
                   genes_found = length(unique(raw_ids$Symbol)),
                   genes_matched_platform = length(unique(ref_ids$Symbol)),
                   duplicate_handling = process_duplicates,
                   sample_count = nrow(expr_data),
                   gene_columns = setdiff(colnames(expr_data), "ID"))

  result_obj <- list(input_params = list(dataset = dataset, gene = gene,
                                         process_duplicates = process_duplicates,
                                         analysis_time = Sys.time()),
                     raw_ids = raw_ids,
                     platform_info = platform_info,
                     ref_ids = ref_ids,
                     expr_data = expr_data,
                     metadata = metadata)
  class(result_obj) <- "cpas_get_expr"
  result_obj
}
