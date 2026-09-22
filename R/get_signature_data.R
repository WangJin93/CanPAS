#' @title Weighted Gene-Signature Score per Sample
#' @description
#' Parses a linear signature string such as
#' \code{"0.3*GAPDH + 0.7*ACTB + 0.5*RPN1"} and computes the weighted sum of
#' expression values for every sample of a dataset.
#' @param signature Character string of the form
#' \code{<weight>*<GENE> +/- <weight>*<GENE> ...} (weights may be negative).
#' @param dataset Dataset accession (GEO or \code{"TCGA-<PROJECT>"}).
#' @param process_duplicates Probe collapsing rule passed to
#' \code{\link{get_expr_data}} (default \code{"max"}).
#' @param allow_missing If \code{FALSE} (default), all signature genes must be
#' measurable on the platform; otherwise the score is computed on the subset of
#' genes that are available (with a message).
#' @return List of class \code{cpas_signature}:
#'   \item{\code{input_params}:}{echoed inputs}
#'   \item{\code{signature_info}:}{\code{original_signature}, \code{genes},
#'     \code{weights} and \code{weight_gene_pairs} (aligned to available genes)}
#'   \item{\code{expr_data}:}{full \code{cpas_get_expr} object}
#'   \item{\code{results_table}:}{data.frame \code{ID}, \code{signature}}
#'   \item{\code{metadata}:}{bookkeeping}
#' @details Weights are applied in the order of the signature string. Genes
#' that cannot be measured are reported; with \code{allow_missing = TRUE} the
#' score uses only measured genes (weights of missing genes dropped).
#' @export
#' @examples
#' \dontrun{
#'    ## A signature is scored inside a real cohort: expression is fetched,
#'    ## probes are collapsed, and missing genes stop the call unless allowed.
#'    s <- get_signature_value("0.5*GAPDH + 0.5*ACTB", "GSE14814")
#'    s$metadata$genes_used; s$metadata$sample_count
#'    s$signature_info$weights         # the parsed weights, in order
#'    head(s$results_table)            # ID + one score per patient
#' 
#'    ## Genes that are not measurable in that cohort: stop, or score the subset.
#'    s2 <- get_signature_value("0.5*GAPDH + 0.5*ACTB", "GSE13507",
#'                               allow_missing = TRUE)
#'    s2$metadata$genes_used; s2$metadata$genes_missing
#' }
get_signature_value <- function(signature, dataset,
                                process_duplicates = "max",
                                allow_missing = FALSE) {
  if (!is.character(signature) || length(signature) != 1L || !nzchar(signature))
    stop("'signature' must be a single non-empty string.")
  if (!grepl("\\*", signature)) stop("Signature must contain '<weight>*<GENE>' terms.")

  # parse "<sign><weight>*<gene>" terms; weights may use scientific notation.
  # Whitespace is allowed everywhere, including between the sign that joins two
  # terms and the weight that follows it, so that the documented form
  # "0.5*GAPDH + 0.5*TNS1" parses as written (the app's parser has always
  # accepted it, and the package's own examples use it).
  num  <- "(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
  pat  <- paste0("[+-]?[[:space:]]*", num, "[[:space:]]*\\*[[:space:]]*",
                 "[A-Za-z][A-Za-z0-9._-]*")
  terms <- regmatches(signature, gregexpr(pat, signature))[[1]]
  if (!length(terms)) stop("Could not parse any '<weight>*<GENE>' term from the signature.")
  # every character of the formula must be consumed by a parsed term: otherwise the
  # score would silently belong to a different formula than the one the user wrote
  residue <- gsub("[[:space:]]", "", gsub(paste0("(", pat, ")"), "", signature, perl = TRUE))
  if (nzchar(residue))
    stop("Unsupported text in the signature: '", residue, "'. Use only '<weight>*<GENE>' ",
         "terms joined by '+' or '-' (for example 0.5*GAPDH + 0.5*TNS1). A hyphen that ",
         "is directly attached to the next weight is read as part of the gene symbol ",
         "(gene symbols may contain '-', e.g. HLA-DRA); write 'GAPDH - 0.3*ACTB' with ",
         "spaces when the hyphen is a minus sign.", call. = FALSE)
  clean <- gsub("[[:space:]]", "", terms)
  weights <- as.numeric(sub("\\*.*$", "", clean))
  genes   <- sub("^.*\\*", "", clean)
  if (anyNA(weights) || any(!is.finite(weights)))
    stop("Could not parse a finite numeric weight in: ", signature)
  dup <- genes[duplicated(genes)]
  if (length(dup)) stop("Duplicate gene in signature: ", paste(unique(dup), collapse = ", "),
                        ". Combine terms or rename first.")

  if (startsWith(dataset, "TCGA-")) {
    tdat <- tcga_get_expr(sub("^TCGA-", "", dataset), genes)
    if (colnames(tdat)[1] == "sample") colnames(tdat)[1] <- "ID"
    expr_obj <- list(expr_data = tdat,
                     metadata = list(genes_matched_platform = length(genes)))
  } else {
    expr_obj <- get_expr_data(dataset, genes, process_duplicates = process_duplicates)
  }
  avail <- intersect(genes, colnames(expr_obj$expr_data))
  missing_genes <- setdiff(genes, avail)
  if (length(missing_genes)) {
    msg <- paste0("Signature gene(s) not measurable in ", dataset, ": ",
                  paste(missing_genes, collapse = ", "))
    if (!allow_missing) stop(msg, ". Use allow_missing = TRUE to score the available subset.")
    message(msg, " - scoring the available subset.")
  }
  if (!length(avail)) stop("None of the signature genes are measurable in ", dataset, ".")

  keep <- genes %in% avail
  w2 <- weights[keep]
  g2 <- genes[keep]
  # align columns with the signature order
  col_ord <- match(g2, colnames(expr_obj$expr_data))
  dat <- expr_obj$expr_data[, c("ID", colnames(expr_obj$expr_data)[col_ord]), drop = FALSE]

  m <- as.matrix(dat[, -1L, drop = FALSE])
  m[] <- apply(m, 2, function(x) suppressWarnings(as.numeric(x)))
  if (anyNA(m)) {
    bad <- colnames(m)[colSums(is.na(m)) > 0]
    warning("Expression contains missing values in gene(s): ",
            paste(bad, collapse = ", "), " - treated as NA (score becomes NA).",
            call. = FALSE)
  }
  score <- drop(m %*% w2)
  results <- data.frame(ID = dat$ID, signature = score, stringsAsFactors = FALSE)

  result_obj <- list(
    input_params = list(signature = signature, dataset = dataset,
                        process_duplicates = process_duplicates,
                        allow_missing = allow_missing,
                        analysis_time = Sys.time()),
    signature_info = list(original_signature = signature,
                          genes = g2, weights = w2,
                          weight_gene_pairs = data.frame(gene = g2, weight = w2)),
    expr_data = expr_obj,
    results_table = results,
    metadata = list(dataset = dataset,
                    signature_length = length(genes),
                    genes_used = length(g2),
                    genes_missing = length(missing_genes),
                    duplicate_handling = process_duplicates,
                    sample_count = nrow(results)))
  class(result_obj) <- "cpas_signature"
  result_obj
}
