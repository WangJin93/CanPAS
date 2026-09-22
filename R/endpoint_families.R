# endpoint_families.R --------------------------------------------------------
# Survival-endpoint families (analysis plan B).
#
#   pooling families (used by meta-analysis / pooled KM):
#     OS  = {OS}
#     DSS = {DSS, CSS, BCSS}          cancer-specific survival
#     DFS = {DFS, RFS, EFS, DFI}      disease-/recurrence-free survival
#     PFS = {PFS, PFI}                progression-free survival
#     MFS = {MFS, DRFS}               metastasis-/distant-recurrence-free
#
#   browsing families (UI grouping): OS, DSS, DFS and PFS (broad, i.e. EP_PFS
#   or EP_MFS).
#
# The mapping is materialised as columns in the packaged `dataset_info`
# table: EndpointFamilies, EP_OS/EP_DSS/EP_DFS/EP_PFS/EP_MFS, EndpointPrimary
# and EndpointDerived (TCGA-derived DFI/PFI). Nothing is changed in the
# database.

.cpas_pooling_families <- c("OS", "DSS", "DFS", "PFS", "MFS")

.cpas_family_tokens <- list(
  OS  = "OS",
  DSS = c("DSS", "CSS", "BCSS"),
  DFS = c("DFS", "RFS", "EFS", "DFI"),
  PFS = c("PFS", "PFI"),
  MFS = c("MFS", "DRFS"))

.cpas_browse_family <- c(OS = "OS", DSS = "DSS", DFS = "DFS",
                         PFS = "PFS", MFS = "PFS")

.cpas_dataset_info <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    e <- new.env(parent = emptyenv())
    utils::data("dataset_info", package = "CanPAS", envir = e)
    cache <<- get("dataset_info", e)
    cache
  }
})

#' @title Pooling family of a survival endpoint
#' @description Maps a raw endpoint token (e.g. \code{"RFS"}) to the pooling
#' family used by meta-analysis and pooled Kaplan-Meier
#' (\code{OS}, \code{DSS}, \code{DFS}, \code{PFS}, \code{MFS}). Family names
#' are returned unchanged, unknown tokens give \code{NA}.
#' @param type Endpoint token or family name.
#' @return Character scalar (family) or \code{NA_character_}.
#' @export
#' @examples
#'   ## Families group the tokens that different cohorts use for the same
#'   ## endpoint, so one analysis can span cohorts.
#'   endpoint_family("RFS")   # "DFS"
#'   endpoint_family("PFI")   # "PFS"
#'   endpoint_family("OS")    # "OS"
#'
#'   ## The catalog records which families each real cohort can be analysed on.
#'   di <- dataset_info[dataset_info$Accession == "GSE13507", ]
#'   di[, c("Accession", "Type", "N", "EndpointFamilies")]
endpoint_family <- function(type) {
  if (is.null(type) || !length(type)) return(NA_character_)
  type <- as.character(type)[1]
  if (type %in% .cpas_pooling_families) return(type)
  hit <- names(.cpas_family_tokens)[vapply(.cpas_family_tokens,
                                           function(v) type %in% v, logical(1))]
  if (length(hit)) hit[1] else NA_character_
}

#' @title Resolve an endpoint family for one cohort
#' @description Returns the concrete endpoint token available for a cohort and
#' endpoint (family name or raw token), using the columns of the packaged
#' \code{dataset_info} catalog (e.g. \code{EP_DFS}).
#' @param dataset Dataset accession (e.g. \code{"GSE31210"}, \code{"TCGA-LGG"}).
#' @param type Family name (\code{"OS"}, \code{"DSS"}, \code{"DFS"},
#' \code{"PFS"}, \code{"MFS"}) or a raw endpoint token.
#' @return Endpoint token (e.g. \code{"RFS"}). A family is resolved to the token
#' the cohort actually provides; a raw token is returned only when it is a known
#' endpoint token, otherwise \code{NA_character_} (which is also returned when the
#' cohort does not provide the family).
#' @export
#' @examples
#'   ## The token a real cohort actually contributes for a family.
#'   endpoint_resolve("GSE31210", "DFS")   # "RFS": that cohort is RFS-based
#'   endpoint_resolve("GSE13507", "OS")    # "OS"
#'   endpoint_resolve("TCGA-LUAD", "PFS")  # "PFI"
#'
#'   ## NA is returned when the cohort has no endpoint of that family.
#'   endpoint_resolve("GSE31210", "DSS")
endpoint_resolve <- function(dataset, type) {
  if (is.null(dataset) || is.null(type)) return(NA_character_)
  dataset <- as.character(dataset)[1]; type <- as.character(type)[1]
  if (!type %in% .cpas_pooling_families) {
    # a raw token is accepted only when it is a real endpoint token
    known <- unique(unlist(.cpas_family_tokens[.cpas_pooling_families],
                           use.names = FALSE))
    return(if (type %in% known) type else NA_character_)
  }
  di <- .cpas_dataset_info()
  col <- paste0("EP_", type)
  if (!col %in% colnames(di) || !dataset %in% di$Accession) return(NA_character_)
  tok <- di[[col]][match(dataset, di$Accession)]
  if (length(tok) && !is.na(tok) && nzchar(tok)) tok else NA_character_
}

# endpoints available for one cohort, labelled by family
.cpas_endpoints_of <- function(dataset) {
  di <- .cpas_dataset_info()
  out <- character(0)
  if (!dataset %in% di$Accession) return(out)
  i <- match(dataset, di$Accession)
  for (f in .cpas_pooling_families) {
    col <- paste0("EP_", f)
    if (!col %in% colnames(di)) next
    tok <- di[[col]][i]
    if (is.na(tok) || !nzchar(tok)) next
    lab <- if (identical(tok, f)) f else sprintf("%s (%s)", f, tok)
    if ("EndpointDerived" %in% colnames(di)) {
      dv <- di$EndpointDerived[i]
      if (!is.na(dv) && grepl(tok, dv, fixed = TRUE)) lab <- paste0(lab, " [derived]")
    }
    out[lab] <- tok
  }
  out
}

#' @title Endpoints available for a cohort, grouped by family
#' @description Returns one row per endpoint family available for a cohort,
#' with the concrete endpoint token, whether the token is a TCGA-derived
#' endpoint, the analysis family and a display label such as
#' \code{"DFS (RFS)"} or \code{"DFS (DFI) [derived]"}.
#' @param dataset Dataset accession (e.g. \code{"GSE31210"}, \code{"TCGA-LGG"}).
#' @param family Optional family filter (one or more of \code{"OS"},
#' \code{"DSS"}, \code{"DFS"}, \code{"PFS"}, \code{"MFS"}).
#' @return data.frame with columns \code{family}, \code{token},
#' \code{derived}, \code{label}.
#' @export
#' @examples
#'   ## What a real cohort can be analysed on, with the token per family.
#'   endpoint_options("TCGA-LGG")
#'   endpoint_options("GSE13507")
#'
#'   ## Only one family at a time, when that is all you need.
#'   endpoint_options("GSE31210", family = "DFS")
endpoint_options <- function(dataset, family = NULL) {
  di <- .cpas_dataset_info()
  dataset <- as.character(dataset)[1]
  out <- data.frame(family = character(0), token = character(0),
                    derived = logical(0), label = character(0),
                    stringsAsFactors = FALSE)
  if (!dataset %in% di$Accession) return(out)
  i <- match(dataset, di$Accession)
  for (f in .cpas_pooling_families) {
    col <- paste0("EP_", f)
    if (!col %in% colnames(di)) next
    tok <- di[[col]][i]
    if (is.na(tok) || !nzchar(tok)) next
    dv <- if ("EndpointDerived" %in% colnames(di)) di$EndpointDerived[i] else NA
    deriv <- !is.na(dv) && grepl(tok, dv, fixed = TRUE)
    lab <- if (identical(tok, f)) f else sprintf("%s (%s)", f, tok)
    if (deriv) lab <- paste0(lab, " [derived]")
    out <- rbind(out, data.frame(family = f, token = tok, derived = deriv,
                                 label = lab, stringsAsFactors = FALSE))
  }
  if (!is.null(family)) out <- out[out$family %in% family, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Resolve a family (or concrete token) against the columns actually present in
# a data.frame; used by every analysis function so that a family such as
# "DFS" is automatically mapped to RFS / DFI / EFS / DFS in the data at hand.
.resolve_type_in_df <- function(df, type) {
  if (is.null(type) || !length(type)) stop("'type' must be provided.")
  type <- as.character(type)[1]
  has <- function(tok) all(c(paste0(tok, "_time"), paste0(tok, "_status")) %in% colnames(df))
  if (has(type)) return(type)
  fam <- endpoint_family(type)
  if (is.na(fam)) stop("Endpoint '", type, "' has no *_time/*_status columns in the data.")
  for (tok in .cpas_family_tokens[[fam]]) if (has(tok)) return(tok)
  stop("Neither family '", fam, "' (",
       paste(.cpas_family_tokens[[fam]], collapse = "/"),
       ") nor token '", type, "' is available in the data.")
}
