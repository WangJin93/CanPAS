# endpoint_semantics.R -------------------------------------------------------
# Endpoint semantics layer (workstream P, spec A1).
#
# The curation side ships one row per cohort x endpoint family in
# data/endpoint_semantics.csv; the package ships a byte-identical copy as
# inst/extdata/endpoint_semantics.csv and reads it here. Every row answers, for
# one cohort and one pooling family, what the endpoint actually means and - the
# column the analysis uses - which pooling class the row belongs to:
#
#   Exact-equivalent   the token IS the family's canonical definition
#                      (OS in OS, DFS in DFS, DSS in DSS, PFS in PFS, MFS in MFS)
#   Clinically-related a different token pooled into the family by the
#                      documented token -> family map (RFS/EFS/DFI in DFS,
#                      CSS/BCSS in DSS, PFI in PFS, DRFS in MFS)
#   Not-poolable       the cohort has an endpoint in this family but it cannot
#                      be pooled
#   Unknown            the deposit does not state what the endpoint means
#   Absent             the cohort has no endpoint of this family at all, so
#                      there is nothing to judge - deliberately NOT a pooling
#                      verdict, and never a candidate for pooling
#
# The companion file is a full grid (one row per cohort x family). Where it
# records an absent cell as "Not-poolable" with no token ("not stated"), the
# accessors report "Absent" instead, because the absence of an endpoint is not a
# pooling verdict; a file that already uses "Absent" is read identically.
#
# When the companion file is not on disk yet, the pooling class is derived
# mechanically from the token -> family map that the package has always used
# (.cpas_family_tokens) - the same rule the spec gives for the exact/related
# split - and labelled as a fallback. Nothing is invented: a family whose token
# is unknown is reported as "Unknown", never guessed.

.cpas_pooling_class_levels <- c("Exact-equivalent", "Clinically-related",
                                "Not-poolable", "Unknown", "Absent")

# Values that mean "no token recorded" in the companion table
.cpas_no_token <- function(x) {
  x <- trimws(as.character(x))
  is.na(x) | !nzchar(x) | tolower(x) %in% c("not stated", "na", "none", "unknown token")
}

.cpas_semantics_columns <- c("Accession", "Type", "Family", "Token", "TokenRole",
                             "SourceField", "EventDefinition", "TimeOrigin",
                             "CensoringRule", "CompetingEvents", "Derived",
                             "PoolingClass", "Evidence", "Note")

# token -> pooling family, derived from the map the package already uses
.cpas_token_family <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    tk <- unlist(.cpas_family_tokens[.cpas_pooling_families], use.names = FALSE)
    fm <- rep(.cpas_pooling_families,
              times = vapply(.cpas_family_tokens[.cpas_pooling_families],
                             length, integer(1)))
    cache <<- stats::setNames(fm, tk)
    cache
  }
})

#' @title Path of the shipped endpoint-semantics companion table
#' @description Returns the path of \code{inst/extdata/endpoint_semantics.csv}
#' inside the installed package (empty string when the file is not shipped).
#' @return Character scalar.
#' @keywords internal
#' @noRd
.cpas_semantics_path <- function() {
  p <- system.file("extdata", "endpoint_semantics.csv", package = "CanPAS")
  if (length(p) != 1L || is.na(p)) "" else p
}

.cpas_semantics_read <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path) || !file.exists(path))
    stop("The endpoint-semantics companion table is not available.\n",
         "  expected at: ", if (is.null(path) || !length(path)) "<none>" else path, "\n",
         "  It is produced by the curation pipeline as ",
         "data/endpoint_semantics.csv and shipped in the package as ",
         "inst/extdata/endpoint_semantics.csv.\n",
         "  Either regenerate it (pipeline step 11_endpoint_families.R) and ",
         "rebuild/reinstall CanPAS, or pass the table yourself through the ",
         "class_table= / path= argument.", call. = FALSE)
  tb <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(tb))
    stop("The endpoint-semantics companion table at ", path,
         " has no rows.", call. = FALSE)
  miss <- setdiff(c("Accession", "Family", "Token", "PoolingClass"), colnames(tb))
  if (length(miss))
    stop("The endpoint-semantics companion table at ", path,
         " is missing the required column(s): ", paste(miss, collapse = ", "),
         ". Expected the frozen schema: ",
         paste(.cpas_semantics_columns, collapse = ", "), ".", call. = FALSE)
  tb$Accession <- as.character(tb$Accession)
  tb$Family <- as.character(tb$Family)
  tb$Token <- as.character(tb$Token)
  tb$PoolingClass <- trimws(as.character(tb$PoolingClass))
  bad <- unique(tb$PoolingClass[!is.na(tb$PoolingClass) &
                                !tb$PoolingClass %in% .cpas_pooling_class_levels])
  if (length(bad))
    stop("The endpoint-semantics companion table at ", path,
         " uses PoolingClass value(s) outside the frozen vocabulary (",
         paste(.cpas_pooling_class_levels, collapse = " / "), "): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  if (any(is.na(tb$PoolingClass) | !nzchar(tb$PoolingClass)))
    stop("The endpoint-semantics companion table at ", path,
         " has blank PoolingClass entries (",
         sum(is.na(tb$PoolingClass) | !nzchar(tb$PoolingClass)),
         " row(s)); every cohort x family cell must carry a class.",
         call. = FALSE)
  dup <- duplicated(paste(tb$Accession, tb$Family, sep = "\r"))
  if (any(dup))
    stop("The endpoint-semantics companion table at ", path,
         " has duplicate cohort x family rows: ",
         paste(unique(paste(tb$Accession, tb$Family, sep = " / ")[dup]),
               collapse = ", "), ".", call. = FALSE)
  tb
}

#' @title Endpoint semantics of the curated cohorts
#' @description Returns the companion table that records, for every cohort and
#' pooling family, what the endpoint actually measures and which pooling class
#' it belongs to. The table is shipped with the package
#' (\code{inst/extdata/endpoint_semantics.csv}) and is a byte-identical copy of
#' the curation pipeline's \code{data/endpoint_semantics.csv}.
#'
#' \code{\link{cpas_meta}} uses it to decide what \code{pooling = "exact"}
#' keeps; \code{\link{endpoint_pooling_class}} is the one-cell accessor.
#' @param accession Optional accession filter (e.g. \code{"GSE31210"}).
#' @param family Optional pooling-family filter (\code{OS}, \code{DSS},
#'   \code{DFS}, \code{PFS}, \code{MFS}).
#' @param path Optional path of the companion CSV; defaults to the shipped copy.
#' @return data.frame with the frozen columns \code{Accession}, \code{Type},
#'   \code{Family}, \code{Token}, \code{TokenRole}, \code{SourceField},
#'   \code{EventDefinition}, \code{TimeOrigin}, \code{CensoringRule},
#'   \code{CompetingEvents}, \code{Derived}, \code{PoolingClass},
#'   \code{Evidence}, \code{Note}. Errors when the shipped table is absent (it
#'   is produced by the curation pipeline, not by the package).
#' @seealso \code{\link{endpoint_pooling_class}}, \code{\link{cpas_meta}}
#' @export
#' @examples
#' \dontrun{
#'   ## One row per cohort x endpoint family, with the pooling class.
#'   es <- endpoint_semantics(family = "DFS")
#'   table(es$PoolingClass)
#'
#'   ## The class of one cell, the value cpas_meta(pooling = "exact") uses.
#'   endpoint_pooling_class("GSE31210", "DFS")   # "Clinically-related" (RFS)
#'   endpoint_pooling_class("GSE13507", "OS")    # "Exact-equivalent"
#'   endpoint_pooling_class("GSE13507", "DFS")   # "Absent" (no DFS endpoint)
#' }
endpoint_semantics <- function(accession = NULL, family = NULL, path = NULL) {
  tb <- .cpas_semantics_read(if (is.null(path)) .cpas_semantics_path() else path)
  if (!is.null(accession)) {
    accession <- as.character(accession)
    tb <- tb[tb$Accession %in% accession, , drop = FALSE]
  }
  if (!is.null(family)) {
    fam <- vapply(as.character(family), function(f) {
      x <- endpoint_family(f)
      if (is.na(x)) as.character(f) else x
    }, character(1), USE.NAMES = FALSE)
    tb <- tb[tb$Family %in% fam, , drop = FALSE]
  }
  rownames(tb) <- NULL
  tb
}

# Derive the pooling class from the token -> family map alone. This is the
# documented rule of the spec and the fallback used while the companion table
# is not on disk; `token` is the token the cohort actually contributes. No token
# means the family is absent for that cohort, which is not a pooling verdict.
.cpas_pooling_class_from_token <- function(token, family) {
  fam <- endpoint_family(family)
  if (is.na(fam)) fam <- as.character(family)[1]
  if (is.na(token) || !length(token) || !nzchar(as.character(token)))
    return("Absent")
  token <- as.character(token)
  if (identical(token, fam)) return("Exact-equivalent")
  tf <- .cpas_token_family()
  if (token %in% names(tf) && identical(unname(tf[[token]]), fam))
    return("Clinically-related")
  if (token %in% names(tf)) return("Not-poolable")   # a token of another family
  "Unknown"
}

# One row per cohort: token, family, class, and where the class came from.
# `table` is the companion table when it is available (or was supplied); the
# token -> family fallback is used for any cohort the table does not cover.
# `tokens` optionally overrides the resolved token (cpas_meta passes the token
# it actually used, which for a cohort outside the catalog can only come from
# the data at hand).
.cpas_pooling_lookup <- function(accessions, families, table = NULL,
                                 path = NULL, tokens = NULL) {
  if (is.null(table)) {
    p <- if (is.null(path)) .cpas_semantics_path() else path
    table <- tryCatch(.cpas_semantics_read(p), error = function(e) NULL)
  }
  accessions <- as.character(accessions)
  families <- vapply(as.character(families), function(f) {
    x <- endpoint_family(f); if (is.na(x)) as.character(f) else x
  }, character(1), USE.NAMES = FALSE)
  if (is.null(tokens)) {
    tokens <- vapply(seq_along(accessions), function(i)
      endpoint_resolve(accessions[i], families[i]), character(1))
    # a cohort that is not in the catalog may still carry the raw token
    unknown <- is.na(tokens)
    tokens[unknown] <- families[unknown]
  }
  tokens <- as.character(tokens)
  cls <- rep(NA_character_, length(accessions))
  src <- rep("token-map fallback", length(accessions))
  role <- ifelse(tokens == families, "primary", "contributing")
  if (!is.null(table)) {
    key <- paste(table$Accession, table$Family, sep = "\r")
    hit <- match(paste(accessions, families, sep = "\r"), key)
    idx <- which(!is.na(hit))
    if (length(idx)) {
      cls[idx] <- as.character(table$PoolingClass[hit[idx]])
      src[idx] <- "companion table"
      tk <- if ("Token" %in% colnames(table)) as.character(table$Token[hit[idx]]) else rep(NA_character_, length(idx))
      tr <- if ("TokenRole" %in% colnames(table)) as.character(table$TokenRole[hit[idx]]) else rep(NA_character_, length(idx))
      has_tok <- !.cpas_no_token(tk)
      # no token recorded -> the family is absent for this cohort, which is not
      # a pooling verdict even where the table spells it "Not-poolable"
      absent <- !has_tok | .cpas_no_token(tr)
      cls[idx[absent]] <- "Absent"
      if (any(has_tok)) tokens[idx[has_tok]] <- tk[has_tok]
      real_role <- !absent & !.cpas_no_token(tr)
      if (any(real_role)) role[idx[real_role]] <- tr[real_role]
      if (any(absent)) { tokens[idx[absent]] <- NA_character_; role[idx[absent]] <- "not stated" }
    }
  }
  miss <- which(is.na(cls) | !nzchar(cls))
  if (length(miss))
    cls[miss] <- vapply(miss, function(i)
      .cpas_pooling_class_from_token(tokens[i], families[i]), character(1))
  data.frame(accession = accessions, family = families, token = tokens,
             token_role = role, pooling_class = cls,
             pooling_class_source = src, stringsAsFactors = FALSE)
}

#' @title Pooling class of one cohort x endpoint-family cell
#' @description Returns the \code{PoolingClass} of one cohort and endpoint
#' family:
#' \describe{
#'   \item{\code{"Exact-equivalent"}}{the token IS the family's canonical
#'     definition (OS in OS, DFS in DFS, ...)}
#'   \item{\code{"Clinically-related"}}{a different token pooled into the family
#'     by the documented token -> family map (RFS/EFS/DFI in DFS, CSS/BCSS in
#'     DSS, PFI in PFS, DRFS in MFS)}
#'   \item{\code{"Not-poolable"}}{the cohort reports an endpoint in this family
#'     that cannot be pooled}
#'   \item{\code{"Unknown"}}{the deposit does not state what the endpoint means}
#'   \item{\code{"Absent"}}{the cohort has no endpoint of this family at all.
#'     This is the absence of an endpoint, \strong{not} a pooling verdict, and
#'     such a cell is never a candidate for pooling. Where the companion table
#'     spells an absent cell \code{"Not-poolable"} with no token
#'     (\code{TokenRole = "not stated"}), this accessor reports
#'     \code{"Absent"} instead; a table that already uses \code{"Absent"} is
#'     read identically.}
#' }
#' The value comes from the shipped companion table
#' (\code{\link{endpoint_semantics}}); when that table is not available the class
#' is derived from the package's documented token -> family map, which is the
#' same rule the table records for the exact/related split.
#' @param accession Cohort accession (e.g. \code{"GSE31210"}).
#' @param family Pooling family (\code{OS}, \code{DSS}, \code{DFS}, \code{PFS},
#'   \code{MFS}) or a raw endpoint token.
#' @param class_table Optional companion table to use instead of the shipped one.
#' @return Character scalar: one of \code{"Exact-equivalent"},
#'   \code{"Clinically-related"}, \code{"Not-poolable"}, \code{"Unknown"} or
#'   \code{"Absent"}.
#' @seealso \code{\link{endpoint_semantics}}, \code{\link{cpas_meta}}
#' @export
#' @examples
#'   ## The token a cohort contributes decides whether it is the canonical
#'   ## definition of the family or a related one pooled into it by design.
#'   endpoint_pooling_class("GSE31210", "DFS")       # "Clinically-related" (RFS)
#'   endpoint_pooling_class("GSE4922_GPL96", "DFS")  # "Exact-equivalent"
#'   endpoint_pooling_class("GSE13507", "MFS")       # "Absent"
endpoint_pooling_class <- function(accession, family, class_table = NULL) {
  if (length(accession) != 1L || length(family) != 1L)
    stop("endpoint_pooling_class() takes one accession and one family.",
         call. = FALSE)
  .cpas_pooling_lookup(accession, family, table = class_table)$pooling_class
}

# TRUE when the companion table is shipped (used by the manifest and the tests)
.cpas_semantics_available <- function() nzchar(.cpas_semantics_path())
