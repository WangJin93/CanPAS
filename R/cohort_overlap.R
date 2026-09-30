# cohort_overlap.R -----------------------------------------------------------
# Shared-patient overlap as a shipped object (workstream P, spec A6).
#
# Sharing patients between two cohorts is the one data property that can make a
# pooled estimate look precise when it is merely double-counted. The curation
# side ships the register as data/suppl/cohort_overlap.csv; the package ships a
# byte-identical copy as inst/extdata/cohort_overlap.csv and reads it here.
#
# Frozen schema: AccessionA, AccessionB, SharedPatients, Basis, Evidence.

.cpas_overlap_columns <- c("AccessionA", "AccessionB", "SharedPatients",
                           "Basis", "Evidence")

.cpas_overlap_path <- function() {
  p <- system.file("extdata", "cohort_overlap.csv", package = "CanPAS")
  if (length(p) != 1L || is.na(p)) "" else p
}

.cpas_overlap_read <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path) || !file.exists(path))
    stop("The cohort-overlap companion table is not available.\n",
         "  expected at: ", if (is.null(path) || !length(path)) "<none>" else path, "\n",
         "  It is produced by the curation pipeline as ",
         "data/suppl/cohort_overlap.csv and shipped in the package as ",
         "inst/extdata/cohort_overlap.csv.\n",
         "  Either regenerate it (pipeline step 26_cohort_overlap.R) and ",
         "rebuild/reinstall CanPAS, or pass the pairs yourself through the ",
         "overlap_table= / path= argument.", call. = FALSE)
  tb <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  miss <- setdiff(c("AccessionA", "AccessionB", "SharedPatients"), colnames(tb))
  if (length(miss))
    stop("The cohort-overlap companion table at ", path,
         " is missing the required column(s): ", paste(miss, collapse = ", "),
         ". Expected the frozen schema: ",
         paste(.cpas_overlap_columns, collapse = ", "), ".", call. = FALSE)
  for (cn in c("Basis", "Evidence"))
    if (!cn %in% colnames(tb)) tb[[cn]] <- NA_character_
  tb <- tb[, .cpas_overlap_columns, drop = FALSE]
  tb$AccessionA <- as.character(tb$AccessionA)
  tb$AccessionB <- as.character(tb$AccessionB)
  tb$SharedPatients <- suppressWarnings(as.integer(tb$SharedPatients))
  tb
}

#' @title Shared-patient overlap register
#' @description Returns the register of cohort pairs that share patients, with
#' the number of shared patients, the basis of the finding and the evidence it
#' was read from. The table is shipped with the package
#' (\code{inst/extdata/cohort_overlap.csv}) and is a byte-identical copy of the
#' curation pipeline's \code{data/suppl/cohort_overlap.csv}.
#'
#' \code{\link{cpas_meta}} consults it through its \code{overlap} argument.
#' @param path Optional path of the register CSV; defaults to the shipped copy.
#' @return data.frame with the frozen columns \code{AccessionA},
#'   \code{AccessionB}, \code{SharedPatients}, \code{Basis}, \code{Evidence}.
#'   Errors when the shipped table is absent (it is produced by the curation
#'   pipeline, not by the package).
#' @seealso \code{\link{cpas_meta}}
#' @export
#' @examples
#' \dontrun{
#'   ov <- cohort_overlap()
#'   head(ov)
#'   ## Which of my cohorts overlap at all?
#'   ov[ov$AccessionA %in% c("GSE11969", "GSE13213") |
#'      ov$AccessionB %in% c("GSE11969", "GSE13213"), ]
#' }
cohort_overlap <- function(path = NULL) {
  tb <- .cpas_overlap_read(if (is.null(path)) .cpas_overlap_path() else path)
  rownames(tb) <- NULL
  tb
}

# Overlap rows restricted to a set of cohorts, as an undirected pair list. The
# return value always carries a data.frame (0 rows when the register is not
# shipped), so callers never have to test for NULL.
.cpas_overlap_pairs_in <- function(cohorts, table = NULL, path = NULL) {
  if (is.null(table)) {
    p <- if (is.null(path)) .cpas_overlap_path() else path
    table <- tryCatch(.cpas_overlap_read(p), error = function(e) NULL)
  }
  if (is.null(table))
    return(list(available = FALSE,
                table = .cpas_empty_df(c("AccessionA", "AccessionB",
                                         "SharedPatients", "Basis", "Evidence"))))
  cohorts <- as.character(cohorts)
  keep <- table$AccessionA %in% cohorts & table$AccessionB %in% cohorts
  out <- table[keep, , drop = FALSE]
  rownames(out) <- NULL
  list(available = TRUE, table = out)
}

# Connected components of the overlap graph over `cohorts` (deterministic)
.cpas_overlap_groups <- function(pairs) {
  if (is.null(pairs) || !nrow(pairs)) return(list())
  nodes <- unique(c(pairs$AccessionA, pairs$AccessionB))
  parent <- stats::setNames(nodes, nodes)
  find <- function(x) { while (parent[[x]] != x) x <- parent[[x]]; x }
  for (i in seq_len(nrow(pairs))) {
    a <- find(pairs$AccessionA[i]); b <- find(pairs$AccessionB[i])
    if (a != b) parent[[a]] <- b
  }
  root <- vapply(nodes, find, character(1))
  grp <- split(nodes, root)
  grp[vapply(grp, length, integer(1)) > 1L]
}

# Deterministic dedupe: within every overlapping group keep the member with the
# larger N, ties broken by accession sort. Returns the kept cohorts and a table
# of what was dropped (and against which kept cohort it overlapped).
.cpas_overlap_dedupe <- function(cohorts, pairs, n_of) {
  grp <- .cpas_overlap_groups(pairs)
  cohorts <- as.character(cohorts)
  drop <- data.frame(cohort = character(0), kept = character(0),
                     shared_patients = integer(0), group = character(0),
                     reason = character(0), stringsAsFactors = FALSE)
  keep <- cohorts
  for (g in grp) {
    members <- sort(intersect(cohorts, g))
    if (length(members) < 2L) next
    n <- vapply(members, function(m) {
      v <- n_of[[m]]; if (is.null(v) || !length(v) || is.na(v)) NA_real_ else as.numeric(v)
    }, numeric(1))
    ord <- order(-n, members)                    # larger N first, ties by accession
    winner <- members[ord][1]
    for (loser in members[ord][-1]) {
      hit <- pairs[(pairs$AccessionA == loser & pairs$AccessionB == winner) |
                   (pairs$AccessionA == winner & pairs$AccessionB == loser), ,
                   drop = FALSE]
      drop <- rbind(drop, data.frame(
        cohort = loser, kept = winner,
        shared_patients = if (nrow(hit)) as.integer(hit$SharedPatients[1]) else NA_integer_,
        group = paste(sort(g), collapse = "|"),
        reason = sprintf("shares patients with %s (kept: larger N, %s vs %s)",
                         winner,
                         ifelse(is.na(n[[winner]]), "unknown", format(n[[winner]])),
                         ifelse(is.na(n[[loser]]), "unknown", format(n[[loser]]))),
        stringsAsFactors = FALSE))
    }
    keep <- setdiff(keep, setdiff(members, winner))
  }
  list(keep = keep, dropped = drop)
}
