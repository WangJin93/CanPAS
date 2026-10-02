# cpas_manifest.R ------------------------------------------------------------
# The analysis manifest (workstream P, spec B1.4).
#
# Primary analyses carry a manifest list element or plotting-object attribute.
# cpas_manifest() retrieves either representation. The manifest answers, in one object, the questions a reviewer asks about a
# reported number: which cohorts, which endpoint token, which pooling classes,
# how the cohorts were chosen, what was dropped and why, which cut-point rule was
# used and how many were searched, whether the search was adjusted for, what the
# proportional-hazards check said, which meta-analysis method produced the pooled
# estimate, tau^2 / I^2 / the prediction intervals, and the versions and time of
# the run. It prints as a human-readable block and converts with as.data.frame()
# to a field/value table.

.cpas_require <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

.cpas_versions <- function() {
  c(R = paste(R.version$major, R.version$minor, sep = "."),
    CanPAS = .cpas_require("CanPAS"),
    survival = .cpas_require("survival"),
    metafor = .cpas_require("metafor"))
}

.cpas_empty_df <- function(cols) {
  d <- as.data.frame(stats::setNames(replicate(length(cols), character(0),
                                               simplify = FALSE), cols),
                     stringsAsFactors = FALSE)
  d
}

.cpas_hash_analyzed_inputs <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 3)
  unname(tools::md5sum(path))
}

# Build a manifest from named fields; anything not supplied is filled in with a
# neutral value so that the field list is identical for every analysis.
.cpas_manifest_new <- function(analysis, cohorts = character(0),
                               n_cohorts = length(cohorts),
                               family = NA_character_, token = NA_character_,
                               tokens = NULL, token_role = NA_character_,
                               accession = NA_character_, raw_token = NA_character_,
                               endpoint_evidence = "not supplied; schema only",
                               pooling = NA_character_,
                               pooling_class = NA_character_,
                               n_input = NA_integer_, n_analyzed = NA_integer_,
                               n_excluded = NA_integer_, events = NA_integer_,
                               marker_requested = NA_character_,
                               fitted_covariates = character(0),
                               marker_definition = NA_character_, covariates_requested = character(0),
                               estimator = NA_character_, inference = NA_character_,
                               dataset_hash = NA_character_, hash_scope = NA_character_,
                               seed = NA_integer_, seed_status = "not recorded; deterministic procedure not asserted",
                               analyzed_data = NULL, input_data = NULL,
                               primary_path = analysis,
                               coverage_status = "schema_only",
                               pooling_classes = NA_character_,
                               pooling_table = NULL,
                               overlap_mode = NA_character_,
                               overlap_pairs = NULL, overlap_dropped = NULL,
                               overlap_source = NA_character_,
                               selection_rule = NA_character_,
                               dropped_rows = NULL, dropped_covariates = NULL,
                               cut_rule = NA_character_,
                               cut_points_searched = NA_integer_,
                               search_adjusted_p = NA_real_,
                               ph_test = NULL,
                               meta_method = NA_character_,
                               auto_repair = NA,
                               pi_method = NA_character_,
                               tau2 = NA_real_, I2 = NA_real_,
                               pi_primary = NULL, pi_primary_rule = NA_character_,
                               pi_alt = NULL, pi_alt_rule = NA_character_,
                               notes = character(0)) {
  if (is.null(tokens)) tokens <- stats::setNames(rep(NA_character_, length(cohorts)), cohorts)
  if (is.null(overlap_pairs))
    overlap_pairs <- .cpas_empty_df(c("AccessionA", "AccessionB", "SharedPatients", "Basis"))
  if (is.null(overlap_dropped))
    overlap_dropped <- .cpas_empty_df(c("cohort", "kept", "shared_patients", "reason"))
  if (is.null(dropped_rows))
    dropped_rows <- .cpas_empty_df(c("cohort", "n_dropped", "reason"))
  if (is.null(dropped_covariates))
    dropped_covariates <- .cpas_empty_df(c("cohort", "variable", "detail", "reason"))
  if (is.null(pooling_table))
    pooling_table <- .cpas_empty_df(c("accession", "family", "token", "token_role",
                                      "pooling_class", "pooling_class_source"))
  if (is.null(pi_primary)) pi_primary <- c(lower = NA_real_, upper = NA_real_)
  if (is.null(pi_alt)) pi_alt <- c(lower = NA_real_, upper = NA_real_)
  if (!is.null(input_data) && is.na(n_input)[1]) n_input <- nrow(input_data)
  if (!is.null(analyzed_data) && is.na(n_analyzed)[1]) n_analyzed <- nrow(analyzed_data)
  if (is.na(n_excluded)[1] && is.finite(n_input)[1] && is.finite(n_analyzed)[1])
    n_excluded <- n_input - n_analyzed
  if (is.na(dataset_hash)[1] && !is.null(analyzed_data)) {
    dataset_hash <- .cpas_hash_analyzed_inputs(analyzed_data)
    hash_scope <- "serialized analyzed_data RDS (not upstream raw dataset)"
  }
  out <- list(
    analysis = analysis,
    cohorts = as.character(cohorts),
    n_cohorts = as.integer(n_cohorts),
    family = as.character(family),
    token = as.character(token),
    tokens = tokens,
    token_role = as.character(token_role),
    accession = as.character(accession),
    raw_token = as.character(raw_token),
    endpoint_evidence = as.character(endpoint_evidence),
    pooling = as.character(pooling),
    pooling_class = as.character(pooling_class),
    pooling_classes = as.character(pooling_classes),
    n_input = as.integer(n_input),
    n_analyzed = as.integer(n_analyzed),
    n_excluded = as.integer(n_excluded),
    events = as.integer(events),
    marker_requested = as.character(marker_requested),
    marker_definition = as.character(marker_definition),
    covariates_requested = as.character(covariates_requested),
    fitted_covariates = as.character(fitted_covariates),
    estimator = as.character(estimator),
    inference = as.character(inference),
    dataset_hash = as.character(dataset_hash),
    hash_scope = as.character(hash_scope),
    seed = as.integer(seed),
    seed_status = as.character(seed_status),
    primary_path = as.character(primary_path),
    coverage_status = as.character(coverage_status),
    pooling_table = pooling_table,
    overlap_mode = as.character(overlap_mode),
    overlap_source = as.character(overlap_source),
    overlap_pairs = overlap_pairs,
    overlap_dropped = overlap_dropped,
    selection_rule = as.character(selection_rule),
    dropped_rows = dropped_rows,
    dropped_covariates = dropped_covariates,
    cut_rule = as.character(cut_rule),
    cut_points_searched = as.integer(cut_points_searched),
    search_adjusted_p = as.numeric(search_adjusted_p),
    ph_test = ph_test,
    meta_method = as.character(meta_method),
    auto_repair = auto_repair,
    pi_method = as.character(pi_method),
    tau2 = as.numeric(tau2),
    I2 = as.numeric(I2),
    pi_primary = pi_primary,
    pi_primary_rule = as.character(pi_primary_rule),
    pi_alt = pi_alt,
    pi_alt_rule = as.character(pi_alt_rule),
    versions = .cpas_versions(),
    timestamp = Sys.time(),
    notes = as.character(notes))
  class(out) <- "cpas_manifest"
  out
}

# One-line rendering of a table field for as.data.frame()/print()
.cpas_manifest_collapse <- function(x, key = NULL, value = "reason", n = 6L) {
  if (is.null(x) || !is.data.frame(x) || !nrow(x)) return("")
  if (!is.null(key) && !key %in% colnames(x)) key <- NULL
  if (!value %in% colnames(x)) value <- colnames(x)[ncol(x)]
  rows <- if (is.null(key)) as.character(x[[value]]) else
    sprintf("%s: %s", as.character(x[[key]]), as.character(x[[value]]))
  rows <- rows[!is.na(rows) & nzchar(rows)]
  if (length(rows) > n)
    rows <- c(rows[seq_len(n)], sprintf("... (+%d more)", length(rows) - n))
  paste(rows, collapse = " | ")
}

.cpas_manifest_num <- function(x, digits = 6) {
  if (is.null(x) || !length(x) || all(is.na(x))) return("NA")
  paste(ifelse(is.na(x), "NA", formatC(as.numeric(x), digits = digits, format = "g")),
        collapse = ", ")
}

# field/value pairs, in a fixed order, for as.data.frame() and print()
.cpas_manifest_pairs <- function(m) {
  pi1 <- m$pi_primary; pi2 <- m$pi_alt
  f <- function(label, value) data.frame(field = label, value = as.character(value),
                                         stringsAsFactors = FALSE)
  toks <- if (length(m$tokens) && any(!is.na(m$tokens)))
    paste(sprintf("%s=%s", names(m$tokens), m$tokens), collapse = ", ") else "NA"
  ph <- if (is.null(m$ph_test)) "not computed" else {
    tb <- m$ph_test$table
    if (is.null(tb) || !"p" %in% colnames(tb)) "not computed" else {
      g <- if ("term" %in% colnames(tb)) which(as.character(tb$term) == "GLOBAL") else integer(0)
      if (length(g))
        sprintf("GLOBAL p = %s (%d cohort x term check(s))",
                format.pval(tb$p[g[1]], digits = 3), max(0L, nrow(tb) - length(g)))
      else sprintf("per-term only (%d check(s), no GLOBAL row)", nrow(tb))
    }
  }
  do.call(rbind, list(
    f("analysis", m$analysis),
    f("cohorts", if (length(m$cohorts)) paste(m$cohorts, collapse = ", ") else "NA"),
    f("n_cohorts", m$n_cohorts),
    f("endpoint_family", m$family),
    f("resolved_token", m$token),
    f("resolved_tokens_per_cohort", toks),
    f("token_role", m$token_role),
    f("accession", m$accession),
    f("raw_token", m$raw_token),
    f("endpoint_evidence", m$endpoint_evidence),
    f("pooling_mode", m$pooling),
    f("pooling_class", m$pooling_class),
    f("n_input", m$n_input),
    f("n_analyzed", m$n_analyzed),
    f("n_excluded", m$n_excluded),
    f("events", m$events),
    f("marker_requested", m$marker_requested),
    f("marker_definition", m$marker_definition),
    f("covariates_requested", if (length(m$covariates_requested)) paste(m$covariates_requested, collapse = ", ") else "NA"),
    f("fitted_covariates", if (length(m$fitted_covariates)) paste(m$fitted_covariates, collapse = ", ") else "NA"),
    f("estimator", m$estimator),
    f("inference", m$inference),
    f("dataset_hash", m$dataset_hash),
    f("hash_scope", m$hash_scope),
    f("seed", m$seed),
    f("seed_status", m$seed_status),
    f("primary_path", m$primary_path),
    f("coverage_status", m$coverage_status),
    f("pooling_classes_pooled",
      if (length(m$pooling_classes) && any(!is.na(m$pooling_classes)))
        paste(m$pooling_classes, collapse = ", ") else "NA"),
    f("pooling_class_source",
      if (nrow(m$pooling_table)) paste(unique(m$pooling_table$pooling_class_source),
                                       collapse = ", ") else "NA"),
    f("overlap_mode", m$overlap_mode),
    f("overlap_source", m$overlap_source),
    f("overlap_pairs", if (nrow(m$overlap_pairs))
      .cpas_manifest_collapse(m$overlap_pairs, key = "AccessionA", value = "AccessionB") else "none"),
    f("overlap_dropped", if (nrow(m$overlap_dropped))
      .cpas_manifest_collapse(m$overlap_dropped, key = "cohort") else "none"),
    f("selection_rule", m$selection_rule),
    f("rows_dropped", if (nrow(m$dropped_rows))
      .cpas_manifest_collapse(m$dropped_rows, key = "cohort") else "none"),
    f("n_rows_dropped_reported", if (nrow(m$dropped_rows))
      sum(suppressWarnings(as.numeric(m$dropped_rows$n_dropped)), na.rm = TRUE) else 0),
    f("covariates_dropped", if (nrow(m$dropped_covariates))
      .cpas_manifest_collapse(m$dropped_covariates, key = "variable") else "none"),
    f("cut_rule", m$cut_rule),
    f("cut_points_searched", m$cut_points_searched),
    f("search_adjusted_p", .cpas_manifest_num(m$search_adjusted_p)),
    f("ph_test_cox_zph", ph),
    f("meta_method", m$meta_method),
    f("auto_repair", if (is.na(m$auto_repair)) "NA" else if (isTRUE(m$auto_repair)) "TRUE" else "FALSE"),
    f("pi_method", if (is.na(m$pi_method)) "NA" else m$pi_method),
    f("tau2", .cpas_manifest_num(m$tau2)),
    f("I2", .cpas_manifest_num(m$I2)),
    f("pi_primary", sprintf("[%s, %s]", .cpas_manifest_num(pi1[["lower"]]),
                            .cpas_manifest_num(pi1[["upper"]]))),
    f("pi_primary_rule", m$pi_primary_rule),
    f("pi_alt", sprintf("[%s, %s]", .cpas_manifest_num(pi2[["lower"]]),
                        .cpas_manifest_num(pi2[["upper"]]))),
    f("pi_alt_rule", m$pi_alt_rule),
    f("versions", paste(sprintf("%s %s", names(m$versions), m$versions), collapse = "; ")),
    f("timestamp", format(m$timestamp, "%Y-%m-%d %H:%M:%S %Z")),
    f("notes", if (length(m$notes)) paste(m$notes, collapse = " | ") else "")))
}

#' @title Analysis manifest of a CanPAS result
#' @description Returns the analysis manifest of a result: the cohorts, the
#' endpoint family and the token each cohort actually contributed, the pooling
#' classes that were pooled, how the cohorts were selected, which rows and
#' covariates were dropped and why, the cut-point rule and how many cut-points
#' were searched, the search-adjusted p-value when one was computed, the
#' proportional-hazards (\code{cox.zph}) result, the meta-analysis method,
#' tau\eqn{^2}, \eqn{I^2}, the prediction interval(s), the R and package
#' versions, and the timestamp.
#'
#' Primary analysis paths return the common schema as a \code{$manifest}
#' element or a plotting-object \code{manifest} attribute; this accessor reads both.
#' Input/analysis counts and serialized-input hashes are recorded where measured.
#' Upstream accession, endpoint evidence and caller seed remain explicitly unknown
#' when they cannot be established from supplied data. Field presence is distinct
#' from complete source provenance; \code{hash_scope} describes the hashed input. The manifest prints as a
#' readable block and converts with \code{as.data.frame()} into a two-column
#' \code{field}/\code{value} table; the structured pieces (per-cohort pooling
#' classes, dropped rows/covariates, overlapping pairs) are data.frames inside
#' the manifest itself.
#' @param x A result object returned by an analysis function
#'   (\code{\link{cpas_meta}}, \code{\link{COX_analysis}},
#'   \code{\link{COX_screen_adjust}}, \code{\link{cpas_km_pooled}}, ...). When
#'   omitted, an empty manifest skeleton with the same field list is returned.
#' @param ... Unused; accepted so that \code{cpas_manifest(result)} and
#'   \code{cpas_manifest()} take the same shape as the other accessors.
#' @return An object of class \code{cpas_manifest} (a named list).
#' @seealso \code{\link{cpas_meta}}
#' @export
#' @examples
#' \dontrun{
#'   m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
#'   cpas_manifest(m)                 # readable block
#'   as.data.frame(cpas_manifest(m))  # field / value table
#'   cpas_manifest(m)$pooling_classes # what was pooled
#' }
cpas_manifest <- function(x = NULL, ...) {
  if (is.null(x))
    return(.cpas_manifest_new("(none)", notes = "empty manifest skeleton"))
  if (inherits(x, "cpas_manifest")) return(x)
  if (is.list(x) && !is.null(x$manifest) && inherits(x$manifest, "cpas_manifest"))
    return(x$manifest)
  am <- attr(x, "manifest", exact = TRUE)
  if (inherits(am, "cpas_manifest")) return(am)
  stop("cpas_manifest() needs a CanPAS analysis result carrying a $manifest ",
       "element, or no argument at all for an empty skeleton; got ",
       paste(class(x), collapse = "/"), ".", call. = FALSE)
}

#' @title Print an analysis manifest
#' @description Prints the manifest as an aligned \code{field: value} block.
#' @param x An object of class \code{cpas_manifest}.
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_manifest <- function(x, ...) {
  if (!inherits(x, "cpas_manifest")) x <- cpas_manifest(x)
  cat("CanPAS analysis manifest\n")
  cat("========================\n")
  df <- .cpas_manifest_pairs(x)
  w <- max(nchar(df$field))
  for (i in seq_len(nrow(df)))
    cat(sprintf("%-*s : %s\n", w, df$field[i], df$value[i]))
  if (nrow(x$pooling_table))
    cat(sprintf("\nper-cohort pooling classes (%d row(s)): $manifest$pooling_table\n",
                nrow(x$pooling_table)))
  invisible(x)
}

#' @title Convert an analysis manifest to a data.frame
#' @description Flattens the manifest into a two-column \code{field}/\code{value}
#' table so it can be written to CSV, compared between runs, or embedded in a
#' report. Table-valued fields are collapsed into one row each; the structured
#' tables remain in the manifest itself (\code{$pooling_table},
#' \code{$dropped_rows}, \code{$dropped_covariates}, \code{$overlap_pairs}).
#' @param x An object of class \code{cpas_manifest}.
#' @param row.names,optional,... Passed to \code{as.data.frame} for compatibility.
#' @return data.frame with columns \code{field} and \code{value}.
#' @export
as.data.frame.cpas_manifest <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!inherits(x, "cpas_manifest")) x <- cpas_manifest(x)
  df <- .cpas_manifest_pairs(x)
  if (!is.null(row.names)) rownames(df) <- row.names
  df
}
