# endpoint_adjudication.R ----------------------------------------------------
# Blinded endpoint adjudication (spec B1).
#
# The endpoint-semantics layer records, for every cohort x family, which token
# the cohort reports, how the event is defined and which pooling class we gave
# it.  That annotation is our own judgement, so it has to be auditable: two (or
# more) raters, blinded to our verdict, judge a sample of cohort-endpoint
# records from the token and the evidence text alone, and the agreement between
# them is reported.
#
# Two functions:
#   endpoint_adjudication()  the blinded rating-sheet generator (stratified
#                            sample by family and by source, stable record ids,
#                            blank per-rater columns)
#   endpoint_agreement()     the agreement calculator (raw agreement and
#                            Cohen's kappa per field, overall pooled agreement,
#                            adjudication rate, the disagreeing records)
#
# Cohen's kappa is implemented here with base R only: the package must not gain
# a dependency for it.  Nothing in this file contacts the network; the rating
# sheet is drawn from the shipped endpoint-semantics table.

# The five judgement fields a rater fills, in the order of the sheet.
.cpas_adjudication_fields <- c("event_definition", "time_origin",
                               "censoring_rule", "competing_events",
                               "pooling_class")

# Landis & Koch (1977) verbal categories for a kappa value.
.cpas_kappa_label <- function(k) {
  if (!length(k) || is.na(k)) return(NA_character_)
  if (k < 0) "poor (less than chance)"
  else if (k <= 0.20) "slight"
  else if (k <= 0.40) "fair"
  else if (k <= 0.60) "moderate"
  else if (k <= 0.80) "substantial"
  else "almost perfect"
}

# Data source of a cohort, from its accession prefix.  The endpoint-semantics
# table carries no source column, so the prefix is the documented proxy:
# GSE -> GEO, TCGA- -> TCGA (UCSC Xena), E- -> EMBL-EBI ArrayExpress, CGGA ->
# CGGA, anything else -> "other (publication / cBioPortal-hosted)".
.cpas_accession_source <- function(accession) {
  a <- as.character(accession)
  out <- rep("other (publication / cBioPortal-hosted)", length(a))
  out[grepl("^GSE", a)] <- "GEO"
  out[grepl("^TCGA", a)] <- "TCGA"
  out[grepl("^E-", a)] <- "EMBL-EBI"
  out[grepl("^CGGA", a)] <- "CGGA"
  out
}

# Cohen's kappa for one pair of nominal ratings.
#
#   po = sum(diagonal) / n          observed agreement
#   pe = sum(row_i * col_i) / n^2   chance agreement
#   kappa = (po - pe) / (1 - pe)
#
# Only pairs where BOTH ratings are present and non-empty take part; the number
# compared is returned so a caller can see how much of the sheet was usable.
# pe = 1 (every rating the same single category) leaves kappa undefined (0/0)
# and is reported as NA with po = 1, not as NaN.
.cpas_cohen_kappa <- function(a, b) {
  ok <- !is.na(a) & !is.na(b) & nzchar(a) & nzchar(b)
  a <- as.character(a[ok]); b <- as.character(b[ok])
  n <- length(a)
  if (n == 0L)
    return(list(kappa = NA_real_, po = NA_real_, pe = NA_real_, n = 0L,
                note = "no record was rated by both raters"))
  lv <- sort(unique(c(a, b)))
  ta <- table(factor(a, levels = lv)); tb <- table(factor(b, levels = lv))
  po <- sum(a == b) / n
  pe <- sum(as.numeric(ta) * as.numeric(tb)) / n^2
  if (is.na(pe) || pe >= 1)
    return(list(kappa = NA_real_, po = po, pe = pe, n = n,
                note = "chance agreement is 1 (a single category); kappa is undefined"))
  list(kappa = (po - pe) / (1 - pe), po = po, pe = pe, n = n, note = NA_character_)
}

# The rater-pair labels of a sheet: rater1/rater2/... in numeric order.
.cpas_adjudication_raters <- function(sheet, fields = .cpas_adjudication_fields) {
  cn <- colnames(sheet)
  ids <- integer(0)
  for (f in fields) {
    m <- regmatches(cn, regexpr(sprintf("^rater([0-9]+)_%s$", f), cn))
    if (length(m)) ids <- c(ids, as.integer(sub(sprintf("_(%s)$", f), "", sub("^rater", "", m))))
  }
  sort(unique(ids))
}

#' @title Blinded endpoint-adjudication rating sheet
#' @description
#' Draws a stratified sample of cohort-endpoint records from the shipped
#' endpoint-semantics table and returns a blinded rating sheet: one row per
#' record with a stable \code{record_id}, the cohort, the endpoint \strong{token}
#' and the \strong{evidence text} only, plus one empty column per judgement field
#' and rater for the raters to fill.
#'
#' The sheet deliberately does \strong{not} carry our own verdict: the
#' \code{Family} assignment, the \code{PoolingClass} and the four semantic fields
#' (\code{EventDefinition}, \code{TimeOrigin}, \code{CensoringRule},
#' \code{CompetingEvents}) are withheld, so a rater judges the endpoint from the
#' token and the evidence text alone and cannot anchor on the annotation under
#' review.
#' @param n Number of records to draw (default 80; capped at the number of
#'   eligible records).
#' @param raters Number of raters the sheet is prepared for (default 2). Any
#'   number >= 1 is accepted; the agreement calculator reports pairwise mean
#'   kappa when more than two raters filled it.
#' @param seed Fixed random seed (default \code{20260101}) making the draw
#'   reproducible. The previous RNG state is restored on exit.
#' @param families Optional character vector of endpoint families to sample from
#'   (default: all five).
#' @param drop_absent Drop records whose token is absent for the family
#'   (\code{default FALSE}; set \code{TRUE} to sample only records that report an
#'   endpoint).
#' @param path Optional file path. When supplied the sheet is also written there
#'   as CSV with empty cells for the blank ratings (and the path is recorded in
#'   \code{attr(x, "csv_path")}).
#' @param table Optional endpoint-semantics table (the frozen schema of
#'   \code{\link{endpoint_semantics}}) used instead of the shipped copy.
#' @return A data frame of class \code{cpas_adjudication_sheet} with columns
#'   \code{record_id}, \code{cohort}, \code{token}, \code{evidence} and, for
#'   every rater \code{r} and field \code{f},
#'   \code{rater<r>_<f>} (\code{NA} until filled). The attributes record the
#'   sampling design (\code{attr(x, "sampling")}), the pooling-class vocabulary
#'   the raters must use (\code{attr(x, "pooling_class_vocabulary")}) and the
#'   seed.
#' @details
#' \strong{Sampling.} Records are stratified by endpoint family and by source
#' (the cohort's accession prefix: GEO, TCGA, EMBL-EBI, CGGA, other), and the
#' \code{n} records are allocated across the strata in proportion to their size
#' (largest-remainder rounding, capped at each stratum's size). A stratum smaller
#' than its share is taken whole and the remainder is redistributed, so exactly
#' \code{n} records are returned whenever the table has that many.
#'
#' \strong{Stable ids.} \code{record_id} is assigned to the whole
#' cohort-endpoint table in catalog order (\code{Accession}, then \code{Family})
#' as \code{REC0001}, \code{REC0002}, ... before any filtering or sampling. A
#' record therefore keeps the same id when \code{n}, \code{seed}, \code{families},
#' \code{drop_absent} or the rater count changes, and two sheets can be compared
#' by id. The consequence is that the ids of a sample are sparse (they are not
#' \code{1..n}).
#'
#' \strong{Blank fields.} Every rater column is \code{NA_character_} so that a
#' missing rating is distinguishable from a rating of \code{""}: the agreement
#' calculator compares only pairs in which both ratings are present and
#' non-empty.
#' @seealso \code{\link{endpoint_agreement}}, \code{\link{endpoint_semantics}},
#'   \code{\link{endpoint_pooling_class}}
#' @export
#' @examples
#' ## Offline: the shipped endpoint-semantics table, no network.
#' sheet <- endpoint_adjudication(n = 6, raters = 2, seed = 42)
#' sheet[, c("record_id", "cohort", "token")]
#' attr(sheet, "pooling_class_vocabulary")
#'
#' ## Fill it in (here: one rater copies the second for the worked example).
#' sheet$rater1_pooling_class <- rep("Exact-equivalent", nrow(sheet))
#' sheet$rater2_pooling_class <- sheet$rater1_pooling_class
#' endpoint_agreement(sheet)$overall_raw_agreement
endpoint_adjudication <- function(n = 80, raters = 2, seed = 20260101,
                                  families = NULL, drop_absent = FALSE,
                                  path = NULL, table = NULL) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1)
    stop("'n' must be a single positive number of records to draw.", call. = FALSE)
  n <- as.integer(n)
  if (!is.numeric(raters) || length(raters) != 1L || is.na(raters) || raters < 1)
    stop("'raters' must be a single positive number of raters.", call. = FALSE)
  raters <- as.integer(raters)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
    stop("'seed' must be a single number.", call. = FALSE)

  tab <- if (is.null(table)) endpoint_semantics() else table
  if (!is.data.frame(tab) ||
      !all(c("Accession", "Family", "Token", "Evidence") %in% colnames(tab)))
    stop("'table' must carry the endpoint-semantics columns ",
         "Accession, Family, Token and Evidence.", call. = FALSE)

  # Stable ids are assigned to the WHOLE cohort-endpoint table, in catalog order,
  # before any filter, so a record keeps its id whatever subset is drawn.
  tab <- tab[order(as.character(tab$Accession), as.character(tab$Family)), , drop = FALSE]
  rownames(tab) <- NULL
  tab$record_id <- sprintf("REC%04d", seq_len(nrow(tab)))

  if (!is.null(families)) {
    tab <- tab[tab$Family %in% as.character(families), , drop = FALSE]
    if (!nrow(tab))
      stop("No record belongs to the requested famil", "ies ",
           paste(families, collapse = ", "), ".", call. = FALSE)
  }
  if (isTRUE(drop_absent)) {
    tok <- as.character(tab$Token)
    tab <- tab[!is.na(tok) & nzchar(tok) & tok != "not stated", , drop = FALSE]
    if (!nrow(tab))
      stop("No record reports an endpoint after 'drop_absent = TRUE'.", call. = FALSE)
  }

  tab$source <- .cpas_accession_source(tab$Accession)
  tab$stratum <- paste(as.character(tab$Family), tab$source, sep = " | ")

  # Hamilton allocation: floor the capped target's quotas, then award each
  # residual seat to a distinct stratum, breaking ties by stratum name.
  sizes <- table(tab$stratum)
  strata <- names(sizes)
  capacity <- as.integer(sizes)
  target <- min(n, nrow(tab))
  quota <- target * capacity / sum(capacity)
  take <- pmin(floor(quota), capacity)
  names(take) <- strata
  remaining <- target - sum(take)
  if (remaining > 0L) {
    room <- which(take < capacity)
    ord <- room[order(-(quota[room] - floor(quota[room])), strata[room])]
    if (remaining > length(ord))
      stop("Cannot allocate the requested sample within stratum capacities.",
           call. = FALSE)
    winners <- ord[seq_len(remaining)]
    take[winners] <- take[winners] + 1L
  }
  if (sum(take) != target || any(take < 0L | take > capacity))
    stop("Invalid sample allocation for the available stratum capacities.",
         call. = FALSE)

  # Deterministic within-stratum draw; the caller's RNG state is restored.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed) %% .Machine$integer.max)
  idx <- unlist(lapply(names(take)[take > 0L], function(s) {
    rows <- which(tab$stratum == s)
    if (length(rows) <= take[[s]]) rows else sort(sample(rows, take[[s]]))
  }), use.names = FALSE)
  out <- tab[idx, , drop = FALSE]
  out <- out[order(out$record_id), , drop = FALSE]
  rownames(out) <- NULL

  res <- data.frame(record_id = out$record_id,
                    cohort = as.character(out$Accession),
                    token = as.character(out$Token),
                    evidence = as.character(out$Evidence),
                    stringsAsFactors = FALSE)
  for (r in seq_len(raters)) {
    for (f in .cpas_adjudication_fields)
      res[[sprintf("rater%d_%s", r, f)]] <- NA_character_
  }
  attr(res, "pooling_class_vocabulary") <- .cpas_pooling_class_levels
  attr(res, "judgement_fields") <- .cpas_adjudication_fields
  attr(res, "seed") <- seed
  attr(res, "sampling") <- list(
    n_requested = n, n_drawn = nrow(res), raters = raters, seed = seed,
    families = if (is.null(families)) unique(as.character(tab$Family)) else as.character(families),
    drop_absent = isTRUE(drop_absent),
    strata = data.frame(stratum = names(sizes),
                        available = as.integer(sizes),
                        drawn = as.integer(take[match(names(sizes), names(take))]),
                        stringsAsFactors = FALSE),
    rule = paste0("stratified by endpoint family and by source (accession prefix); ",
                  "proportional allocation with largest-remainder rounding, capped at ",
                  "the stratum size; the draw is seeded and the ids are assigned ",
                  "before filtering and sampling, so they are stable across n, seed ",
                  "and the family / absent filters"))
  if (!is.null(path)) {
    utils::write.csv(res, path, row.names = FALSE, na = "")
    attr(res, "csv_path") <- path
  }
  class(res) <- c("cpas_adjudication_sheet", "data.frame")
  res
}

#' @title Agreement between endpoint raters
#' @description
#' Takes a filled rating sheet (from \code{\link{endpoint_adjudication}}) and
#' reports how far the raters agree: per-field raw agreement, \strong{Cohen's
#' kappa} per field with its Landis-and-Koch interpretation, the overall pooled
#' agreement, the adjudication rate and the records whose raters disagreed.
#' @param sheet A filled rating sheet: the columns of
#'   \code{\link{endpoint_adjudication}} (\code{record_id}, \code{cohort},
#'   \code{token}, \code{evidence}) plus at least two \code{rater<r>_<field>}
#'   columns per judgement field with at least one non-empty pair.
#' @param fields Judgement fields to score (default: the five the sheet
#'   prepares: \code{event_definition}, \code{time_origin},
#'   \code{censoring_rule}, \code{competing_events}, \code{pooling_class}).
#' @return A list of class \code{cpas_adjudication_agreement}:
#'   \item{\code{n_records}, \code{n_raters}, \code{fields}:}{the size of the
#'     exercise and the fields scored}
#'   \item{\code{per_field}:}{one row per field: \code{field},
#'     \code{n_compared} (record-pairs in which both ratings are present),
#'     \code{n_agree}, \code{raw_agreement}, \code{kappa},
#'     \code{kappa_interpretation}, \code{kappa_basis}, \code{note}}
#'   \item{\code{overall_raw_agreement}:}{pooled over all fields and rater pairs
#'     (agreeing compared pairs / all compared pairs)}
#'   \item{\code{pooled_kappa}:}{the mean of the per-field kappas (the pooled
#'     agreement on the kappa scale)}
#'   \item{\code{adjudication_rate}:}{records with at least one disagreement
#'     divided by the records that carry at least one comparable rating}
#'   \item{\code{disagreements}:}{the disagreeing records, ready for
#'     adjudication: \code{record_id}, \code{cohort}, \code{token},
#'     \code{evidence}, \code{fields_disagreed} and \code{detail} (which rater
#'     said what, field by field)}
#'   \item{\code{pairwise}:}{per field and rater pair, the raw agreement and
#'     kappa; for two raters this is the same number as \code{per_field}}
#'   \item{\code{kappa_basis}, \code{notes}:}{how kappa was aggregated, and any
#'     field whose kappa is undefined}
#' @details
#' \strong{Kappa.} Cohen's kappa is computed from the nominal contingency table
#' of the two raters: \eqn{\kappa = (p_o - p_e) / (1 - p_e)}, with \eqn{p_o} the
#' observed agreement and \eqn{p_e} the agreement expected from the margins.
#' Only rater pairs in which both ratings are present and non-empty are used, so
#' missing ratings reduce \code{n_compared} instead of biasing the estimate, and
#' a field in which every rating is the same single category reports
#' \code{kappa = NA} (chance agreement 1 leaves it undefined) rather than
#' \code{NaN}. With more than two raters every rater pair is scored and the
#' field's kappa is the \strong{mean of the pairwise Cohen's kappas} - kappa is
#' defined for two raters, and this is said explicitly in
#' \code{kappa_basis} and in \code{notes}. For two raters there is exactly one
#' pair and the mean is that pair's kappa.
#' @seealso \code{\link{endpoint_adjudication}}
#' @export
#' @examples
#' ## A tiny synthesised sheet (no real study data): two raters, known table.
#' ## 20 records where both say "A", 5 where rater1 says "A" and rater2 "B",
#' ## 10 the other way round and 15 where both say "B":
#' ##   po = 35/50 = 0.70 ; pe = (25*30 + 25*20)/50^2 = 0.50 ; kappa = 0.40
#' set.seed(1)
#' r1 <- c(rep("A", 25), rep("B", 25))
#' r2 <- c(rep("A", 20), rep("B", 5), rep("A", 10), rep("B", 15))
#' sheet <- data.frame(record_id = sprintf("REC%04d", 1:50),
#'                     cohort = "SYN", token = "OS", evidence = "synthetic",
#'                     rater1_pooling_class = r1, rater2_pooling_class = r2,
#'                     stringsAsFactors = FALSE)
#' ag <- endpoint_agreement(sheet, fields = "pooling_class")
#' ag$per_field[, c("field", "n_compared", "raw_agreement", "kappa",
#'                  "kappa_interpretation")]
endpoint_agreement <- function(sheet, fields = .cpas_adjudication_fields) {
  if (!is.data.frame(sheet))
    stop("'sheet' must be a data frame (a filled rating sheet).", call. = FALSE)
  if (!"record_id" %in% colnames(sheet))
    stop("'sheet' must carry a 'record_id' column (see endpoint_adjudication()).",
         call. = FALSE)
  fields <- as.character(fields)
  raters <- .cpas_adjudication_raters(sheet, fields)
  if (length(raters) < 1L)
    stop("No rater column found. The sheet must carry columns named ",
         "rater<r>_<field> for at least one rater and field.", call. = FALSE)
  if (length(raters) < 2L)
    stop("Agreement needs at least two raters; the sheet has ", length(raters),
         " (rater", paste(raters, collapse = ", rater"), ").", call. = FALSE)
  n <- nrow(sheet)
  pairs <- utils::combn(raters, 2)          # 2 x (R choose 2)
  pair_ids <- apply(pairs, 2, paste, collapse = "-")

  per_field <- list(); pairwise <- list(); notes <- character(0)
  total_num <- 0L; total_den <- 0L
  for (f in fields) {
    cols <- sprintf("rater%d_%s", raters, f)
    if (!all(cols %in% colnames(sheet))) {
      notes <- c(notes, sprintf("field '%s' has no rating columns for all raters; skipped", f))
      next
    }
    num_f <- 0L; den_f <- 0L; kap_f <- numeric(0)
    for (j in seq_along(pair_ids)) {
      a <- as.character(sheet[[sprintf("rater%d_%s", pairs[1, j], f)]])
      b <- as.character(sheet[[sprintf("rater%d_%s", pairs[2, j], f)]])
      kk <- .cpas_cohen_kappa(a, b)
      ok <- !is.na(a) & !is.na(b) & nzchar(a) & nzchar(b)
      num_f <- num_f + sum(a[ok] == b[ok]); den_f <- den_f + sum(ok)
      kap_f <- c(kap_f, kk$kappa)
      pairwise[[length(pairwise) + 1L]] <- data.frame(
        field = f, pair = pair_ids[j], n_compared = kk$n,
        raw_agreement = kk$po, kappa = kk$kappa,
        kappa_interpretation = .cpas_kappa_label(kk$kappa),
        stringsAsFactors = FALSE)
    }
    total_num <- total_num + num_f; total_den <- total_den + den_f
    kap_mean <- if (any(is.finite(kap_f))) mean(kap_f[is.finite(kap_f)]) else NA_real_
    basis <- if (length(raters) > 2L)
      sprintf("mean of the %d pairwise Cohen's kappas over all rater pairs",
              length(pair_ids)) else "Cohen's kappa (two raters)"
    note <- if (!any(is.finite(kap_f)) && den_f > 0L)
      "every usable pair agrees on a single category, so chance agreement is 1 and kappa is undefined" else NA_character_
    if (den_f == 0L) {
      note <- "no record was rated by both members of any rater pair"
      notes <- c(notes, sprintf("field '%s': %s", f, note))
    } else if (!is.na(note)) {
      notes <- c(notes, sprintf("field '%s': %s", f, note))
    }
    per_field[[length(per_field) + 1L]] <- data.frame(
      field = f, n_compared = as.integer(den_f), n_agree = as.integer(num_f),
      raw_agreement = if (den_f > 0L) num_f / den_f else NA_real_,
      kappa = kap_mean,
      kappa_interpretation = .cpas_kappa_label(kap_mean),
      kappa_basis = basis, note = note, stringsAsFactors = FALSE)
  }
  if (!length(per_field))
    stop("None of the requested 'fields' has rating columns for two raters.",
         call. = FALSE)
  per_field <- do.call(rbind, per_field)
  rownames(per_field) <- NULL
  pairwise_df <- do.call(rbind, pairwise)
  rownames(pairwise_df) <- NULL

  # Records whose raters disagree, field by field, ready for adjudication.
  dis <- list()
  for (i in seq_len(n)) {
    bad <- character(0); det <- character(0)
    for (f in fields) {
      cols <- sprintf("rater%d_%s", raters, f)
      if (!all(cols %in% colnames(sheet))) next
      vals <- vapply(cols, function(cc) as.character(sheet[[cc]][i]), character(1))
      have <- !is.na(vals) & nzchar(vals)
      if (sum(have) >= 2L && length(unique(vals[have])) > 1L) {
        bad <- c(bad, f)
        det <- c(det, sprintf("%s: %s", f,
                              paste(sprintf("rater%d=%s", raters[have], vals[have]),
                                    collapse = "; ")))
      }
    }
    if (length(bad)) {
      dis[[length(dis) + 1L]] <- data.frame(
        record_id = as.character(sheet$record_id[i]),
        cohort = if ("cohort" %in% colnames(sheet)) as.character(sheet$cohort[i]) else NA_character_,
        token = if ("token" %in% colnames(sheet)) as.character(sheet$token[i]) else NA_character_,
        evidence = if ("evidence" %in% colnames(sheet)) as.character(sheet$evidence[i]) else NA_character_,
        fields_disagreed = paste(bad, collapse = ", "),
        detail = paste(det, collapse = " || "),
        stringsAsFactors = FALSE)
    }
  }
  dis_df <- if (length(dis)) do.call(rbind, dis) else
    data.frame(record_id = character(0), cohort = character(0),
               token = character(0), evidence = character(0),
               fields_disagreed = character(0), detail = character(0),
               stringsAsFactors = FALSE)
  rownames(dis_df) <- NULL

  # A record counts towards the adjudication rate when at least one field has a
  # usable rating pair; otherwise it carries no comparison to disagree on.
  comparable <- rep(FALSE, n)
  for (f in fields) {
    cols <- sprintf("rater%d_%s", raters, f)
    if (!all(cols %in% colnames(sheet))) next
    for (i in seq_len(n)) {
      vals <- vapply(cols, function(cc) as.character(sheet[[cc]][i]), character(1))
      have <- !is.na(vals) & nzchar(vals)
      if (sum(have) >= 2L) comparable[i] <- TRUE
    }
  }
  n_comparable <- sum(comparable)
  kap_all <- per_field$kappa[is.finite(per_field$kappa)]

  out <- list(
    n_records = n,
    n_raters = length(raters),
    raters = raters,
    fields = per_field$field,
    per_field = per_field,
    overall_raw_agreement = if (total_den > 0L) total_num / total_den else NA_real_,
    pooled_kappa = if (length(kap_all)) mean(kap_all) else NA_real_,
    pooled_kappa_interpretation = .cpas_kappa_label(if (length(kap_all)) mean(kap_all) else NA_real_),
    adjudication_rate = if (n_comparable > 0L) nrow(dis_df) / n_comparable else NA_real_,
    n_comparable_records = as.integer(n_comparable),
    disagreements = dis_df,
    pairwise = pairwise_df,
    kappa_basis = if (length(raters) > 2L)
      sprintf("Cohen's kappa per field is the mean of its %d pairwise kappas (>2 raters)",
              ncol(pairs)) else "Cohen's kappa (nominal, two raters)",
    pooling_class_vocabulary = attr(sheet, "pooling_class_vocabulary"),
    notes = notes)
  class(out) <- "cpas_adjudication_agreement"
  out
}

#' @title Print an endpoint-adjudication agreement result
#' @param x An object of class \code{cpas_adjudication_agreement}.
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_adjudication_agreement <- function(x, ...) {
  cat("CanPAS endpoint adjudication: agreement between", x$n_raters, "rater(s)\n")
  cat(sprintf("records: %d | comparable: %d | adjudication rate: %s\n",
              x$n_records, x$n_comparable_records,
              if (is.na(x$adjudication_rate)) "NA" else
                sprintf("%.1f%% (%d record(s))", 100 * x$adjudication_rate,
                        nrow(x$disagreements))))
  print(x$per_field[, c("field", "n_compared", "raw_agreement", "kappa",
                           "kappa_interpretation")], row.names = FALSE)
  cat(sprintf("overall pooled raw agreement: %s | pooled kappa: %s\n",
              if (is.na(x$overall_raw_agreement)) "NA" else
                sprintf("%.3f", x$overall_raw_agreement),
              if (is.na(x$pooled_kappa)) "NA" else
                sprintf("%.3f (%s)", x$pooled_kappa, x$pooled_kappa_interpretation)))
  cat("kappa basis:", x$kappa_basis, "\n")
  if (length(x$notes)) cat("notes:", paste(x$notes, collapse = " | "), "\n")
  if (nrow(x$disagreements))
    cat(sprintf("%d record(s) need adjudication (see $disagreements).\n",
                nrow(x$disagreements)))
  invisible(x)
}
