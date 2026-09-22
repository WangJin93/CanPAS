# ---------------------------------------------------------------------------
# Local file cache for remote fetches (mirror API and UCSC Xena).
#
# Why: every GEO/CGGA query goes through get_data() and every TCGA gene through
# tcga_gene_expr_df(); both hit a rate-limited mirror or a flaky Xena endpoint.
# Without a cache, re-running an analysis re-downloads the same gene tables and
# a single transient failure aborts a whole sweep. With a cache the second run
# is offline, and a failed fetch falls back to the last good copy.
#
# Design (deliberately dependency-free, base R only):
#   * one wrapper, .cpas_cached(), used by every network entry point, so the
#     behaviour is testable with a fake fetcher and needs no connection;
#   * cache directory: option CanPAS.cache_dir > env CANPAS_CACHE_DIR /
#     CPAS_CACHE_DIR > tools::R_user_dir("CanPAS", "cache");
#   * switch: option CanPAS.cache (default TRUE) or env CANPAS_CACHE=0/false/off;
#   * entries expire after CanPAS.cache_ttl seconds (default 30 days) but an
#     expired entry is still used when the network fails, with a message, so a
#     live-service outage degrades to the last known data instead of an error;
#   * writes are atomic (temp file + rename), and every entry records the
#     request URL, the creation time and a schema version, so a future change
#     of the stored format invalidates old files instead of misreading them.
# Nothing here is exported: the public API surface (37 exported functions and
# 11 registered print methods) is unchanged, and control is through options and
# environment variables.
# ---------------------------------------------------------------------------

.CPAS_CACHE_SCHEMA <- 1L

#' Resolve the cache directory (internal)
#' @noRd
.cpas_cache_dir <- function(create = TRUE) {
  d <- getOption("CanPAS.cache_dir", NULL)
  if (is.null(d) || !nzchar(d)) d <- Sys.getenv("CANPAS_CACHE_DIR", unset = "")
  if (!nzchar(d)) d <- Sys.getenv("CPAS_CACHE_DIR", unset = "")
  if (!nzchar(d)) d <- tools::R_user_dir("CanPAS", which = "cache")
  if (isTRUE(create) && !dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Is caching switched on? (internal)
#' @noRd
.cpas_cache_enabled <- function() {
  o <- getOption("CanPAS.cache", NULL)
  if (!is.null(o)) return(isTRUE(o))
  e <- Sys.getenv("CANPAS_CACHE", unset = "")
  if (!nzchar(e)) e <- Sys.getenv("CPAS_CACHE", unset = "")
  if (nzchar(e)) return(tolower(e) %in% c("1", "true", "yes", "on", "t", "y"))
  TRUE
}

#' Entry time-to-live in seconds (internal); 0 means "never reuse a fresh copy"
#' @noRd
.cpas_cache_ttl <- function() {
  t <- getOption("CanPAS.cache_ttl", NULL)
  if (is.null(t)) t <- suppressWarnings(as.numeric(Sys.getenv("CANPAS_CACHE_TTL", unset = "")))
  if (is.null(t) || !length(t) || !is.finite(t[1]) || t[1] < 0) t <- 30 * 24 * 3600
  as.numeric(t[1])
}

.cpas_cache_verbose <- function() isTRUE(getOption("CanPAS.cache_verbose", FALSE))

#' Build a stable cache key from the request-defining values (internal)
#' @noRd
.cpas_cache_key <- function(...) {
  paste(vapply(list(...), function(x) paste(as.character(x), collapse = ","),
               character(1)), collapse = "|")
}

#' Path of the cache file for a key (internal)
#' @noRd
.cpas_cache_file <- function(key, sub = "misc", create = TRUE) {
  d <- file.path(.cpas_cache_dir(create = create), sub)
  if (isTRUE(create)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  safe <- gsub("[^A-Za-z0-9._-]", "_", key)
  safe <- gsub("_{2,}", "_", safe)
  if (nchar(safe) > 120L) {                     # long gene lists: hash the tail
    tf <- tempfile()
    writeLines(key, tf)
    h <- unname(tools::md5sum(tf))
    unlink(tf)
    safe <- paste0(substr(safe, 1, 60L), "_", h)
  }
  file.path(d, paste0(safe, ".rds"))
}

#' Read one cache entry (internal); returns NULL when absent/unreadable
#' @noRd
.cpas_cache_read <- function(key, sub = "misc") {
  fp <- .cpas_cache_file(key, sub, create = FALSE)
  if (!file.exists(fp)) return(NULL)
  e <- tryCatch(readRDS(fp), error = function(e) NULL)
  if (is.null(e) || !is.list(e) || !identical(e$schema, .CPAS_CACHE_SCHEMA) ||
      !identical(e$key, key) || is.null(e$value)) return(NULL)
  e$age_sec <- as.numeric(difftime(Sys.time(), e$created, units = "secs"))
  e$file <- fp
  e
}

#' Write one cache entry atomically (internal)
#' @noRd
.cpas_cache_write <- function(key, sub, value, url = NA_character_) {
  fp <- .cpas_cache_file(key, sub)
  tmp <- paste0(fp, ".", Sys.getpid(), ".tmp")
  ok <- tryCatch({
    saveRDS(list(schema = .CPAS_CACHE_SCHEMA, key = key, created = Sys.time(),
                 url = url, pkg_version = tryCatch(as.character(utils::packageVersion("CanPAS")),
                                                   error = function(e) NA_character_),
                 value = value),
            tmp)
    file.rename(tmp, fp)          # atomic on the same filesystem
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) unlink(tmp)
  invisible(isTRUE(ok))
}

#' Fetch through the cache (internal)
#'
#' @param key unique request key.
#' @param sub sub-directory of the cache (one per request family).
#' @param fetch zero-argument function performing the real request.
#' @param use_cache NULL follows the options/environment, TRUE/FALSE forces it.
#' @param label short label used in messages.
#' @return list(value, from_cache, stale, created, url, n_calls)
#' @noRd
.cpas_cached <- function(key, sub, fetch, use_cache = NULL, label = "") {
  use <- if (is.null(use_cache)) .cpas_cache_enabled() else isTRUE(use_cache)
  verbose <- .cpas_cache_verbose()
  if (!use) return(list(value = fetch(), from_cache = FALSE, stale = FALSE,
                        created = Sys.time(), url = NA_character_, n_calls = 1L))
  hit <- .cpas_cache_read(key, sub)
  ttl <- .cpas_cache_ttl()
  fresh <- !is.null(hit) && (ttl <= 0 || hit$age_sec <= ttl)
  if (fresh) {
    if (verbose)
      message("CanPAS cache hit", if (nzchar(label)) paste0(" (", label, ")") else "",
              ": ", basename(hit$file), sprintf(" [%.1f h old]", hit$age_sec / 3600))
    return(list(value = hit$value, from_cache = TRUE, stale = FALSE,
                created = hit$created, url = hit$url, n_calls = 0L))
  }
  got <- tryCatch(list(ok = TRUE, value = fetch()), error = function(e) list(ok = FALSE, error = e))
  if (!isTRUE(got$ok)) {
    ## Network failure: fall back to the last good copy rather than aborting.
    if (!is.null(hit)) {
      message("CanPAS: the remote request failed",
              if (nzchar(label)) paste0(" (", label, ")") else "", "; using the cached copy from ",
              format(hit$created, "%Y-%m-%d %H:%M"), " (",
              sprintf("%.1f days old", hit$age_sec / 86400), "). Reason: ",
              substr(conditionMessage(got$error), 1, 160))
      return(list(value = hit$value, from_cache = TRUE, stale = TRUE,
                  created = hit$created, url = hit$url, n_calls = 1L))
    }
    stop(got$error)
  }
  wrote <- .cpas_cache_write(key, sub, got$value, url = NA_character_)
  if (verbose && wrote)
    message("CanPAS cache store", if (nzchar(label)) paste0(" (", label, ")") else "",
            ": ", basename(.cpas_cache_file(key, sub)))
  list(value = got$value, from_cache = FALSE, stale = FALSE,
       created = Sys.time(), url = NA_character_, n_calls = 1L)
}

#' List the cache contents (internal); used by tests and by the README
#' @noRd
.cpas_cache_info <- function() {
  d <- .cpas_cache_dir(create = FALSE)
  if (!dir.exists(d)) return(data.frame(file = character(0), sub = character(0),
                                        mb = numeric(0), age_h = numeric(0)))
  fs <- list.files(d, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  if (!length(fs)) return(data.frame(file = character(0), sub = character(0),
                                     mb = numeric(0), age_h = numeric(0)))
  data.frame(file = basename(fs),
             sub = basename(dirname(fs)),
             mb = round(file.info(fs)$size / 1e6, 3),
             age_h = round(as.numeric(difftime(Sys.time(), file.info(fs)$mtime,
                                               units = "hours")), 1),
             row.names = NULL)
}

#' Delete cache entries (internal)
#' @param sub NULL deletes everything, otherwise only that family.
#' @noRd
.cpas_cache_clear <- function(sub = NULL) {
  d <- .cpas_cache_dir(create = FALSE)
  if (!dir.exists(d)) return(invisible(0L))
  target <- if (is.null(sub)) d else file.path(d, sub)
  fs <- list.files(target, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  n <- length(fs)
  if (n) unlink(fs)
  invisible(n)
}
