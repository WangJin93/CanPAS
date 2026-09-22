#' @title Query the CanPAS GEO API
#' @description
#' Query the CanPAS public API (backed by a MySQL mirror of GEO series-matrix
#' survival / expression tables). The response is returned as a structured
#' object of class \code{cpas_get} together with the exact URL used.
#' @param dataset Mirror table to query: a dataset accession for
#'   \code{action = "expression"} / \code{"surv_data"}, or a platform id
#'   (e.g. \code{"GPL96"}) for \code{action = "gpl"}. The value is sent to the
#'   API as its \code{table} request parameter.
#' expression / survival queries or \code{"GPL570"} for platform (gpl) queries.
#' @param action One of \code{"expression"} (sample-level expression of the
#' requested probes/ids), \code{"gpl"} (probe-to-gene mapping of a platform) or
#' \code{"surv_data"} (survival table of a dataset).
#' @param ids Reference ids (probes for \code{action="expression"}, Entrez gene
#' ids for \code{action="gpl"}). Ignored for \code{action="surv_data"}. May be a
#' character vector; multiple ids are collapsed with \code{","}.
#' @param base_url API endpoint. Defaults to the CanPAS public API.
#' @param timeout_sec Connection timeout in seconds passed to
#' \code{getOption("timeout")} for the duration of the call.
#' @param use_cache Reuse the local file cache. \code{NULL} (default) follows
#'   \code{getOption("CanPAS.cache", TRUE)} and the environment variables
#'   \code{CANPAS_CACHE} / \code{CPAS_CACHE}; \code{TRUE} or \code{FALSE}
#'   forces the choice for this call. Cache files live in
#'   \code{getOption("CanPAS.cache_dir")}, \code{CANPAS_CACHE_DIR} or
#'   \code{tools::R_user_dir("CanPAS", "cache")}, expire after
#'   \code{getOption("CanPAS.cache_ttl")} seconds (30 days by default) and are
#'   used even after expiry when the network request fails, so an offline or
#'   rate-limited mirror degrades to the last known copy instead of an error.
#' @return A list of class \code{cpas_get} with elements:
#'   \item{\code{input_params}:}{inputs echoed back with the request time}
#'   \item{\code{response}:}{decoded JSON response (usually a data.frame)}
#'   \item{\code{url}:}{the constructed request URL}
#'   \item{\code{metadata}:}{api version, response time, result count, and
#'     \code{from_cache} / \code{cache_stale} flags reporting whether the
#'     answer came from the local file cache}
#' @details
#' The public API is rate-limited. Sequential bulk queries (e.g. sweeping many
#' datasets) frequently trigger HTTP 429 responses; retry after a short pause,
#' or leave the local cache on (the default), which answers repeated queries for
#' the same table, platform or probe set from disk.
#' @examples
#' \dontrun{
#'    ## Raw API access behind the prerequisite functions: the mirrored survival
#'    ## table of a real cohort (ID, endpoint columns and clinical covariates).
#'    r <- get_data("GSE14814", "surv_data")
#'    head(r$response)
#' 
#'    ## Sample-level expression for chosen platform ids.
#'    e <- get_data("GSE14814", "expression", ids = c("1053_at", "200801_x_at"))
#'    head(e$response)[, 1:4]
#' }
#' @export
get_data <- function(dataset,
                     action = c("expression", "gpl", "surv_data"),
                     ids = NULL,
                     base_url = "https://www.jingege.wang/bioinformatics/CPAS/api.php",
                     timeout_sec = 120,
                     use_cache = NULL) {
  action <- match.arg(action)
  if (length(dataset) != 1L || is.na(dataset) || !nzchar(dataset))
    stop("'dataset' must be a single non-missing non-empty string.")
  if (!grepl("^[A-Za-z0-9_\\-]+$", dataset))
    stop("'dataset' contains unsupported characters: ", sQuote(dataset))

  input_params <- list(dataset = dataset, action = action, ids = ids,
                       request_time = Sys.time())

  if (action %in% c("expression", "gpl")) {
    if (is.null(ids) || length(ids) == 0L)
      stop("'ids' is required for action = '", action, "'.")
    ids <- as.character(ids)
    if (anyNA(ids) || any(!grepl("^[A-Za-z0-9./,_-]+$", ids)))
      stop("'ids' may only contain letters, digits, '.', '_', '-', '/' and ','.")
    id_str <- paste0(ids, collapse = ",")
    url <- paste0(base_url, "?action=", action, "&table=", dataset, "&ids=", id_str)
  } else {
    url <- paste0(base_url, "?action=", action, "&table=", dataset)
  }

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(1, as.numeric(timeout_sec)))

  .fetch <- function() tryCatch(jsonlite::fromJSON(url),
                  error = function(e) {
                    msg <- conditionMessage(e)
                    if (grepl("429|Too Many Requests", msg, ignore.case = TRUE))
                      stop("CanPAS API rate limit (HTTP 429). Please wait a few ",
                           "seconds and retry, or spread requests over time.", call. = FALSE)
                    if (grepl("cannot open the connection|timed out|timeout",
                              msg, ignore.case = TRUE))
                      stop("CanPAS API request failed (", msg, "). URL: ", url,
                           "\nThe mirror answers HTTP 429 when many cohorts are queried in ",
                           "quick succession; wait a few seconds and retry, or query fewer ",
                           "cohorts at a time. If the host is unreachable, check the network.",
                           call. = FALSE)
                    stop("Failed to query CanPAS API (", msg, "). URL: ", url,
                         call. = FALSE)
                  })

  key <- .cpas_cache_key("api", action, dataset, ids, base_url)
  sub <- switch(action, expression = "api_expression", gpl = "api_gpl",
                surv_data = "api_surv_data", "api")
  got <- .cpas_cached(key, sub, .fetch, use_cache = use_cache,
                      label = paste0(action, " ", dataset))
  res <- got$value
  if (is.null(res)) stop("CanPAS API returned an empty response.")

  metadata <- list(api_version = "v1",
                   response_time = if (isTRUE(got$from_cache)) got$created else Sys.time(),
                   response_format = "json",
                   result_count = if (is.data.frame(res)) nrow(res) else length(res),
                   from_cache = isTRUE(got$from_cache),
                   cache_stale = isTRUE(got$stale))

  result_obj <- list(input_params = input_params,
                     response = res,
                     url = url,
                     metadata = metadata)
  class(result_obj) <- "cpas_get"
  result_obj
}
