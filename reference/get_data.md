# Query the CanPAS GEO API

Query the CanPAS public API (backed by a MySQL mirror of GEO
series-matrix survival / expression tables). The response is returned as
a structured object of class `cpas_get` together with the exact URL
used.

## Usage

``` r
get_data(
  dataset,
  action = c("expression", "gpl", "surv_data"),
  ids = NULL,
  base_url = "https://www.jingege.wang/bioinformatics/CPAS/api.php",
  timeout_sec = 120,
  use_cache = NULL
)
```

## Arguments

- dataset:

  Mirror table to query: a dataset accession for `action = "expression"`
  / `"surv_data"`, or a platform id (e.g. `"GPL96"`) for
  `action = "gpl"`. The value is sent to the API as its `table` request
  parameter, with one substitution: a dash is sent as the underscore the
  mirror uses (`"A5-PCPG"` queries the table `A5_PCPG`), because the API
  interpolates the name into its SQL and cannot serve a dash-spelled
  table. expression / survival queries or `"GPL570"` for platform (gpl)
  queries.

- action:

  One of `"expression"` (sample-level expression of the requested
  probes/ids), `"gpl"` (probe-to-gene mapping of a platform) or
  `"surv_data"` (survival table of a dataset).

- ids:

  Reference ids (probes for `action="expression"`, Entrez gene ids for
  `action="gpl"`). Ignored for `action="surv_data"`. May be a character
  vector; multiple ids are collapsed with `","`.

- base_url:

  API endpoint. Defaults to the CanPAS public API.

- timeout_sec:

  Connection timeout in seconds passed to `getOption("timeout")` for the
  duration of the call.

- use_cache:

  Reuse the local file cache. `NULL` (default) follows
  `getOption("CanPAS.cache", TRUE)` and the environment variables
  `CANPAS_CACHE` / `CPAS_CACHE`; `TRUE` or `FALSE` forces the choice for
  this call. Cache files live in `getOption("CanPAS.cache_dir")`,
  `CANPAS_CACHE_DIR` or `tools::R_user_dir("CanPAS", "cache")`, expire
  after `getOption("CanPAS.cache_ttl")` seconds (30 days by default) and
  are used even after expiry when the network request fails, so an
  offline or rate-limited mirror degrades to the last known copy instead
  of an error.

## Value

A list of class `cpas_get` with elements:

- `input_params`::

  inputs echoed back with the request time

- `response`::

  decoded JSON response (usually a data.frame)

- `url`::

  the constructed request URL

- `metadata`::

  api version, response time, result count, and `from_cache` /
  `cache_stale` flags reporting whether the answer came from the local
  file cache

## Details

The public API is rate-limited. Sequential bulk queries (e.g. sweeping
many datasets) frequently trigger HTTP 429 responses; retry after a
short pause, or leave the local cache on (the default), which answers
repeated queries for the same table, platform or probe set from disk.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## Raw API access behind the prerequisite functions: the mirrored survival
   ## table of a real cohort (ID, endpoint columns and clinical covariates).
   r <- get_data("GSE14814", "surv_data")
   head(r$response)

   ## Sample-level expression for chosen platform ids.
   e <- get_data("GSE14814", "expression", ids = c("1053_at", "200801_x_at"))
   head(e$response)[, 1:4]
} # }
```
