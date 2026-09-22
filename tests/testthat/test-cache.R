## Cache behaviour: no network needed, the fetcher is faked.
## Every test points the cache at a temporary directory and restores the
## options/environment it touches, so the user's real cache is never used.

with_tmp_cache <- function(code, ttl = 30 * 24 * 3600, enabled = TRUE,
                           verbose = FALSE) {
  old <- options(CanPAS.cache_dir = tempfile("cpas_cache_"),
                 CanPAS.cache = enabled,
                 CanPAS.cache_ttl = ttl,
                 CanPAS.cache_verbose = verbose)
  on.exit(options(old), add = TRUE)
  force(code)
}

test_that("the cache directory follows option, then environment, then R_user_dir", {
  with_tmp_cache({
    d <- getOption("CanPAS.cache_dir")
    expect_identical(CanPAS:::.cpas_cache_dir(), d)
    expect_true(dir.exists(d))

    old <- Sys.getenv("CANPAS_CACHE_DIR", unset = NA)
    options(CanPAS.cache_dir = NULL)
    on.exit({
      if (is.na(old)) Sys.unsetenv("CANPAS_CACHE_DIR") else Sys.setenv(CANPAS_CACHE_DIR = old)
    }, add = TRUE)
    Sys.setenv(CANPAS_CACHE_DIR = file.path(tempdir(), "cpas_cache_env"))
    expect_identical(CanPAS:::.cpas_cache_dir(), file.path(tempdir(), "cpas_cache_env"))

    Sys.unsetenv("CANPAS_CACHE_DIR")
    expect_identical(CanPAS:::.cpas_cache_dir(), tools::R_user_dir("CanPAS", "cache"))
  })
})

test_that("caching can be switched off by option or environment", {
  old <- options(CanPAS.cache = FALSE)
  on.exit(options(old), add = TRUE)
  expect_false(CanPAS:::.cpas_cache_enabled())

  options(CanPAS.cache = NULL)
  oldenv <- Sys.getenv("CANPAS_CACHE", unset = NA)
  on.exit({
    if (is.na(oldenv)) Sys.unsetenv("CANPAS_CACHE") else Sys.setenv(CANPAS_CACHE = oldenv)
  }, add = TRUE)
  Sys.setenv(CANPAS_CACHE = "off")
  expect_false(CanPAS:::.cpas_cache_enabled())
  Sys.setenv(CANPAS_CACHE = "1")
  expect_true(CanPAS:::.cpas_cache_enabled())
})

test_that("cache keys map to stable, safe file names", {
  k1 <- CanPAS:::.cpas_cache_key("api", "expression", "GSE31210", c("200801_x_at", "1053_at"))
  k2 <- CanPAS:::.cpas_cache_key("api", "expression", "GSE31210", c("200801_x_at", "1053_at"))
  k3 <- CanPAS:::.cpas_cache_key("api", "expression", "GSE31210", c("1053_at"))
  expect_identical(k1, k2)
  expect_false(identical(k1, k3))

  with_tmp_cache({
    f <- CanPAS:::.cpas_cache_file(k1, "api_expression")
    expect_true(grepl("^[A-Za-z0-9._/-]+$", f))
    expect_identical(basename(f), basename(CanPAS:::.cpas_cache_file(k1, "api_expression")))
  })
  ## a very long gene list is still mapped to a short, unique file name
  long <- CanPAS:::.cpas_cache_key("api", "expression", "GSE31210",
                                   paste0("probe_", 1:400))
  with_tmp_cache({
    f <- CanPAS:::.cpas_cache_file(long, "api_expression")
    expect_lt(nchar(basename(f)), 100L)
    expect_true(grepl("[0-9a-f]{32}\\.rds$", basename(f)))
  })
})

test_that("a fetched value is stored and reused", {
  with_tmp_cache({
    calls <- 0L
    fetch <- function() { calls <<- calls + 1L; data.frame(a = 1:3) }
    k <- CanPAS:::.cpas_cache_key("test", "reuse")
    r1 <- CanPAS:::.cpas_cached(k, "test_family", fetch, label = "unit")
    r2 <- CanPAS:::.cpas_cached(k, "test_family", fetch, label = "unit")
    expect_identical(calls, 1L)
    expect_false(r1$from_cache)
    expect_true(r2$from_cache)
    expect_false(r2$stale)
    expect_identical(r1$value, r2$value)
    expect_identical(r2$n_calls, 0L)
  })
})

test_that("use_cache = FALSE bypasses the cache and writes nothing", {
  with_tmp_cache({
    calls <- 0L
    fetch <- function() { calls <<- calls + 1L; 42L }
    k <- CanPAS:::.cpas_cache_key("test", "bypass")
    CanPAS:::.cpas_cached(k, "test_family", fetch, use_cache = FALSE)
    CanPAS:::.cpas_cached(k, "test_family", fetch, use_cache = FALSE)
    expect_identical(calls, 2L)
    expect_null(CanPAS:::.cpas_cache_read(k, "test_family"))
  })
})

test_that("an expired entry is refetched, and used only when the fetch fails", {
  with_tmp_cache(ttl = 1, {
    k <- CanPAS:::.cpas_cache_key("test", "ttl")
    CanPAS:::.cpas_cached(k, "test_family", function() "v1")
    Sys.sleep(1.2)                      # entry is now older than ttl = 1 s
    calls <- 0L
    r <- CanPAS:::.cpas_cached(k, "test_family", function() { calls <<- calls + 1L; "v2" })
    expect_identical(calls, 1L)         # expired -> refetched
    expect_identical(r$value, "v2")
    expect_false(r$from_cache)

    ## now the network goes down: the same expired entry is better than an error
    Sys.sleep(1.2)
    expect_message(
      r2 <- CanPAS:::.cpas_cached(k, "test_family",
                                  function() stop("cannot open the connection")),
      "cached copy")
    expect_identical(r2$value, "v2")
    expect_true(r2$from_cache)
    expect_true(r2$stale)
  })
})

test_that("a failed fetch without any cached copy still errors", {
  with_tmp_cache({
    k <- CanPAS:::.cpas_cache_key("test", "nofallback")
    expect_error(CanPAS:::.cpas_cached(k, "test_family",
                                       function() stop("boom")),
                 "boom")
  })
})

test_that("a corrupt or foreign entry is ignored rather than trusted", {
  with_tmp_cache({
    k <- CanPAS:::.cpas_cache_key("test", "corrupt")
    fp <- CanPAS:::.cpas_cache_file(k, "test_family")
    writeLines("not an rds file", fp)
    expect_null(CanPAS:::.cpas_cache_read(k, "test_family"))

    ## a valid file under a different key must not be returned for this key
    k2 <- CanPAS:::.cpas_cache_key("test", "other")
    CanPAS:::.cpas_cache_write(k2, "test_family", "value-of-k2")
    expect_null(CanPAS:::.cpas_cache_read(k, "test_family"))

    ## an older schema is ignored
    saveRDS(list(schema = 0L, key = k, created = Sys.time(), value = "old"),
            fp)
    expect_null(CanPAS:::.cpas_cache_read(k, "test_family"))
  })
})

test_that("cache entries can be listed and cleared family by family", {
  with_tmp_cache({
    CanPAS:::.cpas_cache_write(CanPAS:::.cpas_cache_key("test", "a"), "fam_a", 1)
    CanPAS:::.cpas_cache_write(CanPAS:::.cpas_cache_key("test", "b"), "fam_a", 2)
    CanPAS:::.cpas_cache_write(CanPAS:::.cpas_cache_key("test", "c"), "fam_b", 3)
    info <- CanPAS:::.cpas_cache_info()
    expect_identical(nrow(info), 3L)
    expect_setequal(unique(info$sub), c("fam_a", "fam_b"))

    expect_identical(CanPAS:::.cpas_cache_clear("fam_a"), 2L)
    expect_identical(nrow(CanPAS:::.cpas_cache_info()), 1L)
    expect_identical(CanPAS:::.cpas_cache_clear(), 1L)
    expect_identical(nrow(CanPAS:::.cpas_cache_info()), 0L)
  })
})

test_that("get_data() and tcga_gene_expr_df() expose the cache switch", {
  expect_true("use_cache" %in% names(formals(get_data)))
  expect_true("use_cache" %in% names(formals(tcga_gene_expr_df)))
  expect_null(formals(get_data)$use_cache)          # NULL = follow the options

  ## the documented cache controls are real options
  with_tmp_cache({
    expect_identical(CanPAS:::.cpas_cache_dir(), getOption("CanPAS.cache_dir"))
    expect_true(dir.exists(CanPAS:::.cpas_cache_dir()))
  })
})
