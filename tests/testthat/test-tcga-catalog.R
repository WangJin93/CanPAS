# Catalog <-> TCGA helper agreement -------------------------------------------
# Regression test for B1 of the 1.0.0 audit (2026-09-24): `tcga_retained` was a
# hard-coded 15-project vector while the catalog carried 31 TCGA cohorts, so 16
# of them failed in every TCGA helper (tcga_project_dataset / tcga_surv_table /
# tcga_merged / tcga_get_expr / get_expr_data) and canonical_type() returned NA
# for them. The project table is now derived from the shipped catalog, and this
# test runs under R CMD check - unlike the network tests in test-tcga.R, which
# are skip_on_cran() - so the next catalog expansion cannot silently break the
# helpers again. It needs no network, no credentials and no local data root.

test_that("every catalogued TCGA cohort is supported by the TCGA helpers", {
  di   <- CanPAS::dataset_info
  acc  <- as.character(di$Accession)
  tcga <- acc[grepl("^TCGA-", acc)]
  expect_gt(length(tcga), 0L)
  proj <- toupper(sub("^TCGA-", "", tcga))

  # the exported project table agrees with the catalog, cell by cell
  expect_setequal(names(tcga_retained), proj)
  expect_identical(unname(tcga_retained[proj]),
                   as.character(di$Type[match(tcga, acc)]))

  # and every one of them is accepted by the helpers (no network, no data root)
  for (i in seq_along(proj)) {
    p <- proj[i]
    expect_identical(tcga_project_dataset(p),
                     paste0("TCGA.", p, ".sampleMap/HiSeqV2"))
    expect_identical(tcga_project_dataset(tcga[i]),
                     paste0("TCGA.", p, ".sampleMap/HiSeqV2"))
    expect_false(is.na(canonical_type(tcga[i])))
    expect_identical(canonical_type(tcga[i]), as.character(di$Type[match(tcga[i], acc)]))
  }
})

test_that("the built-in fallback TCGA list matches the catalog", {
  # tcga_retained is refreshed from the catalog when the package is loaded; the
  # built-in vector is only the fallback for an unreadable catalog, so it must
  # be kept in step too. This is the assertion that fails on catalog drift.
  di <- CanPAS::dataset_info
  tc <- di[grepl("^TCGA-", as.character(di$Accession)), , drop = FALSE]
  fb <- CanPAS:::.cpas_tcga_retained_builtin
  expect_identical(names(fb), toupper(sub("^TCGA-", "", as.character(tc$Accession))))
  expect_identical(unname(fb), as.character(tc$Type))
})

test_that("the TCGA helpers still reject unknown projects", {
  expect_error(tcga_project_dataset("XXXX"), "Unsupported TCGA project")
  expect_true(is.na(canonical_type("TCGA-XXXX")))
  expect_true(is.na(canonical_type("GSE13507")))   # GEO accession, not a TCGA id
})

# Bundled clinical/survival tables --------------------------------------------
# Regression test for the user-facing failure of 2026-09-30: the two TCGA tables
# lived outside the package and CPAS_DATA_ROOT had no default, so every TCGA
# cohort failed with "Local CanPAS data file not found: /data/tcga/tcga_clinical.rda"
# unless the user configured the variable. The tables now ship in
# inst/extdata/tcga, and the resolution order is
#   explicit override > CPAS_DATA_ROOT > bundled copy > actionable error.
# These tests never need the network and never rely on CPAS_DATA_ROOT being
# pre-set: the helper below sets/clears it (and the override option/binding) per
# test and restores the previous state.

# run `code` with a controlled TCGA environment; restores everything on exit
.with_tcga_env <- function(root = NULL, option_cli = NULL, option_surv = NULL, code) {
  old_env <- Sys.getenv("CPAS_DATA_ROOT", unset = NA)
  old_opt <- options(CanPAS.tcga_clinical_rda = option_cli,
                     CanPAS.tcga_survival_rda = option_surv)
  bindings <- list()
  for (nm in c("TCGA_CLI_RDA", "TCGA_SURV_RDA")) {
    if (exists(nm, envir = globalenv(), inherits = FALSE)) {
      bindings[[nm]] <- get(nm, envir = globalenv())
      rm(list = nm, envir = globalenv())
    }
  }
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("CPAS_DATA_ROOT") else Sys.setenv(CPAS_DATA_ROOT = old_env)
    options(old_opt)
    for (nm in names(bindings)) assign(nm, bindings[[nm]], envir = globalenv())
  }, add = TRUE)
  if (is.null(root)) Sys.unsetenv("CPAS_DATA_ROOT") else Sys.setenv(CPAS_DATA_ROOT = root)
  force(code)
}

.pkg_copy <- function(name) system.file("extdata", "tcga", name, package = "CanPAS")
.checkout <- "/home/Jingle/data/Project/CPAS"   # the author's project checkout

test_that("with CPAS_DATA_ROOT unset the TCGA tables resolve inside the package", {
  .with_tcga_env(root = NULL, code = {
    t <- CanPAS::tcga_local_tables()
    expect_true(t$found)
    expect_identical(names(t$source), c("clinical", "survival"))
    expect_identical(unname(t$source), c("package", "package"))
    expect_length(c(t$clinical, t$survival), 2L)
    expect_true(all(file.exists(c(t$clinical, t$survival))))
    expect_identical(t$clinical, .pkg_copy("tcga_clinical.rda"))
    expect_identical(t$survival, .pkg_copy("tcga_surv.rda"))
    expect_match(t$clinical, "extdata[/\\\\]tcga[/\\\\]tcga_clinical[.]rda$")
    # a byte copy of the two checkout files (116,990 B + 76,102 B): re-saving
    # them with save() would change the size and is not the shipped artefact
    expect_equal(unname(file.size(c(t$clinical, t$survival))), c(116990, 76102))

    # both tables load and keep their documented object names
    e <- new.env()
    load(t$clinical, envir = e); load(t$survival, envir = e)
    expect_true(all(c("tcga_clinical", "tcga_surv") %in% ls(e)))
    expect_true(all(c("sample", "type", "vital_status") %in% colnames(e$tcga_clinical)))
    expect_true(all(c("sample", "OS", "OS.time") %in% colnames(e$tcga_surv)))

    # the clinical table carries exactly the catalog's TCGA projects
    di  <- CanPAS::dataset_info
    acc <- as.character(di$Accession)
    cat_proj <- toupper(sub("^TCGA-", "", acc[grepl("^TCGA-", acc)]))
    expect_setequal(unique(as.character(e$tcga_clinical$type)), cat_proj)
    # every catalogued project has at least one survival row
    m <- merge(e$tcga_clinical[, c("sample", "type")], e$tcga_surv, by = "sample")
    expect_setequal(unique(as.character(m$type)), cat_proj)
  })
})

test_that("with CPAS_DATA_ROOT unset a TCGA cohort can be prepared for analysis", {
  .with_tcga_env(root = NULL, code = {
    # the covariate-bearing survival table the app's TCGA path consumes
    s <- tcga_surv_table("TCGA-LUAD")
    expect_s3_class(s, "data.frame")
    expect_true(all(c("ID", "OS_status", "OS_time", "age", "sex", "stage") %in%
                      colnames(s)))
    expect_equal(nrow(s), 574L)
    expect_true(all(s$OS_status[!is.na(s$OS_status)] %in% c(0, 1)))
    expect_gt(sum(s$OS_status == 1, na.rm = TRUE), 100L)
    # 565 of the 574 LUAD samples carry an OS time; the rest are NA, the
    # documented "endpoint not available for that sample" representation
    expect_true(all(s$OS_time[!is.na(s$OS_time)] >= 0))
    expect_gt(sum(!is.na(s$OS_time) & !is.na(s$OS_status)), 500L)
    # the acceptance form used elsewhere in the package
    expect_identical(tcga_project_dataset("TCGA-LUAD"), "TCGA.LUAD.sampleMap/HiSeqV2")

    # merged data for one cohort, built offline: the survival table plus a
    # synthetic marker (the real expression step needs UCSC Xena, so it is not
    # called here) is analysable without any local data root
    set.seed(1)
    d <- s
    d$TP53 <- stats::rnorm(nrow(d))
    r <- COX_analysis(d, type = "OS", cont_Variates = "TP53", method = "uni")
    expect_s3_class(r, "cpas_COX")
    expect_true(is.finite(r$results_table$HR[1]))
    expect_gt(r$metadata$sample_size, 100L)
  })
})

test_that("CPAS_DATA_ROOT makes the checkout copies win over the bundled ones", {
  skip_if_not(dir.exists(file.path(.checkout, "data", "tcga")),
              "project checkout not available on this host")
  .with_tcga_env(root = .checkout, code = {
    t <- CanPAS::tcga_local_tables()
    expect_true(t$found)
    expect_identical(unname(t$source), c("CPAS_DATA_ROOT", "CPAS_DATA_ROOT"))
    expect_identical(normalizePath(t$clinical),
                     normalizePath(file.path(.checkout, "data/tcga/tcga_clinical.rda")))
    expect_identical(normalizePath(t$survival),
                     normalizePath(file.path(.checkout, "data/tcga/tcga_surv.rda")))
    # the same rows: the bundled copy is byte-identical to the checkout file
    e <- new.env(); load(t$clinical, envir = e)
    expect_equal(nrow(get("tcga_clinical", e)), 12591L)
    expect_equal(nrow(tcga_surv_table("TCGA-LUAD")), 574L)
  })
})

test_that("an explicit override wins over CPAS_DATA_ROOT and the bundle", {
  skip_if_not(dir.exists(file.path(.checkout, "data", "tcga")),
              "project checkout not available on this host")
  cli <- file.path(.checkout, "data/tcga/tcga_clinical.rda")
  .with_tcga_env(root = NULL, option_cli = cli, code = {
    t <- CanPAS::tcga_local_tables()
    expect_identical(t$clinical, cli)
    expect_identical(t$source[["clinical"]], "override")
    # the survival table is not overridden and comes from the package
    expect_identical(t$source[["survival"]], "package")
    expect_true(t$found)
  })
})

test_that("an unresolvable TCGA table stops with both options and the paths tried", {
  fake <- list(clinical = NA_character_, survival = NA_character_,
               source   = c(clinical = NA_character_, survival = NA_character_),
               found    = FALSE,
               tried    = list(clinical = c(override = "/no/such/dir/tcga_clinical.rda",
                                            CPAS_DATA_ROOT = "/no/such/root/data/tcga/tcga_clinical.rda"),
                               survival = c(package = "/no/such/pkg/extdata/tcga/tcga_surv.rda")))
  msg <- tryCatch(CanPAS:::.cpas_need_local(fake), error = function(e) conditionMessage(e))
  expect_match(msg, "Local CanPAS TCGA data file")
  expect_match(msg, "tcga_clinical.rda, tcga_surv.rda")
  expect_match(msg, "reinstall CanPAS")
  expect_match(msg, "CPAS_DATA_ROOT")
  expect_match(msg, "/no/such/dir/tcga_clinical.rda", fixed = TRUE)
  expect_match(msg, "/no/such/root/data/tcga/tcga_clinical.rda", fixed = TRUE)
  # a resolvable set passes without a message
  expect_true(CanPAS:::.cpas_need_local(CanPAS::tcga_local_tables()))
})

test_that("the historical CPAS_DATA_ROOT / TCGA_CLI_RDA bindings still exist", {
  expect_type(CanPAS:::CPAS_DATA_ROOT, "character")
  expect_match(CanPAS:::TCGA_CLI_RDA, "data/tcga/tcga_clinical[.]rda$")
  expect_match(CanPAS:::TCGA_SURV_RDA, "data/tcga/tcga_surv[.]rda$")
})

