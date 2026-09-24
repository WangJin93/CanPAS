# TCGA source for CanPAS.  Design constraints (project convention):
#   * no runtime link to the MPAT project directory;
#   * no per-gene cache files and no full-matrix downloads: expression is
#     fetched per gene from UCSC Xena with the same query used by
#     UCSCXenaShiny (get_data_df), i.e. UCSCXenaShiny:::get_data();
#   * only the clinical / survival tables copied into <CPAS_DATA_ROOT>/data/tcga
#     are used;
#   * TCGA cohorts are exposed with the same "merged" structure as GEO cohorts,
#     so the whole package (COX/KM/meta) can consume them unchanged.

CPAS_DATA_ROOT <- Sys.getenv("CPAS_DATA_ROOT", "")   # no personal default: set it, or the helpers below stop with instructions
TCGA_CLI_RDA  <- file.path(CPAS_DATA_ROOT, "data/tcga", "tcga_clinical.rda")
TCGA_SURV_RDA <- file.path(CPAS_DATA_ROOT, "data/tcga", "tcga_surv.rda")

# Local TCGA tables are shipped outside the package: fail with an actionable message
# instead of an obscure downstream error when the root is not configured.
.cpas_need_local <- function(path) {
  if (!file.exists(path))
    stop("Local CanPAS data file not found: ", path,
         "\nSet it with Sys.setenv(CPAS_DATA_ROOT = \"<path to the CanPAS data root>\"); ",
         "TCGA analyses read <CPAS_DATA_ROOT>/data/tcga/*.rda.",
         call. = FALSE)
  invisible(TRUE)
}

# Supported TCGA projects — derived from the shipped catalog -----------------
# The catalog (data/dataset_info.rda) is the single authority: every row whose
# Accession starts with "TCGA-" is a TCGA cohort the helpers must accept, and
# its Type column is the canonical cancer-type label. Hard-coding a second list
# is what broke 16 of the 31 catalogued projects in 1.0.0, so the list below is
# only the fallback used when the catalog cannot be read; the object actually
# exported is refreshed from the catalog when the package is loaded (see
# .onLoad at the end of this file), and tests/testthat/test-tcga.R fails as soon
# as this fallback and the catalog disagree.
.cpas_tcga_retained_builtin <- c(
  BLCA = "Bladder Cancer", BRCA = "Breast Cancer", CESC = "Cervical Cancer",
  COAD = "Colorectal Cancer", GBM = "Glioma Cancer", LAML = "Leukemia Cancer",
  LGG = "Glioma Cancer", LIHC = "Liver Cancer", LUAD = "Lung Cancer",
  LUSC = "Lung Cancer", OV = "Ovarian Cancer", PAAD = "Pancreatic Cancer",
  PRAD = "Prostate Cancer", READ = "Colorectal Cancer", STAD = "Gastric Cancer",
  KIRC = "Kidney Cancer", THCA = "Thyroid Cancer", HNSC = "Head and Neck Cancer",
  SKCM = "Melanoma", KIRP = "Kidney Cancer", SARC = "Sarcoma",
  ESCA = "Esophageal Cancer", UCEC = "Endometrial Cancer", PCPG = "Pheochromocytoma",
  TGCT = "Testicular Cancer", THYM = "Thymoma", KICH = "Kidney Cancer",
  MESO = "Mesothelioma", UVM = "Uveal Melanoma", ACC = "Adrenocortical Cancer",
  UCS = "Uterine Carcinosarcoma"
)

# The project set actually used by the helpers: read once per session from the
# packaged catalog, with the built-in vector above as the fallback.
.cpas_tcga_projects <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    from_catalog <- tryCatch({
      di  <- .cpas_dataset_info()
      acc <- as.character(di$Accession)
      keep <- grepl("^TCGA-", acc)
      lab <- as.character(di$Type)[keep]
      proj <- toupper(sub("^TCGA-", "", acc[keep]))
      ok <- nzchar(proj) & !is.na(lab) & nzchar(lab)
      v <- lab[ok]
      names(v) <- proj[ok]
      v[!duplicated(names(v))]
    }, error = function(e) NULL)
    cache <<- if (length(from_catalog)) from_catalog else .cpas_tcga_retained_builtin
    cache
  }
})

#' @title TCGA projects available through the TCGA helpers
#' @description Named vector mapping TCGA project abbreviations to the cancer
#' type labels used in CanPAS. It is derived from the shipped catalog (every
#' \code{\link{dataset_info}} row whose \code{Accession} starts with
#' \code{"TCGA-"}; the \code{Type} column supplies the label), so it always
#' agrees with the catalog and every listed project is accepted by
#' \code{\link{tcga_project_dataset}}, \code{\link{tcga_surv_table}},
#' \code{\link{tcga_merged}}, \code{\link{tcga_get_expr}} and
#' \code{\link{canonical_type}}.
#' @format Named character vector, one entry per TCGA cohort in the catalog.
#' @export
tcga_retained <- .cpas_tcga_retained_builtin

#' @title Xena dataset id of a TCGA project
#' @description Builds the UCSC Xena sampleMap/HiSeqV2 dataset id for a TCGA
#' project, e.g. \code{"TCGA.LUAD.sampleMap/HiSeqV2"}.
#' @param dataset TCGA dataset id, either the project abbreviation
#'   (e.g. \code{"LUAD"}) or the accession (\code{"TCGA-LUAD"}).
#' @return A single string.
#' @export
tcga_project_dataset <- function(dataset) {
  project <- toupper(sub("^TCGA-", "", as.character(dataset)[1]))
  supported <- .cpas_tcga_projects()          # catalog-derived, never hard-coded
  if (!project %in% names(supported))
    stop("Unsupported TCGA project: ", project,
         ". Choose one of: ", paste(names(supported), collapse = ", "), ".")
  paste0("TCGA.", project, ".sampleMap/HiSeqV2")
}

#' @title Fetch one gene's TCGA expression from UCSC Xena
#' @description Downloads expression of a single gene for a TCGA project on
#' demand (per-gene query, no full-matrix download) using the internal
#' \code{UCSCXenaShiny:::get_data} helper, and returns a long data.frame of
#' tumour samples only.
#' @param dataset TCGA dataset id (\code{"LUAD"} or \code{"TCGA-LUAD"}).
#' @param gene Single gene symbol.
#' @param max_try Number of attempts before giving up (Xena may be flaky).
#' @param use_cache Reuse the local file cache. \code{NULL} (default) follows
#'   \code{getOption("CanPAS.cache", TRUE)} and \code{CANPAS_CACHE} /
#'   \code{CPAS_CACHE}; \code{TRUE} or \code{FALSE} forces the choice. Xena is
#'   queried one gene at a time and is intermittently unreachable, so a cached
#'   gene table is reused after \code{CanPAS.cache_ttl} seconds only when the
#'   live request fails.
#' @return data.frame with columns \code{sample}, \code{value}.
#' @examples
#' \dontrun{
#'    x <- tcga_gene_expr_df("LUAD", "TP53")   # long table: sample, value
#'    head(x)
#' 
#'    ## The project id and the Xena dataset id it maps to, plus which projects
#'    ## CanPAS retains.
#'    tcga_project_dataset("TCGA-LUAD")
#'    names(tcga_retained)
#' }
#' @export
tcga_gene_expr_df <- function(dataset, gene, max_try = 3, use_cache = NULL) {
  ds <- tcga_project_dataset(dataset)
  .xena <- function() {
    if (!requireNamespace("UCSCXenaShiny", quietly = TRUE))
      stop("Package 'UCSCXenaShiny' is required for TCGA expression queries. ",
           "Install with install.packages('UCSCXenaShiny').")
    gfun <- utils::getFromNamespace("get_data", "UCSCXenaShiny")
    v <- NULL
    for (i in seq_len(max_try)) {
      v <- tryCatch(gfun(ds, gene),
                    error = function(e) { Sys.sleep(3); NULL })
      if (!is.null(v)) break
    }
    if (is.null(v)) stop("Xena query failed for dataset ", ds, " gene ", gene, ".")
    v
  }
  key <- .cpas_cache_key("xena", ds, gene)
  v <- .cpas_cached(key, "xena_gene", .xena, use_cache = use_cache,
                    label = paste0("Xena ", ds, "/", gene))$value
  ok <- grepl("^TCGA-", names(v))
  long <- nchar(names(v)) >= 15L
  ok[long] <- ok[long] & suppressWarnings(as.integer(substr(names(v)[long], 14, 15))) < 10L
  data.frame(sample = names(v)[ok],
             value = suppressWarnings(as.numeric(v[ok])),
             stringsAsFactors = FALSE)
}

#' @title Multi-gene TCGA expression (wide format)
#' @description Fetches each requested gene for a TCGA project from UCSC Xena
#' and merges them into one wide data.frame (rows = samples).
#' @param dataset TCGA dataset id (\code{"LUAD"} or \code{"TCGA-LUAD"}).
#' @param genes Character vector of gene symbols.
#' @return data.frame with first column \code{sample} and one numeric column
#' per requested gene (columns follow the order of \code{genes}).
#' @examples
#' \dontrun{
#'    ## The expression step alone, for a TCGA project (same schema as
#'    ## get_expr_data() for GEO cohorts).
#'    e <- tcga_get_expr("LUAD", c("TP53", "GAPDH"))
#'    colnames(e$expr_data)
#'    head(e$expr_data)
#' }
#' @export
tcga_get_expr <- function(dataset = "LUAD", genes = "TP53") {
  if (missing(genes) || length(genes) == 0L) stop("'genes' is empty.")
  genes <- unique(as.character(genes))
  dfs <- lapply(genes, function(g) {
    d <- tcga_gene_expr_df(dataset, g)
    d$sample <- as.character(d$sample)
    colnames(d)[2] <- g
    d
  })
  m <- Reduce(function(a, b) merge(a, b, by = "sample", all = TRUE), dfs)
  m <- m[order(m$sample), , drop = FALSE]
  rownames(m) <- NULL
  m
}

#' @title TCGA clinical + survival table
#' @description Loads the local TCGA clinical and survival tables (copied into
#' the CanPAS data tree) and returns one row per tumour sample of a project, with
#' CanPAS-compatible columns: \code{ID}, \code{<endpoint>_status} /
#' \code{<endpoint>_time} (years), \code{age}, \code{sex}, \code{stage} and,
#' when present, \code{histology} and \code{grade}.
#' @param dataset TCGA dataset id (\code{"LUAD"} or \code{"TCGA-LUAD"}).
#' @param endpoints Endpoint prefixes to export.
#' @return data.frame (see description).
#' @examples
#' \dontrun{
#'    ## Local clinical table of a project: all endpoint families at once.
#'    s <- tcga_surv_table("LUAD")
#'    colnames(s)
#'    head(s)
#' 
#'    ## Only the endpoints the project really carries.
#'    s2 <- tcga_surv_table("LUAD", endpoints = c("OS", "DSS"))
#'    colnames(s2)
#' }
#' @export
tcga_surv_table <- function(dataset = "LUAD",
                            endpoints = c("OS", "DSS", "DFI", "PFI")) {
  tcga_project_dataset(dataset)          # validates the project abbreviation
  project <- toupper(sub("^TCGA-", "", as.character(dataset)[1]))
.cpas_need_local(TCGA_CLI_RDA); .cpas_need_local(TCGA_SURV_RDA)
  e <- new.env()
  load(TCGA_CLI_RDA, envir = e); load(TCGA_SURV_RDA, envir = e)
  cli <- get("tcga_clinical", envir = e); sv <- get("tcga_surv", envir = e)
  if (!all(c("sample", "type") %in% colnames(cli)))
    stop("tcga_clinical must contain 'sample' and 'type' columns.")
  m <- merge(cli, sv, by = "sample")
  m <- m[m$type == project & !duplicated(m$sample), , drop = FALSE]
  out <- data.frame(ID = as.character(m$sample), stringsAsFactors = FALSE)
  for (ep in endpoints) {
    tcol <- paste0(ep, ".time")
    if (all(c(ep, tcol) %in% colnames(m))) {
      st <- suppressWarnings(as.numeric(m[[ep]]))
      if (!all(st %in% c(0, 1, NA_real_)))
        stop("Endpoint ", ep, " of project ", project, " is not coded 0/1.")
      out[[paste0(ep, "_status")]] <- st
      out[[paste0(ep, "_time")]] <- suppressWarnings(as.numeric(m[[tcol]])) / 365.25
    }
  }
  out$age  <- suppressWarnings(as.numeric(m$age_at_initial_pathologic_diagnosis))
  out$sex  <- tolower(as.character(m$gender))
  out$stage <- as.character(m$ajcc_pathologic_tumor_stage)
  if ("histological_type" %in% colnames(m)) out$histology <- as.character(m$histological_type)
  if ("histological_grade" %in% colnames(m)) out$grade <- as.character(m$histological_grade)
  out
}

#' @title TCGA cohort merged with expression (CanPAS format)
#' @description Combines on-demand Xena expression (\code{\link{tcga_get_expr}})
#' with the local clinical/survival table (\code{\link{tcga_surv_table}}) into a
#' single data.frame whose schema matches GEO merged data: first column
#' \code{ID}, then \code{<type>_time} / \code{<type>_status} and clinical
#' columns, then one numeric column per gene.
#' @param dataset TCGA dataset id (\code{"LUAD"} or \code{"TCGA-LUAD"}).
#' @param genes Gene symbols.
#' @param type Endpoint used downstream (must exist for the project).
#' @return data.frame.
#' @examples
#' \dontrun{
#'    ## TCGA keeps its clinical columns, so covariates are usable directly.
#'    d <- tcga_merged("LUAD", c("TP53", "GAPDH"), type = "OS")
#'    colnames(d)
#' 
#'    r <- COX_analysis(d, type = "OS", cont_Variates = c("TP53", "age"),
#'                       cate_Variates = c("sex", "stage"), method = "uni")
#'    head(r$results_table)
#' }
#' @export
tcga_merged <- function(dataset = "LUAD", genes = c("TP53"), type = "OS") {
  acc <- tcga_project_dataset(dataset)   # validates the project abbreviation
  acc <- sub("^TCGA\\.", "TCGA-", sub("\\.sampleMap.*$", "", acc))
  tok <- endpoint_resolve(acc, type)
  if (is.na(tok)) tok <- as.character(type)[1]
  ex <- tcga_get_expr(dataset, genes)
  sv <- tcga_surv_table(dataset)
  d <- merge(ex, sv, by.x = "sample", by.y = "ID")
  d <- d[!duplicated(d$sample), , drop = FALSE]
  colnames(d)[1] <- "ID"
  sc <- paste0(tok, "_status"); tc <- paste0(tok, "_time")
  if (!all(c(sc, tc) %in% colnames(d)))
    stop("endpoint ", tok, " not in TCGA ", acc)
  attr(d, "endpoint") <- tok
  attr(d, "family") <- endpoint_family(tok)
  d
}

#' @title Canonical CanPAS cancer type of a TCGA project
#' @param dataset TCGA dataset id (\code{"LUAD"} or \code{"TCGA-LUAD"}). A GEO or
#' CGGA accession is not a TCGA id and returns \code{NA}.
#' @return The canonical cancer-type label or \code{NA_character_} when the
#' project is unknown.
#' @examples
#'   ## Real TCGA accessions from the catalog, and the type label CanPAS reports.
#'   canonical_type("TCGA-LUAD")
#'   canonical_type("TCGA-ACC")
#'   canonical_type("TCGA-XXXX")   # unknown project -> NA
#' @export
canonical_type <- function(dataset) {
  project <- toupper(sub("^TCGA-", "", as.character(dataset)[1]))
  supported <- .cpas_tcga_projects()          # catalog-derived, never hard-coded
  if (project %in% names(supported)) unname(supported[[project]]) else NA_character_
}

#' @title Short dataset name of a TCGA project
#' @description Prefixes a TCGA project id with \code{"TCGA-"}. Only the first
#' element of \code{dataset} is used, and the input must be a TCGA id
#' (\code{"LUAD"} or \code{"TCGA-LUAD"}) — a GEO accession such as
#' \code{"GSE13507"} would come back as \code{"TCGA-GSE13507"}, which is why
#' GEO cohorts are labelled by their catalog accessions instead.
#' @param dataset TCGA dataset id (\code{"LUAD"} or \code{"TCGA-LUAD"}).
#' @return A single string, \code{"TCGA-<PROJECT>"} (only the first element of
#' \code{dataset} is used).
#' @examples
#'   ## Short labels used in figures and tables (TCGA ids only).
#'   short_name("LUAD")
#'   short_name("TCGA-LUAD")
#' @export
short_name <- function(dataset) paste0("TCGA-", toupper(sub("^TCGA-", "", as.character(dataset)[1])))

#' @title Unified CanPAS cohort reader (GEO or TCGA)
#' @description Returns one merged data.frame for a cohort regardless of the
#' source: datasets prefixed \code{"TCGA-"} are handled through the TCGA
#' helpers (Xena + local clinical tables); anything else is treated as a GEO
#' accession served by the CanPAS API. The output schema is identical for both
#' sources, so downstream functions of the package work unchanged.
#' @param dataset Dataset identifier, e.g. \code{"GSE14814"} or
#' \code{"TCGA-LUAD"}.
#' @param genes Gene symbols to retrieve.
#' @param type Endpoint family (OS/DSS/DFS/PFS/MFS) or concrete token; families
#'   are resolved to the token available in the cohort and recorded as the
#'   \code{endpoint} attribute of the returned data.frame.
#' @param process_duplicates Probe collapsing rule used for GEO queries.
#' @param clin If \code{TRUE}, clinical covariates supplied for that cohort are
#'   appended: \code{age}, \code{sex}, \code{stage}, \code{T}, \code{N},
#'   \code{M}, \code{grade}, \code{histology}, whichever the cohort carries.
#'   GEO cohorts take them from the mirrored clinical table that
#'   \code{\link{merge_surv_expr}()} returns in \code{raw_surv_data}; TCGA
#'   cohorts already include them. Columns missing for a cohort are simply
#'   absent. The mirror stores these columns as text, so they are converted:
#'   a column whose values all parse as numbers (\code{age}, often \code{T},
#'   \code{N}, \code{M}) becomes numeric, every other one becomes a factor.
#'   Default \code{FALSE} (survival + genes only).
#' @return data.frame: first column \code{ID} followed by survival columns, one
#' numeric column per gene, and (with \code{clin = TRUE}) the clinical
#' covariates available for that cohort.
#' @examples
#' \dontrun{
#'    ## One cohort, one call: the endpoint token is resolved from the catalog,
#'    ## expression is fetched (mirror for GEO, Xena for TCGA), survival merged.
#'    d <- cohort_merged("GSE14814", c("GAPDH", "ACTB"), type = "OS")
#'    attr(d, "family"); attr(d, "endpoint")
#'    colnames(d)
#' 
#'    ## clin = TRUE adds the clinical covariates that cohort really carries, so
#'    ## a multivariable model can use them directly.
#'    d2 <- cohort_merged("GSE13507", "GAPDH", type = "OS", clin = TRUE)
#'    colnames(d2)
#' 
#'    ## A TCGA project goes through the same function and has the same schema.
#'    tcga <- cohort_merged("TCGA-LUAD", c("TP53", "GAPDH"), type = "OS")
#'    colnames(tcga)[1:6]
#' }
#' @export
cohort_merged <- function(dataset, genes = "TP53", type = "OS",
                          process_duplicates = "max", clin = FALSE) {
  tok <- endpoint_resolve(dataset, type)
  if (is.na(tok)) tok <- as.character(type)[1]      # not in catalog / custom data
  if (startsWith(dataset, "TCGA-")) {
    # the TCGA helper already carries the clinical columns of the project
    d <- tcga_merged(sub("^TCGA-", "", dataset), genes, tok)
  } else {
    m <- merge_surv_expr(dataset,
                         get_expr_data(dataset, genes,
                                       process_duplicates = process_duplicates))
    d <- m$merged_data
    if (isTRUE(clin) && !is.null(m$raw_surv_data)) {
      cl <- m$raw_surv_data
      keep <- intersect(c("age", "sex", "stage", "T", "N", "M", "grade", "histology"),
                        colnames(cl))
      keep <- setdiff(keep, colnames(d))            # never duplicate a column
      if ("ID" %in% colnames(cl) && length(keep)) {
        cl <- cl[, c("ID", keep), drop = FALSE]
        # the mirror returns every clinical column as text; restore the evident
        # type so that downstream code sees numbers as numbers ("55.4" -> 55.4)
        # and labels as factors, exactly as the analysis functions expect
        for (v in keep) {
          x <- cl[[v]]
          if (!is.character(x)) next
          num <- suppressWarnings(as.numeric(x))
          cl[[v]] <- if (any(!is.na(num)) && !any(is.na(num) & !is.na(x))) num else factor(x)
        }
        d <- merge(d, cl, by = "ID", all.x = TRUE)
      }
    }
  }
  attr(d, "endpoint") <- tok
  attr(d, "family") <- endpoint_family(tok)
  d
}

# Load-time refresh of the exported TCGA project table ------------------------
# `tcga_retained` is a documented, exported object, so it cannot be a function
# call; rebinding it here from the catalog keeps it equal to the set the
# helpers accept (and to the catalog) for every session, without a second
# hard-coded list. The built-in vector stays as the fallback if the catalog is
# unreadable, which cannot make loading fail.
.onLoad <- function(libname, pkgname) {
  projects <- tryCatch(.cpas_tcga_projects(), error = function(e) NULL)
  if (length(projects)) assign("tcga_retained", projects, envir = asNamespace(pkgname))
  invisible()
}
