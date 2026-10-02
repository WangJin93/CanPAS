# cpas_reproducibility.R -----------------------------------------------------
# The reproducibility index (spec B7) and its pipeline-steps companion.
#
# The package ships the curation pipeline as scripts, and the project keeps the
# registers, the defect reports and the additional-file set outside the package.
# This file is the one place that says, offline and from an installed package,
# what that material is and where each piece lives: the 72 curation steps and
# the consolidated inst/pipeline/ file that provides each one, the exclusion
# register (109 records), the repair / defect records, the frozen-state record,
# and the additional-file layout.
#
# It is an index, not a data dump: the only table it reads is a small shipped CSV
# (a 72-row step inventory and a 14-row material list).  It never reads the
# catalog, the mirror or any cohort table, and it never downloads anything.

#' @title The 72 curation steps of the shipped pipeline
#' @description
#' Returns the step inventory of the curation pipeline that built the CanPAS
#' mirror and catalog: one row per step with its canonical run order, its
#' original script name, what the step does, and which consolidated file under
#' \code{inst/pipeline/} provides it.
#'
#' The inventory is a small CSV shipped with the package
#' (\code{inst/reproducibility/curation_steps.csv}); the consolidated file of
#' every step is resolved with \code{\link[base]{system.file}}, so the result is
#' usable from an installed package.  When the shipped entry point
#' \code{inst/pipeline/run_pipeline.R} is present, each row is cross-checked
#' against its step table as well, so drift between the inventory and the entry
#' point is visible rather than silent.
#' @param package Package to resolve the shipped files in (default
#'   \code{"CanPAS"}); mainly useful for tests.
#' @return A data frame with one row per step (72 rows):
#'   \item{\code{number}}{canonical run order, 1..72}
#'   \item{\code{id}}{the step id used by \code{run_pipeline.R} (e.g. \code{"07b"})}
#'   \item{\code{script}}{the original script name (e.g. \code{"07b_split_tnm.R"})}
#'   \item{\code{what_it_does}}{one-line description of the step}
#'   \item{\code{kind}}{\code{"numbered construction / repair step"} (63 rows) or
#'     \code{"helper / demo / validation script"} (9 rows)}
#'   \item{\code{consolidated_file}}{the \code{inst/pipeline/} file that provides
#'     the step}
#'   \item{\code{installed_path}}{the resolved absolute path of that file
#'     (\code{""} when the package does not ship it)}
#'   \item{\code{available}}{whether \code{installed_path} exists}
#'   \item{\code{matches_entry_point}}{\code{TRUE}/\code{FALSE} when the entry
#'     point could be read and the row agrees with its step table; \code{NA}
#'     when it could not be read}
#' @seealso \code{\link{cpas_reproducibility_index}}
#' @export
#' @examples
#' steps <- cpas_pipeline_steps()
#' nrow(steps)                                  # 72
#' table(steps$kind)
#' head(steps[, c("number", "script", "consolidated_file")])
cpas_pipeline_steps <- function(package = "CanPAS") {
  f <- system.file("reproducibility", "curation_steps.csv", package = package)
  if (!nzchar(f) || !file.exists(f))
    stop("The shipped step inventory (inst/reproducibility/curation_steps.csv) ",
         "was not found in package '", package, "'.", call. = FALSE)
  steps <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("number", "id", "script", "what_it_does", "consolidated_file", "kind")
  if (!all(need %in% colnames(steps)))
    stop("The shipped step inventory is malformed: missing ",
         paste(setdiff(need, colnames(steps)), collapse = ", "), ".", call. = FALSE)
  steps <- steps[, need, drop = FALSE]

  steps$installed_path <- vapply(steps$consolidated_file, function(g)
    system.file("pipeline", g, package = package), character(1))
  steps$available <- nzchar(steps$installed_path) & file.exists(steps$installed_path)

  # Cross-check against the shipped entry point's own step table.
  steps$matches_entry_point <- NA
  rp <- system.file("pipeline", "run_pipeline.R", package = package)
  if (nzchar(rp) && file.exists(rp)) {
    tab <- tryCatch(.cpas_entry_point_steps(rp), error = function(e) NULL)
    if (!is.null(tab)) {
      key <- paste(tab$id, tab$script, tab$group, sep = "\r")
      mine <- paste(steps$id, steps$script, steps$consolidated_file, sep = "\r")
      steps$matches_entry_point <- mine %in% key
    }
  }
  steps
}

# Evaluate only the `steps <- data.frame(...)` expression of the shipped entry
# point, so the check needs no side effects (run_pipeline.R itself sets a working
# directory and checks the data root).
.cpas_entry_point_steps <- function(path) {
  e <- parse(path, keep.source = FALSE)
  expr <- NULL
  for (ex in e) {
    if (is.call(ex) && is.symbol(ex[[1]]) && identical(as.character(ex[[1]]), "<-") &&
        is.symbol(ex[[2]]) && identical(as.character(ex[[2]]), "steps"))
      expr <- ex[[3]]
  }
  if (is.null(expr)) return(NULL)
  eval(expr, envir = new.env(parent = baseenv()))
}

#' @title Index of the reproducibility material
#' @description
#' A tidy data frame indexing the material behind the reported numbers: the 72
#' curation steps of the shipped pipeline (with the consolidated
#' \code{inst/pipeline/} file that provides each one), the exclusion register
#' (109 records), the repair and defect records, the frozen-state record, and the
#' Additional-file layout.
#'
#' The index is assembled from two small CSVs that ship with the package
#' (\code{inst/reproducibility/curation_steps.csv} and
#' \code{inst/reproducibility/materials.csv}) and resolved with
#' \code{\link[base]{system.file}}, so it works offline from an installed
#' package.  It is an index, not a data dump: it carries counts and locations,
#' never the registers' rows, and it never reads the catalog or the mirror.
#'
#' The registers, defect reports and additional files themselves live in the
#' project tree (the reproducibility bundle), not in the package; their
#' \code{location} is reported relative to the project root
#' (\code{CPAS_DATA_ROOT}, default \code{/home/Jingle/data/Project/CPAS}) and is
#' resolved to an absolute existing path in \code{resolved_source} when that root
#' is present on the machine.
#' @param package Package to resolve the shipped files in (default
#'   \code{"CanPAS"}); mainly useful for tests.
#' @param root Optional project root used to resolve the external material
#'   (default \code{CPAS_DATA_ROOT}, else the documented project path).
#' @return A data frame, one row per indexed item, with columns:
#'   \item{\code{section}}{\code{"curation_step"}, \code{"register"},
#'     \code{"repair_record"}, \code{"frozen_state"} or \code{"additional_file"}}
#'   \item{\code{number}}{1..72 for the curation steps, \code{NA} otherwise}
#'   \item{\code{item}}{step id, or the item name for the other sections}
#'   \item{\code{script}}{the original script of a step (\code{NA} otherwise)}
#'   \item{\code{what_it_does}}{what the step or item is}
#'   \item{\code{consolidated_file}}{the \code{inst/pipeline/} file providing a
#'     step (\code{NA} otherwise)}
#'   \item{\code{records}}{the number of records the item carries where that
#'     number is part of the record (109 for the exclusion register, 197 for
#'     Additional file 1, ...); \code{NA} when the item is not a record list}
#'   \item{\code{location}}{where the item lives, relative to the project root}
#'   \item{\code{package_path}}{the path inside the installed package, when the
#'     item ships with the package (\code{NA} otherwise)}
#'   \item{\code{installed_path}, \code{available}}{the resolved
#'     \code{system.file()} path and whether it exists}
#'   \item{\code{resolved_source}}{the absolute path of \code{location} under the
#'     project root when it exists there, else \code{NA}}
#' @seealso \code{\link{cpas_pipeline_steps}}, \code{\link{endpoint_semantics}},
#'   \code{\link{cohort_overlap}}, \code{\link{cpas_manifest}}
#' @export
#' @examples
#' idx <- cpas_reproducibility_index()
#' table(idx$section)
#' subset(idx, section == "register", select = c("item", "records", "location"))
#' subset(idx, section == "curation_step")[1:3, c("number", "script", "consolidated_file")]
cpas_reproducibility_index <- function(package = "CanPAS", root = NULL) {
  steps <- cpas_pipeline_steps(package = package)
  step_rows <- data.frame(
    section = "curation_step",
    number = as.integer(steps$number),
    item = as.character(steps$id),
    script = as.character(steps$script),
    what_it_does = as.character(steps$what_it_does),
    consolidated_file = as.character(steps$consolidated_file),
    records = NA_integer_,
    location = file.path("pipeline", "R", as.character(steps$script)),
    package_path = file.path("pipeline", as.character(steps$consolidated_file)),
    installed_path = as.character(steps$installed_path),
    available = as.logical(steps$available),
    stringsAsFactors = FALSE)

  f <- system.file("reproducibility", "materials.csv", package = package)
  if (!nzchar(f) || !file.exists(f))
    stop("The shipped material index (inst/reproducibility/materials.csv) was ",
         "not found in package '", package, "'.", call. = FALSE)
  mat <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("section", "item", "label", "records", "location", "what_it_is")
  if (!all(need %in% colnames(mat)))
    stop("The shipped material index is malformed: missing ",
         paste(setdiff(need, colnames(mat)), collapse = ", "), ".", call. = FALSE)
  mat_rows <- data.frame(
    section = as.character(mat$section),
    number = NA_integer_,
    item = as.character(mat$item),
    script = NA_character_,
    what_it_does = sprintf("%s -- %s", as.character(mat$label),
                           as.character(mat$what_it_is)),
    consolidated_file = NA_character_,
    records = as.integer(mat$records),
    location = as.character(mat$location),
    package_path = NA_character_,
    installed_path = NA_character_,
    available = NA,
    stringsAsFactors = FALSE)

  idx <- rbind(step_rows, mat_rows)
  rownames(idx) <- NULL

  if (is.null(root)) {
    root <- Sys.getenv("CPAS_DATA_ROOT", unset = "")
    if (!nzchar(root)) root <- "/home/Jingle/data/Project/CPAS"
  }
  root <- path.expand(root)
  resolved <- rep(NA_character_, nrow(idx))
  if (dir.exists(root)) {
    cand <- file.path(root, idx$location)
    ok <- file.exists(cand) | dir.exists(cand)
    resolved[ok] <- cand[ok]
  }
  idx$resolved_source <- resolved
  attr(idx, "n_steps") <- nrow(step_rows)
  attr(idx, "n_materials") <- nrow(mat_rows)
  attr(idx, "steps_source") <- system.file("reproducibility", "curation_steps.csv",
                                           package = package)
  attr(idx, "materials_source") <- f
  attr(idx, "project_root") <- root
  attr(idx, "note") <- paste0("index only: counts and locations, never the ",
                              "registers' rows; no catalog or mirror access")
  idx
}
