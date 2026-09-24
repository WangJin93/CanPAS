#' @title Gene-symbol to Entrez mapping (internal reference)
#' @description Mapping of gene symbols to Entrez gene ids (and Ensembl gene
#' ids when available) used to resolve probes across GEO platforms. Loaded
#' lazily by the package; there is no need to attach it manually.
#' @docType data
#' @format A data.frame with columns \code{Symbol}, \code{Ensembl_gene} and
#' \code{gene_id}.
#' @keywords datasets
"ID_map"

#' @title CanPAS dataset catalog
#' @description One row per cohort in the CanPAS catalog — 179 rows: 143 GEO
#' series, 31 TCGA projects, 3 CGGA glioma cohorts and 2 cBioPortal-hosted
#' studies (`A5-PCPG`, `IMmotion150`) — with the accession, platform (GPL),
#' cancer type and the endpoint / sample-size columns described below. Used by
#' \code{\link{get_expr_data}} to find the platform of a dataset.
#' @docType data
#' @format A data.frame with 26 columns. \code{Accession}, \code{Type},
#' \code{GPL}, \code{N}, \code{SurvivalTypes} (raw endpoint tokens) and the
#' endpoint-family mapping added by the analysis plan B: \code{EndpointFamilies}
#' (browsing families), \code{EP_OS}, \code{EP_DSS}, \code{EP_DFS}, \code{EP_PFS},
#' \code{EP_MFS} (the concrete token available for each pooling family),
#' \code{EndpointPrimary} and \code{EndpointDerived} (TCGA-derived DFI/PFI).
#'
#' Sample-size columns (recomputed from the mirror; see
#' \code{pipeline/R/19_fix_catalog_N.R}):
#' \describe{
#'   \item{\code{N}}{Analysable sample count: samples with both expression data
#'     in the mirror and non-missing \code{time} and \code{status} for the
#'     dataset's \code{EndpointPrimary}. \code{NA} when the dataset has no
#'     expression table or no primary endpoint annotation, i.e. gene-level
#'     analysis is not possible.}
#'   \item{\code{n_expr}, \code{n_surv}}{Samples in the mirror expression table
#'     and rows in the mirror survival table (\code{NA} when absent).}
#'   \item{\code{n_events}}{Events (\code{status = 1}) for
#'     \code{EndpointPrimary}, counted inside the analysable set.}
#'   \item{\code{n_OS}, \code{n_DSS}, \code{n_DFS}, \code{n_PFS}, \code{n_MFS}}{
#'     Analysable sample count per pooling family, resolved to the token each
#'     dataset actually provides (GSE31210 contributes to \code{n_DFS} through
#'     RFS). \code{0} when the family is unavailable.}
#'   \item{\code{expr_in_mirror}}{Whether an expression table exists in the
#'     mirror.}
#' }
#'
#' Two columns describe how a cohort may be used:
#' \describe{
#'   \item{\code{CohortGroup}}{Label of the group of cohorts that share
#'     patients (same study on another platform, or the same series deposited
#'     twice). \code{NA} when the cohort does not overlap any other.
#'     Computed by \code{pipeline/R/26_cohort_overlap.R} from sample titles.}
#'   \item{\code{Note}}{Human-readable flag: the overlapping partners and how
#'     many patients are shared (\code{"OVERLAPS ... do not pool together"}), or
#'     \code{"expression only — no survival table, cannot be analysed"}, or
#'     \code{"no annotated endpoint"}. The App shows it on the Datasets page and
#'     warns on the multi-dataset pages when an overlapping pair is selected.}
#' }
#'
#' Two bookkeeping columns are not used by the analysis functions:
#' \describe{
#'   \item{\code{X}}{Row index carried over from the pre-removal catalog (values
#'     1–195, 179 distinct). Harmless residue; kept so the packaged table stays
#'     cell-identical to \code{data/dataset_info.csv}.}
#'   \item{\code{method}}{Assay / data type recorded for the cohort:
#'     \code{"RNA"} (121), \code{"TCGA-RNAseq"} (31), \code{"RNA-seq"} (8),
#'     \code{"array"} (1), \code{"SRA"} (1), and \code{NA} for the 17
#'     lung-cancer GEO cohorts whose method was not recorded.}
#' }
#' Earlier releases recorded in \code{N} the planned/expression cohort size,
#' which overstated the analysable sample count for many datasets (GSE31210 was
#' listed as 133 while 226 tumours are analysable).
#' @keywords datasets
#' @examples
#'   ## The catalog is the starting point of every real analysis: it records the
#'   ## cohorts mirrored in CanPAS and, per cohort, which endpoint families can be
#'   ## analysed and how many samples are usable.
#'   head(dataset_info[, c("Accession", "Type", "N", "EndpointFamilies")])
#'
#'   ## Cohorts of one cancer type with a DFS-family endpoint, largest first.
#'   lung <- dataset_info[dataset_info$Type == "Lung Cancer" &
#'                         grepl("DFS", dataset_info$EndpointFamilies), ]
#'   lung[order(-lung$N), c("Accession", "N", "EndpointFamilies")]
#'
#'   ## The catalog drives the analysis helpers, e.g. the endpoint token.
#'   endpoint_resolve(lung$Accession[1], "DFS")
"dataset_info"
