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
#' @description One row per cohort in the CanPAS catalog - 197 rows: 147 GEO
#' series, 14 EMBL-EBI cohorts (ArrayExpress/BioStudies), 31 TCGA projects,
#' 3 CGGA glioma cohorts and 2 cBioPortal-hosted studies (\code{A5-PCPG},
#' \code{IMmotion150}) - with the accession, platform (GPL),
#' cancer type and the endpoint / sample-size columns described below. Used by
#' \code{\link{get_expr_data}} to find the platform of a dataset.
#' @docType data
#' @format A data.frame with 32 columns. \code{Accession}, \code{Type},
#' \code{GPL}, \code{N}, \code{SurvivalTypes} (raw endpoint tokens) and the
#' endpoint-family mapping added by the analysis plan B: \code{EndpointFamilies}
#' (browsing families), \code{EP_OS}, \code{EP_DSS}, \code{EP_DFS}, \code{EP_PFS},
#' \code{EP_MFS} (the concrete token available for each pooling family),
#' \code{EndpointPrimary} and \code{EndpointDerived} (TCGA-derived DFI/PFI).
#'
#' Sample-size columns. The delivered local artefacts under \code{data/} are the
#' authoritative layer; the mirror and this packaged object are downstream copies
#' of it, checked against it by \code{pipeline/R/16_verify_catalog_mirror.R}.
#' \code{n_surv} is the row count of the delivered survival table and
#' \code{n_events} the primary endpoint's events among those \code{N} analysable
#' samples, so it is not the event total of the survival table:
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
#'     twice); 33 of the 197 rows carry a label and the rest are \code{NA}.
#'     Computed by \code{pipeline/R/26_cohort_overlap.R} from sample titles.}
#'   \item{\code{Note}}{Free-text provenance and use flag for the cohort,
#'     written by the curation pipeline. 120 of the 197 rows carry no note and
#'     are stored as the literal \code{NA}; the other 77 are ASCII prose. A
#'     note records whichever of the following applies, in the pipeline's own
#'     wording: the platform and the probe-to-gene coverage of the expression
#'     table (\code{GPL...}; \code{probes->genes: ... n/m non-control = p\%});
#'     the sample and event accounting behind \code{N} and \code{n_events},
#'     with paired or replicate designs marked as such because they
#'     double-count patients; how the endpoint was derived (time unit, status
#'     coding, audit or patient-level check); an overlap group
#'     (\code{OVERLAPS ...}, \code{SAME SERIES on another platform - do not
#'     pool together ...}, or \code{TITLE COLLISION (not a shared-patient
#'     group)}); the registration batch date (\code{... batch 2026-09-24});
#'     and \code{[!]} where the note raises a caveat the user should read. The
#'     two rows whose \code{n_convention} is \code{clinical-record}
#'     (\code{GSE325123}, \code{GSE31312}) additionally state the
#'     expression-join-restricted \code{N}/\code{n_events} pair in a
#'     \code{[n_convention=clinical-record: ...]} clause, and
#'     \code{pipeline/R/16_verify_catalog_mirror.R} treats those two rows as
#'     documented convention exceptions, not rule violations. The App shows the
#'     note on the Datasets page and warns on the multi-dataset pages when an
#'     overlapping pair is selected.}
#' }
#'
#' Two bookkeeping columns are not used by the analysis functions:
#' \describe{
#'   \item{\code{X}}{Row index carried over from the pre-removal catalog (values
#'     1-213, 197 distinct). Harmless residue; kept so the packaged table stays
#'     cell-identical to \code{data/dataset_info.csv}.}
#'   \item{\code{method}}{Assay / data type recorded for the cohort:
#'     \code{"RNA"} (136), \code{"TCGA-RNAseq"} (33), \code{"RNA-seq"} (9),
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
