# CanPAS — Cancer Prognosis Analysis Suite

`CanPAS` is an R package for survival analysis of public cancer cohorts.
It serves expression matrices and survival tables from a curated
GEO/CGGA mirror over a public REST API, from EMBL-EBI
ArrayExpress/BioStudies deposits and from TCGA (UCSC Xena expression
plus clinical/survival tables bundled with the package), and exposes one
uniform **merged-data schema** to Kaplan–Meier, Cox, time-dependent ROC,
cross-cohort meta-analysis, pooled Kaplan–Meier and competing-risks
analyses. A bundled Shiny application drives the same exported
functions.

What it adds to the usual single-cohort workflow:

- **Endpoint families with per-cohort token resolution.** Cohorts are
  requested by family (OS, DSS, DFS, PFS, MFS); the token a cohort
  actually reports (`RFS`, `DFS`, `DFI`, …) is resolved per cohort and
  returned with every result, and pooling never crosses families.
- **Five statistical safeguards on by default.** Cut-point rule and
  threshold, events-per-variable warning, proportional-hazards test,
  counted excluded rows, and the resolved endpoint token are printed
  with every result rather than on request.
- **REML by default, with DL and HK beside it.** `cpas_meta(method = )`
  estimates tau^2 by restricted maximum likelihood (`"REML"`, the
  default), DerSimonian-Laird (`"DL"`, the previous default, kept for
  continuity) or Hartung-Knapp-Sidik-Jonkman (`"HK"`); both prediction
  intervals are always returned.
- **Fail-safe rather than silent repair.** `auto_repair = FALSE` is the
  default in the Cox / multivariable paths: a model that would need a
  covariate dropped, a level merged or patients excluded is reported as
  `not estimable` with its reasons and offending terms instead of being
  quietly reduced.
- **An analysis manifest for every result.** `cpas_manifest(result)`
  (also `result$manifest`) records the cohorts, the resolved token of
  each, the pooling classes pooled, what was dropped and why, the
  cut-point rule, the meta-analysis method, tau^2, I^2, both prediction
  intervals, the versions and the timestamp - printable and convertible
  with [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
- **Numerical verification against reference implementations**
  (`survival`, `metafor`, `cmprsk`), including a competing-risks
  variance defect that the comparison exposed and this release corrects.
- **A catalog, not just a downloader**: 197 catalogued cohorts in one
  schema, with sample sizes defined as analysable patients and with
  patient-overlap groups recorded.

> Naming note: the package was originally created under the name “Cancer
> Patient Aftercare System”. The content is prognostic survival analysis
> of public cohorts, not aftercare intervention, so it is now the
> “Cancer Prognosis Analysis Suite” (the acronym is unchanged).

**Status**: version 1.0.0 (first release) · licence MIT · `R CMD check`
0 errors / 0 warnings / 1 note (host library notice) · a tool paper is
in preparation.

**Documentation site**: <https://wangjin93.github.io/CanPAS/> (function
reference, both articles and the 1.0.0 changelog).

## Installation

``` r
# from GitHub
remotes::install_github("WangJin93/CanPAS")

# or from the release tarball
install.packages("CanPAS_1.0.0.tar.gz", repos = NULL, type = "source")
# inside a source checkout: R CMD INSTALL .
```

`Imports` are all on CRAN. Optional packages enable extra features:

| Package                                                     | Enables                                                                                                               |
|-------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| `shiny`, `DT`, `bs4Dash`, `shinyWidgets`, `shinycssloaders` | the bundled application                                                                                               |
| `UCSCXenaShiny`, `UCSCXenaTools`                            | TCGA expression through UCSC Xena                                                                                     |
| `flextable`                                                 | three-line summary table in the Cox page and its Word export (without it the app shows the same table and offers CSV) |
| `maxstat`                                                   | search-adjusted approximate p-value for the cut-point search                                                          |
| `patchwork`                                                 | composing both pooling routes into one figure                                                                         |
| `cmprsk`                                                    | independent Fine–Gray reference in the test suite                                                                     |

``` r
install.packages(c("dplyr", "ggplot2", "ggtext", "gridExtra", "jsonlite",
                   "survival", "survminer", "survivalROC", "flextable"))
install.packages(c("shiny", "DT", "bs4Dash", "shinyWidgets", "shinycssloaders",
                   "flextable", "maxstat", "patchwork", "cmprsk"))
```

## Configuration

| Variable / option                                | Used for                                                                                                                                                                                                                                                                                                                                |
|--------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `CPAS_DATA_ROOT`                                 | **optional override**, needed only to read cohort clinical covariates in the app or to use a project checkout’s `<root>/data/tcga/*.rda` instead of the tables bundled with the package. Unset, TCGA cohorts work out of the box (the clinical/survival tables ship inside CanPAS) and the GEO/CGGA mirror and UCSC Xena work as usual. |
| `CPAS_DB_PASSWORD`                               | the curation pipeline shipped under `inst/pipeline/`, for the steps that write to the mirror (never hard-coded in the package)                                                                                                                                                                                                          |
| `CPAS_PIPELINE_DIR`                              | optional override for the directory holding the consolidated pipeline scripts (defaults to `system.file("pipeline", package = "CanPAS")`)                                                                                                                                                                                               |
| `options(CanPAS.cache_dir=)`, `CANPAS_CACHE_DIR` | where downloaded answers are cached                                                                                                                                                                                                                                                                                                     |
| `options(CanPAS.cache=)`, `CANPAS_CACHE`         | switch the cache off                                                                                                                                                                                                                                                                                                                    |
| `options(CanPAS.cache_ttl=)`                     | cache lifetime in seconds (30 days by default)                                                                                                                                                                                                                                                                                          |

``` r
Sys.setenv(CPAS_DATA_ROOT = "/path/to/CanPAS-data-root")   # optional override; not needed for TCGA
```

GEO and CGGA data are read from the public CanPAS mirror API
(`https://www.jingege.wang/bioinformatics/CPAS/api.php`); it is a live
service and rate-limits bulk sweeps (HTTP 429 → the client backs off and
retries).
[`get_data()`](https://wangjin93.github.io/CanPAS/reference/get_data.md)
accepts a different `base_url` if you host your own mirror.

## Quick start

``` r
library(CanPAS)

## --- one cohort: expression + survival ------------------------------------
e <- get_expr_data("GSE14814", genes = c("GAPDH", "ACTB"))
m <- merge_surv_expr("GSE14814", e)
r <- COX_by_genes(m$merged_data, type = "OS", genes = c("GAPDH", "ACTB"))
r$results_table                       # HR / 95% CI / p / BH-adjusted p

plot_km(m$merged_data, type = "OS", marker = "GAPDH", pval = TRUE)   # High red, Low green
plot_roc(m$merged_data, type = "OS", marker = "GAPDH", predict.time = 3)

## --- endpoint families: OS / DSS / DFS / PFS / MFS -------------------------
endpoint_options("TCGA-LGG")          # families, tokens, derived flags
endpoint_resolve("GSE31210", "DFS")   # -> "RFS": the token that cohort provides
plot_km(m$merged_data, type = "DFS", marker = "GAPDH")

## --- cross-cohort two-stage meta-analysis ----------------------------------
## method = "REML" is the default; "DL" reproduces the earlier estimator and
## "HK" gives the Hartung-Knapp-Sidik-Jonkman interval.
res <- cpas_meta(datasets = c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
res$per_dataset    # per-cohort HR, n, events, resolved token and pooling class
res$pooled         # pooled HR, Q, I2, tau2 and BOTH prediction intervals
res$pooled$pi_lower; res$pooled$pi_alt_lower   # primary and alternative interval
res$pooled$pi_method                           # the rule used ("t" by default)
res$pooled$pi_rule                             # and how that interval is built
## pi_method = "normal" or "HK" selects the other two rules; only the interval
## changes, never the pooled estimate or the heterogeneity statistics
cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS",
          pi_method = "normal")$pooled$pi_lower
res$pooled$pooling_classes                     # the classes actually pooled
res$pooled$overlap_pairs                       # pairs known to share patients
cpas_manifest(res)                             # every setting behind the number
loo_meta(res)      # leave-one-out sensitivity (same pooling method)

## strict mode: only cohorts whose token IS the DFS definition
res_exact <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS",
                       pooling = "exact")
res_exact$per_dataset$pooling_class            # what was kept, and why

## shared patients: warn (default), refuse, or keep the larger cohort
cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS",
          overlap = "dedupe")$excluded          # what dedupe dropped, and why

## --- pooled Kaplan-Meier with the two within-cohort cut rules ---------------
merged <- lapply(c("GSE31210", "GSE37745"), function(a) cohort_merged(a, "GAPDH", type = "RFS"))
names(merged) <- c("GSE31210", "GSE37745")
km <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS", method = "both")
km$cutpoint                                  # the rule and the per-cohort thresholds
km2 <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS",
                      cut = "top_pct", top_pct = 25)   # highest quarter of every cohort
km2$empty_cohorts; km2$skipped_cohorts        # cohorts left out, with the reason

## --- per-cohort expression diagnostics before pooling ----------------------
bd <- batch_diagnostics(c("GSE31210", "GSE37745"), genes = "GAPDH", type = "RFS")
bd$table[, c("dataset", "n", "median", "iqr", "cut_median", "cut_top_pct")]
bd$batch            # how far apart the absolute expression scales are
plot_batch_diagnostics(bd)   # the distributions side by side, cut-points marked

## --- blinded endpoint adjudication -----------------------------------------
## A rating sheet with the cohort, the token and the evidence text only.
sheet <- endpoint_adjudication(n = 12, raters = 2, seed = 20260930)
attr(sheet, "pooling_class_vocabulary")
## ... the raters fill it in; then:
endpoint_agreement(sheet)$per_field[, c("field", "raw_agreement", "kappa")]

## --- reproducibility index (67 steps, registers, additional files) ---------
idx <- cpas_reproducibility_index()
table(idx$section)
subset(idx, section == "curation_step")[, c("number", "script", "consolidated_file")]

## --- competing risks -------------------------------------------------------
cr <- competing_risk_COX(df, time = "OS_time", status = "cause",
                         covariates = c("marker", "age"), etype = 1)
cr$cause_specific; cr$subdistribution
plot_cif(cif_fit(df, time = "OS_time", status = "cause", times = c(1, 3, 5)))

## --- the Shiny application -------------------------------------------------
run_cpas_app()
```

## The catalog

`data(dataset_info)` ships the catalog used by the app and by the paper:
**197 cohorts — 145 GEO, 14 EMBL-EBI (ArrayExpress/BioStudies), 33 TCGA
projects, 3 CGGA and 2 cBioPortal-hosted studies — across 29 cancer
types**, and every row carries a resolved endpoint; together they
contribute **38,981 analysable samples** (median 149 per cohort, range
28–1,210). Two source groups are deliberately kept out of GEO: the
fourteen EMBL-EBI cohorts (`E-MTAB-*`, `E-TABM-*`, `E-MEXP-*`) are
ArrayExpress/BioStudies deposits retrieved from their own SDRF
annotation and processed matrices rather than through GEO’s mirror
actions, and the two cBioPortal-hosted cohorts (A5-PCPG,
Pheochromocytoma; IMmotion150, Kidney Cancer) are studies whose clinical
and expression files are deposited together; counting either as GEO
would give 159 or 147 instead of 145. A small number of cohorts carry
only 28–40 patients and are flagged as such in the catalog `Note` field
(the admission gate was relaxed from \> 50 to \>= 30 patients per cohort
for these additions); one NanoString cohort sits below even that relaxed
gate — `GSE205209` (endometrial cancer, `GPL27956`) was recruited as 29
patients, one of whom has no usable expression array, so it is analysed
at **N = 28** (21 OS events; the delivered survival table has 29 rows)
and was admitted by an explicit author decision recorded in the row’s
`Note`. Sample size means analysable patients: expression data plus a
non-missing time and status for the cohort’s primary endpoint.

### Column semantics

The size columns are not interchangeable, and one rule fixes all of
them:

| Column       | Meaning                                                                                                                                                                                                                                                                                                                                          |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `N`          | **Analysable samples** — present in both the delivered expression and survival tables and carrying a usable *primary* endpoint (status and time non-missing).                                                                                                                                                                                    |
| `n_events`   | Events of the **primary** endpoint **among those `N` samples** — not the event total of the survival table.                                                                                                                                                                                                                                      |
| `n_<family>` | The same rule applied **per family** — samples with a usable endpoint *of that family* — counted **independently of `N`**, so `n_OS`, `n_DSS`, `n_DFS`, `n_PFS` and `n_MFS` may each **exceed** `N`. `NA` when the cohort carries no endpoint in that family.                                                                                    |
| `n_surv`     | Row count of the delivered survival table, **including** rows with no usable endpoint.                                                                                                                                                                                                                                                           |
| `n_expr`     | Sample columns of the delivered expression table; `NA` for the 33 TCGA cohorts, whose expression is fetched on demand from UCSC Xena and whose clinical/survival tables ship inside the package ([`tcga_local_tables()`](https://wangjin93.github.io/CanPAS/reference/tcga_local_tables.md) reports the resolved copies) and are never mirrored. |

`N` and `n_events` therefore describe one endpoint (the primary one),
while `n_<family>` describes each family separately. The independence is
real, not an artefact of a single row: across the 197-row catalog,
`n_<family>` exceeds `N` in **6 cohort-by-family cells** — OS 0, DSS 0,
DFS 1 (`GSE22226_GPL1708` 129 \> 125), PFS 4 (`TCGA-BLCA` 426 \> 425,
`TCGA-LUSC` 543 \> 542, `TCGA-STAD` 445 \> 443, `TCGA-SKCM` 456 \> 455)
and MFS 1 (`GSE45255` 136 \> 134). The verifier therefore does **not**
assert `n_<family> <= N`: these cells are legitimate. The clearest
single-row witness is `GSE205209`: `N = 28`, `n_surv = 29` (the
delivered survival table keeps the recruited patient with no usable
array), `n_OS = 28` and `n_events = 21`, since the 29-row table carries
22 deaths but one of them is that non-analysable patient.

Two rows are deliberately registered at the **audited clinical-record**
level rather than at the join-restricted level, and say so in the
`n_convention` column: `GSE325123` (`N = 105`, 62 events, versus 102/60
join-restricted) and `GSE31312` (`N = 475`, 172 events, versus 470/170).
`19_fix_catalog_N.R` proposed changing these two to the join-restricted
counts; that proposal was reviewed and **rejected**, so do not re-apply
it. The remaining 195 rows carry `n_convention = "join-restricted"`, and
six further columns record the curation decisions: `n_join_dropped` and
`join_drop_reason` (survival rows with a usable primary endpoint that
have no expression row; non-zero for `GSE54460`, `GSE108474` and
`GSE205209`), `admission_gate` and `gate_decision` (`>=30` for 16 rows
and `>50` for 181; `author-decision` only for `GSE205209`), and
`overlap_group` (33 rows carry a shared-patient group).

**Which layer is authoritative:** the delivered local artefacts — the
expression `.rds` under `data/expr/` and the survival tables under
`data/processed/surv/` — are the source of truth. The MySQL mirror and
the packaged `dataset_info` object are downstream copies, and
`pipeline/R/16_verify_catalog_mirror.R` checks both against the local
artefacts (at the 197-row state it reports 0 errors / 0 warnings / 0
info, 197/197 rows fully usable). It raises an **ERROR** when a row is
endpoint-annotated, or has a survival table, yet every `EP_*` column is
`NA`, and a second **ERROR** when the mirror holds an expression table
whose platform/GPL table is missing.

This convention is why a survival table with 522 rows can carry
`n_surv = 522`, `N = 476` and `n_events = 397` at the same time
(GSE108474): 397 is the number of events among the 476 analysable
samples, while 404 is the event count over all 522 rows — a different
quantity, and not the one the catalog records. Patient-overlap groups
are recorded (30 pairs in 14 groups recomputed against this catalog by
`pipeline/R/26_cohort_overlap.R`, with `data(dataset_info)$CohortGroup`
and `$Note` carrying the result). The most recent pair is GSE1379 ×
GSE1378: the same 60 breast patients deposited twice (whole-tissue
sections and microdissected cells), recorded in
`pipeline/ref/cohort_overlap_seed.csv` and in the catalog as
`CohortGroup = GSE1378(+1)`; both series carry their DFS time and status
in the series-matrix `!Sample_description` free text rather than in
`!Sample_characteristics_ch1`. One further title match — the
GSE25066–GSE32918 pair — is a sample-title collision rather than shared
patients, since GSE32918’s titles are panel replicate codes for 172
patients and its genuine duplicate deposit (GSE69051) is not catalogued,
so it is recorded as a note instead of a group
(`pipeline/ref/cohort_overlap_exclude.csv`). The multi-dataset pages
warn when a selection contains two members of one group.

## Endpoint families

A family groups tokens that answer the same clinical question, so
cohorts reporting different token names can be pooled. The token each
cohort contributes is resolved per cohort and always reported.

| Family | Tokens pooled      | Cohorts |
|--------|--------------------|---------|
| OS     | OS                 | 142     |
| DSS    | DSS, CSS, BCSS     | 41      |
| DFS    | DFS, RFS, EFS, DFI | 98      |
| PFS    | PFS, PFI           | 49      |
| MFS    | MFS, DRFS          | 17      |

`DFI` and `PFI` (TCGA) are derived from the original endpoint fields and
are flagged as derived. Pooling inside a family and never across
families is enforced by the package; a mixed-token pool warns, naming
each cohort and its token. A separate browsing vocabulary widens the
progression family to metastasis endpoints, so the Datasets page can
show a PFS cohort count (66) larger than the pooling count (49).

## Statistical safeguards (read before quoting a result)

- **Cut point.** The default is a within-cohort median split; a custom
  rule is a top percentage (entering 25 takes the highest quarter of
  patients) or an absolute log2 threshold. The app’s `Auto` mode
  searches 13 percentile cut points and keeps the largest log-rank
  statistic, so its naive p-value is optimistically biased: with a
  marker generated independently of survival, p \< 0.05 in 23.1–26.7% of
  2,000 replicates per sample size versus 4.45–5.90% for the median
  split, and 2.04–3.53% with the `maxstat` search-adjusted
  approximation. The app prints the number of candidates searched and
  the adjusted p-value when `maxstat` is installed. Pooled analyses
  offer only the two percentile rules; a per-cohort cut-point search is
  deliberately not offered.
- **Sample sizes.** Reported n is the number of rows entering the model;
  rows excluded for a missing endpoint time, status or marker are
  counted and printed.
- **Events per variable.** Multivariable Cox models warn below 10 events
  per variable.
- **Proportional hazards.**
  [`COX_analysis()`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md)
  returns the global `cox.zph` test; p \< 0.05 is flagged.
- **Endpoint tokens and pooling classes.** Every result carries the
  resolved token and the pooling class of that token (`Exact-equivalent`
  when the token *is* the family’s canonical definition,
  `Clinically-related` when a different token is pooled into the family
  by the documented rule, `Not-poolable`, `Unknown`, or `Absent` when
  the cohort has no endpoint in that family at all - the absence of an
  endpoint, not a pooling verdict). Mixed-token pools warn.
  `cpas_meta(pooling = )` chooses between the documented default
  `"family"` (pool Exact-equivalent + Clinically-related, which is what
  makes a cross-token DFS pool of 98 cohorts possible) and the strict
  `"exact"` (Exact-equivalent only, which for DFS keeps the ~20 cohorts
  whose token is literally DFS); either way the classes actually pooled
  are recorded on the result.
- **Prediction intervals.** Three construction rules, selected by
  `cpas_meta(pi_method = )`: `"t"` (the default, the documented and
  cited t(k-2) construction), `"normal"` (the normal approximation) and
  `"HK"` (the Hartung-Knapp-based interval, on the HK-adjusted standard
  error with t(k-1)). The rule actually used is recorded in
  `$pooled$pi_method` and `$pooled$pi_rule`, and the other construction
  is returned beside it as `$pooled$pi_alt_lower/upper` so both can be
  reported. Only the interval changes: the pooled HR, its confidence
  interval and the heterogeneity statistics are identical whatever
  `pi_method` is.
- **Batch effects across cohorts.** Multi-cohort pooling works on
  **within-cohort standardised effect sizes** (a two-stage
  meta-analysis: each cohort is fitted and standardised inside itself
  before pooling), and pooled Kaplan-Meier splits each cohort at its
  **own percentile cut-point**. A constant cross-platform shift in
  expression therefore does not enter the pooled estimate.
  [`batch_diagnostics()`](https://wangjin93.github.io/CanPAS/reference/batch_diagnostics.md)
  reports the per-cohort distribution that both rules are computed from
  (n, median, IQR, min/max, the median-split and percentile cut-points
  used) and
  [`plot_batch_diagnostics()`](https://wangjin93.github.io/CanPAS/reference/plot_batch_diagnostics.md)
  draws the cohorts side by side on one axis.
- **Shared-patient overlap.**
  [`cohort_overlap()`](https://wangjin93.github.io/CanPAS/reference/cohort_overlap.md)
  returns the register of cohorts known to share patients;
  `cpas_meta(overlap = )` can `"warn"` (default: proceed, name the pairs
  and record them), `"refuse"` (stop with an actionable error) or
  `"dedupe"` (keep the larger cohort of every overlapping group,
  deterministically, and record what it dropped).
- **Competing risks.** Only usable when the source data distinguish the
  cause of death; the mirrored GEO tables code non-disease death as
  censored.
- **Not implemented.** Time-varying covariates, landmark/immortal-time
  correction, functional-form modelling, and
  inverse-probability-of-censoring time-dependent ROC. Single-marker
  pages apply no multiplicity control; gene panels are BH-adjusted. (The
  Hartung-Knapp-Sidik-Jonkman variance adjustment is now implemented as
  `cpas_meta(method = "HK")`.)

## The Shiny application

[`run_cpas_app()`](https://wangjin93.github.io/CanPAS/reference/run_cpas_app.md)
starts a ten-page app: Dashboard, Datasets (endpoint family first, then
cancer type, with the token each cohort contributes), Single dataset
analysis (KM, COX, COX by genes), Multi-datasets analysis (COX by
datasets, pooled KM, meta-analysis), Methods and Help. Selection runs
from the endpoint to the cohorts, so a pooling selection cannot lose a
cohort to an endpoint it does not have. The Methods page renders the
data sources, curation steps, harmonisation rules, endpoint-family
definitions, the overlap warning and the caveats.

## Data sources

| Source                        | Expression                                                           | Survival / clinical                                                                                                                    |
|-------------------------------|----------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| GEO (145 cohorts)             | CanPAS MySQL mirror over a public REST API                           | mirror table `<ACC>_surv`                                                                                                              |
| EMBL-EBI (14 cohorts)         | the deposit’s own processed matrix, or CEL files re-processed by RMA | SDRF annotation fields, resolved to patient level                                                                                      |
| CGGA (3 glioma cohorts)       | mirror                                                               | mirror table `CGGA_<ID>_surv`                                                                                                          |
| TCGA (33 projects)            | UCSC Xena, fetched per gene on demand                                | clinical/survival tables bundled with the package (`inst/extdata/tcga/*.rda`); `CPAS_DATA_ROOT` overrides them with a project checkout |
| cBioPortal-hosted (2 cohorts) | mirror (the study’s own RNA-seq matrix)                              | clinical patient files deposited with the study, loaded as local survival tables                                                       |

Cohort data remain the property of the original studies: cite the
GEO/CGGA/TCGA accessions listed in `dataset_info` alongside any result.
The two cBioPortal-hosted studies are counted in their own source
bucket, not as GEO.

## Caching of remote fetches

Queries to the mirror API and to UCSC Xena are cached as `.rds` files,
so repeated analyses of the same cohort, platform or gene do not
download again:

``` r
CanPAS:::.cpas_cache_dir()      # where entries live
CanPAS:::.cpas_cache_info()     # size and age of the cache
CanPAS:::.cpas_cache_clear()    # empty it (internal helpers)
options(CanPAS.cache_dir = "/scratch/cpas_cache",  # move it
        CanPAS.cache_ttl = 7 * 24 * 3600,          # 7 days
        CanPAS.cache = FALSE)                      # or switch it off
get_data("GSE31210", "surv_data", use_cache = FALSE)   # per call
```

Entries expire after `CanPAS.cache_ttl` seconds and are reused after
expiry only when the live request fails, so a flaky endpoint degrades to
the last good copy.
[`get_data()`](https://wangjin93.github.io/CanPAS/reference/get_data.md)
reports `metadata$from_cache` and `metadata$cache_stale`.

### Colour convention in Kaplan–Meier plots

High expression is red, low expression is green in
[`plot_km()`](https://wangjin93.github.io/CanPAS/reference/plot_km.md),
[`plot_cpas_km()`](https://wangjin93.github.io/CanPAS/reference/plot_cpas_km.md)
and the per-cohort grid. Override the session default or a single call:

``` r
options(CanPAS.km_palette = c("#1B7837", "#B2182B"))   # low, high
plot_km(d, type = "OS", marker = "GAPDH", palette = c("steelblue", "firebrick"))
plot_cpas_km_perdataset(km, ncol = 2, draw = FALSE)    # returns a patchwork object
```

## Validation

The estimation paths that are checked against outside references, and
the ones that are not, are listed in the accompanying paper. In short:

- [`COX_by_genes()`](https://wangjin93.github.io/CanPAS/reference/COX_by_genes.md)
  reproduces
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  exactly (wrapper fidelity, 133 estimable cohort–marker combinations).
- [`cpas_meta()`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
  / `meta_pool()` match `metafor::rma(method = "DL")` to 2.84e-14 over
  24 comparable pairs (machine precision under that convention).
- [`competing_risk_COX()`](https://wangjin93.github.io/CanPAS/reference/competing_risk_COX.md)
  matches [`cmprsk::crr()`](https://rdrr.io/pkg/cmprsk/man/crr.html)
  point estimates to 6.35e-07; the standard errors shipped in earlier
  builds were wrong (clustering per `finegray` row instead of per
  patient) and are corrected here — the CanPAS-to-`crr` mean ratio is
  0.9999756 after the fix.
- The pooled Kaplan–Meier log(−log S) delta-method standard error
  matches a hand-computed Greenwood standard error to 1.11e-16.
- `tests/testthat/` holds 119 `test_that` blocks and 505 assertions,
  covering every defect found in the pre-release audit.

## Curation pipeline

The scripted pipeline that built and maintains the curated mirror ships
**inside** the package, as scripts only, under `inst/pipeline/`. It is
the pipeline the paper refers to; it is not a separate download.

- Nine consolidated files, `01_ingest_parse.R` …
  `09_analysis_and_figures.R`, plus the entry point `run_pipeline.R`, a
  `README.md` and `.Renviron.example`. They consolidate the 67
  standalone step scripts the pipeline was originally run from.
- Every original script is embedded **byte-for-byte** inside a
  zero-argument runner (`run_<script>()`), so
  [`source()`](https://rdrr.io/r/base/source.html)-ing a file defines
  functions and has no side effects, and each step still runs in its own
  process. All **491 one-level function definitions** keep an identical
  deparsed body (0 altered; counting rule in `inst/pipeline/README.md`).
- **No data files are added to the package.** The pipeline reads and
  writes an external data root (`CPAS_DATA_ROOT`, by default
  `/home/Jingle/data/Project/CPAS` on the authors’ machine); the
  packaged TCGA tables under `inst/extdata/tcga/` are unrelated to it.
  On a machine without that tree the scripts still document exactly how
  each cohort was curated, but they cannot re-create the mirror without
  the raw archives.
- Locate and run them from R:

``` r
system.file("pipeline", package = "CanPAS")          # the installed script directory
source(system.file("pipeline", "05_clinical_and_scale.R", package = "CanPAS"))
run_07b_split_tnm()                                   # one step, as before
```

``` sh
Rscript "$(Rscript -e 'cat(system.file("pipeline","run_pipeline.R",package="CanPAS"))')" --list
```

## Reproducibility index

[`cpas_reproducibility_index()`](https://wangjin93.github.io/CanPAS/reference/cpas_reproducibility_index.md)
returns a tidy data frame indexing the reproducibility material: the
**67 curation steps** (canonical order, original script name, what the
step does, and the consolidated `inst/pipeline/` file that provides it),
the **exclusion register** (109 records, with its location), the repair
and defect records, the frozen-state record and the Additional-file
layout.
[`cpas_pipeline_steps()`](https://wangjin93.github.io/CanPAS/reference/cpas_pipeline_steps.md)
is the companion that returns just the 67 steps. Both read two small
CSVs shipped in `inst/reproducibility/` and resolve the packaged files
with [`system.file()`](https://rdrr.io/r/base/system.file.html), so they
work offline from an installed package; the registers themselves stay in
the project tree and are reported by count and location, never copied.

``` r
idx <- cpas_reproducibility_index()
table(idx$section)                                  # curation_step / register / ...
subset(idx, section == "register", select = c("item", "records", "location"))
```

## Endpoint adjudication

The endpoint annotation is our own reading of each deposit, so it is
auditable:
[`endpoint_adjudication()`](https://wangjin93.github.io/CanPAS/reference/endpoint_adjudication.md)
draws a stratified sample (by family and by source) of cohort-endpoint
records and returns a **blinded** rating sheet - record id, cohort,
token and evidence text only, with one empty column per judgement field
and rater - and
[`endpoint_agreement()`](https://wangjin93.github.io/CanPAS/reference/endpoint_agreement.md)
scores the filled sheet: per-field raw agreement, Cohen’s kappa with a
Landis-and-Koch label, the overall pooled agreement, the adjudication
rate and the records needing adjudication. Cohen’s kappa is implemented
in the package (base R only, no added dependency); with more than two
raters the per-field kappa is the mean of the pairwise kappas, which the
result states explicitly.

## Repository layout

    R/                  exported functions (48) and internal helpers
    man/                roxygen-generated help pages (66)
    inst/shiny/CanPAS/  the bundled Shiny application (apps/, www/, HELP.md)
    inst/pipeline/      the curation pipeline, shipped as consolidated scripts (no data)
    inst/reproducibility/  the shipped index behind cpas_reproducibility_index() (two CSVs)
    inst/extdata/       endpoint semantics, cohort overlap and the TCGA clinical tables
    data/               dataset_info.rda (the catalog) and ID_map.rda
    tests/testthat/     regression tests
    NEWS.md             change log for 1.0.0

The data root the pipeline operates on is distributed separately from
this package; the curated tables the package serves are already inside
it (`data/dataset_info.rda`) or are served by the public mirror.

## Citation and licence

A tool paper accompanies this release; until it has a DOI, cite the
package version:

> CanPAS: Cancer Prognosis Analysis Suite. R package version 1.0.0.
> <https://github.com/WangJin93/CanPAS>

Licence: **MIT** (see `LICENSE`), so the code can be reused and modified
freely, including in closed-source derivatives. CanPAS is developed by
**Jin Wang**, Soochow University (<Jinwang93@suda.edu.cn>). Note that
some optional features rely on GPL-licensed packages (`survminer`,
`ggtext`, `survivalROC`, `gridExtra`, `flextable`), which keep their own
licences when installed. Cohort data remain the property of the original
studies — cite the GEO/CGGA/TCGA accessions alongside any result.

## Reporting problems

Please open an issue at <https://github.com/WangJin93/CanPAS/issues>
with the output of
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), the call
you made, and — when a number looks wrong — the resolved endpoint token
and sample sizes printed in the result object.
