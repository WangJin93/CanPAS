# CanPAS — Cancer Prognosis Analysis Suite

`CanPAS` is an R package for survival analysis of public cancer cohorts. It serves
expression matrices and survival tables from a curated GEO/CGGA mirror over a public
REST API and from TCGA (UCSC Xena expression plus local clinical tables), and exposes
one uniform **merged-data schema** to Kaplan–Meier, Cox, time-dependent ROC,
cross-cohort meta-analysis, pooled Kaplan–Meier and competing-risks analyses. A bundled
Shiny application drives the same exported functions.

What it adds to the usual single-cohort workflow:

* **Endpoint families with per-cohort token resolution.** Cohorts are requested by
  family (OS, DSS, DFS, PFS, MFS); the token a cohort actually reports (`RFS`, `DFS`,
  `DFI`, …) is resolved per cohort and returned with every result, and pooling never
  crosses families.
* **Five statistical safeguards on by default.** Cut-point rule and threshold,
  events-per-variable warning, proportional-hazards test, counted excluded rows, and
  the resolved endpoint token are printed with every result rather than on request.
* **Numerical verification against reference implementations** (`survival`, `metafor`,
  `cmprsk`), including a competing-risks variance defect that the comparison exposed
  and this release corrects.
* **A catalog, not just a downloader**: 177 catalogued cohorts in one schema, with
  sample sizes defined as analysable patients and with patient-overlap groups recorded.

> Naming note: the package was originally created under the name "Cancer Patient
> Aftercare System". The content is prognostic survival analysis of public cohorts, not
> aftercare intervention, so it is now the "Cancer Prognosis Analysis Suite" (the
> acronym is unchanged).

**Status**: version 1.0.0 (first release) · licence GPL-3 · `R CMD check`
0 errors / 0 warnings / 1 note (host library notice) · a tool paper is in preparation.

## Installation

```r
# from GitHub
remotes::install_github("WangJin93/CanPAS")

# or from the release tarball
install.packages("CanPAS_1.0.0.tar.gz", repos = NULL, type = "source")
# inside a source checkout: R CMD INSTALL .
```

`Imports` are all on CRAN. Optional packages enable extra features:

| Package | Enables |
|---|---|
| `shiny`, `DT`, `bs4Dash`, `shinyWidgets`, `shinycssloaders` | the bundled application |
| `UCSCXenaShiny`, `UCSCXenaTools` | TCGA expression through UCSC Xena |
| `flextable` | three-line summary table in the Cox page and its Word export (without it the app shows the same table and offers CSV) |
| `maxstat` | search-adjusted approximate p-value for the cut-point search |
| `patchwork` | composing both pooling routes into one figure |
| `cmprsk` | independent Fine–Gray reference in the test suite |

```r
install.packages(c("dplyr", "ggplot2", "ggtext", "gridExtra", "jsonlite",
                   "survival", "survminer", "survivalROC", "flextable"))
install.packages(c("shiny", "DT", "bs4Dash", "shinyWidgets", "shinycssloaders",
                   "flextable", "maxstat", "patchwork", "cmprsk"))
```

## Configuration

| Variable / option | Used for |
|---|---|
| `CPAS_DATA_ROOT` | local TCGA clinical/survival tables (`<root>/data/tcga/*.rda`) and cohort clinical covariates in the app. **No default path is assumed**: unset, the GEO/CGGA mirror and UCSC Xena still work, but the TCGA helpers stop with instructions. |
| `CPAS_DB_PASSWORD` | the curation pipeline that maintains the mirror (never hard-coded in the package) |
| `options(CanPAS.cache_dir=)`, `CANPAS_CACHE_DIR` | where downloaded answers are cached |
| `options(CanPAS.cache=)`, `CANPAS_CACHE` | switch the cache off |
| `options(CanPAS.cache_ttl=)` | cache lifetime in seconds (30 days by default) |

```r
Sys.setenv(CPAS_DATA_ROOT = "/path/to/CanPAS-data-root")   # for TCGA and clinical covariates
```

GEO and CGGA data are read from the public CanPAS mirror API
(`https://www.jingege.wang/bioinformatics/CPAS/api.php`); it is a live service and
rate-limits bulk sweeps (HTTP 429 → the client backs off and retries). `get_data()`
accepts a different `base_url` if you host your own mirror.

## Quick start

```r
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
res <- cpas_meta(datasets = c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
res$per_dataset    # per-cohort HR, n, events and the endpoint token used
res$pooled         # pooled HR, Q, I2, tau2 and the prediction interval
loo_meta(res)      # leave-one-out sensitivity

## --- pooled Kaplan-Meier with the two within-cohort cut rules ---------------
merged <- lapply(c("GSE31210", "GSE37745"), function(a) cohort_merged(a, "GAPDH", type = "RFS"))
names(merged) <- c("GSE31210", "GSE37745")
km <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS", method = "both")
km$cutpoint                                  # the rule and the per-cohort thresholds
km2 <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS",
                      cut = "top_pct", top_pct = 25)   # highest quarter of every cohort
km2$empty_cohorts; km2$skipped_cohorts        # cohorts left out, with the reason

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
**177 cohorts — 141 GEO, 31 TCGA projects, 3 CGGA and 2 cBioPortal-hosted studies — across
29 cancer types**, and every row carries a resolved endpoint; together they contribute
**37,052 analysable samples** (median 163 per cohort, range 36–1,210). The two
cBioPortal-hosted cohorts (A5-PCPG, Pheochromocytoma; IMmotion150, Kidney Cancer) are
studies whose clinical and expression files are deposited together, not GEO series, so they
are counted in their own bucket and the GEO count is 141 rather than 143. Sample size means
analysable patients:
expression data plus a non-missing time and status for the cohort's primary endpoint.
Patient-overlap groups are recorded (29 pairs in 13 groups recomputed against this
catalog by `pipeline/R/26_cohort_overlap.R`, with `data(dataset_info)$CohortGroup` and
`$Note` carrying the result). One further title match — the GSE25066–GSE32918 pair — is
a sample-title collision rather than shared patients, since GSE32918's titles are panel
replicate codes for 172 patients and its genuine duplicate deposit (GSE69051) is not
catalogued, so it is recorded as a note instead of a group
(`pipeline/ref/cohort_overlap_exclude.csv`). The multi-dataset pages warn when a
selection contains two members of one group.

## Endpoint families

A family groups tokens that answer the same clinical question, so cohorts reporting
different token names can be pooled. The token each cohort contributes is resolved per
cohort and always reported.

| Family | Tokens pooled | Cohorts |
|---|---|---|
| OS | OS | 125 |
| DSS | DSS, CSS, BCSS | 39 |
| DFS | DFS, RFS, EFS, DFI | 93 |
| PFS | PFS, PFI | 44 |
| MFS | MFS, DRFS | 17 |

`DFI` and `PFI` (TCGA) are derived from the original endpoint fields and are flagged as
derived. Pooling inside a family and never across families is enforced by the package;
a mixed-token pool warns, naming each cohort and its token. A separate browsing
vocabulary widens the progression family to metastasis endpoints, so the Datasets page
can show a PFS cohort count (61) larger than the pooling count (44).

## Statistical safeguards (read before quoting a result)

* **Cut point.** The default is a within-cohort median split; a custom rule is a top
  percentage (entering 25 takes the highest quarter of patients) or an absolute log2
  threshold. The app's `Auto` mode searches 13 percentile cut points and keeps the
  largest log-rank statistic, so its naive p-value is optimistically biased: with a
  marker generated independently of survival, p < 0.05 in 23.1–26.7% of 2,000
  replicates per sample size versus 4.45–5.90% for the median split, and 2.04–3.53%
  with the `maxstat` search-adjusted approximation. The app prints the number of
  candidates searched and the adjusted p-value when `maxstat` is installed. Pooled
  analyses offer only the two percentile rules; a per-cohort cut-point search is
  deliberately not offered.
* **Sample sizes.** Reported n is the number of rows entering the model; rows excluded
  for a missing endpoint time, status or marker are counted and printed.
* **Events per variable.** Multivariable Cox models warn below 10 events per variable.
* **Proportional hazards.** `COX_analysis()` returns the global `cox.zph` test; p < 0.05
  is flagged.
* **Endpoint tokens.** Every result carries the resolved token; mixed-token pools warn.
* **Competing risks.** Only usable when the source data distinguish the cause of death;
  the mirrored GEO tables code non-disease death as censored.
* **Not implemented.** Hartung–Knapp variance adjustment, time-varying covariates,
  landmark/immortal-time correction, functional-form modelling, and inverse-
  probability-of-censoring time-dependent ROC. Single-marker pages apply no
  multiplicity control; gene panels are BH-adjusted.

## The Shiny application

`run_cpas_app()` starts a ten-page app: Dashboard, Datasets (endpoint family first, then
cancer type, with the token each cohort contributes), Single dataset analysis (KM, COX,
COX by genes), Multi-datasets analysis (COX by datasets, pooled KM, meta-analysis),
Methods and Help. Selection runs from the endpoint to the cohorts, so a pooling
selection cannot lose a cohort to an endpoint it does not have. The Methods page renders
the data sources, curation steps, harmonisation rules, endpoint-family definitions, the
overlap warning and the caveats.

## Data sources

| Source | Expression | Survival / clinical |
|---|---|---|
| GEO (140 cohorts) | CanPAS MySQL mirror over a public REST API | mirror table `<ACC>_surv` |
| CGGA (3 glioma cohorts) | mirror | mirror table `CGGA_<ID>_surv` |
| TCGA (31 projects) | UCSC Xena, fetched per gene on demand | local `<CPAS_DATA_ROOT>/data/tcga/*.rda` |

Cohort data remain the property of the original studies: cite the GEO/CGGA/TCGA
accessions listed in `dataset_info` alongside any result.

## Caching of remote fetches

Queries to the mirror API and to UCSC Xena are cached as `.rds` files, so repeated
analyses of the same cohort, platform or gene do not download again:

```r
CanPAS:::.cpas_cache_dir()      # where entries live
CanPAS:::.cpas_cache_info()     # size and age of the cache
CanPAS:::.cpas_cache_clear()    # empty it (internal helpers)
options(CanPAS.cache_dir = "/scratch/cpas_cache",  # move it
        CanPAS.cache_ttl = 7 * 24 * 3600,          # 7 days
        CanPAS.cache = FALSE)                      # or switch it off
get_data("GSE31210", "surv_data", use_cache = FALSE)   # per call
```

Entries expire after `CanPAS.cache_ttl` seconds and are reused after expiry only when
the live request fails, so a flaky endpoint degrades to the last good copy.
`get_data()` reports `metadata$from_cache` and `metadata$cache_stale`.

### Colour convention in Kaplan–Meier plots

High expression is red, low expression is green in `plot_km()`, `plot_cpas_km()` and the
per-cohort grid. Override the session default or a single call:

```r
options(CanPAS.km_palette = c("#1B7837", "#B2182B"))   # low, high
plot_km(d, type = "OS", marker = "GAPDH", palette = c("steelblue", "firebrick"))
plot_cpas_km_perdataset(km, ncol = 2, draw = FALSE)    # returns a patchwork object
```

## Validation

The estimation paths that are checked against outside references, and the ones that are
not, are listed in the accompanying paper. In short:

* `COX_by_genes()` reproduces `survival::coxph()` exactly (wrapper fidelity, 133
  estimable cohort–marker combinations).
* `cpas_meta()` / `meta_pool()` match `metafor::rma(method = "DL")` to 2.84e-14 over 24
  comparable pairs (machine precision under that convention).
* `competing_risk_COX()` matches `cmprsk::crr()` point estimates to 6.35e-07; the
  standard errors shipped in earlier builds were wrong (clustering per `finegray` row
  instead of per patient) and are corrected here — the CanPAS-to-`crr` mean ratio is
  0.9999756 after the fix.
* The pooled Kaplan–Meier log(−log S) delta-method standard error matches a
  hand-computed Greenwood standard error to 1.11e-16.
* `tests/testthat/` holds 116 `test_that` blocks and 492 assertions, covering every
  defect found in the pre-release audit.

## Repository layout

```
R/                  exported functions (37) and internal helpers
man/                roxygen-generated help pages (50)
inst/shiny/CanPAS/  the bundled Shiny application (apps/, www/, HELP.md)
data/               dataset_info.rda (the catalog) and ID_map.rda
tests/testthat/     regression tests
NEWS.md             change log for 1.0.0
```

The curation pipeline that builds and maintains the mirror, and the data root itself,
are distributed separately from this package.

## Citation and licence

A tool paper accompanies this release; until it has a DOI, cite the package version:

> CanPAS: Cancer Prognosis Analysis Suite. R package version 1.0.0.
> https://github.com/WangJin93/CanPAS

Licence: **MIT** (see `LICENSE`), so the code can be reused and modified freely,
including in closed-source derivatives. CanPAS is developed by **Jin Wang**, Soochow
University (<Jinwang93@suda.edu.cn>). Note that some optional features rely on
GPL-licensed packages (`survminer`, `ggtext`, `survivalROC`, `gridExtra`, `flextable`),
which keep their own licences when installed. Cohort data remain the property of the
original studies — cite the GEO/CGGA/TCGA accessions alongside any result.

## Reporting problems

Please open an issue at <https://github.com/WangJin93/CanPAS/issues> with the output of
`sessionInfo()`, the call you made, and — when a number looks wrong — the resolved
endpoint token and sample sizes printed in the result object.
