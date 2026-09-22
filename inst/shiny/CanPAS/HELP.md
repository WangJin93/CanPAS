# CanPAS Shiny App — Help

Launch:

```r
library(CanPAS)
run_cpas_app()            # opens in the default browser
run_cpas_app(port = 3838) # fixed port
```

## Navigation

Two analysis groups, each with three pages (no deeper nesting):

| Page | Purpose |
|---|---|
| **Dashboard** | Introduction and a snapshot of the catalog (cohorts, cancer types, endpoints). |
| **Datasets** | Browse cohorts by cancer type and endpoint family; select a row to make it the shared selection and jump to KM / COX analysis. |
| **Single dataset analysis · KM analysis** | Kaplan-Meier curves for a gene (optional REF_ID probe) or a weighted signature, with median / top percent / absolute threshold / auto cut points. |
| **Single dataset analysis · COX analysis** | Univariable Cox for the marker plus selected clinical covariates, then a multivariable model restricted to the covariates below the univariable p threshold; forest plots, tables, EPV and proportional-hazards diagnostics. Results sit in two tabs: **Forest plots & tables** (both models, marker comparison) and **Print results** (the same two models combined into one three-line table, downloadable as Word or CSV). These forest plots keep the **model order**: sorting them by HR would scatter the levels of a clinical covariate, so each covariate stays together under its own name. The per-gene and multi-dataset pages have one row per gene/cohort and **do** sort by HR. |
| **Single dataset analysis · COX by genes** | One univariable Cox model per gene in ONE dataset, with Benjamini-Hochberg FDR across the genes; forest plot above the table. |
| **Multi-datasets analysis · COX by datasets** | The same univariable Cox (one gene) in several datasets: per-dataset estimates with FDR and a forest plot (no pooling). |
| **Multi-datasets analysis · Pooled KM** | Integrates several datasets in two ways: IPD pooling of the patients (median split inside each dataset) and time-point pooling of S(t) via log(-log S); landmark table, per-dataset curves and per-dataset estimates. |
| **Multi-datasets analysis · Meta-analysis** | Two-stage meta-analysis for one gene/signature (per-SD standardisation, RE/FE pooling, I2, prediction interval, leave-one-out) or a per-gene panel mode (`cpas_meta_panel()`) that pools every gene of a list and adjusts the gene-level FDR. |
| **Help** | This document, rendered inside the app: navigation, page-by-page description, endpoint families, statistical caveats and runnable R examples built from the package's own functions and real cohorts. |
| **Methods** | Data-processing pipeline (produce → clean → mirror), standardisation rules, endpoint-family pooling principles, the family-to-token map and the statistical caveats. |

Every result page shows the **figure above the table**, and all numbers are
printed with **4 decimals** (p-values below 0.0001 are shown as `<0.0001`).

## Shared selection (synchronised across pages)

A single selection — **dataset + survival endpoint** — is shared by all pages:

- selecting a row in **Datasets** immediately updates the shared selection;
- the **KM analysis** and **COX analysis** pages each contain their own
  *Dataset* and *Endpoint* drop-downs, pre-filled from the shared selection;
- changing the dataset on any page updates the others (the endpoint list is
  rebuilt from the endpoints actually available in that cohort);
- every parameter carries the short label shown in the sidebar plus a hover
  explanation behind the small ⓘ next to it; the probe box additionally states
  what the current choice means ("Auto: all probes (6 on this platform) of
  GAPDH are collapsed per sample by the maximum");
- the marker definition (gene / REF_ID probe / signature + probe rule) set on
  the **KM analysis** page is inherited by **COX analysis**.

## Data sources

| Source | Notes | Network |
|---|---|---|
| GEO cohorts | served from the CanPAS MySQL mirror through the public API | required (HTTP 429 = rate limit, retry later) |
| CGGA cohorts (3 glioma cohorts) | served from the same mirror (`CGGA_<ID>_surv`) | required |
| TCGA projects | expression fetched per gene on demand from UCSC Xena; clinical/survival tables stored locally | required for expression |
| TCGA projects (`TCGA-XXXX`) | local clinical/survival tables + on-demand UCSC Xena expression | required for expression |
| Datasets table links | NCBI GEO for GEO cohorts, GDC Portal for TCGA projects | — |

## KM analysis

- **Marker**: a gene (with an auto-populated REF_ID / probe list, default
  "all probes → max") or a weighted signature such as `0.5*GAPDH + 0.5*TNS1`
  with probe collapsing by max / mean / median / min.
- **Cut point**: `50 (top 50% high)` = median split (default), `Auto` = the
  cut-point maximising the log-rank statistic, or `Custom` with a unit of either
  **Top percent (%)** — entering 25 makes the highest quarter of patients High, i.e.
  the threshold is the 75th percentile — or **Absolute threshold** (an expression
  value on the log2 scale). The unit and value boxes appear only for `Custom`.
- Outputs: Kaplan-Meier plot with risk table, downloadable **PDF** (custom
  width/height) and the analysis data table as **CSV**.

## COX analysis

- Marker is inherited from the KM page.
- Clinical covariates are multi-selected (continuous and categorical variables
  are detected automatically); "use all available clinical covariates" is one click.
- Outputs: univariable **and** multivariable forest plots, result tables,
  downloads, plus a marker estimate comparison table.

## Multi-datasets analysis

1. **Meta-analysis** — two-stage per-SD pooling (random or fixed effects) with
   per-cohort estimates, pooled HR / p / I², heterogeneity, leave-one-out
   sensitivity, forest plot and CSV/PDF downloads.
2. **Pooled KM** — IPD pooling, time-point meta pooling or both; pooled and
   stratified log-rank p, landmark survival table (1/3/5 years by default),
   per-cohort KM panels, PDF/CSV downloads.
3. **Gene panel** — univariable Cox for every gene × cohort combination, shown
   as a log2(HR) heat-map (asterisk = p < 0.05) plus a downloadable table.

Data is cached per session, so repeated analyses of the same cohort do not
re-query the API.

## Notes

- Survival time is expressed in years; status is 0 = censored / 1 = event.
- TCGA `DFI`/`PFI` endpoints are additional to the GEO set and appear in the
  endpoint list when available.
- Cross-source analyses (GEO / CGGA / TCGA) must use `cpas_meta()` (per-SD
  standardisation) or within-cohort median splits; raw expression values are not
  comparable across platforms.

## Survival endpoint families

Raw endpoint tokens are grouped into **five pooling families** (used by
meta-analysis and pooled KM) and **four browsing families** (used by the
Datasets filter):

| Pooling family | Member tokens | Typical label |
|---|---|---|
| `OS` | OS | Overall survival |
| `DSS` | DSS, CSS, BCSS | Cancer/disease-specific survival |
| `DFS` | DFS, RFS, EFS, DFI | Disease-/recurrence-free survival |
| `PFS` | PFS, PFI | Progression-free survival |
| `MFS` | MFS, DRFS | Metastasis-/distant-recurrence-free |

Browsing families: `OS`, `DSS`, `DFS` and `PFS` (broad: `EP_PFS` or `EP_MFS`).

Rules of use:

- Single-cohort pages (KM / COX) show the family with the token actually used,
  e.g. `DFS (RFS)` or `DFS (DFI) [derived]`; TCGA-derived DFI/PFI are marked
  `[derived]`.
- `DFS`-family members may be pooled together (event definition is comparable).
  For the progression/metastasis family, pool `PFS + PFI` and `MFS + DRFS`
  **separately** — the endpoint selector lists `PFS (progression)` and
  `MFS (metastasis)` as distinct pooling families.
- Meta-analysis and pooled KM report the endpoint used by each cohort
  (`per_dataset$endpoint`, `dataset_endpoints`), so mixed-member pooling stays
  traceable.
- The raw `SurvivalTypes` column of `dataset_info` is preserved unchanged; the
  family columns (`EndpointFamilies`, `EP_OS … EP_MFS`, `EndpointPrimary`,
  `EndpointDerived`) are added on top. Nothing in the database is modified.

### Family-first analysis

Every analysis entry point accepts an **endpoint family** (or a raw token) and
resolves it to the token actually available in each cohort:

| Function | Behaviour |
|---|---|
| `cohort_merged()`, `tcga_merged()` | family resolved per cohort; the token is stored in the `endpoint` attribute |
| `COX_by_genes()`, `COX_analysis()`, `plot_km()`, `plot_roc()` | `type` may be a family; resolved against the columns present in the data |
| `quick_km()`, `quick_roc()` | family resolved after data retrieval |
| `COX_by_datasets()` | per-dataset resolution; `endpoint` column added to the results |
| `cpas_meta()`, `cpas_km_pooled()` | family pooling; per-cohort endpoint reported (`per_dataset$endpoint`, `dataset_endpoints`) |
| `endpoint_family()`, `endpoint_resolve()`, `endpoint_options()` | mapping / resolution helpers |

The App follows the same rule: endpoint selectors store the **family** and show
the resolved token (`DFS (RFS)`, `DFS (DFI) [derived]`), and every status line
reports the endpoint actually used.

### Mirror coverage and catalog names

The Datasets page lists every cohort of the catalog. Two conventions matter when
a cohort is added:

- **Accession spelling** — catalog and database use an **underscore** before the
  platform (`GSE10885_GPL1390`); the local RDS files use a **dash**
  (`GSE10885-GPL1390.rds`). The API queries the database with the catalog string
  verbatim, so a dash-spelled accession returns an error.
- **Reported n** — KM / COX status lines report the rows actually entering the
  model. Rows without the selected endpoint's time/status, or without a marker
  value, are excluded and the number excluded is printed explicitly.

### Endpoint choice for a new cohort

Before annotating a family, check the censoring structure: a token whose
censored patients all have time 0 (as happens with RFS in `GSE40272`) cannot be
used as a time-to-event endpoint, and a token with (almost) no events (OS in
`GSE40272`: 1 event) cannot support KM/Cox either.

### Catalog conventions

- **Cancer type** labels are title-cased (`Lung Cancer`, `Multiple Myeloma`) and
  are the grouping key of the Multi-datasets page.
- **Source** is derived from the accession: `TCGA-*` → TCGA, `CGGA*` → CGGA,
  everything else → GEO. The Datasets page shows the source and links to the
  matching portal (NCBI GEO / CGGA / GDC).
- **Endpoint columns**: `SurvivalTypes` holds the raw tokens, `EndpointFamilies`
  the families a cohort can be pooled under, and `EP_OS`/`EP_DSS`/`EP_DFS`/`EP_PFS`/
  `EP_MFS` the token the cohort contributes to each family (NA when absent).

### Endpoint families in the interface

The Datasets page filters by family, and every available family is listed
separately with its cohort count:

| Filter | Meaning | Cohorts |
|---|---|---|
| Any | no family filter | 152 |
| OS | overall survival | 98 |
| DSS | disease-specific survival (DSS / CSS / BCSS) | 22 |
| DFS | disease-free survival (DFS / RFS / EFS / DFI) | 80 |
| PFS | progression-free survival (PFS / PFI) | 23 |
| MFS | metastasis-free survival (MFS / DRFS) | 13 |
| PFS or MFS | the broader "no progression / no metastasis" view | 36 |

The **Endpoint families** column shows, per cohort, the family and the token it
contributes — `DFS (RFS), MFS`, `MFS (DRFS)`, `DFS (DFI) [derived]` — i.e. the
same labels used by the endpoint selectors on the analysis pages. The KM, COX
and Multi-datasets pages always run an analysis **per family** and report the
token actually used.

The Dashboard shows the family inventory with cohort counts (and the number of
families) instead of the fine-grained token list.

### Statistical caveats (read before writing a result into a paper)

| Topic | What the app does | What you should do |
|---|---|---|
| **Cut point** | Median (default), custom, or `Auto` = best of 13 percentile cut points by log-rank χ². The naive p-value of an auto cut is optimistically biased: in a simulation where the marker is independent of survival, p < 0.05 occurs in 23.1–26.7% of datasets (auto) versus 4.45–5.90% (median), 2000 replicates at each of n = 50, 100, 200, 400; with the `maxstat` search adjustment it drops to 2.0–3.5%, i.e. the correction is conservative. | Prefer the median/custom rule; when you use `Auto`, quote the search-adjusted p (shown when `maxstat` is installed) or report it as exploratory. |
| **Small cohorts** | Status line shows group sizes, events and median follow-up, and notes < 10 events; multivariable Cox warns below 10 events per variable. | For a 1–3-event cohort, report the curve descriptively only. |
| **Proportional hazards** | `cox.zph` GLOBAL test is included with the multivariable Cox output and flagged when p < 0.05. | Report the test; consider time-stratified or time-varying models if it fails. |
| **Endpoint tokens** | Family-level analysis; the token actually used is always shown, and mixed tokens inside a family raise a warning. | State the per-cohort tokens in the paper (e.g. RFS in GEO, DFI in TCGA). |
| **Sample sizes** | Reported n = rows entering the model; exclusions for missing endpoint/marker are counted. | Use the reported n, not the catalog N. |
| **Not implemented** | Competing risks (Fine-Gray), time-varying covariates, immortal-time/landmark correction, multiple-testing correction for gene panels, prediction intervals. | Handle these outside the app if the question requires them. |

## R examples (real functions, real cohorts)

Everything below runs against the CanPAS mirror with the package's own functions
— no simulated data. These are the same pipelines as in the R help
(`?cohort_merged`, `?COX_analysis`, `?cpas_meta`, `?competing_risk_COX`, ...).
GEO cohorts come from the mirrored database, TCGA projects from UCSC Xena, so
all of them need network access.

### 1. Catalog and the prerequisite chain

```r
library(CanPAS)

## the real catalog: cohorts, cancer types, usable n, endpoint families
head(dataset_info[, c("Accession", "Type", "N", "EndpointFamilies")])

## one cohort, one call: endpoint token resolved, expression fetched
## (probe ids collapsed), survival merged into one analysis-ready table
d <- cohort_merged("GSE14814", c("GAPDH", "ACTB", "TP53"), type = "OS")
colnames(d)                # ID, OS_status, OS_time, GAPDH, ACTB, TP53
attr(d, "family")          # "OS"; attr(d, "endpoint") gives the token used

## the two steps explicitly, when you need the pieces
expr <- get_expr_data("GSE14814", c("GAPDH", "ACTB", "TP53"), process_duplicates = "max")
m    <- merge_surv_expr("GSE14814", expr)
d2   <- m$merged_data                   # what cohort_merged() returns
cl   <- m$raw_surv_data                 # + the cohort's clinical columns

## clinical covariates, for a cohort that has them
dc <- cohort_merged("GSE13507", c("GAPDH", "ACTB"), type = "OS", clin = TRUE)   # covariates
colnames(dc)               # ..., age, sex, T, N, M, grade, histology
```

### 2. Single-dataset analysis (what the KM / COX pages do)

```r
## KM with median split, log-rank p and risk table
plot_km(d, type = "OS", marker = "GAPDH", pval = TRUE)

## time-dependent ROC
plot_roc(d, type = "OS", marker = "GAPDH", predict.time = 3)$auc

## univariable Cox for the marker plus real covariates
r <- COX_analysis(dc, type = "OS", cont_Variates = c("GAPDH", "age"),
                  cate_Variates = c("grade", "N"), method = "uni")
head(r$results_table)

## multivariable Cox: covariates that cannot be co-estimated are located,
## removed and reported (never silently, never with meaningless HRs)
r2 <- COX_analysis(dc, type = "OS", cont_Variates = c("GAPDH", "age"),
                   cate_Variates = c("grade", "N"), method = "multi")
r2$metadata$dropped_covariates   # variable, detail, reason
r2$metadata$reduced              # TRUE when only one covariate survived
forest_plot(r2$results_table)

## one univariable model per gene, FDR across genes
COX_by_genes(d, type = "OS", genes = c("GAPDH", "ACTB", "TP53"))$results_table
```

### 3. Multi-dataset analysis

```r
## the same gene in several cohorts (no pooling)
COX_by_datasets(c("GSE14814", "GSE31210"), gene = "GAPDH", type = "OS")$results_table

## pooled KM: IPD pooling and time-point pooling of S(t)
merged <- list(GSE31210 = cohort_merged("GSE31210", "GAPDH", type = "RFS"),
               GSE37745 = cohort_merged("GSE37745", "GAPDH", type = "RFS"))
km <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS", method = "both")
km$dataset_endpoints; km$logrank_p
plot_cpas_km(km)
plot_cpas_km_perdataset(km, ncol = 2)

## two-stage meta-analysis of one gene, with leave-one-out sensitivity
mm <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
mm$pooled
plot_meta_forest(mm)
loo_meta(mm)

## per-gene meta panel (same as the Meta-analysis page, panel mode)
cpas_meta_panel(c("GSE31210", "GSE37745"), genes = c("TP53", "GAPDH"), type = "DFS")$table
```

### 4. Signatures and TCGA projects

```r
## weighted signature, scored inside a real cohort
sig <- get_signature_value("0.5*GAPDH + 0.5*ACTB", "GSE14814")
sig$metadata$genes_used; sig$signature_info$weights
head(sig$results_table)          # ID + one score per patient

## TCGA goes through the same reader and keeps its clinical columns
tcga <- cohort_merged("TCGA-LUAD", c("TP53", "GAPDH"), type = "OS")
colnames(tcga)
COX_analysis(tcga, type = "OS", cont_Variates = c("TP53", "age"),
             cate_Variates = c("sex", "stage"), method = "uni")$results_table

## the TCGA pieces individually
tcga_gene_expr_df("LUAD", "TP53")        # Xena, per gene, long table
tcga_surv_table("LUAD")                  # local clinical table, all endpoints
tcga_merged("LUAD", c("TP53", "GAPDH"), type = "OS")
```

### 5. Competing risks (package functions, not wired into the app)

```r
## GSE14814 carries OS and DSS, so both causes are real:
## 1 = cancer death (DSS event), 2 = death from other causes
d  <- cohort_merged("GSE14814", "GAPDH", type = "OS", clin = TRUE)
cr <- data.frame(time   = d$OS_time,
                 status = ifelse(d$OS_status == 0, 0,
                                 ifelse(d$DSS_status == 1, 1, 2)),
                 age    = d$age,          # numeric -> one coefficient
                 sex    = factor(d$sex))  # factor  -> one coefficient per level

cif_fit(cr, time = "time", status = "status", times = c(1, 3, 5))$table
plot_cif(cif_fit(cr, time = "time", status = "status"))

## stratified by a real categorical column
ci_sex <- cif_fit(cr, time = "time", status = "status", group = "sex")
ci_sex$group_levels; ci_sex$group_sizes

## 'group' must be categorical: passing the continuous age directly is refused
## (it would create one stratum per distinct value); group_cut asks for a
## documented split and returns the cut point used
ci_age <- cif_fit(cr, time = "time", status = "status",
                  group = "age", group_cut = "median")
ci_age$cut_points; ci_age$group_levels
cc <- competing_risk_COX(cr, time = "time", status = "status",
                         covariates = c("age", "sex"))
cc$cause_specific; cc$subdistribution      # cause-specific HR / Fine-Gray sHR
cc$covariate_types; cc$diagnostics         # typing used, and estimability per model
```

### Choosing an endpoint first

Every analysis starts from the **endpoint family** (OS, DSS, DFS, PFS, MFS), then
from the cohorts that carry it:

1. **Datasets** page: pick the endpoint family and the table lists only the cohorts
   that have it, with a column giving the token each one contributes there
   (`RFS`, `DFI`, `PFI`, ...). Narrow by cancer type afterwards if you want.
2. **COX by datasets / Pooled KM / Meta-analysis**: Cancer type -> **Endpoint** ->
   Datasets. The dataset list is restricted to cohorts that carry the chosen family,
   so a pooling run never silently loses a cohort because it lacks that endpoint.
   If a cohort still arrives from the shared selection without the family, the page
   says so instead of dropping it silently.
3. **Single-dataset pages** (KM, COX, COX by genes): keep the reverse order —
   pick the cohort, then the endpoint, because there you are looking at what that
   one cohort offers (`DFS (RFS)`, `PFS (PFI) [derived]`, ...).

The same rule decides both the filter and the analysis: a cohort counts as carrying a
family when the catalog says so, or when its raw endpoint token maps to that family.

### Cut point on the pooled KM page

Two rules, each applied inside every cohort separately, so the groups stay comparable
across platforms:

* **50 (top 50% high)** — patients above their own cohort's median are High. This is
  the default.
* **Top percent (%)** — the highest *x*% of that cohort by expression are High, the
  rest Low: entering **25** makes the top quarter of the cohort High. The threshold is
  that cohort's (100 − *x*)th percentile, so **25** uses the cohort's 75th percentile.
  Where several patients tie exactly at the threshold slightly more than *x*% can end
  up in High; the actual threshold and group sizes are printed with the result.

If a rule leaves one side empty in a cohort — which can still happen with a top-percent
split when the marker is constant in that cohort — the cohort is reported in the result
and left out of the pool rather than quietly changing it. Cohorts dropped for
insufficient data (no resolved endpoint, missing columns, fewer than 10 complete
time/status/marker rows, or a constant marker) are listed in the status line under
"Skipped", so the pool never shrinks without a printed reason. The `cpas_km_pooled()`
script interface additionally accepts `cut = "custom"` with one absolute log2 threshold
applied in every cohort; it is not offered on this page.

A per-cohort search for the best log-rank cut point is **not** offered on this page:
that search inflates the p-value (see the cut-point simulation in the paper), and
running it inside every cohort before pooling compounds the inflation. Use the
single-dataset KM page for the search-adjusted cut point.

### What "both" shows on the pooled KM page

`Pooling = Both (IPD + time-point meta)` runs the two routes and draws them **in one
figure**:

* **(A) IPD** — all patients of the selected cohorts stacked into one table and one
  Kaplan-Meier curve, with the number-at-risk table underneath and both the
  unstratified and the cohort-stratified log-rank p-values in the status line.
* **(B) Time-point meta** — inside each cohort, S(t) is estimated by Kaplan-Meier at
  the 1 / 3 / 5-year landmarks and pooled across cohorts by inverse variance (RE or
  FE); the figure shows the pooled curve with its confidence band and marks each
  landmark estimate with its confidence interval, and the table below lists S(t),
  the number of contributing cohorts `k`, I² and the heterogeneity p per landmark.

Compare the two: if the IPD curve and the landmark points disagree, the cohorts
differ in baseline risk, which the meta route quantifies as I². `IPD` alone and
`Time-point meta` alone draw only their own panel (and then the meta-pooling
estimator / landmark controls are hidden, because nothing would read them).
