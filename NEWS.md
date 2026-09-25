# CanPAS 1.0.0

First public release. CanPAS is a curated cross-archive cancer prognosis resource
(GEO mirror, CGGA, TCGA), a scripted curation pipeline and an R package with a
bundled Shiny application; this section documents the state of that first release.

## Catalog growth to 196 cohorts: GSE1379 restored, TCGA-CHOL and TCGA-DLBC added (2026-09-25)

The catalog now holds **196 cohorts (38,953 analysable samples)** across 29 cancer types:
GEO 144, EMBL-EBI 14, TCGA 33, CGGA 3 and cBioPortal-hosted 2. Family coverage is now
OS 141, DSS 41, DFS 98, PFS 49 and MFS 17; the cohort-by-family audit table is
196 x 5 = 980 cells. No existing cohort, no cohort data file and no mirror table was
modified.

* **GSE1379 restored as a breast-cancer DFS cohort.** 60 patients and 28 DFS events
  (32 censored) on `GPL1223`, whole-tissue sections. The series has **no**
  `!Sample_characteristics_ch1` rows at all: time and status live in the series-matrix
  `!Sample_description` free text (`Clinical information: ...;DFS=<months>;Status=recur|non-recur`),
  and the survival table was parsed from there
  (`pipeline/R/39_build_gse1379_surv.R`, platform map `42_upload_gpl1223.R`). The
  earlier record that the series had "no clinical annotation" was wrong and has been
  corrected in the exclusion register.
* **The same-patient set is now recorded.** GSE1378 is the microdissected-cell version of
  the **same 60 patients** (identical 60/60 case ids). It stays uncatalogued, and the pair
  is registered in `pipeline/ref/cohort_overlap_seed.csv` and in the catalog as
  `CohortGroup = GSE1378(+1)` with an explanatory `Note`, so pooling counts those patients
  once. `pipeline/R/26_cohort_overlap.R` now reports **30 pairs in 14 groups (33 cohorts)**
  where it previously reported 29 in 13 (32); the app's Methods page and `HELP.md` were
  updated to match.
* **TCGA-CHOL (45 patients, 23 OS events) and TCGA-DLBC (47 patients, 9 OS events) added.**
  Both are Xena-backed (`expr_in_mirror = FALSE`): expression is fetched per gene on demand
  and the clinical/survival tables are local, so neither adds a mirror table
  (`pipeline/R/40_add_tcga_chol_dlbc.R`, `41_register_final3_rows.R`). TCGA-CHOL is
  registered under the existing **Liver Cancer** type, the same vocabulary already used for
  the earlier ICC cohort E-MTAB-6389, so no new cancer type is introduced; TCGA-DLBC joins
  the existing **Lymphoma** type. Each carries four endpoints (CHOL: OS 45/23, DSS 43/20,
  DFI 32/12, PFI 45/24; DLBC: OS 47/9, DSS 47/4, DFI 27/4, PFI 47/12). **TCGA-DLBC is
  thin** — only 9 OS events — and its OS time distribution contains a `time = 0` row, which
  the methods pages discuss alongside other zero-follow-up rows.
* **`16_verify_catalog_mirror.R` gained an ERROR check.** A row that is endpoint-annotated,
  or has a survival table, while every `EP_*` column is `NA` silently disappears from every
  family and from family-paired pooling. That condition is now an **ERROR** (the earlier
  `SurvivalTypes`-based "rows without endpoint annotation" count cannot see it). At the
  196-row state the script reports **0 errors / 0 warnings / 0 info**, 196/196 rows fully
  usable.
* **Additional file 2 refreshed.** The delivered exclusion list was a 2026-09-15 snapshot of
  35 records while the authoritative register `pipeline/out/excluded.csv` had grown to
  **109 records** over 106 accessions (74 added since that snapshot, four reasons corrected
  this round: GSE1378, GSE1379, GSE21501, GSE35629). The file is now a byte-identical copy of
  the register, the three superseded-append pairs (`GSE33630`, `GSE60542`, `GSE205209`) are
  documented, and the paper's Table 1 Block C count and its cross-reference were corrected
  (109, pointing at Additional file 2 rather than Additional file 3).
* **App and Help text refreshed with the new counts**, the new overlap pair and the new
  curation steps 39-42 (54 numbered steps, 63 R files); Figure 2 and Figure 9 were re-rendered
  from the 196-row build and the Word/PDF export re-run.

## Catalog growth to 193 cohorts: fourteen EMBL-EBI cohorts added (2026-09-25)

The catalog now holds **193 cohorts (38,801 analysable samples)** in **five** source
buckets instead of four: GEO 143, **EMBL-EBI 14**, TCGA 31, CGGA 3 and cBioPortal-hosted 2.
No existing cohort, no cohort data file and no mirror table was modified.

* **A new retrieval route: EMBL-EBI ArrayExpress/BioStudies.** Fourteen deposits
  (`E-MTAB-*`, `E-TABM-*`, `E-MEXP-*`) are retrieved from their own SDRF annotation and
  their deposited processed matrix — or, where none was deposited, from the raw CEL files
  re-processed here by RMA — instead of through the GEO mirror. They are their own source
  bucket in `dataset_info`, in the app and in the paper; folding them into GEO would report
  157 GEO cohorts rather than 143. Only one needed a platform map built from the GEO
  platform SOFT because no Bioconductor annotation package exists for it (`GPL16686`,
  `GPL17585`); measured probe coverage is recorded per cohort in the catalog `Note`.
* **Patient-level registration for row-level deposits.** Several deposits describe one
  patient in more than one row (a tumour and a normal column, a second non-expression
  assay, or `GeoMx` ROIs). These are collapsed to one row per patient before `N` is
  computed, and each collapse is recorded in the `Note` field together with the row-level
  figure it replaces.
* **Two small cohorts admitted under the relaxed gate.** `E-MTAB-1719` (mesothelioma,
  N = 34) and `E-MEXP-2780` (pancreatic, N = 30) were admitted under the gate relaxed from
  > 50 to >= 30 patients and are flagged `small cohort` in `Note`, as the earlier relaxed-gate
  additions are. The catalog now holds seven cohorts at 30–40 analyzable patients.
* **`16_verify_catalog_mirror.R` was hardened.** An earlier EMBL-EBI batch registered eleven
  rows whose `SurvivalTypes` was populated but whose `EndpointFamilies`, `EndpointPrimary` and
  every `EP_*` column were `NA` — invisible to the "rows without endpoint annotation" count
  (which reads `SurvivalTypes`), so those cohorts were silently absent from every family and
  from family-paired pooling. The check now raises an **ERROR** for "endpoint columns all NA
  while the row is endpoint-annotated or has a survival table", and it is evaluated locally
  as well as against the mirror. At the 193-row state the script reports
  **0 errors / 0 warnings / 0 info** with 193/193 rows fully usable.
* **Six-row arithmetic audit resolved with no row changed.** A sample audit flagged six rows
  where the catalog's arithmetic appeared to disagree with the local tables (`TCGA-MESO`,
  `TCGA-UCEC` `n_surv`; `CGGA_693`, `CGGA_301`, `CGGA_325`, `GSE108474` `n_events`). Re-reading
  the delivered artefacts showed every catalog value already correct under one convention, now
  written up in `README.md` § *Column semantics*: `n_surv` is the row count of the delivered
  survival table, and `n_events` counts the primary endpoint's events **among the `N`
  analysable samples**, not among all rows of the survival table. The audit's alternative
  figures were a different quantity in each case (the raw TCGA project sample count before the
  build drops endpoint-uninformative samples; the all-rows event count). Verified against every
  non-TCGA row that has both a delivered expression and a delivered survival table carrying its
  primary endpoint (112 rows): the convention reproduces `N` and `n_events` for 110, the two
  exceptions being documented patient/clinical-level audit cohorts (`GSE31312`, `GSE325123`).
  No catalog cell, and therefore no catalog hash, changed.
* **Documentation and the app follow the five buckets.** `README.md`, `HELP.md`, the packaged
  `dataset_info` help page and the Shiny application's Datasets page, Methods page (source
  table and pipeline step list) and source links all distinguish the EMBL-EBI bucket and
  link it to ArrayExpress/BioStudies. The catalog CSV and the packaged `dataset_info.rda` are
  cell-identical (0 of 5,018 cells differ).

## Fixes from the full package audit of 2026-09-24 (two blocking defects)

A full package audit against the 179-row catalog found two blocking defects; both are fixed
here and re-verified. No catalog row, no cohort data and no mirror table was changed.

* **Every catalogued TCGA project now works (16 of the 31 failed before).**
  `tcga_project_dataset()`, `tcga_surv_table()`, `tcga_merged()`, `tcga_get_expr()` and
  `canonical_type()` validated their input against a hard-coded 15-project vector
  (`tcga_retained`) while the catalog had grown to 31 TCGA cohorts, so ACC, ESCA, HNSC,
  KICH, KIRC, KIRP, MESO, PCPG, SARC, SKCM, TGCT, THCA, THYM, UCEC, UCS and UVM raised
  `Unsupported TCGA project` in every TCGA helper, `canonical_type()` returned `NA` for
  them, and the Shiny app (which builds its cohort lists from the catalog and routes every
  `TCGA-*` accession through those helpers) could not analyse them. The supported set is
  now **derived from the shipped catalog** — every `dataset_info` row whose `Accession`
  starts with `TCGA-`, with the `Type` column as the label — the exported `tcga_retained`
  object is refreshed from the catalog when the package is loaded (the previous vector
  survives only as an internal fallback), and the new
  `tests/testthat/test-tcga-catalog.R` fails as soon as the fallback and the catalog
  disagree. It needs no network and is **not** `skip_on_cran()`, unlike the four network
  tests in `test-tcga.R` that let this drift through unnoticed. Verified for all 31
  cohorts: `tcga_project_dataset()`, `tcga_surv_table()`,
  `tcga_merged(<project>, "TP53")`, `get_expr_data("TCGA-<project>", "TP53")` and
  `canonical_type()` all succeed (31/31, 2026-09-24).
* **An accession whose catalog spelling contains a dash is now readable (`A5-PCPG`).**
  `get_data()` sent the accession to the API verbatim, but the mirror names every table with
  an underscore and the API interpolates the name into its SQL: a dash-spelled table cannot
  be served at all (measured — an existing `X-Y` table answers HTTP 500 while its `X_Y`
  twin answers 200), so `get_data("A5-PCPG", "surv_data")`, `get_expr_data("A5-PCPG", …)`
  and `cohort_merged("A5-PCPG", …)` all failed with HTTP 500, and the cohort was
  unreachable in the app. `get_data()` now sends a dash as the underscore the mirror uses
  (`A5-PCPG` → table `A5_PCPG`) while the catalog accession stays authoritative and is
  echoed back unchanged; this is a no-op for the other 178 accessions (no other non-TCGA
  accession contains a dash). `HELP.md` documents the rule, and
  `pipeline/R/16_verify_catalog_mirror.R` gained an explicit check that no mirror table
  name contains a dash, so a future hyphenated accession cannot be uploaded into a table
  the API can never read. Verified: `A5-PCPG` survival (77 rows), expression (77 × 2) and
  `cohort_merged("A5-PCPG", "TP53", "OS")` all return data; `16_verify_catalog_mirror.R`
  reports 0 errors / 0 warnings / 0 info for all 179 cohorts.

Documentation and app-text corrections in the same round: `dataset_info.Rd` now documents
the `X` and `method` columns and no longer claims "one row per GEO dataset" (the catalog is
143 GEO + 31 TCGA + 3 CGGA + 2 cBioPortal-hosted); `cpas_km_pooled()` is documented as
returning a plain `list` (it has no class attribute); the `canonical_type()` and
`short_name()` examples use TCGA accessions and `short_name()` is documented as using only
the first element instead of being "vectorised"; the app's Methods page reports the current
overlap register (29 pairs in 13 groups), the current re-verification count (159 cohort
tables, 156 clean), the honest four source buckets (GEO 143 / TCGA 31 / CGGA 3 /
cBioPortal-hosted 2), the analysable-sample definition of `N`, and the full 45-step
pipeline inventory; `HELP.md` carries the updated family counts (OS 128, DSS 39, DFS 94,
PFS 45, MFS 17, "PFS or MFS" 62, any 179) and the fourth source bucket.

## Catalog expanded to 179 rows: two small cohorts admitted under a relaxed gate

Two small cohorts were added under an admission gate the author relaxed from **> 50 to
>= 30 patients per cohort**; both are flagged as such in the catalog `Note` field
(`small cohort: N=34; gate relaxed to >=30 per author` and the N=35 analogue). They are
the only catalog rows below 50 patients.

* **GSE76019** (Adrenocortical Cancer, N = 34, 12 EFS events). Expression from GEO
  GSE76019 on GPL13158 (already log-scale) and clinical from the same series'
  `!Sample_characteristics_ch1` (`efs.time` in years, `efs.event`). Only EFS is
  registered: the series carries no OS or DSS. The `EFS` token maps to the DFS family
  through the existing rule in `11_endpoint_families.R`, which was not changed. This is
  the paediatric COG ARAR0332 cohort (PMID 27307598).
* **GSE76039** (Thyroid Cancer, N = 35, 29 OS deaths). One study split across two
  deposits: expression from GEO GSE76039 on GPL570, and clinical from the same study's
  cBioPortal archive `thyroid_mskcc_2016` (poorly differentiated / anaplastic thyroid
  cancer, MSK, JCI 2016, PMID 26878173), joined on `!Sample_title` against
  `SAMPLE_ID` and cross-checked against `OTHER_SAMPLE_ID` (37/37 matched both ways).
  This is the same expression-plus-publication-table pattern as the earlier GSE3218
  build.

The catalog now holds **179 cohorts across 29 cancer types** (143 GEO, 31 TCGA,
3 CGGA and 2 cBioPortal-hosted studies, the last counted in their own bucket rather
than as GEO), endpoint families **OS 128 / DSS 39 / DFS 94 / PFS 45 / MFS 17** (counts
of non-NA `EP_*`; DFS now carries two `EFS` tokens) and **37,121 analysable samples**
(median 161 per cohort, range 34–1,210). `16_verify_catalog_mirror.R` reports
0 errors / 0 warnings / 0 info; `data(dataset_info)` was rebuilt and verified
cell-by-cell against the CSV (179 x 26 cells, 0 differences), and
`19_fix_catalog_N.R` was not run. The three still-open cancer types are unchanged:
Endometrial Cancer (best route 29 patients, one short of the gate), Uterine
Carcinosarcoma (time and status only inside KM figures) and Thymoma (no time-and-status
pair anywhere). The patient-overlap register is unchanged by this batch: recomputing it
with `pipeline/R/26_cohort_overlap.R` (dry run) against the 179-row catalog still gives
**29 pairs in 13 groups** among 32 cohorts.

## Catalog expanded to 177 rows: three cohorts closing three cancer types

Three cohorts built through the supplementary/bespoke route were added, each the first
independent (non-TCGA) cohort for a cancer type that until now had only a TCGA project:
**A5-PCPG** (Pheochromocytoma, N = 70, 19 OS events; cBioPortal clinical plus the same
study's RNA-seq CPM, GPL24676), **IMmotion150** (Kidney Cancer, N = 263, 164 PFS events;
cBioPortal clinical plus the same study's TPM, GPL24676) and **GSE3218** (Testicular Cancer,
N = 74, 27 OS events; survival from the study's publication table via Europe PMC, expression
from GEO GSE3218 on GPL96). The catalog now holds **177 cohorts across 29 cancer types**
(141 GEO, 31 TCGA, 3 CGGA and **2 cBioPortal-hosted studies, which are not GEO series**),
endpoint families **OS 127 / DSS 39 / DFS 93 / PFS 45 / MFS 17** (counts of non-NA `EP_*`)
and **37,052 analysable samples** (median 163 per cohort, range 36-1,210).
`16_verify_catalog_mirror.R` reports 0 errors / 0 warnings / 0 info; `data(dataset_info)`
was rebuilt and verified cell-by-cell against the CSV (177 x 26 cells, 0 differences), and
`19_fix_catalog_N.R` was not run. The two cBioPortal cohorts are counted in their own source
bucket rather than as GEO, so the GEO count is 141 and not 143: their clinical and
expression files are deposited together by the same study rather than through GEO.

**Mirror change beyond the cohort uploads.** The DB table `GPL24676` was stale (38,594
rows, built before the local probe map was extended), which left the two new cBioPortal
cohorts only about 21-30 % queryable through the app's `ID_map -> <GPL> -> expr` path. An
additive append of **26,137 probe-to-gene rows** (no row updated and none deleted) brings
the table to 123,646 rows, matching `data/processed/gpl/GPL24676.rds`; coverage is now
A5-PCPG 68.9 % and IMmotion150 64.8 %. Rollback artefact:
`pipeline/backup/gpl_db_pre_suppl_20260924/GPL24676_db_rows.csv.gz`.

Recorded rather than smoothed over:

* **IMmotion150**'s PFS is measured from the start of ICI treatment (the trial endpoint),
  matching the GSE159067 / GSE162520 convention, not from diagnosis.
* **A5-PCPG**'s longest OS is 456 months on a censored patient, taken verbatim from the
  source file, and the cBioPortal study's `cancerTypeId` is mislabelled `hnsc` although the
  cohort is the A5 Consortium PPGL series.
* **GSE3218**'s survival comes from the PMC4666461 S1 publication table (74 of its 108
  patients land in GSE3218; the 34-patient GSE10783 validation arm is below the 50-patient
  threshold and is not catalogued). GEO carries no clinical fields for GSE3218, and the
  table's 2-year DFS / 5-year DSS columns are milestone binaries, so they were deliberately
  not registered as DFS or DSS.

The patient-overlap register is unchanged by this batch: recomputing it with
`pipeline/R/26_cohort_overlap.R` (dry run) against the 177-row catalog still gives
**29 pairs in 13 groups** among 32 cohorts.

## Catalog expanded with 19 GEO cohorts (155 -> 174 rows)

Nineteen GEO cohorts built through the normal pipeline now carry expression, survival and
GPL tables in the mirror: four head and neck (GSE65858, GSE117973, GSE27020, GSE159067),
five melanoma (GSE65904, GSE198430, GSE198431, GSE22153 and the GeoMx cohort GSE325123),
three sarcoma (GSE71118, GSE30929, GSE271517), four lymphoma (GSE31312, GSE32918,
GSE23501, GSE248835), one uveal melanoma (GSE22138), one mesothelioma (GSE183088) and
GSE162520. The catalog now holds **174 cohorts across 29 cancer types** (140 GEO, 31 TCGA,
3 CGGA; `Lymphoma` is the new cancer type), endpoint families **OS 125 / DSS 39 / DFS 93 /
PFS 44 / MFS 17** (counts of non-NA `EP_*`), and **36,645 analysable samples** (median 163
per cohort, range 36-1,210). Every row carries a resolved endpoint and
`16_verify_catalog_mirror.R` reports 0 errors / 0 warnings / 0 info; `data(dataset_info)`
was rebuilt and verified cell-by-cell against the CSV.

Three points are recorded rather than smoothed over:

* **GSE162520** was labelled "Head and Neck" in the candidate list, but its GEO title, its
  overall design ("92 patients with surgically treated NSCLC") and every per-patient
  diagnosis are non-small cell lung cancer, so it is catalogued as **Lung Cancer**.
* **GSE31312** (lymphoma) is catalogued at the clinical-table size **N = 475 / 172 events**;
  joined to the GEO samples by depository id only **470 patients and 170 events** have an
  expression profile (5 clinical ids have no GEO sample; 28 further GEO samples have no
  clinical row). The GEO-side 470/170 is written into that row's `Note`.
* **GSE325123** (melanoma, GeoMx DSP) is catalogued at the **patient level, N = 105 /
  62 events** (the `included in_analysis` ROIs aggregated per patient); **102 of those
  patients have an OS time and 60 of them died**, and the 102/60 count is the GEO-side
  analysable figure in that row's `Note`.

The patient-overlap register behind `CohortGroup` and `Note` was rewritten for the 174-row
catalog with `pipeline/R/26_cohort_overlap.R --write`: **29 pairs in 13 groups** share
patients (8 detected by sample title, 14 by the same-series/different-platform rule, 7 by
the verified GEO seed list). One further title match — **GSE25066–GSE32918** — is a
sample-title collision rather than a shared patient set: GSE32918's titles are panel
replicate codes (`1`, `1_Rep1`, …; 77 of 249 titles) for 172 patients, its genuine
duplicate deposit is GSE69051 (deliberately not catalogued) and the two cohorts are of
different cancer types. It is therefore recorded as a `Note` instead of a group, listed in
`pipeline/ref/cohort_overlap_exclude.csv`, and the script preserves the provenance notes
already present in the catalog instead of overwriting the column.

The offline GEO screen behind the batch was repaired in
`pipeline/R/36_screen_geo_candidates.R` and `pipeline/R/37_geo_expansion_screen.R`: the
survival-keyword hint no longer acts as a hard gate (it had silently dropped real cohorts
such as GSE304059 and GSE22138), time/event roles are decided by the decoded token and the
values rather than by the key text, the event column is chosen with the maximum-overlap
time/event pair instead of the first decodable 0/1 column, and the header splitter was
vectorised. On the re-scan the passing candidates went from 23 to 37 with no cohort
dropped, and three parse defects that had reported GSE183088 as 1/1, GSE23501 with the
wrong `status` and GSE30929/GSE71118 with no endpoint were corrected.

## Licence change to MIT, `flextable` made optional, author metadata

* The package is now licensed **MIT** (`License: MIT + file LICENSE`; the full MIT text is in
  `LICENSE.md`, next to the two-field `LICENSE` stub R expects), replacing GPL-3. The motivation is reuse: MIT allows modification and
  closed-source derivatives, while GPL-3 already allowed modification but required
  derivatives to stay GPL-3. `COPYING` (the GPL-3 text) was removed.
* `Authors@R` now names the author and maintainer (**Jin Wang**,
  <Jinwang93@suda.edu.cn>), replacing the placeholder team entry.
* `flextable` moved from `Imports` to `Suggests`: the three-line summary table of the
  Cox page is still a `flextable` when the package is installed, and falls back to the
  plain data frame with the same numbers otherwise (the app renders it as a table and
  offers the CSV download; the Word download explains that `flextable` is needed).
* README licence/dependency sections updated.

## Catalog expanded with three supplementary GEO cohorts (152 -> 155 rows)

Three cohorts whose survival data sits outside the series matrix were built and uploaded
through the normal pipeline: **GSE108474** (REMBRANDT glioma, GPL570; time and event both
from the clinical supplement joined by subject id at 550/550 = 100%; OS, N = 476, 397 events;
eight all-NA `*_duplicate*` expression columns dropped, 550 -> 542; zero overlap with
CGGA_301/325/693 so it enters standalone), **GSE53625** (esophageal squamous cell carcinoma,
GPL18109; 358 samples, OS, 212 events) and **GSE102238** (pancreatic ductal adenocarcinoma,
GPL19072; 100 samples, OS, 60 events). The two Agilent platforms carried no gene annotation in
GEO, so probe maps were obtained with `AnnoProbe::idmap(type = "pipe")` (73.0% and 67.3% of
probes mapped to Entrez ids) in the usual `ID_REF` + `ENTREZ_GENE_ID` layout. Survival,
expression and GPL tables are in the mirror and `16_verify_catalog_mirror.R` reports 0 errors.
The catalog now holds 155 cohorts across 28 cancer types (121 GEO, 31 TCGA,
3 CGGA), families OS 116 / DSS 38 / DFS 91 / PFS 39 / MFS 14, and
33,814 analysable samples.

## Catalog expanded with 16 TCGA projects (136 -> 152 rows)

Sixteen TCGA projects that were present in the local clinical/survival tables but had never
been catalogued are now first-class cohorts, all with more than 50 analysable samples:
KIRC (603), THCA (571), HNSC (563), SKCM (456), KIRP (320), SARC (264), ESCA (195),
UCEC (193), PCPG (185), TGCT (137), THYM (120), KICH (89), MESO (86), UVM (79), ACC (77)
and UCS (57). This brings the catalog to 152 cohorts across 28 cancer types
(118 GEO, 31 TCGA, 3 CGGA), with endpoint families OS 113, DSS 38,
DFS 91, PFS 39, MFS 14; two projects were held back because they fall below the
threshold (CHOL 45, DLBC 47 OS samples). The per-project tables were generated with the
exact layout and clinical conventions of the existing TCGA tables - the generator reproduces
all fifteen pre-existing TCGA tables cell-by-cell (0 differences across 13 columns) - and
`data(dataset_info)` was rebuilt and verified against the CSV cell-by-cell (0 differences).
Note that TCGA cohorts are served from `data/tcga/*.rda` plus on-demand expression, so the
mirror still holds no TCGA tables by design.

## Catalog: endpoint-less cohorts removed entirely (152 -> 136 rows)

The thirteen rows listed below were removed first, and the three remaining unannotated rows
(GSE40272_GPL15971, GSE40272_GPL15972, GSE40272_GPL9497 - the prostate platform rows whose
DFS event counts are 2, 2 and 3) were removed in the same pass, so the catalog now contains
**no cohort without an endpoint annotation**. The tool therefore ships 136 cohorts
(118 GEO, 15 TCGA projects, 3 CGGA) across the same 14 cancer types, with the endpoint
families (OS 97, DSS 22, DFS 78, PFS 23, MFS 14), the 28,890 analysable samples and the
136 endpoint-annotated cohorts unchanged; the overlap structure is 26 pairs in 13 groups.
GSE40272_GPL15973 (DFS, 40 analysable samples, 10 events) remains as the prostate cohort.

## Catalog: 13 endpoint-less cohorts removed

Thirteen catalogue rows that carried no endpoint annotation were dropped from the catalog
(152 -> 139 rows; GEO 134 -> 121, TCGA 15 and CGGA 3 unchanged; 14 cancer types unchanged):
GSE10885_GPL5325, GSE10885_GPL885, GSE1378, GSE1379, GSE20624_GPL5325, GSE20624_GPL885,
GSE20624_GPL887, GSE22226_GPL4133, GSE35629_GPL1390, GSE35629_GPL5325, GSE35629_GPL887,
GSE32062_GPL570 and GSE14520_GPL571. They are either reference-array submissions with no
clinical annotation at all, or platform subsets of series whose patients are already
represented by other catalog rows (GSE10885/GSE20624/GSE35629 are the same breast
tissue-bank series deposited repeatedly, and GSE1378/GSE1379 are the same 60 patients
measured twice). Nothing analysis-facing changes: the endpoint families and their cohort
counts (OS 97, DSS 22, DFS 78, PFS 23, MFS 14), the 136 cohorts with a resolved endpoint,
the 28,890 analysable samples and the overlap structure among the remaining rows
(32 pairs in 14 groups) are unaffected except for the removal of pairs that involved these
rows. The three remaining unannotated rows are the GSE40272 prostate platform rows, whose
DFS event counts (2, 2 and 3) stay below the five-event annotation rule.

## Publication preparation: repository URL, licence file, no personal default path

* `DESCRIPTION` now carries `URL` and `BugReports` pointing at the public repository
  (`https://github.com/WangJin93/CanPAS`) and a rewritten `Description` that names the
  endpoint-family/token resolution, the default-visible safeguards and the bundled app.
* `COPYING` (the full GPL-3 text) was shipped at that point; it was removed in the
  licence change below, when the package moved to MIT.
* The `CPAS_DATA_ROOT` fallback no longer hard-codes a developer's home directory: it
  defaults to an empty value, and the TCGA helpers stop with instructions naming the
  variable and the expected `<root>/data/tcga/*.rda` files. Set the variable before
  TCGA or clinical-covariate analyses; GEO/CGGA mirror access and UCSC Xena gene
  fetching are unaffected.
* `README.md` rewritten for the public repository (install from GitHub, configuration,
  catalog and family counts, safeguards, app pages, cache, K-M colour convention,
  validation summary, layout, citation).
* `.gitignore` added for R session artefacts.

## Cut points: "top percent" now means the top x% of a cohort, and it is available when pooling

Two corrections to how high/low groups are formed, both visible in the app:

* **`Top percent (%)` was inverted.** On the single-dataset KM page, entering 25
  split the cohort at its 25th percentile, which put the *top 75%* of patients in
  High. The threshold is now the (100 - x)th percentile, so 25 means the highest
  quarter of patients, matching the label. The value box is restricted to 1-99.
* **The pooled KM page gained the same rule.** `cpas_km_pooled()` takes
  `cut = "top_pct"` with `top_pct` (1-99): each cohort is split at its own
  (100 - top_pct)th percentile, so `top_pct = 25` makes the highest quarter of every
  cohort High and the rest Low. The page exposes it as `Top percent (%)` next to the
  existing `50 (top 50% high)` median split, and the result records the rule, the
  per-cohort thresholds (`cohort_thresholds`) and the group sizes. `top_pct = 50`
  reproduces the median split exactly.
* **The absolute-threshold rule moved to the script interface.** `cut = "custom"`
  with `cut_value` (one log2 value applied in every cohort) is unchanged for scripts
  but is no longer offered on the pooled page, where the two percentile rules are the
  ones a reader can reproduce from the reported numbers. A per-cohort search for the
  best cut point remains deliberately unavailable (see the cut-point simulation).

Cohorts are no longer dropped silently from a pool: a cohort with no resolved
endpoint, a missing column, fewer than 10 complete (time, status, marker) rows or a
constant marker is returned in `skipped_cohorts` / `skipped_reasons` and printed in
the app status line, alongside the existing `empty_cohorts` for thresholds that empty
one side of a cohort.

## Multivariable Cox: a model that had to be reduced is returned and flagged

`COX_analysis(method = "multi")` no longer refuses the analysis when the
automatic repair cannot keep two covariates. It now searches down to a single
covariate, and if that is all that can be estimated the model is fitted,
returned, and flagged:

* `metadata$final_covariates` - the covariates actually in the model;
* `metadata$reduced` - `TRUE` when fewer than `min_covariates` covariates
  survived the repair;
* `metadata$reduced_note` - a sentence naming the survivors and the reasons the
  others could not be co-estimated, stating that the fit is not an adjusted
  model;
* an R warning carrying the same information.

`min_covariates` therefore marks the point below which a model stops being a
multivariable (adjusted) model; it is no longer a hard stop. Only a model with
no estimable covariate at all is refused, and that error now lists both the
covariates already removed and the reason the final fit failed, instead of
reporting only the last convergence message.

## App

* **`COX_screen_adjust()`** returns the multivariate metadata as
  `multi_metadata` and states in a message when that model had to be reduced.

* **COX analysis page.** The multivariable step is now isolated from the
  univariable step: a multivariable model that cannot be estimated no longer
  discards the univariable forest plot and table. The panel shows a red alert
  with the estimator reasons, the multivariable figure shows the reasons
  instead of an empty frame, and a reduced model is labelled on the figure
  itself ("REDUCED model: k covariate(s) - not adjusted for confounding") and
  in the status line.


* **Brand logo.** `inst/shiny/CanPAS/www/logo.jpg` is shown in the top-left brand
  block of the header (46 px tall, aspect preserved, served through the
  `cpas-assets` resource path so it resolves from any working directory). The
  block is painted `#f3f5f5` to match the logo's own background, and the block
  collapses to the logo alone when the sidebar is collapsed.
* **Dashboard.** The per-source count boxes ("GEO cohorts (curated mirror)",
  "CGGA cohorts (glioma)", "TCGA projects (on-demand Xena)") were removed; the
  "Database at a glance" card now reports cohorts in catalog, cancer types and
  survival endpoint families only.


* **Figure blocks.** On all six analysis pages the figure is now followed
  immediately by its own controls: width and height in pixels, and two download
  buttons (PNG, PDF). Changing the size resizes the figure in the app at once
  (the analysis result is already cached, so only the drawing repeats), and the
  download keeps exactly that size — PNG at the requested pixel size, PDF at the
  same physical size at 96 px per inch. Previously the size inputs sat below the
  result table, governed only the PDF, and the figure on screen was a fixed
  460–620 px. Pages with two figures (COX analysis: univariable and
  multivariable; Pooled KM: pooled curve and per-dataset grid) get one block per
  figure.
* Table downloads (CSV) stay below their tables, unchanged.
* Implemented in `inst/shiny/CanPAS/apps/core.R` (`.fig_dims()`, `ui_fig_block()`,
  `server_fig_block()`), so every page shares one implementation; the drawing
  function is used for the on-screen plot and for both downloads, which is what
  keeps the three from diverging. Regression tests in
  `tests/testthat/test-app-figure-block.R` pin the size contract (PNG pixels are
  read back from the file header) and check that no page keeps the old PDF panel.


* **TCGA cohorts now work in the package-level readers.** `get_expr_data()`
  fetched expression from the mirror for every accession, but the mirror holds no
  TCGA expression table, so `get_expr_data("TCGA-LUAD", ...)`,
  `merge_surv_expr()`, `COX_by_datasets()`, `cpas_meta()` and `cpas_meta_panel()`
  failed with `action=gpl&table=TCGA_HiSeqV2` HTTP 500 (the App was unaffected: it
  has its own TCGA wrapper). `get_expr_data()` now routes `TCGA-*` accessions to
  UCSC Xena (`tcga_gene_expr_df`, per gene, log2 values, no probe resolution;
  genes Xena does not carry are reported and skipped) and `merge_surv_expr()` reads
  the local TCGA clinical table instead of the mirror. Regression tests live in
  `tests/testthat/test-tcga.R` (skip without network or `UCSCXenaShiny`).
* **Endpoint annotation now follows a minimum-events rule** (`>= 5`, the
  `min_events` default of the pooling functions): GSE40272_GPL15971/GPL15972/
  GPL9497 (2, 2, 3 DFS events) lost their annotation while GSE40272_GPL15973
  (10 events) kept it.
* **Two cohorts that had no survival table were built and added**
  (`pipeline/R/25_fix_endpoint_annotations.R`): GSE5327 (58 breast cancers,
  metastasis-free survival, 11 events — the fields were present in the stamped
  pheno but no table had ever been generated) and GSE7849 (78 breast cancers,
  DFS, 14 events, parsed from the raw series matrix where the fields are stored
  as `Key = Value` pairs that the original parser dropped).
* **Cohort-patient overlap is now detected and reported**
  (`pipeline/R/26_cohort_overlap.R`). Comparing sample titles across cohorts finds
  15 pairs of cohorts in 10 groups that share patients — the same study on two
  platforms (GSE3494_GPL96/GPL97: 179 patients; GSE37642_GPL96/GPL97: 422;
  GSE9782_GPL96/GPL97: 264; GSE17536/GSE17537/GSE17538: 177 and 55) and whole
  series deposited twice (GSE2990 vs GSE6532_GPL96: 189; GSE11969 vs GSE13213: 87;
  GSE10885_GPL1390 vs GSE20624_GPL1390: 91). Pooling such a pair counts patients
  more than once, which narrows confidence intervals and biases I², Q and the
  prediction interval. The catalog gained `CohortGroup` and `Note` columns, the
  Datasets page shows the note, and the three multi-dataset pages display a
  warning for the current selection.
* **Datasets page** shows `N (analysable)` and `Events` (events of the primary
  endpoint), with `—` for cohorts that have no expression table, and a note
  explaining both columns.
* **Methods page** documents the expression scale, the conversion and the fact
  that hazard ratios are per 1 unit of the stored (log2) values.


* **Cancer type selector added to the three single-dataset pages** (KM analysis,
  COX analysis, COX by genes), matching the multi-dataset pages. The type lists
  every cancer type with its number of datasets (a single-dataset question may
  target a type with only one cohort), the dataset selector is filtered to that
  type, and the shared cohort is preselected when it belongs to the chosen type.
* All six analysis pages now start from the same selection path:
  Cancer type -> Dataset(s) -> Endpoint family.


* **COX analysis page: the inherited marker block was removed.** The page used to
  display "Marker (inherited from KM analysis)" and silently reuse the marker
  chosen on the KM page. It now has its own marker controls (gene + optional
  REF_ID probe, or a signature formula with the probe-collapsing rule), so the
  single-dataset Cox analysis is self-contained; the KM page no longer writes the
  shared marker state.
* **Cancer type selector added to the three multi-dataset pages** (COX by datasets,
  Pooled KM, Meta-analysis). The cancer type drives the dataset list (types with
  at least two datasets are offered), the shared cohort is preselected when it
  belongs to that type, and the endpoint family choices are derived from the
  selected datasets as before.


* The analysis result is frozen at run time: the status line, figure subtitle and
  download file names use the cohort/endpoint/marker snapshot, so changing the
  shared selection can no longer label one cohort's data with another's name.
* KM and COX pages no longer raise errors before the first click; failed runs
  clear the figure/table and refuse downloads instead of serving the previous
  result; the PDF device is always closed.
* The auto cut-point reports the number of candidate cut points searched, the
  naive p-value and (when `maxstat` is installed) the search-adjusted p-value,
  and the threshold is drawn in the figure.
* KM/COX status lines report events, group sizes, median follow-up, the rows
  excluded before analysis, the probe-collapsing rule actually used, the
  complete-case/events-per-variable accounting and the proportional-hazards test.
* The cache is session-scoped; the Multi-datasets page inherits the shared cohort
  selection; failures are surfaced as notifications; a new "Statistical caveats"
  card documents the cut-point bias, small-cohort and competing-risk limitations.

## Multivariable Cox: locate and drop the offending covariate instead of refusing the model

A multivariable model that cannot be estimated (complete separation, collinearity
or a constant covariate) used to fail as a whole:

```
Multivariate Cox model is not estimable: Ran out of iterations and did not
converge. Complete separation, collinear or constant covariates produce
meaningless coefficients; remove or combine the offending variable.
```

The user was left with no model, no plot and no indication of *which* covariate
was responsible. `COX_analysis(method = "multi")` now handles this in three
steps and reports every decision in `metadata$dropped_covariates`:

1. **Screening** (`.cpas_COX_screen()`): a continuous covariate that is constant
   or nearly constant (>= 95% of patients share one value) is dropped; a
   categorical level with fewer than 5 patients, with no events, or with events
   for every patient is handled before fitting, since those levels are the usual
   cause of separation.
2. **Level repair**: such levels are merged into one `Other` level (keeping their
   patients). If `Other` remains perfectly separated, the patients in those
   levels are excluded instead and the covariate is kept with its remaining
   levels; if fewer than two levels would remain, the covariate is dropped.
3. **Iterative removal**: whatever still makes the joint model non-estimable is
   removed one covariate at a time (the covariate behind a non-finite
   coefficient, an extreme coefficient or a large standard error). When the
   offender cannot be isolated — typically a collinear pair — each removal is
   tried in turn and the first estimable subset is kept.

The model is fitted on what remains, so the result table and forest plot are
produced as usual, and the reason for every removal, merge or exclusion is
recorded and printed. The safeguard is unchanged in substance: no meaningless
coefficient is published. The error is still raised when nothing can be
estimated, with the accumulated reasons instead of a generic message.

Example (GSE13507, marker + age + grade + M + N; 165 complete cases, 69 events):
`M` is dropped because its `M1` level holds 7 patients and all 7 had an event,
and `N` is kept with its `N3`/`NX` single-patient levels excluded, so the model
is estimated on marker, age, grade and `N` with 163 patients. Both decisions are
printed in the application status line and in a note above the forest plot.

New tests: `tests/testthat/test-multivariable-screening.R` (6 cases: estimable
model untouched, constant covariate, single-patient level, all-events level,
collinear pair, separation of the marker itself).

## Documentation

* `COX_analysis()` documentation regenerated, plus `dataset_info`,
  `get_expr_data` and `merge_surv_expr`, which had lost the TCGA-branch wording.
  `R CMD check` no longer reports codoc mismatches.

## Data (mirror + catalog)

* **Expression values brought to one scale: `log2(value + 1)`.** The GEO
  expression tables were stored exactly as deposited, so 28 datasets held linear
  intensities (platform-dependent units, e.g. MAS5 values with SD in the
  10^3-10^5 range) while 117 held log2 values (`pipeline/R/18_log2_transform.R`,
  `20_log2_transform_db.R`, `21_reupload_log2_tables.R`). Consequences of the old
  mix: hazard ratios per raw unit were unreadable at 4 decimals (HR = 1.0000-1.0002)
  and incomparable across cohorts, and the proportional-hazards model assumed
  log-hazard linear in raw intensity for some cohorts and in log2 intensity for
  others. The affected datasets are listed in
  `pipeline/out/REPORT_log2_transform.md`; the pre-conversion files are kept in
  `pipeline/backup/expr_pre_log2/`.
* **Hazard ratios fitted before the harmonisation are not comparable with new ones** for those
  28 datasets. Kaplan-Meier curves, log-rank tests and p-value ordering can also
  shift, because log2 is a nonlinear transformation. CGGA_301/325/693 were already
  on a log scale (`10_parse_cgga.R` declares this per dataset) and were not touched.
* **Survival-data defects repaired** (`pipeline/R/23_fix_surv_ids_and_endpoints.R`):
  the survival tables of GSE4573 (129 samples, 68 OS events) and GSE3494_GPL96/97
  (179 samples, 36 DSS events each) had lost their sample ids — `row_names` were
  row numbers, so they shared no sample with the expression tables and the
  cohorts read as unanalysable. Row alignment was verified against
  `data/pheno/<ACC>.rds` (`geo_accession`, ages matching row by row) before the
  ids were restored. GSE48075's `OS_status` column was entirely empty because the
  source field is `os censor` (uncensored/censored) rather than a parsed status;
  it is now derived from the source (73 analysable, 45 events) with the time
  column unchanged. GSE70768's primary endpoint moved from OS (1 event in 57
  records) to RFS (29 events in 41 records), following the existing rule that an
  endpoint is annotated only when the data support it.
* **No read-only tables remain in the mirror.** 11 of 293 tables (all lung
  expression tables) rejected writes with `Table 'X' is read only`; they were
  rebuilt as writable copies with `CREATE ... LIKE`, `INSERT ... SELECT`,
  `CHECKSUM` verification and an atomic `RENAME TABLE`
  (`pipeline/R/24_make_tables_writable.R`). Table names, row counts and contents
  are unchanged (checksums identical).
* **`dataset_info` sample-size columns corrected** (`pipeline/R/19_fix_catalog_N.R`).
  `N` is now the analysable sample count (expression AND non-missing time/status for
  the cohort's primary endpoint) instead of the legacy planned/expression size —
  GSE31210 was listed as 133 while 226 tumours are analysable, GSE14520_GPL3921 as
  445 while 221 are. Added `n_expr`, `n_surv`, `n_events`, per-family counts
  (`n_OS` ... `n_MFS`) and `expr_in_mirror`; 34 rows are now `NA` because they have
  no expression table (or no endpoint annotation) and cannot be analysed at gene
  level. The mirror upload gate (`06_upload_db.R`) now uses `n_expr`, and
  `16_verify_catalog_mirror.R` checks the catalog against the mirror on every run.

## Note on pooled analysis

`cpas_meta()` / `cpas_meta_panel()` standardise the marker within each dataset
before pooling, so pooled estimates were already scale-free. Single-dataset
hazard ratios printed per raw unit should still be reported per SD (or per log2
unit) when cohorts come from different platforms.

## Fixed: the meta-analysis forest plot could not be drawn

`plot_meta_forest()` still referenced the per-dataset column by its original name
(`table`). After the rename that symbol resolved to base `table()`, so the plot
failed with "cannot coerce type 'closure' to vector of type 'character'" - the
Multi-datasets page could not draw a forest plot at all. The plot now uses the
`dataset` column, and a regression test renders it from a `cpas_meta()` object.

Also in this release:

* `forest_plot()` accepts the `gene` column returned by `COX_by_genes()` (it used
  to require `Variates`, which is why the gene-panel forest plot failed), and
  `forest_plot()`, `plot_meta_forest()` and `plot_meta_panel()` take a `digits`
  argument (default 4) for the printed estimates.

## Shiny app restructured (figure above table, 4 decimals everywhere)

| Group | Pages |
|---|---|
| - | Dashboard, Datasets |
| **Single dataset analysis** | KM analysis · COX analysis (univariable + multivariable restricted to the significant covariates) · COX by genes (one model per gene, BH-FDR, forest plot) |
| **Multi-datasets analysis** | COX by datasets (same gene across datasets, forest plot) · Pooled KM (IPD pooling and/or time-point pooling of S(t), landmark table, per-dataset curves and estimates) · Meta-analysis (one gene/signature, or a per-gene panel with FDR) |
| - | Methods |

* every result page shows the figure above the table and downloads below it;
* all numbers are printed with 4 decimals (`.dtx4()` rounds server-side, p-values
  below 0.0001 are shown as `<0.0001`);
* the COX page gained a threshold input and now fits the multivariable model only
  with the covariates whose univariable p is below it, reporting the selection,
  the events per variable, the proportional-hazards test and skipped variables;
* the old single "Multi-datasets analysis" page (meta + pooled KM + cohort x gene
  grid) was split into the three pages above; the cohort x gene grid moved into
  the per-gene panel mode of the Meta-analysis page (`cpas_meta_panel()`).

## New: per-gene meta-analysis panel

`cpas_meta_panel(datasets, genes, type)` runs the two-stage meta-analysis of
`cpas_meta()` for every gene of a panel over the same datasets and returns one
table of pooled estimates with Benjamini-Hochberg FDR across the genes:

* retrieval happens **once per dataset for all genes** (a panel of G genes over D
  datasets costs D retrievals instead of G x D), and expression can also be
  supplied pre-merged through `merged =`;
* the table carries `k`, `datasets_used`, `total_n`, `total_events`, pooled
  `HR`/`lower`/`upper`, `p`, `P_adj`, `I2`, `tau2`, the prediction interval and
  the heterogeneity p per gene;
* `per_dataset` keeps the long table (gene x dataset with the endpoint token,
  n, events, HR, CI, p) and `endpoints` the token map, so a supplementary table
  can state which endpoint each dataset contributed;
* `errors` records, per gene, why a dataset was excluded (gene missing on the
  platform, constant marker, too few events), `fetch_errors` records retrieval
  failures, and `notes` collects the mixed-token warnings;
* the FDR family contains only genes that were actually pooled (k >= 2); a
  single-dataset row keeps `P_adj = NA` and is counted in
  `metadata$genes_single_dataset`;
* `plot_meta_panel()` draws the pooled hazard ratios as a forest plot with point
  size by number of pooled datasets, colour by FDR significance and the per-gene
  I2 in the labels; `print.cpas_meta_panel()` summarises coverage.

`cpas_meta()` also reports the per-dataset reasons when no dataset can be pooled
("no dataset produced estimable results. Reasons: ...") instead of a bare message.

## Breaking: `quick_uni_multi_COX()` renamed to `COX_screen_adjust()`

The old name described neither the input nor the model: the function screens every
candidate covariate univariably, keeps those below `p.threshold`, fits ONE
multivariable model with them and formats both sets of estimates as a flextable.
The new name states the workflow (screen, then adjust) and the documentation now
spells out that the selection is exploratory and that the adjusted p-values are
not corrected for it. The returned list is unchanged.

## Breaking: the last three naming inconsistencies

* `get_data(dataset, action, ids)` — the argument was `table`; the documentation
  now states that the value is the mirror table behind the action (a dataset
  accession for `"expression"` / `"surv_data"`, a platform id such as `"GPL96"`
  for `"gpl"`).
* `forest_plot(COX_out, ...)` — the argument was `cox_out`.
* Returned structures use the dataset vocabulary:
  * `cpas_meta()` returns `per_dataset` (was `per_cohort`) and that table's first
    column is `dataset` (was `table`);
  * `cpas_km_pooled()` returns `df$dataset`, `datasets` and `dataset_endpoints`
    (were `df$cohort`, `cohorts`, `cohort_endpoints`);
  * `plot_cpas_km_percohort()` is now `plot_cpas_km_perdataset()`;
  * `print.cpas_meta()` prints "Per-dataset:".
  Prose keeps the word cohort where it describes the study unit
  ("one cohort per dataset" is not enforced); only identifiers changed.

Update calls with:

```sh
sed -i 's/per_cohort/per_dataset/g; s/cohort_endpoints/dataset_endpoints/g; \
        s/plot_cpas_km_percohort/plot_cpas_km_perdataset/g; s/cox_out/COX_out/g' *.R
```

## Breaking: one vocabulary for the study axis (`dataset`) and one spelling of COX

Follows the earlier rename; the old names are **removed** (no aliases, no
deprecation warnings), since the package is released only in this version.

Dataset arguments (were `table =`, `acc =`, `project =`):

| Function | New argument | Was |
|---|---|---|
| `get_expr_data(dataset, genes)` | `dataset` | `table` |
| `merge_surv_expr(dataset, expr_data)` | `dataset` | `table` |
| `get_signature_value(signature, dataset)` | `dataset` | `table` |
| `quick_km(dataset, genes)`, `quick_roc(dataset, genes)` | `dataset` | `table` |
| `endpoint_resolve(dataset, type)`, `endpoint_options(dataset)` | `dataset` | `acc` |
| `cohort_merged(dataset, genes, type)` | `dataset` | `acc` |
| `tcga_project_dataset(dataset)`, `tcga_gene_expr_df(dataset, gene)`, `tcga_get_expr(dataset, genes)`, `tcga_surv_merged` helpers | `dataset` | `project` |

The TCGA helpers now accept both forms of the id (`"LUAD"` and `"TCGA-LUAD"`).
`get_data(table = , action = )` keeps `table`, because that low-level accessor
addresses any mirror table (a dataset, a `GPL` platform or a survival table), not
a dataset specifically.

COX spelling in identifiers: `COX_analysis()` (unchanged), `COX_by_genes()`
(was `cox_by_genes()`), `COX_by_datasets()` (was `cox_by_datasets()`),
`competing_risk_COX()` (was `competing_risk_cox()`), `quick_uni_multi_COX()`
(was `quick_uni_multi_cox()`); classes `cpas_COX`, `cpas_COX_by_genes`,
`cpas_COX_by_datasets` with the matching `print()` methods. `survival::coxph()`
and `survival::cox.zph()` are untouched.

Complete rename of the earlier internal names (both steps together):

```sh
grep -rl 'univariable_cox\|multisets_COX' --include='*.R' . | \
  xargs sed -i 's/univariable_cox(/COX_by_genes(/g; s/multisets_COX(/COX_by_datasets(/g; \
                s/cox_by_genes(/COX_by_genes(/g; s/cox_by_datasets(/COX_by_datasets(/g; \
                s/competing_risk_cox(/competing_risk_COX(/g; s/quick_uni_multi_cox(/quick_uni_multi_COX(/g; \
                s/markers =/genes =/g; s/tables =/datasets =/g; s/cohorts =/datasets =/g'
```

## Breaking: two functions renamed to name the axis they vary

The package is not released yet, so the old names are **removed** (no aliases, no
deprecation warnings).

| Purpose | New name | Removed name |
|---|---|---|
| One dataset, many genes: one univariable Cox model per gene, Benjamini-Hochberg FDR across the genes | `COX_by_genes(df, type, genes)` | `univariable_cox()` |
| Many datasets, one gene: one univariable Cox model per dataset (no pooling; use `cpas_meta()` to pool) | `COX_by_datasets(datasets, gene, type)` | `multisets_COX()` |

Consequences of the rename:

* arguments follow the same vocabulary: `genes =` (was `markers =`) and
  `datasets =` (was `tables =` / `cohorts =`);
* `cpas_meta()` now takes its study axis through `datasets =`;
* the result classes are `cpas_COX_by_genes` and `cpas_COX_by_datasets`
  (`print()` methods renamed accordingly), and the `COX_by_genes()` results table
  column `ID` is now `gene`, with `metadata$failed_genes` (was `failed_markers`);
* the Shiny gene panel and `HELP.md` use the new names.

Rename all call sites with, for example:

```sh
grep -rl 'univariable_cox\|multisets_COX' --include='*.R' . | \
  xargs sed -i 's/univariable_cox(/COX_by_genes(/g; s/multisets_COX(/COX_by_datasets(/g; s/markers =/genes =/g'
```

## Fine-Gray subdistribution standard errors were too small

`competing_risk_COX()` fitted the weighted expansion produced by
`survival::finegray()` without a usable cluster term: the patient id was attached
after the expansion (and `finegray()` keeps only the variables named in its
formula, so the cluster argument silently evaluated to `NULL`). The reported
standard error was therefore the naive model-based one.

* Measured deviation from `cmprsk::crr()` before the fix: **19.6%, 21.0%, 23.1%,
  19.6% (user's scenarios: 14.7-20.6%) - i.e. the SE was 15-23% too small.**
* After the fix (the id is carried through the expansion formula and used as
  `cluster`): deviation **0.01-0.07%** across five simulated scenarios, with
  coefficients unchanged and still matching `crr` to ~1e-6.
* A regression test now requires the SE to agree with `cmprsk::crr()` within 2%.

## Complete separation was still published as a hazard ratio

`coxph()` does not fail on complete separation: it warns ("Loglik converged
before variable 1; coefficient may be infinite") and returns a finite but
meaningless estimate (e.g. HR 1.78e+09 with a 95% CI of [0, Inf]), so the
previous guard (finite coefficient and SE, `fit$fail == 0`) did not trigger and
`failed_markers` stayed empty.

* New internal helper `.cpas_COX_diag()` captures the warnings emitted during the
  fit and applies numeric red flags (non-finite coefficient/SE, SE <= 0,
  non-convergence, |log HR| > 20, SE > 10 on the log-HR scale) together with
  `fit$fail`.
* `univariable_cox()` now returns an `NA` row, records the marker in
  `failed_markers` and adds `metadata$failure_reasons` with the reason
  (separation, non-convergence, ...).
* `COX_analysis()` skips the offending covariate in `method = "uni"` and refuses
  the model in `method = "multi"` with an explanatory error; the skipped reasons
  are kept in `metadata$skipped_variables`.
* `cpas_meta()` excludes a cohort whose marker effect is not estimable and
  records the reason in `errors`, instead of pooling an HR of 1e+09.
* `coxph()` fits now keep `model = TRUE, x = TRUE` so that `survival::cox.zph()`
  can rebuild the data outside the fitting helper.

## New analysis features

* **Competing risks** — `cif_fit()` computes Aalen-Johansen cumulative incidence
  functions (with pointwise confidence intervals and optional group split),
  `plot_cif()` draws them, and `competing_risk_COX()` reports the cause-specific
  Cox model and the Fine-Gray subdistribution hazard model side by side.
  Implemented with `survival::finegray()` + weighted `coxph()` (cluster-robust
  variance); the point estimates were validated against `cmprsk::crr()`.
* **Multiple-testing correction** — `univariable_cox()` and `multisets_COX()` now
  return a Benjamini-Hochberg `P_adj` column (plus `P_adj_text` and
  `metadata$n_significant_fdr005` for the former). The Shiny gene panel shows the
  FDR, sorts by it, and lists the skipped cohort x gene combinations.
* **Meta-analysis reporting** — `cpas_meta()`/`meta_pool()` now return a 95%
  prediction interval for a new cohort (`pi_lower`, `pi_upper`; defined for
  k >= 3), `print.cpas_meta()` shows it and the forest plot annotates it.

## Statistical guards (previously silent failure modes)

* A single estimable cohort no longer produces `NaN` pooled estimates; the
  heterogeneity statistics are `NA` and the cohort's own estimate is returned.
* Signatures reject unsupported text (`log(0.5*GAPDH)+0.5*ACTB`, `1/0*GAPDH`) and
  parse scientific notation (`1e-3*GAPDH` is 0.001, not -3).
* `univariable_cox()`, `plot_km()`, `COX_analysis()` and `cpas_km_pooled()`
  validate the status coding (0/1) and follow-up times (>= 0).
* Non-convergent, separated or zero-variance fits are reported as failures
  (`failed_markers`, `skipped_variables`) instead of finite hazard ratios with
  95% CI [0, Inf].
* `COX_analysis()` warns when events per variable is below 10, skips a covariate
  that cannot be estimated instead of aborting, keeps non-syntactic covariate
  names, and returns a `cox.zph` proportional-hazards test.
* `plot_roc()` refuses horizons at which the time-dependent AUC is not estimable
  (beyond the longest follow-up, no event before it, or no control after it).
* `merge_surv_expr()` rejects duplicated expression sample IDs and reports zero
  overlap instead of returning an empty 0 x 0 data frame.
* `cpas_km_pooled()` errors when no cohort contributes an event, handles
  landmarks beyond the follow-up without producing `NaN`, and warns when cohorts
  inside one family contribute different endpoint tokens.
* `loo_meta()` works with two cohorts (and errors clearly with one).

## Documentation and API

* `cohorts` is the preferred argument name for `cpas_meta()`/`multisets_COX()`
  (`tables` still works) and `genes` for `get_expr_data()` (`gene` still works).
* Seven `@examples` blocks that were silently dropped by roxygen2 (unclosed
  `\dontrun{`) are restored; `cpas_km_pooled()` has a documented `@return`;
  `run_cpas_app()` documents the packages it really needs.
* README.md added; credentials are read from the environment
  (`CPAS_DB_PASSWORD`, `CPAS_DATA_ROOT`) — no hard-coded database password
  anywhere in the repository.

## Data corrections and pipeline

* `pipeline/R/16_verify_catalog_mirror.R` checks the catalog against the mirror
  (dash-spelled accessions, duplicated accessions, annotated-but-unavailable rows,
  orphan tables, survival status coding) and exits non-zero on inconsistencies, so
  it can run in CI.
* `pipeline/run_pipeline.R` is a driver for the whole pipeline
  (`--list`, `--check`, `--from/--to`, `--only`, `--dry-run`) with prerequisite
  checks and a run log.
* **GSE86166 (breast, n=366) survival status recoded 1/2 -> 0/1**
  (`pipeline/R/17_recod_gse86166_status.R`). The source codes `vital status`
  1 = alive / 2 = dead and `recurrence` 1 = recurrence / 2 = none, so the mirrored
  OS/RFS status columns were inverted relative to the package convention. Evidence:
  the associated publication (PMID 28629479 / PMC5477261) reports "71.9% were
  alive", i.e. 263 of 366 patients, exactly the count of `vital status = 1`, and
  in the data 284 of 284 patients with `recurrence = 2` have RFS time equal to OS
  time (censored at last follow-up) while 69 of 71 with `recurrence = 1` have
  RFS < OS. The cohort now yields 103 OS events (28.1%) and 71 RFS events.

## Endpoint first, then cohorts

The analysis flow now starts from the endpoint rather than from the cohort list,
so the app can only offer cohorts that can actually answer the chosen endpoint.

* **Datasets page**: the endpoint family filter sits above the cancer type. Once a
  single family is selected the table gains a column showing the token each cohort
  contributes at that endpoint (e.g. `RFS` for GSE31210 at DFS, `DFI` for TCGA-LUAD),
  so "which cohorts can answer this question" is readable at a glance instead of
  being inferred from the "Endpoint families" column.
* **COX by datasets / Pooled KM / Meta-analysis**: the order is now
  Cancer type -> **Endpoint family** -> Datasets. The dataset list contains only
  cohorts that carry the chosen family (for Bladder cancer, DFS offers the two
  DFS cohorts while OS offers all four), and the family list is labelled with the
  number of cohorts that carry it. Before this, cohorts were chosen first and any
  cohort without the family was dropped silently when the family was picked.
* **Shared selection**: choosing a cohort keeps the endpoint family you are working
  in whenever that cohort supports it; it now only falls back to the cohort's first
  family when it does not. Previously the family was reset on every cohort change,
  so a family chosen on the Datasets page was silently replaced.
* The user-facing "can this cohort answer this family" rule is defined once
  (`.family_hit()` in the app's `core.R`) and uses the same `EP_*` catalog columns
  the analysis functions resolve tokens from, so the filter and the analysis cannot
  disagree; cohorts that reach a page anyway (through the shared selection) are
  reported in a notification instead of being dropped quietly.

## Kaplan-Meier group colours

`plot_km()`, `plot_cpas_km()` and `plot_cpas_km_perdataset()` now draw the
**high-expression group in red and the low-expression group in green**. Before
this change they used survminer's default hue palette, which gives its first
colour to the first factor level; since the levels are `c("Low", "High")`, the
low group came out in the warm colour and the high group in the cool one.

* Change the two colours with `options(CanPAS.km_palette = c(low, high))`, or for
  a single call with `plot_km(..., palette = c("blue", "orange"))` — an explicit
  `palette =` always wins over the option.
* The pooled KM page now offers exactly two cut rules: **50 (top 50% high)** — each
  cohort split at its own median — and **Custom**, one absolute threshold applied in
  every cohort (above = High, at or below = Low), which is comparable across cohorts
  because every mirrored matrix is log2. A cohort the threshold leaves one-sided is
  reported (`empty_cohorts`) and excluded from the pool; a latent bug meant such a
  cohort used to reach the final check instead (`min(table(group)) < 1` can never be
  true), which is fixed. The status line prints the rule, the threshold and any
  excluded cohort.
* `plot_cpas_km()` gained `which = c("auto", "ipd", "meta")`. The default keeps the
  documented behaviour (the IPD curve whenever that route was run), but a
  `method = "both"` result holds both routes, and the Shiny pooled-KM page now draws
  them **together** in one figure — previously it silently showed only the IPD curve
  because a single call returns that whenever it exists. The meta curve now also
  marks the requested landmark estimates with their confidence intervals, and the
  pooled figure block's default canvas grew to 680x780 to hold both panels.
* `plot_cpas_km_perdataset()` gained `draw = FALSE`, which returns a `patchwork`
  object instead of drawing immediately. The `grid.arrange()` return value cannot
  be composed reliably (patchwork drops it silently), so figures that place this
  grid next to other panels should use `draw = FALSE`.
* Fixed a silent panel loss in the multi-dataset figure scripts: a
  `ggsurvplot` object has no `+` method, so `plot_cpas_km(km) + ggtitle(...)`
  returned `NULL` and the pooled Kaplan-Meier panel disappeared from Fig. 6-8.
  The scripts now add the title to `$plot` and compose the risk table explicitly.

## Remote fetches are cached on disk

`get_data()` (the mirror API: expression, probe-to-gene and survival tables) and
`tcga_gene_expr_df()` (UCSC Xena, one gene at a time) now keep their answers in a
local file cache, so a second analysis of the same cohort or gene is offline and a
transient outage no longer aborts a sweep.

* Default location: `tools::R_user_dir("CanPAS", "cache")`; move it with
  `options(CanPAS.cache_dir = "...")` or the environment variable
  `CANPAS_CACHE_DIR` (also `CPAS_CACHE_DIR`).
* Switch it off per R session with `options(CanPAS.cache = FALSE)`, per call with
  `get_data(..., use_cache = FALSE)` / `tcga_gene_expr_df(..., use_cache = FALSE)`,
  or per process with `CANPAS_CACHE=0`.
* Entries expire after `options(CanPAS.cache_ttl = <seconds>)` (30 days by default).
  An expired entry is reused **only** when the live request fails, with a message
  naming the entry's age, so a rate-limited or unreachable endpoint degrades to the
  last known copy instead of an error.
* Every entry stores the request key, the creation time and a schema version;
  entries that fail to parse, belong to another key or carry an older schema are
  ignored and refetched. Writes are atomic (temp file + rename), so an interrupted
  run cannot leave a half-written entry behind.
* `get_data()` reports `metadata$from_cache` and `metadata$cache_stale`, so a caller
  can tell a fresh download from a cache answer.
* Inspection and removal (internal helpers, no export added):
  `CanPAS:::.cpas_cache_info()` lists entries with size and age,
  `CanPAS:::.cpas_cache_clear(sub = NULL)` removes them —
  or delete the cache directory.

## Known limitations (unchanged)

Competing-risk adjustment is available only when the source data distinguish the
cause of death (the mirrored GEO survival tables code non-disease death as
censored); no time-varying covariates, landmark/immortal-time correction,
prediction-interval graphics for KM, or FDR-based multiplicity control for the
KM/COX pages (single-marker analyses).
