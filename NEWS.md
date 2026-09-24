# CanPAS 1.0.0

First public release. CanPAS is a curated cross-archive cancer prognosis resource
(GEO mirror, CGGA, TCGA), a scripted curation pipeline and an R package with a
bundled Shiny application; this section documents the state of that first release.

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
