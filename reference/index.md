# Package index

## Endpoints and endpoint semantics

The pooling families, the token each cohort resolves to, and the pooling
class of every cohort x family cell.

- [`endpoint_family()`](https://wangjin93.github.io/CanPAS/reference/endpoint_family.md)
  : Pooling family of a survival endpoint
- [`endpoint_resolve()`](https://wangjin93.github.io/CanPAS/reference/endpoint_resolve.md)
  : Resolve an endpoint family for one cohort
- [`endpoint_options()`](https://wangjin93.github.io/CanPAS/reference/endpoint_options.md)
  : Endpoints available for a cohort, grouped by family
- [`endpoint_semantics()`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md)
  : Endpoint semantics of the curated cohorts
- [`endpoint_pooling_class()`](https://wangjin93.github.io/CanPAS/reference/endpoint_pooling_class.md)
  : Pooling class of one cohort x endpoint-family cell

## Endpoint adjudication

A blinded rating sheet for the endpoint annotation, and the agreement
calculator for a filled sheet.

- [`endpoint_adjudication()`](https://wangjin93.github.io/CanPAS/reference/endpoint_adjudication.md)
  : Blinded endpoint-adjudication rating sheet
- [`endpoint_agreement()`](https://wangjin93.github.io/CanPAS/reference/endpoint_agreement.md)
  : Agreement between endpoint raters

## Data access

The uniform “merged data” schema and the local/TCGA sources behind it.

- [`get_data()`](https://wangjin93.github.io/CanPAS/reference/get_data.md)
  : Query the CanPAS GEO API
- [`get_expr_data()`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)
  : Retrieve Expression Data for Genes of a Dataset
- [`merge_surv_expr()`](https://wangjin93.github.io/CanPAS/reference/merge_surv_expr.md)
  : Merge Survival and Expression Data of a Dataset
- [`cohort_merged()`](https://wangjin93.github.io/CanPAS/reference/cohort_merged.md)
  : Unified CanPAS cohort reader (GEO or TCGA)
- [`get_signature_value()`](https://wangjin93.github.io/CanPAS/reference/get_signature_value.md)
  : Weighted Gene-Signature Score per Sample
- [`tcga_get_expr()`](https://wangjin93.github.io/CanPAS/reference/tcga_get_expr.md)
  : Multi-gene TCGA expression (wide format)
- [`tcga_merged()`](https://wangjin93.github.io/CanPAS/reference/tcga_merged.md)
  : TCGA cohort merged with expression (CanPAS format)
- [`tcga_surv_table()`](https://wangjin93.github.io/CanPAS/reference/tcga_surv_table.md)
  : TCGA clinical + survival table
- [`tcga_local_tables()`](https://wangjin93.github.io/CanPAS/reference/tcga_local_tables.md)
  : Local TCGA clinical/survival tables used by the TCGA helpers
- [`tcga_project_dataset()`](https://wangjin93.github.io/CanPAS/reference/tcga_project_dataset.md)
  : Xena dataset id of a TCGA project
- [`tcga_retained`](https://wangjin93.github.io/CanPAS/reference/tcga_retained.md)
  : TCGA projects available through the TCGA helpers
- [`tcga_gene_expr_df()`](https://wangjin93.github.io/CanPAS/reference/tcga_gene_expr_df.md)
  : Fetch one gene's TCGA expression from UCSC Xena
- [`canonical_type()`](https://wangjin93.github.io/CanPAS/reference/canonical_type.md)
  : Canonical CanPAS cancer type of a TCGA project
- [`short_name()`](https://wangjin93.github.io/CanPAS/reference/short_name.md)
  : Short dataset name of a TCGA project

## Survival analysis

- [`plot_km()`](https://wangjin93.github.io/CanPAS/reference/plot_km.md)
  : Kaplan-Meier Survival Plot by Marker
- [`quick_km()`](https://wangjin93.github.io/CanPAS/reference/quick_km.md)
  : Quick Kaplan-Meier Plots for a Dataset
- [`plot_roc()`](https://wangjin93.github.io/CanPAS/reference/plot_roc.md)
  : Time-dependent ROC Curve for a Marker
- [`quick_roc()`](https://wangjin93.github.io/CanPAS/reference/quick_roc.md)
  : Quick Time-dependent ROC Plots for a Dataset
- [`COX_analysis()`](https://wangjin93.github.io/CanPAS/reference/COX_analysis.md)
  : Univariate or Multivariate Cox Regression
- [`COX_by_genes()`](https://wangjin93.github.io/CanPAS/reference/COX_by_genes.md)
  : Cox Regression for Many Genes (Panel Screening)
- [`COX_by_datasets()`](https://wangjin93.github.io/CanPAS/reference/COX_by_datasets.md)
  : Cox Regression for Many Datasets (Single Gene)
- [`COX_screen_adjust()`](https://wangjin93.github.io/CanPAS/reference/COX_screen_adjust.md)
  : Univariate -\> Multivariate Cox Analysis Workflow
- [`forest_plot()`](https://wangjin93.github.io/CanPAS/reference/forest_plot.md)
  : Forest Plot of Cox Regression Results
- [`competing_risk_COX()`](https://wangjin93.github.io/CanPAS/reference/competing_risk_COX.md)
  : Cause-specific and subdistribution (Fine-Gray) Cox models
- [`cif_fit()`](https://wangjin93.github.io/CanPAS/reference/cif_fit.md)
  : Cumulative incidence functions for competing risks
- [`plot_cif()`](https://wangjin93.github.io/CanPAS/reference/plot_cif.md)
  : Plot cumulative incidence curves

## Cross-cohort meta-analysis

Two-stage pooling with method = REML / DL / HK / FE, pooling = family /
exact and overlap = warn / refuse / dedupe.

- [`cpas_meta()`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
  : Two-stage cross-cohort integrative (meta) analysis
- [`cpas_meta_panel()`](https://wangjin93.github.io/CanPAS/reference/cpas_meta_panel.md)
  : Per-gene meta-analysis panel across datasets
- [`plot_meta_forest()`](https://wangjin93.github.io/CanPAS/reference/plot_meta_forest.md)
  : Forest plot of an integrative (meta) analysis
- [`plot_meta_panel()`](https://wangjin93.github.io/CanPAS/reference/plot_meta_panel.md)
  : Plot a per-gene meta-analysis panel
- [`loo_meta()`](https://wangjin93.github.io/CanPAS/reference/loo_meta.md)
  : Sensitivity analysis: leave-one-out
- [`cpas_km_pooled()`](https://wangjin93.github.io/CanPAS/reference/cpas_km_pooled.md)
  : Pooled Kaplan-Meier (marker split inside every cohort)
- [`plot_cpas_km()`](https://wangjin93.github.io/CanPAS/reference/plot_cpas_km.md)
  : Plot pooled / meta Kaplan-Meier result
- [`plot_cpas_km_perdataset()`](https://wangjin93.github.io/CanPAS/reference/plot_cpas_km_perdataset.md)
  : Grid of per-cohort Kaplan-Meier panels

## Shared-patient overlap

- [`cohort_overlap()`](https://wangjin93.github.io/CanPAS/reference/cohort_overlap.md)
  : Shared-patient overlap register

## Batch diagnostics

The per-cohort expression distribution and the within-cohort cut-points,
before anything is pooled.

- [`batch_diagnostics()`](https://wangjin93.github.io/CanPAS/reference/batch_diagnostics.md)
  : Per-cohort expression diagnostics for a multi-cohort analysis
- [`plot_batch_diagnostics()`](https://wangjin93.github.io/CanPAS/reference/plot_batch_diagnostics.md)
  : Plot per-cohort expression distributions

## Reproducibility

The per-result analysis manifest, and the index of the curation steps,
registers and additional files.

- [`cpas_manifest()`](https://wangjin93.github.io/CanPAS/reference/cpas_manifest.md)
  : Analysis manifest of a CanPAS result
- [`cpas_reproducibility_index()`](https://wangjin93.github.io/CanPAS/reference/cpas_reproducibility_index.md)
  : Index of the reproducibility material
- [`cpas_pipeline_steps()`](https://wangjin93.github.io/CanPAS/reference/cpas_pipeline_steps.md)
  : The 67 curation steps of the shipped pipeline

## Package and application

- [`run_cpas_app()`](https://wangjin93.github.io/CanPAS/reference/run_cpas_app.md)
  : Launch the CanPAS Shiny application

## Datasets

- [`dataset_info`](https://wangjin93.github.io/CanPAS/reference/dataset_info.md)
  : CanPAS dataset catalog
- [`ID_map`](https://wangjin93.github.io/CanPAS/reference/ID_map.md) :
  Gene-symbol to Entrez mapping (internal reference)

## Print and conversion methods

The S3 print methods of the result objects, and the manifest’s
as.data.frame().

- [`print(`*`<cpas_COX>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_COX.md)
  : Print Method for cpas_COX Class
- [`print(`*`<cpas_COX_by_datasets>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_COX_by_datasets.md)
  : Print Method for cpas_COX_by_datasets Class
- [`print(`*`<cpas_COX_by_genes>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_COX_by_genes.md)
  : Print Method for cpas_COX_by_genes Class
- [`print(`*`<cpas_adjudication_agreement>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_adjudication_agreement.md)
  : Print an endpoint-adjudication agreement result
- [`print(`*`<cpas_batch_diagnostics>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_batch_diagnostics.md)
  : Print a batch-diagnostics result
- [`print(`*`<cpas_cif>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_cif.md)
  : Print method for cpas_cif objects
- [`print(`*`<cpas_competing>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_competing.md)
  : Print method for cpas_competing objects
- [`print(`*`<cpas_get>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_get.md)
  : Print Method for cpas_get Class
- [`print(`*`<cpas_get_expr>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_get_expr.md)
  : Print Method for cpas_get_expr Class
- [`print(`*`<cpas_manifest>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_manifest.md)
  : Print an analysis manifest
- [`print(`*`<cpas_merge>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_merge.md)
  : Print Method for cpas_merge Class
- [`print(`*`<cpas_meta>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_meta.md)
  : Print method for cpas_meta objects
- [`print(`*`<cpas_meta_panel>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_meta_panel.md)
  : Print method for cpas_meta_panel objects
- [`print(`*`<cpas_not_estimable>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_not_estimable.md)
  : Print Method for cpas_not_estimable Class
- [`print(`*`<cpas_signature>`*`)`](https://wangjin93.github.io/CanPAS/reference/print.cpas_signature.md)
  : Print Method for cpas_signature Class
- [`as.data.frame(`*`<cpas_manifest>`*`)`](https://wangjin93.github.io/CanPAS/reference/as.data.frame.cpas_manifest.md)
  : Convert an analysis manifest to a data.frame
