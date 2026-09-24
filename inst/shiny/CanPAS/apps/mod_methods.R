# Methods module: data processing pipeline & endpoint pooling --------------
# Documentation page. Static tables describe the pipeline actually used in
# pipeline/R/*; the family table is filled live from the catalog so it always
# matches the shipped data.

.methods_card <- function(title, ..., icon = NULL) {
  bs4Dash::bs4Card(
    title = title, width = 12, status = "primary", solidHeader = FALSE,
    collapsible = TRUE, collapsed = FALSE, ...
  )
}

.methods_tbl <- function(df, cols = NULL) {
  if (!is.null(cols)) df <- df[, cols, drop = FALSE]
  hdr <- tags$tr(lapply(names(df), function(n) tags$th(n)))
  rows <- lapply(seq_len(nrow(df)), function(i)
    tags$tr(lapply(df[i, , drop = FALSE], function(v) tags$td(HTML(as.character(v))))))
  tags$table(class = "table table-sm table-striped",
             tags$thead(hdr), tags$tbody(rows))
}

ui_mod_methods <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(12,
      h2("Methods: data processing and endpoint pooling"),
      h5("How the cohorts, survival tables and endpoint families in this app were produced"),
      hr()
    ),
    column(12,
      .methods_card(
        "1 · Data sources",
        uiOutput(ns("sources")),
        p(class = "note", "Expression values are read from the MySQL mirror at analysis
          time; TCGA expression is fetched on demand from UCSC Xena and TCGA clinical
          tables are stored locally. Survival tables always come from the mirror
          (GEO, CGGA, cBioPortal-hosted studies) or the local TCGA table (TCGA)."),
        p(class = "note", tags$b("Expression scale: "), "all mirrored GEO and CGGA
          matrices are on the log2 scale. 28 datasets were deposited by GEO as linear
          intensities (platform-dependent units, e.g. MAS5 values in the 10^3-10^5
          range) and were converted with ", tags$code("log2(value + 1)"),
          " (scripts 18 / 20 / 21); CGGA RNA-seq matrices come from the source as
          log2(RSEM + 1). Analysis therefore models log-hazard as a linear function of
          log2 expression for every cohort. Hazard ratios are still reported ",
          tags$b("per 1 unit of the stored values"), " (i.e. per 1 log2 unit, roughly a
          doubling of expression) and are not comparable with estimates computed before
          this conversion.")
      ),
      .methods_card(
        "2 · Processing pipeline (produce → clean → mirror)",
        p("Each step is a script under ", tags$code("pipeline/R/"), "; the outputs are the
          files the package ships and the tables the mirror serves."),
        uiOutput(ns("pipeline"))
      ),
      .methods_card(
        "3 · Cleaning and standardisation rules",
        uiOutput(ns("clean")),
        p(class = "note", "Standardised cohorts were re-verified one by one
          (159 cohort tables: 156 clean, 3 with one flagged column each; report: ",
          tags$code("REPORT_clinical_normalization_verify.md"),
          "). Tokens that could not be mapped to a standard level are listed in
          ", tags$code("clinical_std_stage_unmapped.csv"), " rather than silently coerced.")
      ),
      .methods_card(
        "4 · Endpoint families: pooling principles",
        tags$ol(
          tags$li(tags$b("A family groups tokens that answer the same clinical question."),
                  " Cohorts that report RFS and cohorts that report DFI can therefore be
                  analysed together as DFS, instead of being scattered over four
                  separate endpoint labels."),
          tags$li(tags$b("Pooling happens inside a family, never across families."),
                  " OS is never mixed with DSS, and DFS is never mixed with PFS."),
          tags$li(tags$b("The token is resolved per cohort and always reported."),
                  " An analysis is requested by family; the resolver picks the token each
                  cohort actually provides. Every result table, KM risk table and status
                  line states the token used (e.g. 'DFS (RFS)', 'DFS (DFI) [derived]')."),
          tags$li(tags$b("Derived endpoints are flagged."),
                  " TCGA DFI/PFI are derived endpoints; they are marked ", tags$code("[derived]"),
                  " in the endpoint selector and recorded in the ", tags$code("EndpointDerived"),
                  " catalog column."),
          tags$li(tags$b("An endpoint is annotated only when the data support it."),
                  " A token needs events and real follow-up: in GSE40272 the RFS censored
                  patients all have time 0 (no follow-up time), so RFS was not annotated,
                  and OS has 1 event out of 84, so OS was not annotated either; the
                  per-platform DFS tables were annotated instead."),
          tags$li(tags$b("Heterogeneity is reported, not hidden."),
                  " Cross-cohort estimates are pooled with a random-effects model by
                  default and shown together with the per-cohort estimates, the I2
                  statistic and a leave-one-out table."),
          tags$li(tags$b("Rows excluded before analysis are counted."),
                  " Samples without the selected endpoint's time/status, or without a marker
                  value, are removed before fitting and the number removed is printed in
                  the status line, so the reported n is the n actually entering the model.")
        ),
        uiOutput(ns("famtable"))
      ),
      .methods_card(
        "5 · Statistical caveats and what the app reports",
        tags$ol(
          tags$li(tags$b("Cut point."), " The default is the cohort median (a pre-specified rule).
                  The ", tags$code("Auto"), " option searches 13 percentile cut points and keeps the
                  one with the largest log-rank statistic, so its naive p-value is biased upwards.
                  A simulation with a marker that is ", tags$b("independent"), " of survival gives
                  p < 0.05 in 23.1-26.7% of datasets for the auto rule versus 4.45-5.90% for the
                  median split (2000 replicates at each of n = 50, 100, 200 and 400). When the ",
                  tags$code("maxstat"), " package is installed the app also prints the
                  search-adjusted p-value (Lau92, falling back to exactGauss); with that adjustment
                  the same simulation gives 2.0-3.5%, i.e. the correction is conservative rather
                  than exact. Otherwise the app states that the p-value is unadjusted and
                  exploratory. The cut point itself is drawn in the figure, not only in the side
                  panel."),
          tags$li(tags$b("Small cohorts."), " Every KM/COX status line now reports the group sizes,
                  the number of events and the median follow-up, and flags a cohort with fewer than
                  10 events. Multivariable models warn when events per variable is below 10."),
          tags$li(tags$b("Proportional hazards."), " The multivariable Cox output includes a
                  ", tags$code("cox.zph"), " GLOBAL test; a p below 0.05 is printed together with a
                  caution that the HRs should be interpreted with care."),
          tags$li(tags$b("Endpoint tokens."), " Analysis is requested by family and the token each
                  cohort contributes is always reported. When cohorts inside one family contribute
                  different tokens (for example RFS and DFI under DFS) the pooling function warns and
                  the per-cohort token stays in the result table."),
          tags$li(tags$b("Sample sizes."), " The reported n is the number of rows that actually
                  entered the model; rows removed for a missing endpoint time/status or a missing
                  marker are counted and printed."),
          tags$li(tags$b("Landmarks beyond follow-up."), " A landmark later than a cohort's longest
                  follow-up carries that cohort's last observed S(t) forward, and the number of
                  contributing cohorts (k) is reported for every landmark."),
          tags$li(tags$b("Not implemented."), " Competing-risks models (Fine-Gray), time-varying
                  covariates, landmark/immortal-time corrections, multiple-testing correction for
                  gene panels and prediction-interval reporting in the meta-analysis are not
                  provided; treat panel-level p-values as screening results."))
      ),
      .methods_card(
        "5b · Cohorts that share patients",
        p("Some cohorts in the catalog are not independent: a study measured on two
          platforms appears as two accessions, and a few series were deposited twice
          under different GSE numbers. Comparing the sample titles across all cohorts
          (", tags$code("26_cohort_overlap.R"), ") finds 29 such pairs (32 cohorts) in
          13 groups,
          e.g. GSE3494_GPL96/GPL97 (179 shared patients), GSE37642_GPL96/GPL97 (422),
          GSE9782_GPL96/GPL97 (264), GSE17536 and GSE17537 against GSE17538_GPL570
          (177 and 55), GSE2990 against GSE6532_GPL96 (189), GSE11969 against
          GSE13213 (87) and GSE10885_GPL1390 against GSE20624_GPL1390 (91)."),
        p(class = "note", "Pooling two cohorts of the same group counts those patients
          more than once, which narrows the confidence interval and biases I², Q and
          the prediction interval — the reported precision is then partly an artefact
          of double counting. The catalog records the group in ", tags$code("CohortGroup"),
          " and a readable flag in ", tags$code("Note"), "; the Datasets page shows the
          note and every multi-dataset page warns when the current selection contains
          two members of one group. Keep one cohort per group, or state the overlap
          explicitly.")
      ),
      .methods_card(
        "6 · Catalog columns and how the app uses them",
        uiOutput(ns("columns"))
      )
    )
  )
}

server_mod_methods <- function(id, dataset_info) {
  moduleServer(id, function(input, output, session) {
    di <- .norm_catalog(dataset_info)

    output$sources <- renderUI({
      src <- .cohort_source(di$Accession)
      d <- data.frame(
        Source = c("GEO", "CGGA", "TCGA", "cBioPortal-hosted (non-GEO)"),
        Cohorts = c(sum(src == "GEO"), sum(src == "CGGA"), sum(src == "TCGA"),
                    sum(src == "cBioPortal")),
        Expression = c("MySQL mirror (gene-level queries over the API)",
                       "MySQL mirror",
                       "UCSC Xena, fetched per gene on demand",
                       "MySQL mirror (the study's own RNA-seq matrix: log2 CPM for A5-PCPG, log2 TPM for IMmotion150)"),
        Survival = c("MySQL mirror (&lt;ACC&gt;_surv)",
                     "MySQL mirror (CGGA_&lt;ID&gt;_surv)",
                     "Local clinical/survival table",
                     "MySQL mirror (&lt;ACC&gt;_surv); clinical patient files deposited with the study"),
        `Cohort ids` = c("GSE…, GSE…_GPL…",
                         paste(di$Accession[src == "CGGA"], collapse = ", "),
                         paste0("TCGA-", c("BLCA", "BRCA", "…"), collapse = ", "),
                         paste(di$Accession[src == "cBioPortal"], collapse = ", ")),
        check.names = FALSE, stringsAsFactors = FALSE)
      .methods_tbl(d)
    })

    output$pipeline <- renderUI({
      d <- data.frame(
        Step = c("01 parse", "02 / 12 platform map", "03 survival table",
                 "04 QC", "05 catalog plan", "06 mirror upload",
                 "07 clinical standardisation", "08 verification",
                 "09 survival sync", "10 CGGA", "11 endpoint families",
                 "13 small cohorts", "14 / 15 per-platform split",
                 "16 catalog ↔ mirror check", "17 status recoding",
                 "18 / 20 expression scale", "19 catalog sample sizes",
                 "22 log2 impact report", "23 / 24 data repairs",
                 "25 / 26 annotation & overlap", "27 catalog notes",
                 "35–38 TCGA / GEO expansion",
                 "90–99 supplementary & expansion builds",
                 "100–102 relaxed gate & platform names"),
        Script = c("01_parse_gse.R", "02_gpl_map.R, 12_gpl_map_symbol.R",
                   "03_surv_table.R", "04_qc_report.R", "05_dataset_plan.R",
                   "06_upload_db.R", "07_standardize_clinical.R, 07b, 07c",
                   "08_verify_normalization.R", "09_update_db_surv.R",
                   "10_parse_cgga.R", "11_endpoint_families.R",
                   "13_upload_small_cohorts.R", "14_split_gse40272_surv.R, 15_split_gse40272_expr.R",
                   "16_verify_catalog_mirror.R", "17_recod_gse86166_status.R",
                   "18_log2_transform.R, 20_log2_transform_db.R, 21_reupload_log2_tables.R",
                   "19_fix_catalog_N.R", "22_report_log2_impact.R",
                   "23_fix_surv_ids_and_endpoints.R, 24_make_tables_writable.R",
                   "25_fix_endpoint_annotations.R, 26_cohort_overlap.R",
                   "27_catalog_notes.R",
                   "35_add_tcga_cohorts.R, 36_screen_geo_candidates.R, 37_geo_expansion_screen.R, 38_geo_expansion_annotation_check.R",
                   "90_build_gse108474_suppl.R, 91_complete_gse14520_surv.R, 92_build_geo_expansion_expr_pheno.R, 93_build_geo_expansion_bespoke.R, 94_update_catalog_geo_expansion.R, 95_build_suppl_expansion.R, 96_build_suppl_expansion.R, 98_update_catalog_suppl_expansion.R, 99_extend_gpl_db.R",
                   "100_build_relaxed_gate.R, 101_update_catalog_relaxed_gate.R, 102_extend_gpl4133_agilent_name.R"),
        What = c(
          "GEO series matrix -> expression table (ID_REF + one column per sample)",
          "Platform annotation -> probe to Entrez map; symbol-based mapping when the platform has no Entrez column",
          "Clinical/survival table per cohort (endpoint time in years, status 0/1)",
          "Sample-level QC report",
          "Cohort plan and catalog (accession, cancer type, platform, method)",
          "Expression and platform tables into the MySQL mirror",
          "Stage / T / N / M, sex, grade and age brought to one vocabulary",
          "Re-check of every standardised cohort (per-column token audit)",
          "Standardised survival tables synchronised to the mirror",
          "CGGA glioma cohorts parsed and uploaded (3 cohorts)",
          "Endpoint tokens grouped into families; catalog annotated per cohort",
          "Cohorts below the historical N > 79 gate added to the mirror",
          "Studies stored only at GSE level split per platform (sample lists from the GEO series matrix)",
          "Standing consistency gate between the catalog and the mirror: accession spelling, orphan tables, sample sizes and 0/1 status coding",
          "Endpoint status column of GSE86166 recoded from its published coding to 0/1",
          "Expression values brought to one scale: log2(value + 1) for the 28 datasets deposited as linear intensities (files were stored as deposited by GEO)",
          "Catalog sample sizes recomputed from the mirror: N = analysable samples, plus n_expr / n_surv / n_events and per-family counts",
          "Before/after report of the expression-scale conversion (which cohorts changed, by how much)",
          "Survival-table repairs: restored lost sample ids (GSE4573, GSE3494), derived an empty status column from its source field (GSE48075) and re-annotated an endpoint with 1 event (GSE70768); read-only MyISAM tables rebuilt as writable copies with checksum verification",
          "Endpoint annotation rule (>= 5 events) enforced per cohort; two cohorts with no survival table built from their source files (GSE5327 MFS, GSE7849 DFS); cohort-patient overlap detected from sample titles and recorded in the catalog (CohortGroup/Note)",
          "CohortGroup / Note bookkeeping written back into the catalog without touching the data tables",
          "16 TCGA cohorts added (catalog TCGA rows 15 -> 31); candidate GEO series screened and their annotations checked before inclusion",
          "Supplementary-file and expansion builds for individual cohorts (GSE108474, GSE14520, the GEO expansion sets); GPL24676 extended additively for the two cBioPortal-hosted studies",
          "Two small cohorts built and catalogued under the relaxed gate (N >= 30); GPL4133 Agilent platform name extended"),
        check.names = FALSE, stringsAsFactors = FALSE)
      tagList(
        .methods_tbl(d),
        p(class = "note", tags$code("pipeline/R/"), " holds ", tags$b("45 numbered steps"),
          " (01–27, 35–38, 90–102; the table above lists every one of them) plus 9 helper,
          demo and validation scripts (", tags$code("batch_integrate.R"),
          ", ", tags$code("batch_integrate2.R"), ", ", tags$code("demo_meta_lung.R"),
          ", ", tags$code("demo_tcga_integration.R"), ", ", tags$code("demo_tcga_ondemand.R"),
          ", ", tags$code("demo_unified_reader.R"), ", ", tags$code("test_cpas_dataset.R"),
          ", ", tags$code("test_cpas_GSE44001.R"), ", ", tags$code("validate_cpas.R"),
          "), i.e. 54 R files.")
      )
    })

    output$clean <- renderUI({
      d <- data.frame(
        Field = c("Endpoint time", "Endpoint status", "Stage", "T", "N", "M",
                  "Sex", "Grade", "Age", "Missing tokens"),
        Rule = c("Converted to years; patients without follow-up time are excluded from the analysis and counted",
                 "1 = event, 0 = censored/alive; other codings remapped",
                 "0, I, II, III, IV",
                 "T0, T1, T2, T3, T4, Tis, TX",
                 "N0, N1, N2, N3, NX",
                 "M0, M1, MX",
                 "male / female",
                 "G1–G4 (numeric scores 5–10 kept as published)",
                 "Numeric",
                 "-, NA, unknown, not reported → NA (never silently kept as a level)"),
        check.names = FALSE, stringsAsFactors = FALSE)
      .methods_tbl(d)
    })

    output$famtable <- renderUI({
      fams <- list(OS = "OS", DSS = c("DSS", "CSS", "BCSS"),
                   DFS = c("DFS", "RFS", "EFS", "DFI"),
                   PFS = c("PFS", "PFI"), MFS = c("MFS", "DRFS"))
      rows <- lapply(names(fams), function(f) {
        col <- paste0("EP_", f)
        tok <- if (col %in% names(di)) as.character(di[[col]]) else rep(NA_character_, nrow(di))
        used <- sort(unique(tok[!is.na(tok)]))
        der  <- if ("EndpointDerived" %in% names(di))
          as.character(di$EndpointDerived) else rep(NA_character_, nrow(di))
        n_der <- sum(!is.na(tok) & !is.na(der) &
                       mapply(function(t, d) grepl(t, d, fixed = TRUE), tok, der))
        data.frame(Family = f,
                   `Tokens pooled` = paste(fams[[f]], collapse = ", "),
                   `Token used in this catalog` = if (length(used)) paste(used, collapse = ", ") else "—",
                   Cohorts = sum(!is.na(tok)),
                   `Of which derived` = n_der,
                   check.names = FALSE, stringsAsFactors = FALSE)
      })
      also <- data.frame(
        Family = "PFS (browser view)",
        `Tokens pooled` = "PFS, PFI + MFS, DRFS",
        `Token used in this catalog` = "as above",
        Cohorts = sum(!is.na(di$EP_PFS) | !is.na(di$EP_MFS)),
        `Of which derived` = sum(grepl("PFI", as.character(di$EndpointDerived), fixed = TRUE)),
        check.names = FALSE, stringsAsFactors = FALSE)
      d <- do.call(rbind, c(rows, list(also)))
      tagList(
        h5("Family to token mapping", class = "section-title"),
        .methods_tbl(d),
        p(class = "note", "The Datasets page lists every family separately (OS, DSS,
          DFS, PFS, MFS) with its cohort count and additionally offers the broader
          'PFS or MFS' view, which is the union of the two; the analysis pages keep
          PFS and MFS as separate families and never pool them together.")
      )
    })

    output$columns <- renderUI({
      d <- data.frame(
        Column = c("Accession", "Type", "SurvivalTypes", "EndpointFamilies",
                   "EP_OS, EP_DSS, EP_DFS, EP_PFS, EP_MFS", "EndpointPrimary",
                   "EndpointDerived", "GPL", "N", "method"),
        Meaning = c("Cohort identifier as used by the mirror and the API",
                    "Cancer type (grouping key of the Multi-datasets page)",
                    "Raw endpoint tokens available for the cohort",
                    "Families the cohort can be pooled under",
                    "The token this cohort contributes to each family (NA when absent)",
                    "Token used as the default endpoint on the analysis pages",
                    "Tokens that are derived rather than directly reported",
                    "Platform of the expression table",
                    "Analysable samples: expression and the EndpointPrimary time/status are both present (so N differs from n_expr for cohorts whose endpoint covers fewer samples)",
                    "Assay / data type"),
        UsedBy = c("every page", "Datasets filter, Multi-datasets grouping",
                   "Datasets table", "Datasets filter, all family selectors",
                   "family resolution (KM, COX, COX by genes, COX by datasets, Pooled KM, Meta-analysis)",
                   "preselected endpoint", "derived-endpoint flag",
                   "Datasets table, REF_ID lookup", "Datasets table",
                   "Datasets table"),
        check.names = FALSE, stringsAsFactors = FALSE)
      .methods_tbl(d)
    })
  })
}
