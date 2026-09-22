# The analysis pages keep their parameter labels short and put the explanation
# into a hover popup next to the label (the browser's own title tooltip, so it
# cannot be clipped by the sidebar's scroll container). The probe box also says
# in plain words what the current selection does, in particular that the
# default takes the maximum over all probes of the gene.

labels_ready <- function() {
  for (p in c("shiny", "bs4Dash", "DT")) if (!requireNamespace(p, quietly = TRUE)) return(FALSE)
  f <- system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS")
  nzchar(f) && file.exists(f)
}

label_env <- function() {
  suppressMessages({ library(shiny); library(bs4Dash); library(DT) })
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  env <- new.env(parent = globalenv())
  for (f in c("core.R", "mod_welcome.R", "mod_datasets.R", "mod_km.R", "mod_cox.R",
              "mod_cox_by_genes.R", "mod_cox_by_datasets.R", "mod_pooled_km.R",
              "mod_meta.R", "mod_methods.R", "mod_help.R"))
    sys.source(file.path(ad, f), envir = env)
  env
}

test_that(".ref_hint explains all-probes, single-probe and TCGA cases", {
  skip_if_not(labels_ready(), "shiny app files unavailable")
  env <- label_env()
  h <- get(".ref_hint", envir = env)
  all_probes <- h("GSE14814", "GAPDH", "__all__", 6L)
  expect_match(all_probes, "all probes")
  expect_match(all_probes, "maximum")
  expect_match(all_probes, "6 on this platform")
  one <- h("GSE14814", "GAPDH", "212581_x_at", 6L)
  expect_match(one, "212581_x_at")
  expect_match(one, "without collapsing")
  tcga <- h("TCGA-LUAD", "TP53", "__all__", 1L)
  expect_match(tcga, "no probe concept")
  expect_match(tcga, "Xena")
  expect_null(h(NULL, "GAPDH", "__all__"))
  expect_null(h("GSE14814", "", "__all__"))
})

test_that("parameter labels are short and carry a hover explanation", {
  skip_if_not(labels_ready(), "shiny app files unavailable")
  env <- label_env()
  pages <- c("ui_mod_km", "ui_mod_cox", "ui_mod_cox_by_genes", "ui_mod_cox_by_datasets",
             "ui_mod_pooled_km", "ui_mod_meta")
  for (p in pages) {
    h <- paste(as.character(get(p, envir = env)("x")), collapse = "")
    expect_gte(lengths(regmatches(h, gregexpr("title=", h))), 3L)  # explanations present
    # the old sentence-style labels are gone from the sidebar
    for (old in c("continuous / categorical are detected automatically",
                  "keep covariates with", "Covariates adjusted inside each dataset",
                  "Gene symbols (comma separated)", "Minimum events per dataset",
                  "Datasets (>= 2, multi-select)", "Cut point (group definition)"))
      expect_false(grepl(old, h, fixed = TRUE), info = paste(p, old))
  }
  # and the explanation text really is in the tooltip of the probe control
  cox <- paste(as.character(get("ui_mod_cox", envir = env)("cox")), collapse = "")
  expect_match(cox, "Auto: all probes", fixed = TRUE)
  expect_match(cox, "collapsed per sample by the maximum", fixed = TRUE)
})

test_that("the probe hint follows the selection in the COX module", {
  skip_if_not(labels_ready(), "shiny app files unavailable")
  env <- label_env()
  env$.ref_probes <- function(acc, gene) c("p1", "p2", "p3")
  rv <- shiny::reactiveValues(sel = list(acc = "GSE14814", type = "Test",
                                        family = "OS", ep = "OS"))
  testServer(get("server_mod_cox", envir = env),
             args = list(id = "cox", rv = rv, dataset_info = CanPAS::dataset_info), {
    session$setInputs(kind = "gene", gene = "GAPDH", ref = "__all__", rule = "max",
                      show_all = FALSE, clin = character(0), p_thr = 0.05)
    txt <- gsub("<[^>]*>", "", paste(as.character(output$ref_hint), collapse = ""))
    expect_match(txt, "all probes")
    expect_match(txt, "maximum")
    expect_match(txt, "3 on this platform")
    session$setInputs(ref = "p2")
    txt2 <- gsub("<[^>]*>", "", paste(as.character(output$ref_hint), collapse = ""))
    expect_match(txt2, "p2")
    expect_match(txt2, "without collapsing")
  })
})

# --- Endpoint-first selection: the family decides which cohorts are offered ---
# The Datasets page and the three multi-dataset pages list the endpoint family
# before the cohorts, and the cohort list is restricted to cohorts that carry
# that family. These tests pin the helpers behind that filter, which are shared
# by the UI and (through the same EP_* columns) by the analysis path, so the two
# cannot disagree about which cohorts can answer an endpoint.

endpoint_ready <- function() {
  f <- system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS")
  nzchar(f) && file.exists(f)
}

ep_env <- function() {
  env <- new.env(parent = globalenv())
  sys.source(system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS"), envir = env)
  env$catalog <- env$.norm_catalog(CanPAS::dataset_info)
  env
}

test_that("cohorts are filtered to those carrying the chosen endpoint family", {
  skip_if_not(endpoint_ready(), "shiny app files unavailable")
  e <- ep_env(); cat_d <- e$catalog
  os <- e$.cohorts_with_family(cat_d, "Lung Cancer", "OS")
  dfs <- e$.cohorts_with_family(cat_d, "Lung Cancer", "DFS")
  expect_gt(length(os), 0)
  expect_true(all(os %in% cat_d$Accession))
  ## the DFS list must be a strict subset of the lung cohorts and every member
  ## must really carry DFS in the catalog annotation
  expect_true(all(dfs %in% cat_d$Accession[cat_d$Type == "Lung Cancer"]))
  expect_true(all(vapply(dfs, function(a) {
    o <- e$.endpoint_opts(cat_d, a); "DFS" %in% o$family }, logical(1))))
  expect_lt(length(dfs), length(e$.peer_datasets(cat_d, "Lung Cancer")))
  ## 'PFS or MFS' is the union
  expect_setequal(e$.cohorts_with_family(cat_d, "Lung Cancer", "__pfsmfs__"),
                  union(e$.cohorts_with_family(cat_d, "Lung Cancer", "PFS"),
                        e$.cohorts_with_family(cat_d, "Lung Cancer", "MFS")))
})

test_that("the family list counts the cohorts that carry each endpoint", {
  skip_if_not(endpoint_ready(), "shiny app files unavailable")
  e <- ep_env()
  ch <- e$.ep_choices_for_type(e$catalog, "Lung Cancer")
  expect_true("OS" %in% unname(ch))
  labs <- names(ch)
  expect_true(any(grepl("^OS \\([0-9]+\\)$", labs)))
  ## the count in the label equals the number of cohorts returned for that family
  i <- which(unname(ch) == "OS")
  n_lab <- as.integer(sub("^OS \\(([0-9]+)\\)$", "\\1", labs[i]))
  expect_identical(n_lab, length(e$.cohorts_with_family(e$catalog, "Lung Cancer", "OS")))
  ## a type with no endpoint at all still yields a usable selector
  expect_true(length(e$.ep_choices_for_type(e$catalog, "Nonsense Type")) >= 1)
})

test_that("changing family drops non-qualifying cohorts and reports them", {
  skip_if_not(endpoint_ready(), "shiny app files unavailable")
  e <- ep_env(); cat_d <- e$catalog
  all_lung <- e$.peer_datasets(cat_d, "Lung Cancer")
  r <- e$.refresh_cohorts(cat_d, "Lung Cancer", "DFS", selected = all_lung)
  expect_true(all(r$selected %in% r$choices))
  expect_true(all(r$choices %in% e$.cohorts_with_family(cat_d, "Lung Cancer", "DFS")))
  expect_setequal(r$dropped, setdiff(all_lung, r$choices))   # nothing vanishes silently
  ## an empty selection falls back to the first cohorts that qualify
  r2 <- e$.refresh_cohorts(cat_d, "Lung Cancer", "DFS", selected = character(0))
  expect_gt(length(r2$selected), 0)
})

test_that("the token column reports the token a cohort contributes", {
  skip_if_not(endpoint_ready(), "shiny app files unavailable")
  e <- ep_env(); cat_d <- e$catalog
  ## GSE31210 reports RFS, which belongs to the DFS family
  expect_identical(e$.family_token_label(cat_d, "GSE31210", "DFS"), "RFS")
  expect_identical(e$.family_token_label(cat_d, "TCGA-LUAD", "DFS"), "DFI")
  ## a family the cohort does not have yields an empty cell
  expect_identical(e$.family_token_label(cat_d, "GSE31210", "MFS"), "")
  ## OS is its own token
  expect_identical(e$.family_token_label(cat_d, "TCGA-LIHC", "OS"), "OS")
})

test_that("selecting a cohort keeps the family the user is working in", {
  skip_if_not(endpoint_ready(), "shiny app files unavailable")
  e <- ep_env()
  ## the app passes a reactiveValues object, so the helper's writes have to be
  ## observed inside isolate() (reactive values cannot be read outside a
  ## reactive context, and a plain list would not see the writes at all)
  rv <- shiny::reactiveValues(sel = list(acc = NULL, family = NULL, ep = NULL))
  shiny::isolate({
    ## user works at DFS in GSE31210 ...
    e$.set_shared_selection(rv, e$catalog, "GSE31210")
    rv$sel$family <- "DFS"
    ## ... then picks another cohort that also has DFS: the family must survive
    e$.set_shared_selection(rv, e$catalog, "GSE37745")
    expect_identical(rv$sel$family, "DFS")
    expect_identical(rv$sel$ep, "RFS")
    ## a cohort that does not have DFS falls back to one of its own families
    e$.set_shared_selection(rv, e$catalog, "TCGA-LUAD")
    expect_true(rv$sel$family %in% e$.endpoint_opts(e$catalog, "TCGA-LUAD")$family)
  })
})
