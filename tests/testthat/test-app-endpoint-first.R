# Endpoint-first selection on the three single-dataset pages ---------------
# The KM / COX / COX by genes pages choose the endpoint FAMILY first, exactly
# like the Datasets page and the multi-dataset pages: the cancer-type and the
# dataset lists then contain only the cohorts that actually carry that family,
# so an analysis can never be started on a cohort that lacks the endpoint. The
# concrete token a cohort contributes is still resolved per cohort.
#
# Two layers are pinned here:
#   * the shared core.R helpers (the same ones the UI and the analysis path
#     use), checked directly against the packaged `dataset_info` catalog;
#   * the three module servers, driven with shiny::testServer(). The mock
#     session does not apply update*Input() messages to its inputs, so the
#     choices are captured by stubbing updateSelectInput()/updateSelectizeInput()
#     in the module environment and the browser's acknowledgement is emulated
#     with session$setInputs() - i.e. exactly what the real client does.
#
# Every check is printed as a PASS/FAIL table by the last test in this file.

ep_first_ready <- function() {
  for (p in c("shiny", "bs4Dash", "DT")) if (!requireNamespace(p, quietly = TRUE)) return(FALSE)
  f <- system.file("shiny", "CanPAS", "apps", "core.R", package = "CanPAS")
  nzchar(f) && file.exists(f)
}

ep_first_env <- function() {
  suppressMessages({ library(shiny); library(bs4Dash); library(DT) })
  ad <- system.file("shiny", "CanPAS", "apps", package = "CanPAS")
  env <- new.env(parent = globalenv())
  for (f in c("core.R", "mod_welcome.R", "mod_datasets.R", "mod_km.R", "mod_cox.R",
              "mod_cox_by_genes.R", "mod_cox_by_datasets.R", "mod_pooled_km.R",
              "mod_meta.R", "mod_methods.R", "mod_help.R"))
    sys.source(file.path(ad, f), envir = env)
  env$catalog <- env$.norm_catalog(CanPAS::dataset_info)
  env
}

# ---- PASS/FAIL table ------------------------------------------------------
EP_FIRST <- new.env(parent = emptyenv())
EP_FIRST$rows <- data.frame(check = character(0), result = character(0),
                            detail = character(0), stringsAsFactors = FALSE)

ep_check <- function(id, ok, detail = "") {
  ok <- isTRUE(ok)
  EP_FIRST$rows <- rbind(EP_FIRST$rows,
                         data.frame(check = id, result = if (ok) "PASS" else "FAIL",
                                    detail = detail, stringsAsFactors = FALSE))
  expect_true(ok, info = paste0(id, if (nzchar(detail)) paste0(" [", detail, "]") else ""))
  invisible(ok)
}

ep_report <- function() {
  r <- EP_FIRST$rows
  if (!nrow(r)) return(invisible(NULL))
  w <- max(nchar(r$check))
  cat("\n==== endpoint-first single-dataset pages: PASS/FAIL ====\n")
  for (i in seq_len(nrow(r)))
    cat(sprintf("  [%s] %-*s %s\n", r$result[i], w, r$check[i], r$detail[i]))
  cat(sprintf("  ---- %d PASS / %d FAIL ----\n",
              sum(r$result == "PASS"), sum(r$result == "FAIL")))
  invisible(r)
}

ep_first_offered <- function(ds) unlist(unname(ds), use.names = FALSE)
ep_first_families <- c("OS", "DFS", "PFS", "MFS")

ep_first_token <- function(e, acc, fam) {
  o <- e$.endpoint_opts(e$catalog, acc)
  i <- match(fam, o$family)
  if (is.na(i)) NA_character_ else o$token[i]
}

# ---- the sidebar order: Endpoint, then Cancer type, then Dataset ----------
test_that("the single-dataset sidebars list Endpoint before Cancer type and Dataset", {
  skip_if_not(ep_first_ready(), "shiny app files available")
  env <- ep_first_env()
  for (ui in c("ui_mod_km", "ui_mod_cox", "ui_mod_cox_by_genes")) {
    h <- paste(as.character(env[[ui]]("x")), collapse = "")
    pos <- vapply(c("ep", "ctype", "ds"),
                  function(k) regexpr(sprintf('id="x-%s"', k), h, fixed = TRUE)[1], numeric(1))
    ep_check(sprintf("%s: Endpoint -> Cancer type -> Dataset in the sidebar", ui),
             all(pos > 0) && pos[["ep"]] < pos[["ctype"]] && pos[["ctype"]] < pos[["ds"]],
             paste(names(pos), pos, collapse = " "))
  }
  km <- paste(as.character(env$ui_mod_km("x")), collapse = "")
  ep_check("the Endpoint control says it is chosen first and restricts the dataset list",
           grepl("Chosen first", km, fixed = TRUE) &&
             grepl("cancer type and dataset lists below then offer only the cohorts that carry it",
                   km, fixed = TRUE))
  ep_check("the Dataset control says the list is restricted to carrying cohorts",
           grepl("Only the cohorts that carry the endpoint family chosen above", km, fixed = TRUE))
})

# ---- the family filter behind the three pages -----------------------------
test_that("the endpoint selector offers the families the catalog carries", {
  skip_if_not(ep_first_ready(), "shiny app files available")
  e <- ep_first_env(); cat_d <- e$catalog
  ch <- e$.family_choices(cat_d)
  n <- e$.family_counts(cat_d)
  ep_check("family choices are the catalog families in canonical order",
           identical(unname(ch), e$.family_order), paste(unname(ch), collapse = ","))
  ep_check("each family label carries its cohort count",
           identical(names(ch), sprintf("%s (%d)", unname(ch), n[unname(ch)])),
           paste(names(ch), collapse = " | "))
  ep_check("an empty catalog still yields a usable family selector",
           length(e$.family_choices(e$catalog[0, , drop = FALSE])) >= 1)
})

test_that("only cohorts carrying the chosen family are offered, per family", {
  skip_if_not(ep_first_ready(), "shiny app files available")
  e <- ep_first_env(); cat_d <- e$catalog
  for (fam in ep_first_families) {
    ds <- e$.dataset_choices_with_family(cat_d, "__all__", fam)
    accs <- ep_first_offered(ds)
    carry <- vapply(accs, function(a) fam %in% e$.endpoint_opts(cat_d, a)$family, logical(1))
    ep_check(sprintf("family %s: every offered dataset carries it", fam),
             all(carry),
             sprintf("%d offered, %d without", length(accs), sum(!carry)))
    ep_check(sprintf("family %s: the offer is exactly the carrying cohorts", fam),
             setequal(accs, e$.cohorts_with_family(cat_d, "__all__", fam)) &&
               length(accs) == length(e$.cohorts_with_family(cat_d, "__all__", fam)),
             sprintf("%d offered vs %d carrying", length(accs),
                     length(e$.cohorts_with_family(cat_d, "__all__", fam))))
    ## the choices are still grouped by cancer type, as .dataset_choices() does
    grp_types <- vapply(names(ds), function(g)
      identical(unname(ds[[g]]), accs[as.character(cat_d$Type[match(accs, cat_d$Accession)]) == g]),
      logical(1))
    ep_check(sprintf("family %s: dataset choices stay grouped by cancer type", fam),
             all(grp_types) && setequal(names(ds),
               unique(as.character(cat_d$Type[match(accs, cat_d$Accession)]))),
             paste(names(ds), collapse = ","))
  }
})

test_that("a chosen cancer type narrows the dataset list to its carrying cohorts", {
  skip_if_not(ep_first_ready(), "shiny app files available")
  e <- ep_first_env(); cat_d <- e$catalog
  dfs_lung <- e$.cohorts_with_family(cat_d, "Lung Cancer", "DFS")
  offered <- ep_first_offered(e$.dataset_choices_with_family(cat_d, "Lung Cancer", "DFS"))
  ep_check("DFS + Lung Cancer offers exactly the Lung DFS cohorts",
           setequal(offered, dfs_lung) && length(offered) == length(dfs_lung),
           paste(offered, collapse = ","))
  ep_check("DFS + Lung Cancer offers nothing outside that type",
           all(unique(cat_d$Type[match(offered, cat_d$Accession)]) == "Lung Cancer"))
  ep_check("DFS + Lung Cancer is strictly smaller than all Lung cohorts",
           length(offered) < length(e$.peer_datasets(cat_d, "Lung Cancer")))
  ep_check("All types + DFS offers the DFS cohorts of every type",
           setequal(ep_first_offered(e$.dataset_choices_with_family(cat_d, "__all__", "DFS")),
                    e$.cohorts_with_family(cat_d, "__all__", "DFS")))
})

test_that("cancer types without a cohort carrying the family are not offered", {
  skip_if_not(ep_first_ready(), "shiny app files available")
  e <- ep_first_env(); cat_d <- e$catalog
  all_types <- sort(unique(as.character(cat_d$Type)))
  for (fam in ep_first_families) {
    ch <- e$.type_choices_with_family(cat_d, fam)
    vals <- unname(ch)
    carrying <- all_types[vapply(all_types, function(t)
      length(e$.cohorts_with_family(cat_d, t, fam)) > 0, logical(1))]
    ep_check(sprintf("family %s: All types comes first", fam),
             identical(vals[1], "__all__"))
    ep_check(sprintf("family %s: offers exactly the types with a carrying cohort", fam),
             setequal(setdiff(vals, "__all__"), carrying),
             paste(setdiff(vals, "__all__"), collapse = ","))
    ep_check(sprintf("family %s: each offered type is labelled with its carrying count", fam),
             identical(names(ch)[-1],
                       sprintf("%s (%d)", carrying,
                               vapply(carrying, function(t)
                                 length(e$.cohorts_with_family(cat_d, t, fam)), integer(1)))))
  }
  ## the concrete "not carried" case quoted in the report: no lung cohort has MFS
  mfs <- e$.type_choices_with_family(cat_d, "MFS")
  ep_check("MFS: Lung Cancer is not offered (no lung cohort carries MFS)",
           !("Lung Cancer" %in% unname(mfs)) && "Lung Cancer" %in% all_types,
           paste(unname(mfs), collapse = ","))
})

test_that("the shared selection picks a family the cohort actually carries", {
  skip_if_not(ep_first_ready(), "shiny app files available")
  e <- ep_first_env(); cat_d <- e$catalog
  ep_check("a family carried by the cohort is kept",
           identical(e$.family_for_cohort(cat_d, "GSE31210", "DFS"), "DFS"))
  ep_check("a family the cohort lacks falls back to its primary family",
           identical(e$.family_for_cohort(cat_d, "GSE31210", "MFS"),
                     e$.endpoint_opts(cat_d, "GSE31210")$family[1]) &&
             identical(e$.family_for_cohort(cat_d, "TCGA-LIHC", "MFS"), "OS"),
           paste(e$.family_for_cohort(cat_d, "GSE31210", "MFS"),
                 e$.family_for_cohort(cat_d, "TCGA-LIHC", "MFS")))
  ep_check("no cohort (NULL) yields no family",
           is.na(e$.family_for_cohort(cat_d, NULL, "DFS")))
  ep_check("the concrete token of the selected cohort/family is resolved",
           identical(e$.resolve_family(cat_d, "GSE31210", "DFS"), "RFS") &&
             identical(e$.resolve_family(cat_d, "TCGA-LUAD", "DFS"), "DFI") &&
             identical(e$.resolve_family(cat_d, "GSE14814", "OS"), "OS"),
           sprintf("GSE31210/DFS=%s TCGA-LUAD/DFS=%s GSE14814/OS=%s",
                   e$.resolve_family(cat_d, "GSE31210", "DFS"),
                   e$.resolve_family(cat_d, "TCGA-LUAD", "DFS"),
                   e$.resolve_family(cat_d, "GSE14814", "OS")))
})

# ---- the three module servers --------------------------------------------
# Stub the two update functions so the choices the server sends can be read
# back, and emulate the browser by feeding the values into the mock session.
ep_first_capture <- function(env, cap) {
  env$updateSelectInput <- function(session, inputId, label = NULL, choices = NULL,
                                    selected = NULL, ...) {
    if (identical(inputId, "ep")) cap$ep <- list(choices = choices, selected = selected)
    if (identical(inputId, "ctype")) cap$ctype <- list(choices = choices, selected = selected)
    if (!is.null(selected)) cap$pending[[inputId]] <- selected
    invisible(NULL)
  }
  env$updateSelectizeInput <- function(session, inputId, label = NULL, choices = NULL,
                                       selected = NULL, ...) {
    if (identical(inputId, "ds")) cap$ds <- list(choices = choices, selected = selected)
    if (!is.null(selected)) cap$pending[[inputId]] <- selected
    invisible(NULL)
  }
  ## the probe lookup would hit the mirror; the selection logic does not need it
  env$.ref_probes <- function(acc, gene) c("probe1", "probe2")
  invisible(env)
}

ep_first_flow <- function(session, rv, e, cap, tag) {
  cat_d <- e$catalog
  apply_pending <- function() {
    for (i in 1:6) {                       # the browser applies updates until quiet
      p <- cap$pending; cap$pending <- list()
      if (!length(p)) break
      do.call(session$setInputs, p)
    }
  }
  inp <- function() session$input
  offered <- function() ep_first_offered(cap$ds$choices)
  carries <- function(accs, fam)
    vapply(accs, function(a) fam %in% e$.endpoint_opts(cat_d, a)$family, logical(1))

  session$flushReact(); apply_pending()
  ep_check(sprintf("%s: first render offers the family selector first", tag),
           identical(inp()$ep, "OS") && identical(unname(cap$ep$choices), e$.family_order),
           paste0("ep=", inp()$ep, " choices=", paste(unname(cap$ep$choices), collapse = ",")))
  ep_check(sprintf("%s: first render defaults the cancer type to All types", tag),
           identical(inp()$ctype, "__all__"), paste0("ctype=", inp()$ctype))
  ep_check(sprintf("%s: first render offers only OS cohorts", tag),
           length(offered()) > 0 && all(carries(offered(), "OS")),
           sprintf("%d offered", length(offered())))

  ## user picks DFS: the type list and the dataset list follow
  session$setInputs(ep = "DFS"); apply_pending()
  ep_check(sprintf("%s: DFS -> every offered dataset carries DFS", tag),
           length(offered()) > 0 && all(carries(offered(), "DFS")),
           sprintf("%d offered, %d without", length(offered()), sum(!carries(offered(), "DFS"))))
  ep_check(sprintf("%s: DFS -> type and dataset lists are the DFS ones", tag),
           identical(unname(cap$ctype$choices), unname(e$.type_choices_with_family(cat_d, "DFS"))) &&
             setequal(offered(), e$.cohorts_with_family(cat_d, "__all__", "DFS")))

  ## narrow to Lung Cancer: exactly the lung cohorts that carry DFS
  session$setInputs(ctype = "Lung Cancer"); apply_pending()
  expected <- e$.cohorts_with_family(cat_d, "Lung Cancer", "DFS")
  ep_check(sprintf("%s: DFS + Lung Cancer -> exactly the Lung DFS cohorts", tag),
           setequal(offered(), expected) && length(offered()) == length(expected),
           paste(offered(), collapse = ","))

  ## selecting a cohort sets the shared selection and its token
  session$setInputs(ds = "GSE31210"); apply_pending()
  ep_check(sprintf("%s: selecting GSE31210 at DFS resolves RFS", tag),
           identical(rv$sel$acc, "GSE31210") && identical(rv$sel$family, "DFS") &&
             identical(rv$sel$ep, "RFS"),
           paste(rv$sel$acc, rv$sel$family, rv$sel$ep))

  ## family switch to one this cohort does not carry: never left without the endpoint
  session$setInputs(ep = "MFS"); apply_pending()
  ep_check(sprintf("%s: MFS -> every offered dataset carries MFS", tag),
           length(offered()) > 0 && all(carries(offered(), "MFS")),
           sprintf("%d offered, %d without", length(offered()), sum(!carries(offered(), "MFS"))))
  ep_check(sprintf("%s: MFS -> Lung Cancer is not offered as a type", tag),
           !("Lung Cancer" %in% unname(cap$ctype$choices)),
           paste(unname(cap$ctype$choices), collapse = ","))
  ep_check(sprintf("%s: MFS -> the shared cohort carries MFS and has its token", tag),
           identical(rv$sel$family, "MFS") && rv$sel$acc %in% offered() &&
             identical(rv$sel$ep, ep_first_token(e, rv$sel$acc, "MFS")),
           paste(rv$sel$acc, rv$sel$family, rv$sel$ep))

  ## PFS follows the same rule (OS/PFS/MFS are all required to behave alike)
  session$setInputs(ep = "PFS"); apply_pending()
  ep_check(sprintf("%s: PFS -> every offered dataset carries PFS", tag),
           length(offered()) > 0 && all(carries(offered(), "PFS")) &&
             setequal(offered(), e$.cohorts_with_family(cat_d, "__all__", "PFS")),
           sprintf("%d offered", length(offered())))

  ## a cohort pushed from the shared selection (Datasets row / jump button):
  ## TCGA-LIHC has no MFS, so the endpoint control moves to one of its families
  session$setInputs(ep = "MFS"); apply_pending()
  shiny::isolate(e$.set_shared_selection(rv, cat_d, "TCGA-LIHC"))
  session$flushReact(); apply_pending()
  ep_check(sprintf("%s: incoming TCGA-LIHC (no MFS) moves the endpoint to its family", tag),
           identical(inp()$ep, e$.family_for_cohort(cat_d, "TCGA-LIHC", "MFS")) &&
             identical(inp()$ds, "TCGA-LIHC") && identical(rv$sel$acc, "TCGA-LIHC") &&
             inp()$ep %in% e$.endpoint_opts(cat_d, "TCGA-LIHC")$family,
           paste(inp()$ep, inp()$ds, rv$sel$ep))

  ## ... and when the incoming cohort does carry the family, the family is kept
  session$setInputs(ep = "DFS"); apply_pending()
  shiny::isolate(e$.set_shared_selection(rv, cat_d, "TCGA-LUAD"))
  session$flushReact(); apply_pending()
  ep_check(sprintf("%s: incoming TCGA-LUAD keeps DFS and selects the cohort", tag),
           identical(inp()$ep, "DFS") && identical(inp()$ds, "TCGA-LUAD") &&
             identical(rv$sel$ep, "DFI"),
           paste(inp()$ep, inp()$ds, rv$sel$ep))
}

test_that("the three single-dataset servers are endpoint-first", {
  skip_if_not(ep_first_ready(), "shiny app files available")
  env <- ep_first_env()
  specs <- list(
    list(fun = env$server_mod_km, id = "km"),
    list(fun = env$server_mod_cox, id = "cox"),
    list(fun = env$server_mod_cox_by_genes, id = "genes"))
  for (spec in specs) {
    cap <- new.env(parent = emptyenv()); cap$pending <- list()
    ep_first_capture(env, cap)
    rv <- shiny::reactiveValues(sel = list(acc = NULL, type = NULL, family = NULL, ep = NULL))
    eval(bquote(shiny::testServer(.(spec$fun),
      args = list(id = .(spec$id), rv = rv, dataset_info = CanPAS::dataset_info), {
        ep_first_flow(session, rv, ep_first_env(), cap, .(spec$id))
      })))
  }
})

test_that("endpoint-first PASS/FAIL table", {
  skip_if_not(ep_first_ready(), "shiny app files unavailable")
  expect_true(nrow(EP_FIRST$rows) >= 70)   # every check above really ran
  ep_report()
})
