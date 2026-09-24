# ------------------------------------------------------------------
# CanPAS app - pure analysis helpers (no UI). Sourced by app.R
# ------------------------------------------------------------------

.STD_EP <- c("OS", "DSS", "DFS", "RFS", "PFS", "MFS", "DRFS", "EFS", "DFI", "PFI")

.pkg_data2 <- function(nm) {                 # load package lazy dataset
  e <- new.env(parent = emptyenv())
  utils::data(list = nm, package = "CanPAS", envir = e)
  get(nm, envir = e)
}

# ensure the catalog has the SurvivalTypes column (older data may lack it)
.norm_catalog <- function(di) {
  if (is.null(di)) return(di)
  di <- as.data.frame(di, stringsAsFactors = FALSE)
  if (!"SurvivalTypes" %in% colnames(di)) di$SurvivalTypes <- NA_character_
  # sample-size columns: N = analysable samples (expression + primary endpoint
  # time/status non-missing), n_events = events of the primary endpoint.
  if (!"n_expr" %in% colnames(di))   di$n_expr <- NA_integer_
  if (!"n_surv" %in% colnames(di))   di$n_surv <- NA_integer_
  if (!"n_events" %in% colnames(di)) di$n_events <- NA_integer_
  if (!"expr_in_mirror" %in% colnames(di)) di$expr_in_mirror <- !is.na(di$n_expr)
  # cohort-overlap bookkeeping (part of this release): present so pages can rely on it
  if (!"CohortGroup" %in% colnames(di)) di$CohortGroup <- NA_character_
  if (!"Note" %in% colnames(di)) di$Note <- NA_character_
  di
}

# NA-safe display for count columns
.dash <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "\u2014"
  x
}

# safe comma-joined survival-type token list
.surv_tokens <- function(di) {
  sv <- di$SurvivalTypes
  if (is.null(sv)) return(character(0))
  sv <- as.character(sv)
  tok <- unique(trimws(unlist(strsplit(sv[!is.na(sv) & nzchar(sv)], ","))))
  tok[nzchar(tok) & tok %in% c("OS", "DSS", "DFS", "RFS", "PFS",
                               "MFS", "DRFS", "EFS", "DFI", "PFI")]
}

.has_ep <- function(df, ep) {
  all(c(paste0(ep, "_time"), paste0(ep, "_status")) %in% colnames(df))
}

.surv_tbl <- function(acc) {
  if (startsWith(acc, "TCGA-")) {          # TCGA: local clinical/survival tables
    return(tcga_surv_table(sub("^TCGA-", "", acc)))
  }
  r <- get_data(acc, "surv_data")
  df <- r$response
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
    stop("No survival data returned for cohort ", acc, ".")
  colnames(df)[1] <- "ID"
  df$ID <- as.character(df$ID)
  df
}

.mk_surv <- function(acc, ep) {
  df <- .surv_tbl(acc)
  if (!.has_ep(df, ep))
    stop("Cohort ", acc, " has no ", ep, " survival data (available: ",
         paste(intersect(.STD_EP,
                         unique(sub("_(time|status)$", "",
                                    colnames(df)[grepl("_(time|status)$",
                                                       colnames(df))]))),
               collapse = ", "), "）。")
  out <- data.frame(ID = df$ID, stringsAsFactors = FALSE)
  out[[paste0(ep, "_time")]]   <- suppressWarnings(as.numeric(df[[paste0(ep, "_time")]]))
  out[[paste0(ep, "_status")]] <- suppressWarnings(as.numeric(df[[paste0(ep, "_status")]]))
  out
}

# join clinical columns available in the local surv mirrors (<root>/data/processed/surv)
.loc_clin <- function(acc) {
  root <- Sys.getenv("CPAS_DATA_ROOT", "")
  dir  <- file.path(root, "data", "processed", "surv")
  if (!dir.exists(dir)) {
    if (!isTRUE(getOption("cpas.clin_warned"))) {
      message("Clinical covariates are unavailable: local survival directory not found (",
              dir, "). Set CPAS_DATA_ROOT to the CanPAS data root to enable them.")
      options(cpas.clin_warned = TRUE)
    }
    return(NULL)
  }
  pat <- paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", acc), "(_[^_]+)?_surv\\.rds$")
  fs  <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!length(fs)) return(NULL)
  s <- readRDS(fs[1])
  if (is.null(s) || !nrow(s)) return(NULL)
  keep <- intersect(c("age", "sex", "stage", "T", "N", "M",
                      "grade", "histology"), colnames(s))
  if (!length(keep)) return(NULL)
  cl <- s[, keep, drop = FALSE]
  cl$ID <- rownames(s)
  cl
}

# reference probes (REF_ID) for a gene on the dataset's platform
.ref_probes <- function(acc, gene) {
  if (startsWith(acc, "TCGA-")) return(gene)   # TCGA has no probe concept
  ID_map <- .pkg_data2("ID_map")
  di     <- .pkg_data2("dataset_info")
  row0   <- di[di$Accession == acc, , drop = FALSE]
  if (!nrow(row0)) stop("Dataset not found in the catalog: ", acc)
  gpl <- row0$GPL[1]
  gids <- unique(as.character(ID_map$gene_id[ID_map$Symbol == gene]))
  gids <- gids[!is.na(gids)]
  if (!length(gids)) stop("Gene not found in ID_map: ", gene)
  ref <- get_data(gpl, "gpl", ids = gids)
  if (is.null(ref$response) || !nrow(ref$response))
    stop("No probe for gene ", gene, " on platform ", gpl, ".")
  sort(unique(as.character(ref$response$row_names)))
}

# probe-level expression -> data.frame(ID, marker)
.expr_by_probe <- function(acc, probe) {
  if (startsWith(acc, "TCGA-")) {          # the "probe" is a gene symbol for TCGA
    d <- tcga_get_expr(sub("^TCGA-", "", acc), probe)
    colnames(d)[1] <- "ID"; colnames(d)[2] <- "marker"
    d$marker <- suppressWarnings(as.numeric(d$marker))
    return(d)
  }
  r <- get_data(acc, "expression", ids = probe)
  ex <- r$response
  if (is.null(ex) || !is.data.frame(ex) || !nrow(ex))
    stop("No expression data for probe/gene ", probe, " in this cohort.")
  samples <- setdiff(colnames(ex), colnames(ex)[1])
  vals <- suppressWarnings(as.numeric(unlist(ex[1, samples], use.names = FALSE)))
  data.frame(ID = samples, marker = vals, stringsAsFactors = FALSE)
}

# gene aggregated marker (rule: max/mean/median/min)
.gene_marker <- function(acc, gene, rule = "max") {
  if (startsWith(acc, "TCGA-")) {
    d <- tcga_get_expr(sub("^TCGA-", "", acc), gene)
    colnames(d)[1] <- "ID"; colnames(d)[2] <- "marker"
    d$marker <- suppressWarnings(as.numeric(d$marker))
    return(d)
  }
  e <- get_expr_data(acc, gene, process_duplicates = rule)
  m <- e$expr_data
  if (!gene %in% colnames(m)) stop("Gene ", gene, " is not mapped on this cohort platform.")
  data.frame(ID = m$ID,
             marker = suppressWarnings(as.numeric(m[[gene]])),
             stringsAsFactors = FALSE)
}

# signature (weighted linear) marker
.sig_marker <- function(acc, formula_txt, rule = "max", allow_missing = FALSE) {
  pat <- "[+-]?[[:space:]]*[0-9]*\\.?[0-9]+[[:space:]]*\\*[[:space:]]*[A-Za-z][A-Za-z0-9._-]*"
  terms <- regmatches(formula_txt, gregexpr(pat, formula_txt))[[1]]
  if (!length(terms)) stop("Could not parse the signature formula, e.g. 0.5*GAPDH + 0.5*TNS1")
  clean <- gsub("[[:space:]]", "", terms)
  w <- as.numeric(sub("\\*.*$", "", clean)); g <- sub("^.*\\*", "", clean)
  if (anyNA(w)) stop("Failed to parse weights in: ", formula_txt)
  dup <- g[duplicated(g)]
  if (length(dup)) stop("Duplicated gene in signature: ", paste(unique(dup), collapse = ", "))
  if (startsWith(acc, "TCGA-")) {
    expr <- tcga_get_expr(sub("^TCGA-", "", acc), g)
    colnames(expr)[1] <- "ID"
  } else {
    e <- get_expr_data(acc, g, process_duplicates = rule)
    expr <- e$expr_data
  }
  avail <- intersect(g, colnames(expr))
  miss  <- setdiff(g, avail)
  if (length(miss)) {
    msg <- paste0("Signature gene(s) not measurable: ", paste(miss, collapse = ", "))
    if (!allow_missing) stop(msg, " (enable 'allow missing genes' to score the subset)")
    message(msg, "; scoring the measurable subset.")
  }
  if (!length(avail)) stop("None of the signature genes are measurable.")
  keep <- g %in% avail
  w2 <- w[keep]; g2 <- g[keep]
  m <- as.matrix(expr[, g2, drop = FALSE])
  m[] <- apply(m, 2, function(x) suppressWarnings(as.numeric(x)))
  score <- drop(m %*% w2)
  data.frame(ID = expr$ID, marker = score, stringsAsFactors = FALSE)
}

# unified marker data.frame: survival (ep) + marker + optional clinical
.build_marker_df <- function(acc, ep, kind = c("gene", "signature"),
                             gene = NULL, ref = NULL,
                             formula_txt = NULL,
                             rule = "max", allow_missing = FALSE) {
  kind <- match.arg(kind)
  marker_df <- if (kind == "gene") {
    if (!is.null(ref) && nzchar(ref) && ref != "__all__") {
      .expr_by_probe(acc, ref)
    } else {
      .gene_marker(acc, gene, rule = rule)
    }
  } else {
    .sig_marker(acc, formula_txt, rule = rule, allow_missing = allow_missing)
  }
  surv <- .mk_surv(acc, ep)
  d <- merge(surv, marker_df, by = "ID")
  cl <- .loc_clin(acc)
  if (!is.null(cl)) {
    need <- setdiff(intersect(colnames(cl), colnames(d)), "ID")
    # only add clinical columns missing from the merged table (avoid .x/.y clash)
    cl <- cl[, c("ID", setdiff(colnames(cl), c("ID", need))), drop = FALSE]
  }
  if (!is.null(cl) && ncol(cl) > 1L) d <- merge(d, cl, by = "ID", all.x = TRUE)
  # Rows that the survival model would silently drop (missing endpoint time or
  # status) or that cannot be grouped (missing marker) are removed here, so the
  # reported n / events match what plot_km() and COX_analysis() actually use.
  tc <- paste0(ep, "_time"); sc <- paste0(ep, "_status")
  n_drop_surv <- if (all(c(tc, sc) %in% names(d)))
    sum(is.na(d[[tc]]) | is.na(d[[sc]])) else 0L
  if (n_drop_surv) d <- d[!is.na(d[[tc]]) & !is.na(d[[sc]]), , drop = FALSE]
  n_drop_marker <- sum(is.na(d$marker))
  if (n_drop_marker) d <- d[!is.na(d$marker), , drop = FALSE]
  attr(d, "dropped_na") <- c(surv = n_drop_surv, marker = n_drop_marker)
  d
}

# best cut-point by maximising the log-rank chi-square
.cut_auto <- function(df, probs = seq(0.2, 0.8, by = 0.05), min_grp = 5) {
  ep <- sub("_time$", "",
            colnames(df)[grepl("_time$", colnames(df))][1])
  d <- data.frame(time = as.numeric(df[[paste0(ep, "_time")]]),
                  status = as.numeric(df[[paste0(ep, "_status")]]),
                  marker = as.numeric(df$marker))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 20) stop("Not enough complete observations to search for an optimal cut-point.")
  cand <- unique(stats::quantile(d$marker, probs, na.rm = TRUE, type = 7))
  cand <- cand[is.finite(cand) & cand > min(d$marker) & cand < max(d$marker)]
  best <- NULL; n_cand <- 0L
  for (cp in cand) {
    hi <- sum(d$marker > cp); lo <- sum(d$marker <= cp)
    if (hi < min_grp || lo < min_grp) next
    g <- factor(ifelse(d$marker > cp, "High", "Low"),
                levels = c("Low", "High"))
    sd <- survival::survdiff(survival::Surv(time, status) ~ g, data = d)
    ch <- sd$chisq
    n_cand <- n_cand + 1L
    if (is.null(best) || ch > best$ch)
      best <- list(cut = cp, ch = ch,
                   p = 1 - stats::pchisq(ch, length(sd$n) - 1))
  }
  # the naive p-value is the maximum of n_cand log-rank statistics: it must be
  # reported as unadjusted. When 'maxstat' is installed a search-adjusted
  # p-value is added.
  # Search-adjusted p-value for a maximally selected log-rank statistic.
  # maxstat's approximations are not all reliable (Lau94 can exceed 1, HL can be
  # NaN), so the first method that returns a value inside (naive p, 1] is used;
  # in practice this is Lau92, with exactGauss as fallback.
  adj <- NA_real_; adj_method <- NA_character_
  if (!is.null(best) && requireNamespace("maxstat", quietly = TRUE)) {
    for (pm in c("Lau92", "exactGauss")) {
      ms <- tryCatch(suppressWarnings(
        maxstat::maxstat.test(survival::Surv(time, status) ~ marker, data = d,
                              smethod = "LogRank", pmethod = pm)),
        error = function(e) NULL)
      if (is.null(ms)) next
      pv <- suppressWarnings(as.numeric(ms$p.value))
      if (is.finite(pv) && pv >= best$p && pv <= 1) { adj <- pv; adj_method <- pm; break }
    }
  }
  if (is.null(best))
    list(cut = stats::median(d$marker), ch = NA_real_, p = NA_real_,
         n_cand = n_cand, adjusted_p = NA_real_, adjusted_method = NA_character_)
  else list(cut = best$cut, ch = best$ch, p = best$p, n_cand = n_cand,
            adjusted_p = adj, adjusted_method = adj_method)
}

# ---- session-level fetch cache (avoid repeated queries) ----
.memo <- new.env(parent = emptyenv())
.cached <- function(key, fn) {
  if (exists(key, envir = .memo, inherits = FALSE)) return(get(key, envir = .memo))
  v <- fn()
  assign(key, v, envir = .memo)
  v
}

# ---- endpoints common to several cohorts (standardised set) ----
.common_endpoints <- function(cohorts,
                              prefer = c("OS", "DSS", "DFS", "RFS", "PFS",
                                         "MFS", "DRFS", "EFS", "DFI", "PFI")) {
  if (!length(cohorts)) return(character(0))
  per <- lapply(cohorts, function(acc) {
    tryCatch({
      df <- .surv_tbl(acc)
      unique(sub("_(time|status)$", "",
                 colnames(df)[grepl("_(time|status)$", colnames(df))]))
    }, error = function(e) character(0))
  })
  common <- Reduce(intersect, per)
  c(intersect(prefer, common), setdiff(common, prefer))
}

# ---- clinical covariates common to several cohorts (for adjusted meta) ----
.common_clinical <- function(cohorts) {
  if (!length(cohorts)) return(character(0))
  per <- lapply(cohorts, function(acc) {
    cl <- tryCatch(.loc_clin(acc), error = function(e) NULL)
    if (is.null(cl)) return(character(0))
    setdiff(colnames(cl), "ID")
  })
  Reduce(intersect, per)
}

# ---- shared selection & selector choices (used by Datasets/KM/COX pages) ----
.dataset_choices <- function(di) {
  di <- .norm_catalog(di)
  grp <- split(di$Accession, di$Type)
  lapply(grp, function(x) stats::setNames(x, x))
}

.endpoint_opts <- function(di, acc) {
  di <- .norm_catalog(di)
  out <- data.frame(family = character(0), token = character(0),
                    derived = logical(0), label = character(0),
                    stringsAsFactors = FALSE)
  i <- match(acc, di$Accession)
  if (is.na(i)) return(out)
  for (f in c("OS", "DSS", "DFS", "PFS", "MFS")) {
    col <- paste0("EP_", f)
    if (!col %in% colnames(di)) next
    tok <- di[[col]][i]
    if (is.na(tok) || !nzchar(tok)) next
    dv <- if ("EndpointDerived" %in% colnames(di)) di$EndpointDerived[i] else NA
    deriv <- !is.na(dv) && grepl(tok, dv, fixed = TRUE)
    lab <- if (identical(tok, f)) f else sprintf("%s (%s)", f, tok)
    if (deriv) lab <- paste0(lab, " [derived]")
    out <- rbind(out, data.frame(family = f, token = tok, derived = deriv,
                                 label = lab, stringsAsFactors = FALSE))
  }
  rownames(out) <- NULL
  out
}

# ---- family-aware cohort selection (shared by the Datasets and multi pages) --
# One definition of "this cohort can answer this family", used by the UI filter
# and by the analysis path, so the two can never disagree.
.family_hit <- function(d, fam) {
  if (identical(fam, "__pfsmfs__"))
    return(.family_hit(d, "PFS") | .family_hit(d, "MFS"))
  hit <- .family_mask(d, fam)
  todo <- which(!hit)
  if (length(todo) && "SurvivalTypes" %in% colnames(d)) {
    hit[todo] <- vapply(as.character(d$SurvivalTypes)[todo], function(s) {
      if (is.na(s) || !nzchar(s)) return(FALSE)
      ts <- trimws(strsplit(s, ",")[[1]]); ts <- ts[nzchar(ts)]
      if (!length(ts)) return(FALSE)
      fam %in% vapply(ts, CanPAS::endpoint_family, character(1))
    }, logical(1))
  }
  hit
}

# cohorts of one cancer type that carry a family (sorted accessions)
.cohorts_with_family <- function(catalog, type, fam) {
  if (is.null(type) || !nzchar(type) || is.null(fam) || !nzchar(fam)) return(character(0))
  d <- catalog[catalog$Type == type, , drop = FALSE]
  if (!nrow(d)) return(character(0))
  ## "any" keeps every cohort; every other value (including the combined
  ## "PFS or MFS") is resolved through .family_hit(), which expands the combined
  ## value itself — skipping the filter here once returned the whole catalog
  if (!identical(fam, "__all__"))
    d <- d[.family_hit(d, fam), , drop = FALSE]
  sort(as.character(d$Accession))
}

# families offered for one cancer type, labelled with their cohort count
.ep_choices_for_type <- function(catalog, type) {
  if (is.null(type) || !nzchar(type)) return(stats::setNames("OS", "OS"))
  n <- vapply(.family_order, function(f) length(.cohorts_with_family(catalog, type, f)),
              integer(1))
  fams <- .family_order[n > 0]
  out <- stats::setNames(fams, sprintf("%s (%d)", fams, n[fams]))
  if (all(c("PFS", "MFS") %in% fams))
    out <- c(out, stats::setNames("__pfsmfs__",
              sprintf("PFS or MFS (%d)", n[["PFS"]] + n[["MFS"]])))
  if (!length(out)) out <- stats::setNames("OS", "OS")
  out
}

# the token one cohort contributes to a family ("" when it has none)
.family_token_label <- function(di, acc, fam) {
  o <- .endpoint_opts(di, acc)
  if (!nrow(o)) return("")
  if (identical(fam, "__pfsmfs__")) {
    tok <- o$token[o$family %in% c("PFS", "MFS")]
    return(if (length(tok)) paste(tok, collapse = " or ") else "")
  }
  if (identical(fam, "__all__")) return(paste(o$label, collapse = ", "))
  i <- match(fam, o$family)
  if (is.na(i)) return("")
  o$token[i]
}

# after a family change: keep the cohorts that still qualify and report the
# ones that no longer do, instead of dropping them silently
.refresh_cohorts <- function(catalog, type, fam, selected = character(0), n = 6) {
  ch <- .cohorts_with_family(catalog, type, fam)
  sel <- intersect(selected %||% character(0), ch)
  dropped <- setdiff(selected %||% character(0), ch)
  if (!length(sel)) sel <- head(ch, n)
  list(choices = ch, selected = unique(sel), dropped = dropped)
}


# endpoint selector: label = "DFS (RFS)"; value = FAMILY (resolved at run time)
.endpoint_choices <- function(di, acc) {
  o <- .endpoint_opts(di, acc)
  if (!nrow(o)) return(c(OS = "OS"))
  stats::setNames(o$family, o$label)
}

# pooling families available in EVERY one of the given cohorts
.common_families <- function(di, cohorts) {
  di <- .norm_catalog(di)
  if (!length(cohorts)) return(character(0))
  fams <- c("OS", "DSS", "DFS", "PFS", "MFS")
  ok <- vapply(fams, function(f) {
    col <- paste0("EP_", f)
    if (!col %in% colnames(di)) return(FALSE)
    all(vapply(cohorts, function(a) {
      i <- match(a, di$Accession)
      !is.na(i) && !is.na(di[[col]][i]) && nzchar(di[[col]][i])
    }, logical(1)))
  }, logical(1))
  fams[ok]
}

# concrete endpoint token used by one cohort for a pooling family
.resolve_family <- function(di, acc, family) {
  di <- .norm_catalog(di)
  col <- paste0("EP_", family)
  i <- match(acc, di$Accession)
  if (is.na(i) || !col %in% colnames(di)) return(NA_character_)
  tok <- di[[col]][i]
  if (is.na(tok) || !nzchar(tok)) NA_character_ else tok
}

# "DFS (RFS)" style label
.family_label <- function(family, token) {
  if (is.na(token)) family
  else if (identical(as.character(family), as.character(token))) family
  else sprintf("%s (%s)", family, token)
}

.set_shared_selection <- function(rv, di, acc) {
  di <- .norm_catalog(di)
  o <- .endpoint_opts(di, acc)
  ## Keep the family the user is working in whenever the new cohort supports it;
  ## only fall back to that cohort's first family when it does not. Before this,
  ## selecting a cohort reset the family, so a family chosen on the Datasets page
  ## was silently replaced by the cohort's first family.
  keep <- !is.null(rv$sel$family) && nrow(o) && rv$sel$family %in% o$family
  rv$sel$acc <- acc
  rv$sel$type <- di$Type[di$Accession == acc][1]
  rv$sel$family <- if (keep) rv$sel$family else if (nrow(o)) o$family[1] else "OS"
  row <- o[o$family == rv$sel$family, , drop = FALSE]
  rv$sel$ep <- if (nrow(row)) row$token[1] else NA_character_
  invisible(rv)
}

# "RFS (DFS family)" style status label
.endpoint_status_label <- function(df_family, token) {
  if (is.null(df_family) || is.na(df_family)) return(as.character(token))
  if (is.null(token) || is.na(token)) return(as.character(df_family))
  sprintf("%s (%s family)", token, df_family)
}

# Source of a cohort: GEO / EMBL-EBI / CGGA / TCGA / cBioPortal-hosted --------
# The accession prefix decides the bucket. The catalog holds FIVE groups:
# GSE* (143), E-* (14: ArrayExpress/BioStudies deposits such as E-MTAB-1727),
# TCGA-* (31), CGGA* (3) and two studies that are none of these because their
# clinical (and expression) files were deposited with the study on cBioPortal,
# not as a GEO series (A5-PCPG, IMmotion150). Folding the EMBL-EBI cohorts into
# GEO, or the cBioPortal pair into either, would misreport both the count and
# the source link, so each is its own bucket.
.cohort_source <- function(acc) {
  acc <- as.character(acc)
  ifelse(startsWith(acc, "TCGA-"), "TCGA",
         ifelse(startsWith(acc, "CGGA"), "CGGA",
                ifelse(startsWith(acc, "GSE"), "GEO",
                       ifelse(startsWith(acc, "E-"), "EMBL-EBI", "cBioPortal"))))
}

# EMBL-EBI cohorts live in ArrayExpress/BioStudies; the accession is enough to
# build the BioStudies landing page (E-MTAB-1727 -> S-BSST or the ArrayExpress
# experiment page, which redirects).
.cohort_link_embl <- function(acc) {
  sprintf('<a href="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/%s" target="_blank">ArrayExpress</a>',
          as.character(acc))
}

# cBioPortal study URL of a cohort that is neither GEO, TCGA nor CGGA. The study
# id is read from the catalog's Note when it quotes one ("cBioPortal <study>"),
# with the cBioPortal data page as the fallback.
.cohort_link_cbioportal <- function(acc) {
  acc <- as.character(acc)
  di <- tryCatch(.pkg_data2("dataset_info"), error = function(e) NULL)
  note <- if (!is.null(di) && all(c("Accession", "Note") %in% colnames(di)))
    as.character(di$Note[match(acc, as.character(di$Accession))]) else rep(NA_character_, length(acc))
  note[is.na(note)] <- ""
  id <- sub("^.*cBioPortal[[:space:]]+([A-Za-z0-9_.-]+).*$", "\\1", note)
  ok <- nzchar(id) & id != note
  ifelse(ok,
         sprintf('<a href="https://www.cbioportal.org/study/summary?id=%s" target="_blank">cBioPortal</a>', id),
         '<a href="https://www.cbioportal.org/datasets" target="_blank">cBioPortal</a>')
}

.cohort_link <- function(acc) {
  acc <- as.character(acc)
  src <- .cohort_source(acc)
  ifelse(src == "TCGA",
         sprintf('<a href="https://portal.gdc.cancer.gov/projects/%s" target="_blank">GDC</a>', acc),
         ifelse(src == "CGGA",
                sprintf('<a href="http://www.cgga.org.cn/" target="_blank">CGGA</a>'),
                ifelse(src == "GEO",
                       sprintf('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" target="_blank">NCBI</a>', acc),
                       ifelse(src == "EMBL-EBI",
                              .cohort_link_embl(acc),
                              .cohort_link_cbioportal(acc)))))
}

# Families a cohort provides, in the canonical family order --------------
.family_order <- c("OS", "DSS", "DFS", "PFS", "MFS")

.families_of <- function(di, acc) {
  o <- .endpoint_opts(di, acc)
  o$family
}

# One-line family description used in tables (labels as in the selectors):
#   "DFS (RFS), MFS"  /  "DFS (DFI) [derived]"
.family_display <- function(di, acc) {
  o <- .endpoint_opts(di, acc)
  if (!nrow(o)) return("")
  paste(o$label, collapse = ", ")
}

# Catalog-wide family availability with cohort counts --------------------
.family_counts <- function(di) {
  di <- .norm_catalog(di)
  out <- integer(0)
  for (f in .family_order) {
    col <- paste0("EP_", f)
    if (!col %in% colnames(di)) { out[[f]] <- 0L; next }
    v <- as.character(di[[col]])
    out[[f]] <- sum(!is.na(v) & nzchar(v))
  }
  out
}

# true for rows that provide the family, using the EP_* columns when present
.family_mask <- function(d, fam) {
  col <- paste0("EP_", fam)
  if (col %in% colnames(d)) {
    v <- as.character(d[[col]])
    !is.na(v) & nzchar(v)
  } else {
    ef <- if ("EndpointFamilies" %in% colnames(d)) as.character(d$EndpointFamilies) else rep(NA_character_, nrow(d))
    !is.na(ef) & grepl(fam, ef, fixed = TRUE)
  }
}

# ---- presentation helpers: every number is shown with 4 decimals ----------
.fmt4 <- function(x, digits = 4) {
  if (is.null(x) || !length(x)) return(character(0))
  ifelse(is.na(suppressWarnings(as.numeric(x))), "",
         formatC(suppressWarnings(as.numeric(x)), format = "f", digits = digits))
}

# p-values keep 4 decimals, but a value below the display resolution is shown as
# "<0.0001" instead of a misleading 0.0000
.fmtp4 <- function(p, digits = 4) {
  if (is.null(p) || !length(p)) return(character(0))
  p <- suppressWarnings(as.numeric(p))
  thr <- 10^(-digits)
  ifelse(is.na(p), "",
         ifelse(p < thr, paste0("<", formatC(thr, format = "f", digits = digits)),
                formatC(p, format = "f", digits = digits)))
}

# numeric columns to 4 decimals, p-value columns to 4-decimal text
.dtx4 <- function(d, digits = 4, pageLength = 12,
                  p_pattern = "^(p|P|Pvalue|P_adj|p_het|p_heterogeneity|p_stouffer|P_adj_text)$",
                  skip = c("N", "n", "events", "k", "df", "total_n", "total_events")) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  if (!nrow(d)) return(DT::datatable(d, rownames = FALSE,
                                     options = list(dom = "t", pageLength = pageLength)))
  pcols <- intersect(grep(p_pattern, colnames(d), value = TRUE), colnames(d))
  for (cn in pcols) d[[cn]] <- .fmtp4(d[[cn]], digits)
  num <- setdiff(colnames(d)[vapply(d, is.numeric, logical(1))], c(pcols, skip))
  dt <- DT::datatable(.round4(d), rownames = FALSE, filter = "top",
                      options = list(pageLength = pageLength, scrollX = TRUE, dom = "ftp"))
  if (length(num)) dt <- DT::formatRound(dt, columns = num, digits = digits)
  dt
}

.round4 <- function(d, digits = 4) {
  for (cn in colnames(d)) if (is.numeric(d[[cn]])) d[[cn]] <- round(d[[cn]], digits)
  d
}

# ---- shared retrieval for the multi-dataset pages -------------------------
# Builds, per dataset, the merged data frame carrying a `marker` column, with a
# session-scoped cache. Returns list(merged, tokens, failed, family).
.multi_bundle <- function(session, catalog, accs, family, spec, min_n = 2) {
  if (length(accs) < min_n)
    stop(sprintf("Please select at least %d dataset(s).", min_n))
  tokens <- stats::setNames(vapply(accs, function(a)
    .resolve_family(catalog, a, family), character(1)), accs)
  usable <- accs[!is.na(tokens[accs])]
  errs <- stats::setNames(
    rep(sprintf("no %s endpoint", family), length(setdiff(accs, usable))),
    setdiff(accs, usable))
  merged <- lapply(usable, function(acc) {
    tok <- tokens[[acc]]
    key <- paste(session$token, acc, tok, spec$kind, spec$gene %||% "",
                 spec$sig %||% "", spec$rule, spec$allow, sep = "|")
    tryCatch(.cached(key, function() .build_marker_df(
      acc = acc, ep = tok, kind = spec$kind, gene = spec$gene, ref = NULL,
      formula_txt = spec$sig, rule = spec$rule, allow_missing = spec$allow)),
      error = function(e) e)
  })
  names(merged) <- usable
  bad <- vapply(merged, inherits, logical(1), "error")
  if (any(bad))
    errs <- c(errs, vapply(merged[bad], function(e) conditionMessage(e), character(1)))
  merged <- merged[!bad]
  if (length(merged) < min_n)
    stop(sprintf("Fewer than %d dataset(s) could be retrieved.", min_n),
         if (length(errs)) paste0(" Failures: ",
           paste(names(errs), errs, sep = ": ", collapse = "; ")) else "")
  list(merged = merged, tokens = tokens[names(merged)], failed = errs, family = family)
}

# marker specification shared by the multi-dataset pages
.marker_spec <- function(input) {
  kind <- input$mkind %||% "gene"
  list(kind = kind,
       gene = if (kind == "gene") (input$mgene %||% "GAPDH") else NULL,
       sig  = if (kind == "signature") (input$msig %||% "0.5*GAPDH + 0.5*TNS1") else NULL,
       rule = if (kind == "gene") (input$mrule_g %||% "max") else (input$mrule %||% "max"),
       allow = isTRUE(input$mallow))
}

# endpoint choices labelled with the token each cohort contributes
.endpoint_labels <- function(catalog, accs, fams) {
  vapply(fams, function(f) {
    toks <- unique(vapply(accs, function(a) .resolve_family(catalog, a, f), character(1)))
    sprintf("%s (%s)", f, paste(stats::na.omit(toks), collapse = ", "))
  }, character(1))
}


# ---- figure block: size controls under the figure, downloads at that size ----
# Every analysis page renders its figure through these three helpers, so that
#   (1) the width/height controls sit directly under the figure,
#   (2) changing them resizes the figure in the app immediately (debounced),
#   (3) the downloaded file has exactly the displayed size: PNG in pixels,
#       PDF at the same physical size (96 px per inch).
.fig_dims <- function(input, key, w0, h0) {
  w <- suppressWarnings(as.numeric(input[[paste0(key, "_width")]]))
  h <- suppressWarnings(as.numeric(input[[paste0(key, "_height")]]))
  # the inputs do not exist yet on the first render, so length-0 values are
  # normal here and must fall back to the defaults
  if (length(w) != 1L || !is.finite(w) || w < 200) w <- w0
  if (length(h) != 1L || !is.finite(h) || h < 150) h <- h0
  list(w = as.integer(round(w)), h = as.integer(round(h)))
}

ui_fig_block <- function(ns, key, w = 600, h = 520) {
  tagList(
    uiOutput(ns(paste0(key, "_box"))),
    br(),
    fluidRow(
      bs4Dash::column(3, numericInput(ns(paste0(key, "_width")), "Figure width (px)",
                                      value = w, min = 200, max = 3000, step = 50)),
      bs4Dash::column(3, numericInput(ns(paste0(key, "_height")), "Figure height (px)",
                                      value = h, min = 150, max = 3000, step = 20)),
      bs4Dash::column(3, br(), downloadButton(ns(paste0(key, "_png")),
                                              "Download PNG", class = "btn-success")),
      bs4Dash::column(3, br(), downloadButton(ns(paste0(key, "_pdf")),
                                              "Download PDF", class = "btn-success"))
    ),
    p(class = "note", "The figure resizes as you change these values, and the
      download keeps that size: PNG at exactly this pixel size, PDF at the same
      physical size (96 px per inch).")
  )
}

# draw(): draws the figure to the current device (used for the on-screen plot
#          and for both downloads, so the three can never diverge)
# ready(): TRUE when there is a completed analysis to download
server_fig_block <- function(input, output, session, key, draw, ready = NULL,
                             w0 = 600, h0 = 520, filename = NULL) {
  ns <- session$ns
  # no debounce: the analysis result is already cached in the calling module, so
  # re-drawing on every change is cheap and the figure follows the controls live
  dims <- shiny::reactive(.fig_dims(input, key, w0, h0))

  output[[key]] <- shiny::renderPlot({ draw() })

  output[[paste0(key, "_box")]] <- shiny::renderUI({
    d <- dims()
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns(key), width = paste0(d$w, "px"), height = paste0(d$h, "px")))
  })

  ok <- function() {
    if (is.null(ready)) return(TRUE)
    isTRUE(tryCatch(ready(), error = function(e) FALSE))
  }
  base_name <- function() {
    b <- tryCatch(if (is.null(filename)) paste0("figure_", key) else filename(),
                  error = function(e) paste0("figure_", key))
    if (is.null(b) || !nzchar(b)) paste0("figure_", key) else b
  }
  mk <- function(ext, open) shiny::downloadHandler(
    filename = function() paste0(base_name(), ext),
    content = function(file) {
      if (!ok()) {
        shiny::showNotification(
          "No completed analysis to download - run the analysis first.",
          type = "warning")
        return(invisible(NULL))
      }
      d <- dims()
      open(file, d)
      on.exit(grDevices::dev.off(), add = TRUE)   # never leak a graphics device
      draw()
    })
  output[[paste0(key, "_png")]] <- mk(".png", function(file, d)
    grDevices::png(file, width = d$w, height = d$h, units = "px", res = 96))
  output[[paste0(key, "_pdf")]] <- mk(".pdf", function(file, d)
    grDevices::pdf(file, width = d$w / 96, height = d$h / 96))
  invisible(dims)
}

# ---- cohort overlap warning ------------------------------------------------
# Cohorts in the same CohortGroup share patients (same study on two platforms,
# or the same series deposited twice). Pooling them would count patients more
# than once, so the multi-dataset pages warn before the analysis is run.
.overlap_groups <- function(catalog, accs) {
  if (!"CohortGroup" %in% colnames(catalog)) return(character(0))
  accs <- unique(as.character(accs))
  g <- as.character(catalog$CohortGroup[match(accs, as.character(catalog$Accession))])
  g <- g[!is.na(g)]
  unique(g[duplicated(g)])
}

.overlap_warning <- function(catalog, accs) {
  bad <- .overlap_groups(catalog, accs)
  if (!length(bad)) return(NULL)
  if (!"CohortGroup" %in% colnames(catalog)) return(NULL)
  accs <- unique(as.character(accs))
  det <- vapply(bad, function(g) {
    memb <- accs[!is.na(match(accs, as.character(catalog$Accession))) &
                   as.character(catalog$CohortGroup[match(accs, as.character(catalog$Accession))]) == g]
    sprintf("%s: %s", g, paste(memb, collapse = ", "))
  }, character(1))
  shinydashboard::box(
    status = "warning", width = 12, solidHeader = FALSE,
    title = "Overlapping cohorts selected",
    p(tags$b("These datasets share patients "),
      "(same study on another platform, or the same series deposited twice). ",
      "Pooling them counts those patients more than once, which narrows the ",
      "confidence interval and biases I², Q and the prediction interval. ",
      "Keep one cohort per group, or report the overlap explicitly."),
    tags$ul(lapply(det, tags$li)))
}

# ---- cancer type -> datasets --------------------------------------------
# Only types with at least two datasets can be analysed across datasets (a
# single-dataset question belongs to the Single dataset pages).
.multi_types <- function(catalog) {
  tab <- table(catalog$Type)
  sort(names(tab)[tab >= 2])
}

.peer_datasets <- function(catalog, type) {
  if (is.null(type) || !nzchar(type)) return(character(0))
  sort(as.character(catalog$Accession[catalog$Type == type]))
}

# Fill "cancer type" and "datasets" selectors; keeps the shared cohort selected.
# Fill "cancer type" -> "endpoint family" -> "datasets" on the multi-dataset
# pages. The family comes first, so the cohort list only offers cohorts that can
# answer the chosen endpoint; the shared cohort is preselected when it qualifies.
.init_type_family_cohorts <- function(session, catalog, shared = NULL, n = 6) {
  tp <- .multi_types(catalog)
  if (!length(tp)) return(invisible(NULL))
  st <- if (!is.null(shared)) catalog$Type[catalog$Accession == shared][1] else NA_character_
  sel_type <- if (!is.na(st) && st %in% tp) st else tp[1]
  updateSelectInput(session, "mtype", choices = tp, selected = sel_type)
  ch_ep <- .ep_choices_for_type(catalog, sel_type)
  shared_fam <- if (!is.null(shared) && !is.na(shared))
    .endpoint_opts(catalog, shared)$family else character(0)
  sel_fam <- intersect(shared_fam, unname(ch_ep))[1]
  if (length(sel_fam) != 1L || is.na(sel_fam)) sel_fam <- unname(ch_ep)[1]
  updateSelectInput(session, "mep", choices = ch_ep, selected = sel_fam)
  r <- .refresh_cohorts(catalog, sel_type, sel_fam, selected = shared, n = n)
  updateSelectizeInput(session, "mcohorts", choices = r$choices, selected = r$selected)
  invisible(r)
}

# ---- cancer type -> dataset (single-dataset pages) ----------------------
# Every type is offered here (a single-dataset question may target a type with
# only one cohort), with the number of datasets in the label.
.single_type_choices <- function(catalog) {
  tab <- table(catalog$Type)
  nm <- sort(names(tab))
  stats::setNames(nm, sprintf("%s (%d)", nm, as.integer(tab[nm])))
}

.dataset_choices_of <- function(catalog, type) {
  if (is.null(type) || !nzchar(type)) return(.dataset_choices(catalog))
  d <- catalog[catalog$Type == type, , drop = FALSE]
  .dataset_choices(d)
}

# Fill "cancer type" + "dataset" on a single-dataset page, keeping the shared
# cohort selected when it belongs to the chosen type.
.init_single_selectors <- function(session, catalog, shared = NULL) {
  ch <- .single_type_choices(catalog)
  if (!length(ch)) return(invisible(NULL))
  st <- if (!is.null(shared)) catalog$Type[catalog$Accession == shared][1] else NA_character_
  sel_type <- if (!is.na(st) && st %in% unname(ch)) st else unname(ch)[1]
  updateSelectInput(session, "ctype", choices = ch, selected = sel_type)
  ds <- .dataset_choices_of(catalog, sel_type)
  accs <- unlist(unname(ds), use.names = FALSE)
  updateSelectizeInput(session, "ds", choices = ds,
                       selected = if (!is.null(shared) && shared %in% accs) shared else accs[1],
                       server = TRUE)
  invisible(sel_type)
}

# ---- help document ------------------------------------------------------
# The app ships HELP.md. It is rendered as HTML when a markdown renderer is
# available (commonmark comes with shiny, markdown is the alternative) and
# falls back to escaped pre-formatted text, so the Help page always works.
.help_md_path <- function() {
  cand <- c(file.path(getwd(), "HELP.md"),
            system.file("shiny", "CanPAS", "HELP.md", package = "CanPAS"))
  cand <- cand[file.exists(cand)]
  if (!length(cand)) return(NULL)
  cand[1]
}

.md_to_html <- function(path) {
  if (is.null(path) || !file.exists(path))
    return(tags$p("Help document not found in this installation."))
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (requireNamespace("commonmark", quietly = TRUE))
    return(HTML(commonmark::markdown_html(txt, extensions = TRUE)))
  if (requireNamespace("markdown", quietly = TRUE))
    return(shiny::includeMarkdown(path))
  esc <- gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", txt)))
  HTML(paste0("<pre class='cpas-help-pre'>", esc, "</pre>"))
}

# ---- short labels with hover explanations --------------------------------
# Sidebar labels stay short; the sentence that used to sit in the label moves
# into a hover popup next to it. The popup is the browser's own title tooltip,
# so it is never clipped by the sidebar's scroll container and needs no JS.
.tip <- function(text) {
  tags$span(class = "cpas-tip", title = text,
            style = "cursor:help; color:#0C4B72; margin-left:4px;",
            shiny::icon("info-circle"))
}

# A label for a form control: the text plus its hover explanation.
.lab <- function(label, text = NULL) {
  if (is.null(text) || !nzchar(text)) return(label)
  tagList(label, .tip(text))
}

# What the marker actually uses, given the REF_ID / probe box:
#  * "__all__" (or empty) -> every probe of the gene, collapsed per sample by
#    the chosen rule (max in gene mode);
#  * a specific probe id  -> that probe alone, no collapsing;
#  * a TCGA cohort        -> no probes exist, the gene-level log2 value is used.
.ref_hint <- function(acc, gene, ref, probe_count = NA_integer_) {
  if (is.null(acc) || !nzchar(acc %||% "")) return(NULL)
  if (!nzchar(gene %||% "")) return(NULL)
  if (startsWith(acc, "TCGA-"))
    return(sprintf(paste0("TCGA cohort: no probe concept - the gene-level log2 value of %s ",
                          "(UCSC Xena) is used, so nothing is collapsed."), gene))
  if (is.null(ref) || !nzchar(ref) || identical(ref, "__all__"))
    return(sprintf(paste0("Auto: all probes%s of %s are collapsed per sample by the ",
                          "maximum (rule = max). Pick a single probe below to use it as it is."),
                   if (!is.na(probe_count) && probe_count > 0)
                     sprintf(" (%d on this platform)", probe_count) else "", gene))
  sprintf("Probe %s only: its value is used as it is, without collapsing.", ref)
}
