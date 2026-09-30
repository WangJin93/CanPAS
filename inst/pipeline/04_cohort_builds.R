# ===========================================================================
# CanPAS curation pipeline -- 04_cohort_builds
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 95_build_suppl_expansion.R
#   2. 96_build_suppl_expansion.R
#   3. 100_build_relaxed_gate.R
#   4. 104_build_embl_cohorts.R
#   5. 106_build_embl_step2.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "04_cohort_builds.R", package = "CanPAS"))
#   run_95_build_suppl_expansion()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================

# --- shared helper (glue) ---------------------------------------------------
# 95_build_suppl_expansion.R calls scale_verdict(), which is only defined by
# 96_build_suppl_expansion.R (byte-identical copies also exist in 100/104/106).
# The definition below is copied verbatim from 96 and exposed at file level so
# that run_95_build_suppl_expansion() resolves it; each runner keeps its own
# byte-identical local copy of the definition inside its verbatim block.
# ----------------------------------------------------------------------------
scale_verdict <- function(m) {
  v <- as.matrix(m); storage.mode(v) <- "double"; v <- v[is.finite(v)]
  if (!length(v)) return("empty")
  if (min(v) < 0) return("already_log")
  if (median(v) > 0 && max(v) / median(v) <= 50) return("already_log")
  if (median(v) <= 0 && max(v) <= 25) return("already_log")
  "linear_needs_log2"
}


# ---------------------------------------------------------------------------
# run_95_build_suppl_expansion()  <-  verbatim pipeline/R/95_build_suppl_expansion.R
# ---------------------------------------------------------------------------
run_95_build_suppl_expansion <- function() {
# 95_build_suppl_expansion.R ------------------------------------------------
# 补充表批次（supplement / publication-table route）：为 catalog 仍缺口的癌种建队列。
#
# 本批次的两个队列，生存数据**不在任何 GEO series matrix** 里，而来自
# 「配套发表的临床表 / cBioPortal 策展的 publication table」：
#
#   IMmotion150 (Kidney Cancer, PFS)
#     McDermott et al. Nat Med 2018 (PMID 29867230) —— IMmotion150 随机 II 期试验
#     (atezolizumab ± bevacizumab vs sunitinib, 263 例初治转移性 ccRCC)。
#     表达 = CRI iAtlas 协调重处理的 bulk RNA-seq (gene-level TPM)；
#     生存 = 同一研究配套的患者级临床表 (PFS_MONTHS / PFS_STATUS)。
#     cBioPortal study id = rcc_iatlas_immotion150_2018（非 TCGA）。
#
#   A5-PCPG (Pheochromocytoma / Paraganglioma, OS)
#     A5 Consortium, Nature 2025 (PMID 40097403) —— SDHB 胚系突变 PPGL 队列，
#     94 个肿瘤 (79 例患者) WGS + 91 个肿瘤 WTS。
#     表达 = 同研究 WTS 的 gene-level CPM；
#     生存 = 配套患者级临床表 (OS_MONTHS / OS_STATUS)。
#     cBioPortal study id = hnsc_a5consortium_2025（注意其 studyId 前缀是 hnsc，
#     但 cancer type / oncotree 全是 PGNG=Paraganglioma；非 TCGA）。
#
# 与既有 pipeline 的接口（与 90_build_gse108474_suppl.R / 92 / 93 同构）：
#   data/expr/<ACC>.rds              ID_REF + 每样本一列（gene symbol 行，未折叠）
#   data/pheno/<ACC>.rds             样本级临床表（行名 = 样本 id）
#   data/processed/surv/<ACC>_surv.rds
#                                    <TOK>_status(0/1) + <TOK>_time(年) + 临床列
#   data/processed/gpl/RNAseq_gene_PLAT.rds
#                                    合成平台表：gene symbol -> ENTREZ（两个队列共用；
#                                    与 CGGA_*_PLAT 一样，命名以 _PLAT 结尾，
#                                    16_verify 的孤儿表检查按此约定排除）
#
# 关键口径：
#   * 样本去重 —— 表达表保留全部样本（IMmotion150 263 / A5 91），surv 表按
#     **患者**去重（A5 91 样本 -> 77 患者，保留优先 primary 的样本），与
#     GSE32918 (249 arrays -> 172 patients) 的处理一致。
#   * 时间单位 —— 原始为月，落库一律换算为**年**（与 03_surv_table.R 约定一致）。
#   * log2 —— 两个矩阵都是**线性**尺度（TPM / CPM：median≈0、max 很大），
#     按 18_log2_transform.R 的判定规则（D2/D3/V 投票）先判后再套 log2(x+1)；
#     本脚本**不**运行 18 的 --apply（避免二次转化既有队列），只对新文件套用，
#     并在 pipeline/out/suppl_expansion_log2.csv 记录转化前后统计。
#   * 本脚本**不**写 data/dataset_info.csv / .rda（由 96 负责）。
#
# 用法: Rscript pipeline/R/95_build_suppl_expansion.R [ROOT]
# ----------------------------------------------------------------------------
# --- explicit output root (A10) ------------------------------------------------
# Derived outputs are written under an explicit, env-overridable root
# (CPAS_OUT_ROOT, default <CPAS_DATA_ROOT>/pipeline/out).  See pipeline/OUTPUT_LAYOUT.md.
# Uses the shared helper when it is present and falls back to identical local
# definitions otherwise (so a consolidated / installed copy is self-contained).
# With the default root every path below resolves exactly where it did before:
# this makes the location explicit and overridable, it changes no computation.
.cpas_helper <- file.path(Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS"),
                          "pipeline/R/00_output_root.R")
if (file.exists(.cpas_helper)) source(.cpas_helper)
if (!exists("cpas_out", mode = "function")) {
  .cpas_root <- function() {
    r <- Sys.getenv("CPAS_DATA_ROOT", unset = "")
    if (!nzchar(r)) r <- "/home/Jingle/data/Project/CPAS"
    path.expand(r)
  }
  cpas_root <- .cpas_root
  cpas_out_root <- function() {
    o <- Sys.getenv("CPAS_OUT_ROOT", unset = "")
    if (nzchar(o)) path.expand(o) else file.path(.cpas_root(), "pipeline", "out")
  }
  cpas_out <- function(...) file.path(cpas_out_root(), ...)
  cpas_data <- function(...) file.path(.cpas_root(), "data", ...)
  cpas_suppl <- function(...) file.path(.cpas_root(), "data", "suppl", ...)
}

suppressMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)

SUPPL <- file.path(ROOT, "data/suppl")
OUT   <- cpas_out_root()
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GPL_NAME <- "RNAseq_gene_PLAT"

# ---- cBioPortal 临床表解析（前 4 行以 '#' 开头的元数据，随后一行短名表头）----
read_cbio_clin <- function(path) {
  x <- data.table::fread(path, sep = "\t", header = FALSE, quote = "",
                         na.strings = c("", "NA"), data.table = FALSE,
                         colClasses = "character")
  keep <- !startsWith(x[[1]], "#")
  x <- x[keep, , drop = FALSE]
  hdr <- as.character(x[1, ])
  d <- x[-1, , drop = FALSE]
  colnames(d) <- hdr
  rownames(d) <- NULL
  d[] <- lapply(d, function(z) { z <- trimws(z); z[z %in% c("", "NA")] <- NA; z })
  d
}

# ---- 尺度判定（逐字复刻 18_log2_transform.R 的 survey()/classify()）-------
#   D1 min<0 / D2 median>0 & max/median<=50 / D3 median<=0 & max<=25 -> already_log
#   V1..V4 四条独立证据，>=3 票 -> already_log，否则 needs_log2
classify_scale <- function(v) {
  v <- v[is.finite(v)]
  nz <- v[v != 0]
  q <- quantile(v, c(.25, .5, .75, .99), na.rm = TRUE, names = FALSE)
  d <- list(min = min(v), q25 = q[1], median = q[2], q75 = q[3], p99 = q[4], max = max(v),
            pct_zero = round(100 * mean(v == 0), 1),
            pct_neg  = round(100 * mean(v < 0), 2),
            pct_nonint_nonzero = round(100 * (if (length(nz)) mean(abs(nz - round(nz)) >= 1e-9) else NA_real_), 3))
  ratio <- if (is.finite(d$median) && d$median > 0) d$max / d$median else Inf
  votes <- c(V1 = isTRUE(d$min >= 0 && d$min <= 1),
             V2 = isTRUE(!is.na(d$pct_nonint_nonzero) && d$pct_nonint_nonzero >= 90),
             V3 = isTRUE(d$p99 <= 20),
             V4 = isTRUE(d$max <= 25))
  nv <- sum(votes)
  if (d$min < 0) return(list(verdict = "already_log", evidence = "D1", votes = nv, stats = d))
  if (d$median > 0 && ratio <= 50) return(list(verdict = "already_log", evidence = "D2", votes = nv, stats = d))
  if (d$median <= 0 && d$max <= 25) return(list(verdict = "already_log", evidence = "D3", votes = nv, stats = d))
  if (nv >= 3) return(list(verdict = "already_log", evidence = sprintf("V%d/4", nv), votes = nv, stats = d))
  list(verdict = "needs_log2", evidence = sprintf("V%d/4", nv), votes = nv, stats = d)
}

# ---- 读表达矩阵（第一列 gene symbol），折叠重复 symbol（取均值）-----------
load_expr <- function(path, value_label) {
  dt <- data.table::fread(path, sep = "\t", header = TRUE, quote = "",
                          na.strings = c("", "NA"), data.table = TRUE)
  sym <- as.character(dt[[1]])
  m  <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(m) <- "double"
  colnames(m) <- colnames(dt)[-1]
  n_dup <- sum(duplicated(sym))
  if (n_dup) {
    # gene-level 矩阵：重复 symbol 不是重复探针，按均值折叠，保证 ID_REF 唯一
    # （镜像表 row_names 必须唯一；既有 92 的 GSE248835 是保留后缀的做法，
    #   那是因为它的行标识本身带探针后缀，而这里是纯 gene symbol）
    ok <- !is.na(sym) & nzchar(sym)
    m <- m[ok, , drop = FALSE]; sym <- sym[ok]
    rs  <- rowsum(m, group = sym)                 # 重复 symbol 取均值
    cnt <- as.vector(table(sym)[rownames(rs)])
    m   <- rs / cnt
    rownames(m) <- rownames(rs)
  } else {
    rownames(m) <- sym
  }
  message(sprintf("  %s: %d x %d (duplicated symbols collapsed: %d)",
                  basename(path), nrow(m), ncol(m), n_dup))
  list(mat = m, n_dup = n_dup, label = value_label)
}

# ---- 探针/基因 -> ENTREZ（与 92_build_geo_expansion_expr_pheno.R 同构）----
symbol2entrez <- function(sym) {
  sym <- unique(as.character(sym))
  out <- setNames(rep(NA_character_, length(sym)), sym)
  db <- tryCatch(getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"), error = function(e) NULL)
  if (is.null(db)) return(out)
  m <- tryCatch(suppressWarnings(AnnotationDbi::select(db, keys = sym, columns = "ENTREZID",
                                                      keytype = "SYMBOL")),
                error = function(e) NULL)
  if (is.null(m)) return(out)
  m <- m[!is.na(m$ENTREZID), , drop = FALSE]
  if (nrow(m)) {
    agg <- tapply(m$ENTREZID, m$SYMBOL, function(z) paste(sort(unique(z)), collapse = " /// "))
    out[names(agg)] <- as.character(agg)
  }
  out
}

manifest <- list()
log_rows <- list()
summ_rows <- list()

# ===========================================================================
# 队列 1: IMmotion150 —— Kidney Cancer, PFS
# ===========================================================================
message("== IMmotion150 (Kidney Cancer, PFS)")
acc <- "IMmotion150"
clin_p <- read_cbio_clin(file.path(SUPPL, "IMmotion150_data_clinical_patient.txt.gz"))
clin_s <- read_cbio_clin(file.path(SUPPL, "IMmotion150_data_clinical_sample.txt.gz"))
stopifnot("PATIENT_ID" %in% colnames(clin_p), "SAMPLE_ID" %in% colnames(clin_s))

ex <- load_expr(file.path(SUPPL, "IMmotion150_data_mrna_seq_tpm.txt.gz"), "TPM")
mat <- ex$mat
# 表达列 = 样本 id，与临床 sample 表的 SAMPLE_ID 对齐（本研究中 1 患者 = 1 样本）
n_match <- mean(colnames(mat) %in% clin_s$SAMPLE_ID)
message(sprintf("  join key: expr colname == clinical SAMPLE_ID  (%d/%d = %.1f%%)",
                sum(colnames(mat) %in% clin_s$SAMPLE_ID), ncol(mat), 100 * n_match))
stopifnot(n_match > 0.99)
sm <- clin_s[match(colnames(mat), clin_s$SAMPLE_ID), , drop = FALSE]
stopifnot(!anyDuplicated(sm$PATIENT_ID), !anyDuplicated(sm$SAMPLE_ID))

pt <- clin_p[match(sm$PATIENT_ID, clin_p$PATIENT_ID), , drop = FALSE]
status_raw <- pt[["PFS_STATUS"]]; time_raw <- suppressWarnings(as.numeric(pt[["PFS_MONTHS"]]))
pfs_status <- ifelse(grepl("^1", status_raw), 1L,
              ifelse(grepl("^0", status_raw), 0L, NA_integer_))
ok <- !is.na(pfs_status) & !is.na(time_raw)
message(sprintf("  PFS usable %d/%d, events %d (time range %.2f-%.2f months)",
                sum(ok), nrow(pt), sum(pfs_status[ok] == 1),
                min(time_raw[ok]), max(time_raw[ok])))

# log2 判定（线性 TPM -> log2(x+1)）
sv <- classify_scale(as.vector(mat)); vd <- sv$verdict
message("  scale verdict: ", vd, " (", sv$evidence, ")")
if (vd == "needs_log2") {
  before <- range(mat, na.rm = TRUE)
  mat <- log2(mat + 1)
  log_rows[[length(log_rows) + 1L]] <- data.frame(
    dataset = acc, verdict = vd, evidence = sv$evidence,
    min_before = before[1], median_before = sv$stats$median, max_before = before[2],
    min_after = min(mat, na.rm = TRUE), max_after = max(mat, na.rm = TRUE),
    action = "log2(x+1)", stringsAsFactors = FALSE)
  message("  applied log2(x+1)")
} else {
  log_rows[[length(log_rows) + 1L]] <- data.frame(
    dataset = acc, verdict = vd, evidence = sv$evidence,
    min_before = sv$stats$min, median_before = sv$stats$median, max_before = sv$stats$max,
    min_after = NA_real_, max_after = NA_real_, action = "none", stringsAsFactors = FALSE)
}

expr_df <- data.frame(ID_REF = rownames(mat), mat, check.names = FALSE)
saveRDS(expr_df, file.path(ROOT, "data/expr", paste0(acc, ".rds")))

pheno <- data.frame(sample_id = sm$SAMPLE_ID, patient_id = sm$PATIENT_ID,
                    sample_type = sm$SAMPLE_TYPE, metastasized = sm$METASTASIZED,
                    cancer_type_detailed = sm$CANCER_TYPE_DETAILED,
                    immune_subtype = sm$IMMUNE_SUBTYPE,
                    row.names = sm$SAMPLE_ID, stringsAsFactors = FALSE)
saveRDS(pheno, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))

stage <- ifelse(is.na(pt[["CLINICAL_STAGE"]]), NA_character_,
                toupper(trimws(pt[["CLINICAL_STAGE"]])))
surv <- data.frame(
  PFS_status = pfs_status,
  PFS_time   = time_raw / 12,                       # 月 -> 年
  sex        = NA_character_,
  age        = NA_real_,
  histology  = "Renal Cell Carcinoma",
  stage      = stage,
  clinical_stage = pt[["CLINICAL_STAGE"]],
  clinical_benefit = pt[["CLINICAL_BENEFIT"]],
  ici_rx     = pt[["ICI_RX"]],
  ici_target = pt[["ICI_TARGET"]],
  progression = pt[["PROGRESSION"]],
  responder  = pt[["RESPONDER"]],
  sample_type = sm$SAMPLE_TYPE,
  metastasized = sm$METASTASIZED,
  immune_subtype = sm$IMMUNE_SUBTYPE,
  row.names = sm$SAMPLE_ID, stringsAsFactors = FALSE)
saveRDS(surv, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))
summ_rows[[length(summ_rows) + 1L]] <- data.frame(
  acc = acc, type = "Kidney Cancer", token = "PFS",
  n_expr_samples = ncol(mat), n_surv_rows = nrow(surv), N = sum(ok),
  n_events = sum(pfs_status[ok] == 1), join_key = "expr colname == clinical SAMPLE_ID (1 sample/patient)",
  source = "cBioPortal rcc_iatlas_immotion150_2018 (McDermott Nat Med 2018, PMID 29867230)",
  stringsAsFactors = FALSE)

# ===========================================================================
# 队列 2: A5-PCPG —— Pheochromocytoma / Paraganglioma, OS
# ===========================================================================
message("== A5-PCPG (Pheochromocytoma/Paraganglioma, OS)")
acc2 <- "A5-PCPG"
clin_p2 <- read_cbio_clin(file.path(SUPPL, "A5PCPG_data_clinical_patient.txt.gz"))
clin_s2 <- read_cbio_clin(file.path(SUPPL, "A5PCPG_data_clinical_sample.txt.gz"))

ex2 <- load_expr(file.path(SUPPL, "A5PCPG_data_mrna_seq_cpm.txt.gz"), "CPM")
mat2 <- ex2$mat
n_match2 <- mean(colnames(mat2) %in% clin_s2$SAMPLE_ID)
message(sprintf("  join key: expr colname == clinical SAMPLE_ID  (%d/%d = %.1f%%)",
                sum(colnames(mat2) %in% clin_s2$SAMPLE_ID), ncol(mat2), 100 * n_match2))
stopifnot(n_match2 > 0.99)
sm2 <- clin_s2[match(colnames(mat2), clin_s2$SAMPLE_ID), , drop = FALSE]

# 患者级去重：优先 primary 样本，其次 sample id 最小者（每个患者保留 1 个样本）
is_primary <- grepl("Primary", sm2$SAMPLE_TYPE, ignore.case = TRUE)
ord <- order(sm2$PATIENT_ID, !is_primary, sm2$SAMPLE_ID)
keep_idx <- ord[!duplicated(sm2$PATIENT_ID[ord])]
keep <- sort(keep_idx)
message(sprintf("  patient dedup: %d samples -> %d patients (%d with >1 sample); kept primary: %d",
                nrow(sm2), length(unique(sm2$PATIENT_ID)),
                sum(table(sm2$PATIENT_ID) > 1), sum(is_primary[keep])))
kept_sm <- sm2[keep, , drop = FALSE]

pt2 <- clin_p2[match(kept_sm$PATIENT_ID, clin_p2$PATIENT_ID), , drop = FALSE]
os_status <- ifelse(grepl("^1", pt2$OS_STATUS), 1L,
             ifelse(grepl("^0", pt2$OS_STATUS), 0L, NA_integer_))
os_time_m <- suppressWarnings(as.numeric(pt2$OS_MONTHS))
ok2 <- !is.na(os_status) & !is.na(os_time_m)
message(sprintf("  OS usable %d/%d patients, deaths %d (time range %.2f-%.2f months)",
                sum(ok2), length(keep), sum(os_status[ok2] == 1),
                min(os_time_m[ok2]), max(os_time_m[ok2])))

vd2 <- scale_verdict(as.vector(mat2))
message("  scale verdict: ", vd2)
if (vd2 == "needs_log2") {
  before <- range(mat2, na.rm = TRUE)
  mat2 <- log2(mat2 + 1)
  log_rows[[length(log_rows) + 1L]] <- data.frame(
    dataset = acc2, verdict = vd2, min_before = before[1], max_before = before[2],
    min_after = min(mat2, na.rm = TRUE), max_after = max(mat2, na.rm = TRUE),
    action = "log2(x+1)", stringsAsFactors = FALSE)
} else {
  log_rows[[length(log_rows) + 1L]] <- data.frame(
    dataset = acc2, verdict = vd2, min_before = min(mat2), max_before = max(mat2),
    min_after = NA_real_, max_after = NA_real_, action = "none", stringsAsFactors = FALSE)
}

# 表达表保留全部 91 个样本（与 GSE32918 的 249 arrays 一致），surv 表只保留
# 患者级去重后的样本行 -> 镜像 join 后不会把同一患者重复计入。
expr2_df <- data.frame(ID_REF = rownames(mat2), mat2, check.names = FALSE)
saveRDS(expr2_df, file.path(ROOT, "data/expr", paste0(acc2, ".rds")))

pheno2 <- data.frame(sample_id = sm2$SAMPLE_ID, patient_id = sm2$PATIENT_ID,
                     sample_type = sm2$SAMPLE_TYPE,
                     tumor_location = sm2$TUMOR_LOCATION,
                     oncotree_code = sm2$ONCOTREE_CODE,
                     row.names = sm2$SAMPLE_ID, stringsAsFactors = FALSE)
saveRDS(pheno2, file.path(ROOT, "data/pheno", paste0(acc2, ".rds")))

sex2 <- ifelse(is.na(kept_sm$SEX), NA_character_,
        ifelse(grepl("^m", kept_sm$SEX, ignore.case = TRUE), "male", "female"))
surv2 <- data.frame(
  OS_status = os_status,
  OS_time   = os_time_m / 12,                        # 月 -> 年
  sex       = sex2,
  age       = NA_real_,
  histology = ifelse(is.na(kept_sm$ONCOTREE_CODE), "Paraganglioma", kept_sm$ONCOTREE_CODE),
  stage     = NA_character_,
  sample_type = kept_sm$SAMPLE_TYPE,
  location_of_primary_pcpg = kept_sm$LOCATION_OF_PRIMARY_PCPG,
  tumor_location = kept_sm$TUMOR_LOCATION,
  catecholamine_profile = kept_sm$CATECHOLAMINE_PROFILE,
  ki67_staining = suppressWarnings(as.numeric(kept_sm$KI67_STAINING)),
  sdhb_staining = kept_sm$SDHB_STAINING,
  developed_metastatic_disease = pt2$DEVELOPED_METASTATIC_DISEASE,
  disease_burden_at_presentation = pt2$DISEASE_BURDEN_AT_PRESENTATION,
  methylation_cluster = kept_sm$METHYLATION,
  row.names = kept_sm$SAMPLE_ID, stringsAsFactors = FALSE)
saveRDS(surv2, file.path(ROOT, "data/processed/surv", paste0(acc2, "_surv.rds")))
summ_rows[[length(summ_rows) + 1L]] <- data.frame(
  acc = acc2, type = "Pheochromocytoma", token = "OS",
  n_expr_samples = ncol(mat2), n_surv_rows = nrow(surv2), N = sum(ok2),
  n_events = sum(os_status[ok2] == 1),
  join_key = "expr colname == clinical SAMPLE_ID; surv rows = 1 sample/patient (primary preferred)",
  source = "cBioPortal hnsc_a5consortium_2025 (A5 Consortium, Nature 2025, PMID 40097403)",
  stringsAsFactors = FALSE)

# ===========================================================================
# 合成平台注释表（gene-level，两个队列共用）
# ===========================================================================
genes <- unique(c(rownames(mat), rownames(mat2)))
dest <- file.path(ROOT, "data/processed/gpl", paste0(GPL_NAME, ".rds"))
existing <- if (file.exists(dest)) readRDS(dest) else NULL
add <- setdiff(genes, if (is.null(existing)) character(0) else rownames(existing))
gpl_tab <- data.frame(gene_id = symbol2entrez(add), row.names = add, stringsAsFactors = FALSE)
if (!is.null(existing)) gpl_tab <- rbind(existing, gpl_tab)
gpl_tab <- gpl_tab[!duplicated(rownames(gpl_tab)), , drop = FALSE]
saveRDS(gpl_tab, dest)
message(sprintf("== %s: %d rows, %d with ENTREZ (%.1f%%) -> %s",
                GPL_NAME, nrow(gpl_tab), sum(!is.na(gpl_tab$gene_id)),
                100 * mean(!is.na(gpl_tab$gene_id)), dest))

# ---- 记录 -----------------------------------------------------------------
logdf <- do.call(rbind, log_rows)
write.csv(logdf, file.path(OUT, "suppl_expansion_log2.csv"), row.names = FALSE)
summ <- do.call(rbind, summ_rows)
write.csv(summ, file.path(OUT, "suppl_expansion_build.csv"), row.names = FALSE)
cat("\n=== 95 build summary ===\n"); print(summ, row.names = FALSE)
cat("\n=== log2 decisions ===\n"); print(logdf, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# run_96_build_suppl_expansion()  <-  verbatim pipeline/R/96_build_suppl_expansion.R
# ---------------------------------------------------------------------------
run_96_build_suppl_expansion <- function() {
# 96_build_suppl_expansion.R --------------------------------------------------
# 补充文件 / 发表表格路线的队列构建（2026-09-24 第二轮）。
# 矩阵级路线已穷尽（37 号筛选：8 个 gap 癌种 0 个队列通过 >50 time+status 闸门），
# 本脚本构建「临床来自 cBioPortal / 独立队列临床文件 + 表达来自独立表达矩阵」的队列。
#
# 两个队列（数据源见 pipeline/out/suppl_expansion_sources.csv）：
#
#   1) A5-PCPG   Pheochromocytoma（嗜铬细胞瘤/副神经节瘤）
#      cBioPortal study hnsc_a5consortium_2025
#        = "Hereditary SDHB-Mutant Pheochromocytomas and Paragangliomas
#           (A5 Consortium, Nature Comm 2025)"（cancerTypeId 被错标为 hnsc，
#           队列本体是 PPGL；已用 clinical-attributes 的 PCPG 专有字段核对）
#      临床 data_clinical_patient.txt  : OS_MONTHS/OS_STATUS(1:DECEASED / 0:LIVING)
#      + data_clinical_sample.txt      : SAMPLE_TYPE
#      表达 data_mrna_seq_cpm.txt      : Hugo_Symbol × SAMPLE_ID，已是 log2 CPM
#                                        （min<0 -> 不再做 log2）
#      患者级去重：每位患者保留 1 个样本（优先 Primary；否则样本名最小者）
#      终点：OS（>50 患者）。患者级 70 例可用；DFS 仅 35 例可用 -> <50，不登记。
#
#   2) IMmotion150  Kidney Cancer（独立于 TCGA-KIRC，Roche/Genentech II 期试验）
#      cBioPortal study rcc_iatlas_immotion150_2018 (n=263，1 患者 1 样本)
#      临床 data_clinical_patient.txt : PFS_MONTHS/PFS_STATUS(1:Progressed /
#                                        0:Not_Progressed)，自治疗开始计
#      表达 data_mrna_seq_tpm.txt     : Hugo_Symbol × PATIENT_ID
#      终点：PFS 263/164。口径与本批既有 ICI 治疗队列一致
#            （GSE159067 / GSE162520 也是治疗起始 PFS）。
#
# 产物（沿用 92/93 的格式约定）：
#   data/expr/<ACC>.rds            ID_REF + 每样本一列（基因级，不折叠）
#   data/pheno/<ACC>.rds           样本级临床表
#   data/processed/surv/<ACC>_surv.rds  <TOK>_status(0/1) + <TOK>_time(年) + 临床列
#   data/processed/gpl/GPL24676.rds    按需把基因符号并入既有 GPL 注释表
#   pipeline/out/suppl_expansion_build_log.csv
#
# 用法: Rscript pipeline/R/96_build_suppl_expansion.R [ROOT]
# ----------------------------------------------------------------------------
# --- explicit output root (A10) ------------------------------------------------
# Derived outputs are written under an explicit, env-overridable root
# (CPAS_OUT_ROOT, default <CPAS_DATA_ROOT>/pipeline/out).  See pipeline/OUTPUT_LAYOUT.md.
# Uses the shared helper when it is present and falls back to identical local
# definitions otherwise (so a consolidated / installed copy is self-contained).
# With the default root every path below resolves exactly where it did before:
# this makes the location explicit and overridable, it changes no computation.
.cpas_helper <- file.path(Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS"),
                          "pipeline/R/00_output_root.R")
if (file.exists(.cpas_helper)) source(.cpas_helper)
if (!exists("cpas_out", mode = "function")) {
  .cpas_root <- function() {
    r <- Sys.getenv("CPAS_DATA_ROOT", unset = "")
    if (!nzchar(r)) r <- "/home/Jingle/data/Project/CPAS"
    path.expand(r)
  }
  cpas_root <- .cpas_root
  cpas_out_root <- function() {
    o <- Sys.getenv("CPAS_OUT_ROOT", unset = "")
    if (nzchar(o)) path.expand(o) else file.path(.cpas_root(), "pipeline", "out")
  }
  cpas_out <- function(...) file.path(cpas_out_root(), ...)
  cpas_data <- function(...) file.path(.cpas_root(), "data", ...)
  cpas_suppl <- function(...) file.path(.cpas_root(), "data", "suppl", ...)
}

suppressMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
SUPPL <- file.path(ROOT, "data/suppl")
stopifnot(dir.exists(SUPPL))
dir.create(file.path(ROOT, "data/expr"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/pheno"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/processed/surv"), showWarnings = FALSE, recursive = TRUE)

read_clin <- function(f) {
  # cBioPortal 临床文件：前 4 行是 '#' 注释（描述/类型/优先级），第 5 行才是表头。
  # fread 不会自动跳过 '#' 行（实测），所以先 readLines 过滤再解析。
  ln <- readLines(gzfile(file.path(SUPPL, f)), warn = FALSE)
  ln <- ln[!grepl("^#", ln)]
  x <- fread(text = paste(ln, collapse = "\n"), sep = "\t", header = TRUE,
             na.strings = c("NA", "NaN", ""), check.names = FALSE)
  as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
}

read_expr <- function(f) {
  x <- fread(cmd = paste("zcat", shQuote(file.path(SUPPL, f))), sep = "\t",
             header = TRUE, check.names = FALSE)
  ids <- as.character(x[[1]])
  m <- as.data.frame(x[, -1, with = FALSE], check.names = FALSE)
  m <- as.data.frame(lapply(m, function(v) suppressWarnings(as.numeric(v))),
                     check.names = FALSE)
  # 基因符号级矩阵：同名符号（如 SNORD 家族 / 假基因）按行均值折叠，
  # 否则 06_upload_db.R 的 column_to_rownames("ID_REF") 会因重复行名报错。
  if (anyDuplicated(ids)) {
    n0 <- length(ids)
    mm <- as.matrix(m); storage.mode(mm) <- "double"
    s <- rowsum(mm, group = ids, reorder = TRUE, na.rm = TRUE)
    c <- rowsum(matrix(as.numeric(!is.na(mm)), nrow = nrow(mm)), group = ids,
                reorder = TRUE, na.rm = TRUE)
    c[c == 0] <- 1
    m <- as.data.frame(s / c, check.names = FALSE)
    ids <- rownames(m)
    message("  collapsed ", n0 - length(ids), " duplicate gene symbols by row-mean -> ",
            length(ids), " rows")
  }
  rownames(m) <- ids
  m
}


# GEO series matrix 的轻量解析（表头 + 表达体），不依赖 Biobase
read_gsm_matrix <- function(path) {
  ln <- readLines(gzfile(path), warn = FALSE)
  beg <- grep("^!series_matrix_table_begin", ln); end <- grep("^!series_matrix_table_end", ln)
  if (length(beg) != 1 || length(end) != 1) stop("series matrix table markers not found in ", path)
  hl <- ln[grepl("^!Sample_", ln)]
  keys <- sub("\\t.*$", "", sub("^!Sample_", "", hl))
  vals <- lapply(hl, function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
  names(vals) <- keys
  body <- ln[(beg + 1L):(end - 1L)]
  dt <- fread(text = paste(body, collapse = "\n"), sep = "\t", header = TRUE,
              check.names = FALSE, na.strings = c("", "NA", "null"))
  ids <- as.character(dt[[1]])
  m <- as.data.frame(dt[, -1, with = FALSE], check.names = FALSE)
  m <- as.data.frame(lapply(m, function(v) suppressWarnings(as.numeric(v))), check.names = FALSE)
  rownames(m) <- ids
  pd <- data.frame(row.names = vals$geo_accession, stringsAsFactors = FALSE)
  for (k in names(vals)) {
    v <- vals[[k]]
    pd[[k]] <- if (length(v) == nrow(pd)) v else rep(NA_character_, nrow(pd))
  }
  list(mat = m, pdata = pd)
}

# 与 92_build_geo_expansion_expr_pheno.R 完全相同的「是否已是 log」判定
scale_verdict <- function(m) {
  v <- as.matrix(m); storage.mode(v) <- "double"; v <- v[is.finite(v)]
  if (!length(v)) return("empty")
  if (min(v) < 0) return("already_log")
  if (median(v) > 0 && max(v) / median(v) <= 50) return("already_log")
  if (median(v) <= 0 && max(v) <= 25) return("already_log")
  "linear_needs_log2"
}

# 基因符号 -> ENTREZ（并入既有 GPL 注释表；与 92 的 build_gpl_map 同口径）
extend_gpl <- function(gpl, ids) {
  dest <- file.path(ROOT, "data/processed/gpl", paste0(gpl, ".rds"))
  existing <- if (file.exists(dest)) readRDS(dest) else
    data.frame(gene_id = character(0))
  miss <- setdiff(unique(as.character(ids)), rownames(existing))
  if (!length(miss)) { message("  ", gpl, ": nothing to add (", nrow(existing), " rows)"); return(existing) }
  db <- tryCatch(getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"), error = function(e) NULL)
  gid <- setNames(rep(NA_character_, length(miss)), miss)
  if (!is.null(db)) {
    m <- tryCatch(suppressWarnings(AnnotationDbi::select(db, keys = miss,
             columns = "ENTREZID", keytype = "SYMBOL")), error = function(e) NULL)
    if (!is.null(m)) {
      m <- m[!is.na(m$ENTREZID), , drop = FALSE]
      if (nrow(m)) {
        agg <- tapply(m$ENTREZID, m$SYMBOL, function(z) paste(sort(unique(z)), collapse = " /// "))
        gid[names(agg)] <- as.character(agg)
      }
    }
  }
  # 注意：不能用 data.frame(..., row.names = miss) + rbind —— 实测 rbind 会把
  # 新块的字符行名丢掉换成 1..n（GPL24676 曾因此被写坏 59409 行，已修复）。
  # 这里显式拼接行名，并强制 gene_id 为字符。
  add_gid <- unname(gid)
  rn <- c(rownames(existing), names(gid))
  out <- data.frame(gene_id = c(as.character(existing$gene_id), add_gid),
                    stringsAsFactors = FALSE)
  attr(out, "row.names") <- rn
  out <- out[!duplicated(rownames(out)), , drop = FALSE]
  # 防御：调用点若误传了 data.frame 的自增行名（1..n），这里必须报错而不是静默写坏表
  if (!all(unique(as.character(ids)) %in% rownames(out)))
    stop("extend_gpl: not all requested ids became row names (ids arg looks wrong)")
  saveRDS(out, dest)
  message("  ", gpl, ": +", length(miss), " symbol rows -> ", nrow(out),
          " (", sum(!is.na(out$gene_id)), " with ENTREZ)")
  out
}

sanitize <- function(nm) {
  nm <- tolower(gsub("[^A-Za-z0-9]+", "_", nm))
  nm <- gsub("^_+|_+$", "", nm)
  make.unique(nm)
}

log_rows <- list()
add_log <- function(...) log_rows[[length(log_rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# =============================================================== 1) A5-PCPG ==
message("=== A5-PCPG (Pheochromocytoma) ===")
acc <- "A5-PCPG"; gpl <- "GPL24676"; tok <- "OS"
pat <- read_clin("A5PCPG_data_clinical_patient.txt.gz")
smp <- read_clin("A5PCPG_data_clinical_sample.txt.gz")
ex  <- read_expr("A5PCPG_data_mrna_seq_cpm.txt.gz")
message("  patients=", nrow(pat), " samples=", nrow(smp), " expr=", paste(dim(ex), collapse = "x"))
req <- c("PATIENT_ID","SEX","OS_MONTHS","OS_STATUS","DFS_STATUS","DFS_MONTHS",
         "DEVELOPED_METASTATIC_DISEASE","LOCATION_OF_METASTASIS",
         "LARGEST_PRIMARY_DIMENSIONS_(CM)","HISTORY_HYPERTENSION","FAMILY_HISTORY_PPGL",
         "POST_DIAGNOSIS_FOLLOW_UP_MONTHS")
miss <- setdiff(req, names(pat)); if (length(miss)) stop("A5 patient cols missing: ", paste(miss, collapse=","))
s2p <- setNames(smp$PATIENT_ID, smp$SAMPLE_ID)
s2t <- setNames(smp$SAMPLE_TYPE, smp$SAMPLE_ID)
cols <- colnames(ex)
stopifnot(all(cols %in% names(s2p)))
pid <- unname(s2p[cols])
# 患者级去重：优先 Primary，其次样本名最小
pri <- ifelse(grepl("Primary", unname(s2t[cols])), 0L, 1L)
ord  <- order(pid, pri, cols)
keep <- ord[!duplicated(pid[ord])]
ex1  <- ex[, keep, drop = FALSE]
keep_s <- colnames(ex1); keep_p <- unname(s2p[keep_s])
message("  dedup: ", length(cols), " expr samples -> ", length(keep_s), " patients")
ex_out <- data.frame(ID_REF = rownames(ex1), ex1, check.names = FALSE)
rownames(ex_out) <- NULL
saveRDS(ex_out, file.path(ROOT, "data/expr", paste0(acc, ".rds")))

p <- pat[match(keep_p, pat$PATIENT_ID), , drop = FALSE]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
os_time_m <- num(p$OS_MONTHS)
os_st <- ifelse(grepl("^1", p$OS_STATUS), 1, ifelse(grepl("^0", p$OS_STATUS), 0, NA))
sv <- data.frame(
  OS_status = os_st,
  OS_time   = os_time_m / 12,
  sex       = tolower(as.character(p$SEX)),
  sample_type = unname(s2t[keep_s]),
  metastatic_disease = as.character(p$DEVELOPED_METASTATIC_DISEASE),
  location_of_metastasis = as.character(p$LOCATION_OF_METASTASIS),
  largest_primary_dim_cm = num(p$`LARGEST_PRIMARY_DIMENSIONS_(CM)`),
  history_hypertension = as.character(p$HISTORY_HYPERTENSION),
  family_history_ppgl = as.character(p$FAMILY_HISTORY_PPGL),
  post_diagnosis_follow_up_months = num(p$POST_DIAGNOSIS_FOLLOW_UP_MONTHS),
  row.names = keep_s, stringsAsFactors = FALSE, check.names = FALSE)
# 列名消毒（DB 表列名不允许括号/空格）
cn <- sanitize(colnames(sv)); sv <- as.data.frame(sv, check.names = FALSE); colnames(sv) <- cn
stopifnot("os_status" %in% colnames(sv), "os_time" %in% colnames(sv))
saveRDS(sv, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))
ph <- data.frame(patient_id = keep_p, sample_type = unname(s2t[keep_s]),
                 row.names = keep_s, stringsAsFactors = FALSE)
saveRDS(ph, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
invisible(extend_gpl(gpl, ex_out$ID_REF))
n_ok <- sum(!is.na(sv$os_status) & !is.na(sv$os_time)); ev <- sum(sv$os_status == 1, na.rm = TRUE)
dfs_st <- ifelse(grepl("^1", p$DFS_STATUS), 1, ifelse(grepl("^0", p$DFS_STATUS), 0, NA))
n_dfs <- sum(!is.na(dfs_st) & !is.na(num(p$DFS_MONTHS)))
cat(sprintf("  A5-PCPG: %d samples | OS usable %d (events %d) | DFS usable %d (below 50, not registered) | scale=%s\n",
            nrow(sv), n_ok, ev, n_dfs, scale_verdict(ex1)))
add_log(acc = acc, type = "Pheochromocytoma", gpl = gpl, token = tok, n_expr = nrow(sv),
        N_patients = n_ok, n_events = ev, scale = scale_verdict(ex1),
        detail = sprintf("A5 Consortium PPGL (cBioPortal hnsc_a5consortium_2025); %d expr samples -> %d patients (primary preferred); DFS usable only %d (<50) so OS only",
                         length(cols), length(keep_s), n_dfs))

# ============================================================= 2) IMmotion150 ==
message("=== IMmotion150 (Kidney Cancer) ===")
acc <- "IMmotion150"; gpl <- "GPL24676"; tok <- "PFS"
pat <- read_clin("IMmotion150_data_clinical_patient.txt.gz")
smp <- read_clin("IMmotion150_data_clinical_sample.txt.gz")
ex  <- read_expr("IMmotion150_data_mrna_seq_tpm.txt.gz")
message("  patients=", nrow(pat), " samples=", nrow(smp), " expr=", paste(dim(ex), collapse = "x"))
s2p <- setNames(smp$PATIENT_ID, smp$SAMPLE_ID)
cols <- colnames(ex)
stopifnot(all(cols %in% names(s2p)))
pid <- unname(s2p[cols])
stopifnot(!anyDuplicated(pid))                      # 1 患者 1 样本
p <- pat[match(pid, pat$PATIENT_ID), , drop = FALSE]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
pfs_m <- num(p$PFS_MONTHS)
pfs_st <- ifelse(grepl("^1", p$PFS_STATUS), 1, ifelse(grepl("^0", p$PFS_STATUS), 0, NA))
sv <- data.frame(
  PFS_status = pfs_st,
  PFS_time   = pfs_m / 12,
  pfs_from_treatment_months = pfs_m,
  clinical_stage = as.character(p$CLINICAL_STAGE),
  ici_rx = as.character(p$ICI_RX),
  ici_target = as.character(p$ICI_TARGET),
  responder = as.character(p$RESPONDER),
  row.names = cols, stringsAsFactors = FALSE, check.names = FALSE)
cn <- sanitize(colnames(sv)); colnames(sv) <- cn
saveRDS(sv, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))
ph <- data.frame(patient_id = pid,
                 sample_type = as.character(smp$SAMPLE_TYPE[match(cols, smp$SAMPLE_ID)]),
                 cancer_type_detailed = as.character(smp$CANCER_TYPE_DETAILED[match(cols, smp$SAMPLE_ID)]),
                 immune_subtype = as.character(smp$IMMUNE_SUBTYPE[match(cols, smp$SAMPLE_ID)]),
                 row.names = cols, stringsAsFactors = FALSE)
saveRDS(ph, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
ex_out <- data.frame(ID_REF = rownames(ex), ex, check.names = FALSE); rownames(ex_out) <- NULL
saveRDS(ex_out, file.path(ROOT, "data/expr", paste0(acc, ".rds")))
invisible(extend_gpl(gpl, ex_out$ID_REF))
n_ok <- sum(!is.na(sv$pfs_status) & !is.na(sv$pfs_time)); ev <- sum(sv$pfs_status == 1, na.rm = TRUE)
cat(sprintf("  IMmotion150: %d samples | PFS usable %d (events %d) | scale=%s\n",
            nrow(sv), n_ok, ev, scale_verdict(ex)))
add_log(acc = acc, type = "Kidney Cancer", gpl = gpl, token = tok, n_expr = nrow(sv),
        N_patients = n_ok, n_events = ev, scale = scale_verdict(ex),
        detail = "IMmotion150 RCC trial (cBioPortal rcc_iatlas_immotion150_2018), 263 patients 1:1; PFS from start of ICI treatment (same convention as GSE159067/GSE162520)")


# =============================================================== 3) GSE3218 ===
# Testicular Cancer（成人男性生殖细胞瘤；独立于 TCGA-TGCT）
#   表达 : GEO series matrix GSE3218-GPL96（107 个肿瘤样本，RMA，已是 log2）
#   临床 : 发表表格 —— PMC4666461 (PLoS One 2015, PMID 26624623) S1 Table
#          "clinical features for the patients included in this study"
#          = data/suppl/GSE3218_PMC4666461_S1_Table_clinical_features.xlsx
#          列: F/u time (1st tx to last f/u) in yrs  +  Surv status at last f/u (1=dead, 2=alive)
#          连接键: 表格 Sample（如 052B）== GEO !Sample_title 里的 [0-9]{3}[A-Z] 词元
#          108 行临床中 74 行落在 GSE3218（全部 Expression Training == Yes）；
#          另 34 行是验证集 GSE10783（患者级 34 例 < 50，单独登记为不足门槛，见 excluded.csv）。
#   ★ 不做 2y DFS / 5y DSS：那是里程碑（milestone）二分指标，不是 time-to-event。
message("=== GSE3218 (Testicular Cancer) ===")
acc <- "GSE3218"; gpl <- "GPL96"; tok <- "OS"
gm <- read_gsm_matrix(file.path(ROOT, "data/raw", "GSE3218-GPL96_series_matrix.txt.gz"))
ex <- gm$mat                                     # 行=探针, 列=GSM
gsm <- colnames(ex); titles <- gm$pdata$title
sid <- stringr::str_extract(titles, "[0-9]{3}[A-Z]")
message("  samples=", ncol(ex), " probes=", nrow(ex), " parsed ids=", sum(!is.na(sid)))
cl <- as.data.frame(readxl::read_excel(file.path(SUPPL, "GSE3218_PMC4666461_S1_Table_clinical_features.xlsx"),
                                       sheet = 1), check.names = FALSE)
fu <- suppressWarnings(as.numeric(cl[["F/u time (1st tx to last f/u) in yrs"]])); names(fu) <- cl$Sample
ss <- suppressWarnings(as.numeric(cl[["Surv status at last f/u (1=dead, 2= alive)"]])); names(ss) <- cl$Sample
risk <- cl[["IGCCCG Risk group (1=good, 2=int, 3=poor)"]]; names(risk) <- cl$Sample
hist2 <- cl[["Histology2"]]; names(hist2) <- cl$Sample
site1 <- cl[["Site1"]]; names(site1) <- cl$Sample
resp  <- cl[["Best response7"]]; names(resp) <- cl$Sample
# 每位患者只在一个样本上填生存（若同一 id 在矩阵里出现 >1 次）
first <- !duplicated(sid) & !is.na(sid)
os_t <- ifelse(first, unname(fu[sid]), NA_real_)
os_s <- ifelse(first & !is.na(unname(ss[sid])), ifelse(unname(ss[sid]) == 1, 1, 0), NA_real_)
sv <- data.frame(
  OS_status  = os_s,
  OS_time    = os_t,
  site       = ifelse(first, unname(site1[sid]), NA_character_),
  histology  = ifelse(first, unname(hist2[sid]), NA_character_),
  igcccg_risk_group = ifelse(first, unname(risk[sid]), NA_character_),
  best_response = ifelse(first, unname(resp[sid]), NA_character_),
  row.names = gsm, stringsAsFactors = FALSE)
colnames(sv) <- sanitize(colnames(sv))
saveRDS(sv, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))
ex_out <- data.frame(ID_REF = rownames(ex), ex, check.names = FALSE); rownames(ex_out) <- NULL
saveRDS(ex_out, file.path(ROOT, "data/expr", paste0(acc, ".rds")))
ph <- data.frame(patient_id = sid, title = titles, row.names = gsm, stringsAsFactors = FALSE)
saveRDS(ph, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
g96 <- readRDS(file.path(ROOT, "data/processed/gpl", paste0(gpl, ".rds")))
n_in <- mean(as.character(ex_out$ID_REF) %in% rownames(g96))
n_ok <- sum(!is.na(sv$os_status) & !is.na(sv$os_time)); ev <- sum(sv$os_status == 1, na.rm = TRUE)
cat(sprintf("  GSE3218: %d samples (%d patients) | OS usable %d (events %d) | probes in GPL96 map %.1f%% | scale=%s\n",
            nrow(sv), length(unique(na.omit(sid))), n_ok, ev, 100 * n_in, scale_verdict(ex)))
add_log(acc = acc, type = "Testicular Cancer", gpl = gpl, token = tok, n_expr = nrow(sv),
        N_patients = n_ok, n_events = ev, scale = scale_verdict(ex),
        detail = "GSE3218-GPL96 adult male GCT (RMA); clinical from PMC4666461 (PLoS One 2015 PMID 26624623) S1 Table, join on [0-9]{3}[A-Z] token of !Sample_title; the other 34 clinical rows are the GSE10783 validation arm (<50 patients, not registered). 2y DFS / 5y DSS are milestone binaries -> not used")

res <- do.call(rbind, log_rows)
write.csv(res, cpas_out("suppl_expansion_build_log.csv"), row.names = FALSE)
cat("\n=== 96 build summary ===\n"); print(res, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# run_100_build_relaxed_gate()  <-  verbatim pipeline/R/100_build_relaxed_gate.R
# ---------------------------------------------------------------------------
run_100_build_relaxed_gate <- function() {
# 100_build_relaxed_gate.R ---------------------------------------------------
# 松弛闸门批次（>=30 患者带可用 time+status 对；旧闸门 >50）的队列构建。
# 沿用 92/93/96 的产物约定：
#   data/expr/<ACC>.rds                 ID_REF + 每样本一列（探针级）
#   data/pheno/<ACC>.rds                样本级临床
#   data/processed/surv/<ACC>_surv.rds  <TOK>_status(0/1) + <TOK>_time(年) + 临床列
#
# 本批第一个队列：
#   GSE76019  Adrenocortical Cancer（儿童肾上腺皮质癌，COG ARAR0332，PMID 27307598）
#     表达 : GEO series matrix GSE76019-GPL13158（54715 探针 × 34 样本，已是 log 尺度）
#     临床 : 同系列 !Sample_characteristics_ch1 —— histology=ACC / Stage /
#            efs.time（年；由整天数换算，例 1.48665297741273 y = 543 d）/
#            efs.event（0/1）
#     患者级 : 34/34 有可用 (efs.time, efs.event) 对，12 个事件
#     终点 : EFS —— 11_endpoint_families.R 的 FAMILY_OF 已把 EFS 归入 DFS 家族
#     ★ 34 < 旧闸门 50，按作者 2026-09-24 决定以「松弛闸门 >=30」登记，catalog Note 必须写明
#
# 用法: Rscript pipeline/R/100_build_relaxed_gate.R [ROOT]
# ----------------------------------------------------------------------------
# --- explicit output root (A10) ------------------------------------------------
# Derived outputs are written under an explicit, env-overridable root
# (CPAS_OUT_ROOT, default <CPAS_DATA_ROOT>/pipeline/out).  See pipeline/OUTPUT_LAYOUT.md.
# Uses the shared helper when it is present and falls back to identical local
# definitions otherwise (so a consolidated / installed copy is self-contained).
# With the default root every path below resolves exactly where it did before:
# this makes the location explicit and overridable, it changes no computation.
.cpas_helper <- file.path(Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS"),
                          "pipeline/R/00_output_root.R")
if (file.exists(.cpas_helper)) source(.cpas_helper)
if (!exists("cpas_out", mode = "function")) {
  .cpas_root <- function() {
    r <- Sys.getenv("CPAS_DATA_ROOT", unset = "")
    if (!nzchar(r)) r <- "/home/Jingle/data/Project/CPAS"
    path.expand(r)
  }
  cpas_root <- .cpas_root
  cpas_out_root <- function() {
    o <- Sys.getenv("CPAS_OUT_ROOT", unset = "")
    if (nzchar(o)) path.expand(o) else file.path(.cpas_root(), "pipeline", "out")
  }
  cpas_out <- function(...) file.path(cpas_out_root(), ...)
  cpas_data <- function(...) file.path(.cpas_root(), "data", ...)
  cpas_suppl <- function(...) file.path(.cpas_root(), "data", "suppl", ...)
}

suppressMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
stopifnot(dir.exists(file.path(ROOT, "data/raw")))
dir.create(file.path(ROOT, "data/expr"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/pheno"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/processed/surv"), showWarnings = FALSE, recursive = TRUE)

# cBioPortal 临床文件读取（与 96 同实现：前 4 行是 '#' 注释，之后才是表头）
read_clin <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE)
  ln <- ln[!grepl("^#", ln)]
  x <- data.table::fread(text = paste(ln, collapse = "\n"), sep = "\t", header = TRUE,
                         na.strings = c("NA", "NaN", ""), check.names = FALSE)
  as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
}

# --- 与 96 完全相同的工具函数 ------------------------------------------------
sanitize <- function(nm) {
  nm <- tolower(gsub("[^A-Za-z0-9]+", "_", nm))
  nm <- gsub("^_+|_+$", "", nm)
  make.unique(nm)
}
scale_verdict <- function(m) {
  v <- as.matrix(m); storage.mode(v) <- "double"; v <- v[is.finite(v)]
  if (!length(v)) return("empty")
  if (min(v) < 0) return("already_log")
  if (median(v) > 0 && max(v) / median(v) <= 50) return("already_log")
  if (median(v) <= 0 && max(v) <= 25) return("already_log")
  "linear_needs_log2"
}
# GEO series matrix 轻量解析（表头 + 表达体），与 96 同实现
read_gsm_matrix <- function(path) {
  ln <- readLines(gzfile(path), warn = FALSE)
  beg <- grep("^!series_matrix_table_begin", ln); end <- grep("^!series_matrix_table_end", ln)
  if (length(beg) != 1 || length(end) != 1) stop("series matrix table markers not found in ", path)
  hl <- ln[grepl("^!Sample_", ln)]
  keys <- sub("\\t.*$", "", sub("^!Sample_", "", hl))
  vals <- lapply(hl, function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
  names(vals) <- keys
  body <- ln[(beg + 1L):(end - 1L)]
  dt <- fread(text = paste(body, collapse = "\n"), sep = "\t", header = TRUE,
              check.names = FALSE, na.strings = c("", "NA", "null"))
  ids <- as.character(dt[[1]])
  m <- as.data.frame(dt[, -1, with = FALSE], check.names = FALSE)
  m <- as.data.frame(lapply(m, function(v) suppressWarnings(as.numeric(v))), check.names = FALSE)
  rownames(m) <- ids
  pd <- list()
  for (k in unique(keys)) {
    rows <- vals[keys == k]
    for (r in rows) if (length(r) == ncol(m)) pd[[length(pd) + 1L]] <- r
  }
  # pdata：所有「每个样本恰好一个值」的 !Sample_ 字段（title / supplementary_file / platform_id ...）
  pd <- list(); n <- ncol(m)
  for (k in unique(keys)) {
    rows <- vals[keys == k]
    if (length(rows) == 1 && length(rows[[1]]) == n) pd[[k]] <- rows[[1]]
  }
  list(mat = m, chars = pd_chars(vals, keys),
       pdata = data.frame(pd, row.names = vals$geo_accession,
                          stringsAsFactors = FALSE, check.names = FALSE))
}
# 供 read_gsm_matrix 用：仅取 characteristics 行（可能一 key 多行）
pd_chars <- function(vals, keys) {
  out <- list()
  for (k in unique(keys)) if (grepl("^characteristics", k)) for (r in vals[keys == k]) out[[length(out) + 1L]] <- r
  out
}
# 把 characteristics 行（"key: value"）整理成 sample x key 的 data.frame
chars_df <- function(gm) {
  n <- ncol(gm$mat)
  keys <- character(0); store <- list()
  for (r in gm$chars) {
    if (length(r) != n) next
    k <- sub(":.*$", "", r[1])
    if (!grepl(":", r[1])) next
    if (!k %in% keys) { keys <- c(keys, k); store[[k]] <- rep(NA_character_, n) }
    store[[k]] <- ifelse(grepl(":", r), sub("^[^:]*:\\s*", "", r), store[[k]])
  }
  d <- as.data.frame(store, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(d) <- colnames(gm$mat)
  d
}

log_rows <- list()
add_log <- function(...) log_rows[[length(log_rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# ============================================================ 1) GSE76019 ====
message("=== GSE76019 (Adrenocortical Cancer) ===")
acc <- "GSE76019"; gpl <- "GPL13158"; tok <- "EFS"
gm <- read_gsm_matrix(file.path(ROOT, "data/raw", "GSE76019_series_matrix.txt.gz"))
ex <- gm$mat
ch <- chars_df(gm)
gsm <- colnames(ex)
message("  samples=", ncol(ex), " probes=", nrow(ex), " char keys=", paste(colnames(ch), collapse = ","))
stopifnot(all(c("histology", "efs.time", "efs.event") %in% colnames(ch)))
# 只保留 ACC（本系列 34 例全为 ACC，防御性过滤）
is_acc <- ch$histology == "ACC"
if (!all(is_acc)) message("  non-ACC samples dropped: ", sum(!is_acc))
ex <- ex[, is_acc, drop = FALSE]; ch <- ch[is_acc, , drop = FALSE]; gsm <- gsm[is_acc]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
efs_t <- num(ch$efs.time); efs_s <- num(ch$efs.event)
n_ok <- sum(!is.na(efs_t) & !is.na(efs_s)); ev <- sum(efs_s == 1, na.rm = TRUE)
cat(sprintf("  GSE76019: %d samples | EFS usable %d (events %d) | time range %.3f-%.3f y | scale=%s\n",
            ncol(ex), n_ok, ev, min(efs_t, na.rm = TRUE), max(efs_t, na.rm = TRUE), scale_verdict(ex)))
stopifnot(n_ok >= 30)                                  # 松弛闸门硬约束
stopifnot(!anyNA(gsm))

sv <- data.frame(EFS_status = efs_s, EFS_time = efs_t,
                 stage = as.character(ch$Stage), histology = as.character(ch$histology),
                 row.names = gsm, stringsAsFactors = FALSE, check.names = FALSE)
colnames(sv) <- sanitize(colnames(sv))
stopifnot("efs_status" %in% colnames(sv), "efs_time" %in% colnames(sv))
saveRDS(sv, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))

ex_out <- data.frame(ID_REF = rownames(ex), ex, check.names = FALSE); rownames(ex_out) <- NULL
stopifnot(!anyDuplicated(ex_out$ID_REF))
saveRDS(ex_out, file.path(ROOT, "data/expr", paste0(acc, ".rds")))

ph <- data.frame(patient_id = sub(".*COG patient ", "COG", gm$pdata$title[is_acc]),
                 title = gm$pdata$title[is_acc], stage = as.character(ch$Stage),
                 row.names = gsm, stringsAsFactors = FALSE)
saveRDS(ph, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))

g <- readRDS(file.path(ROOT, "data/processed/gpl", paste0(gpl, ".rds")))
n_in <- mean(ex_out$ID_REF %in% rownames(g))
cat(sprintf("  probes in %s map: %.1f%%\n", gpl, 100 * n_in))

add_log(acc = acc, type = "Adrenocortical Cancer", gpl = gpl, token = tok,
        n_expr = ncol(ex), N_patients = n_ok, n_events = ev, scale = scale_verdict(ex),
        detail = sprintf("GSE76019 pediatric ACC (COG ARAR0332, PMID 27307598); %d/%d usable EFS pairs; gate relaxed to >=30 (N=%d < 50)",
                         n_ok, ncol(ex), n_ok))

# ============================================================ 2) GSE76039 ====
# Thyroid Cancer（低分化/未分化甲状腺癌 PDTC/ATC）
#   表达 : GEO GSE76039-GPL570（54675 探针 × 37 样本；该系列是 MSK JCI 2016 研究的
#          表达子集，GEO 里没有生存字段）
#   临床 : cBioPortal study thyroid_mskcc_2016（Poorly-Differentiated and Anaplastic
#          Thyroid Cancers, MSK, JCI 2016）——同一研究在 cBioPortal 的临床存档，
#          117 患者带 OS_MONTHS/OS_STATUS
#   连接键: GSE76039 !Sample_title（s_JF_thy_NNN_P）== cBioPortal SAMPLE_ID；
#          并用 CEL 文件名里的 6 位病理号（如 105028T）== cBioPortal OTHER_SAMPLE_ID 交叉核对
#          -> 37/37 全部匹配
#   患者级: 37 例中 35 例有可用 (OS_MONTHS, OS_STATUS)，29 个事件
#   ★ 35 < 旧闸门 50 -> 按松弛闸门 >=30 登记，Note 必须写明
message("=== GSE76039 (Thyroid Cancer) ===")
acc <- "GSE76039"; gpl <- "GPL570"; tok <- "OS"
gm <- read_gsm_matrix(file.path(ROOT, "data/raw", "GSE76039_series_matrix.txt.gz"))
ex <- gm$mat
gsm <- colnames(ex); ttl <- gm$pdata$title
stopifnot(length(ttl) == ncol(ex))
clin_p <- read_clin(file.path(ROOT, "data/suppl/thyroid_mskcc_2016_data_clinical_patient.txt.gz"))
clin_s <- read_clin(file.path(ROOT, "data/suppl/thyroid_mskcc_2016_data_clinical_sample.txt.gz"))
message("  expr=", paste(dim(ex), collapse = "x"), " | cBio patients=", nrow(clin_p),
        " samples=", nrow(clin_s))
stopifnot(all(c("PATIENT_ID", "OS_STATUS", "OS_MONTHS") %in% names(clin_p)))
stopifnot(all(c("SAMPLE_ID", "PATIENT_ID", "OTHER_SAMPLE_ID") %in% names(clin_s)))
# 连接键 1：title == SAMPLE_ID
m_s <- match(ttl, clin_s$SAMPLE_ID)
# 连接键 2（交叉核对）：CEL 文件名里的 <digits>T == OTHER_SAMPLE_ID
cel_tok <- sub(".*_(\\d+T)_.*", "\\1", gm$pdata$supplementary_file)
ok2 <- !is.na(m_s) & !is.na(cel_tok) & cel_tok == clin_s$OTHER_SAMPLE_ID[m_s]
cat(sprintf("  join title==SAMPLE_ID: %d/%d ; cross-check CEL token==OTHER_SAMPLE_ID: %d\n",
            sum(!is.na(m_s)), length(ttl), sum(ok2)))
stopifnot(all(!is.na(m_s)), sum(ok2) >= 30)
pid <- clin_s$PATIENT_ID[m_s]
p <- clin_p[match(pid, clin_p$PATIENT_ID), , drop = FALSE]
os_m <- suppressWarnings(as.numeric(as.character(p$OS_MONTHS)))
os_s <- ifelse(grepl("^1", p$OS_STATUS), 1, ifelse(grepl("^0", p$OS_STATUS), 0, NA))
n_ok <- sum(!is.na(os_m) & !is.na(os_s)); ev <- sum(os_s == 1, na.rm = TRUE)
cat(sprintf("  GSE76039: %d samples | OS usable %d (events %d) | scale=%s\n",
            ncol(ex), n_ok, ev, scale_verdict(ex)))
stopifnot(n_ok >= 30)
sv <- data.frame(OS_status = os_s, OS_time = os_m / 12,
                 age = suppressWarnings(as.numeric(as.character(p$AGE))),
                 sex = tolower(as.character(p$SEX)),
                 pdtc_definition = as.character(p$PDTC_DEFINITION),
                 m_stage = as.character(p$M_STAGE),
                 path_n_stage = as.character(p$PATH_N_STAGE),
                 path_t_stage = as.character(p$PATH_T_STAGE),
                 row.names = gsm, stringsAsFactors = FALSE, check.names = FALSE)
colnames(sv) <- sanitize(colnames(sv))
saveRDS(sv, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))
ex_out <- data.frame(ID_REF = rownames(ex), ex, check.names = FALSE); rownames(ex_out) <- NULL
stopifnot(!anyDuplicated(ex_out$ID_REF))
saveRDS(ex_out, file.path(ROOT, "data/expr", paste0(acc, ".rds")))
ph <- data.frame(patient_id = pid, title = ttl, sample_id = clin_s$SAMPLE_ID[m_s],
                 row.names = gsm, stringsAsFactors = FALSE)
saveRDS(ph, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
g <- readRDS(file.path(ROOT, "data/processed/gpl", paste0(gpl, ".rds")))
cat(sprintf("  probes in %s map: %.1f%%\n", gpl, 100 * mean(ex_out$ID_REF %in% rownames(g))))
add_log(acc = acc, type = "Thyroid Cancer", gpl = gpl, token = tok,
        n_expr = ncol(ex), N_patients = n_ok, n_events = ev, scale = scale_verdict(ex),
        detail = sprintf("GSE76039 expression (MSK PDTC/ATC JCI 2016) + cBioPortal thyroid_mskcc_2016 clinical; join title==SAMPLE_ID, cross-checked by CEL token==OTHER_SAMPLE_ID (37/37); %d/%d usable OS; gate relaxed to >=30",
                         n_ok, ncol(ex)))

saveRDS(log_rows, cpas_out("relaxed_gate_build_log.rds"))
lg <- do.call(rbind, log_rows)
write.csv(lg, cpas_out("relaxed_gate_build_log.csv"), row.names = FALSE)
cat("\n--- build log ---\n"); print(lg[, c("acc", "type", "token", "n_expr", "N_patients", "n_events", "scale")], row.names = FALSE)
cat("OK\n")
}

# ---------------------------------------------------------------------------
# run_104_build_embl_cohorts()  <-  verbatim pipeline/R/104_build_embl_cohorts.R
# ---------------------------------------------------------------------------
run_104_build_embl_cohorts <- function() {
# 104_build_embl_cohorts.R ----------------------------------------------------
# Build the EMBL-EBI (ArrayExpress / BioStudies) cohorts that passed the
# author's gate and the GEO-independence check (embl_overlap_recheck.md).
#
# Products, per cohort (identical conventions to 90/92/93/100):
#   data/expr/<ACC>.rds                 ID_REF + one column per sample (probe level,
#                                       no probe collapsing - same as GSE39582/GSE108474)
#   data/processed/surv/<ACC>_surv.rds  <TOK>_status(0/1) + <TOK>_time(YEARS) + clinical
#   data/pheno/<ACC>.rds                patient-level clinical (sample rows = patients)
#
# Patient-level only: every table is de-duplicated to distinct patients.
# log2 is applied ONLY when the matrix is genuinely linear; the verdict is printed
# and recorded. `18_log2_transform.R --apply` is never invoked.
#
# Usage: Rscript pipeline/R/104_build_embl_cohorts.R <ACC|all> [ROOT]
# ----------------------------------------------------------------------------
suppressMessages({ library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
WANT <- if (length(args) >= 1) args[1] else "all"
ROOT <- if (length(args) >= 2) args[2] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
SUP  <- file.path(ROOT, "data/suppl/embl_build")
dir.create(file.path(ROOT, "data/expr"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/processed/surv"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/pheno"), showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ---------- shared helpers ---------------------------------------------------
read_sdrf <- function(acc) {
  for (p in c(file.path(ROOT, "data/suppl/embi_scan/AUDIT/sdrf", paste0(acc, ".sdrf.txt")),
              file.path(ROOT, "data/suppl/embi_scan/downloads", paste0(acc, ".sdrf.txt")),
              file.path(SUP, acc, paste0(acc, ".sdrf.txt")))) {
    if (file.exists(p)) {
      x <- read.delim(p, sep = "\t", header = TRUE, check.names = FALSE, quote = "",
                      stringsAsFactors = FALSE, colClasses = "character")
      names(x) <- trimws(names(x))
      # SDRFs repeat column names (e.g. several `Unit[time unit]` columns); make them
      # unique so `x[[name]]` cannot silently return the wrong duplicate.
      names(x) <- make.unique(names(x), sep = "__")
      return(x)
    }
  }
  stop(acc, ": SDRF not found")
}

# unnamed helper: pick the first column whose name matches any pattern
pick_col <- function(df, pats, required = TRUE) {
  nm <- names(df)
  for (p in pats) {
    i <- grep(p, nm, ignore.case = TRUE)
    if (length(i)) return(nm[i[1]])
  }
  if (required) stop("no column matching: ", paste(pats, collapse = " | "),
                     " | have: ", paste(nm, collapse = ", "))
  NA_character_
}

# strict first-number extraction (a naive gsub keeps 'e' from words like "years");
# also resolves the European decimal comma used by E-MTAB-4032 ("55,7" = 55.7) and
# space thousands separators ("1 694" = 1694) without breaking "2.349 years".
num <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  eu <- grepl("^[+-]?[0-9]+,[0-9]{1,2}$", x)
  x[eu] <- sub(",", ".", x[eu])
  x <- gsub("(?<=[0-9])[ \u00a0](?=[0-9]{3}(\\D|$))", "", x, perl = TRUE)
  out <- suppressWarnings(as.numeric(x))
  pat <- "[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?"
  bad <- is.na(out) & !is.na(x) & nzchar(x) & x != "NA"
  if (any(bad)) {
    xx <- x[bad]
    pos <- regexpr(pat, xx)
    mm <- rep(NA_character_, length(xx))
    okm <- pos > 0
    if (any(okm)) mm[okm] <- regmatches(xx, pos)
    out[bad] <- suppressWarnings(as.numeric(mm))
  }
  out
}

# scale verdict - same rule as 18_log2_transform.R / 92_build's scale_verdict()
scale_verdict <- function(m) {
  v <- as.matrix(m); storage.mode(v) <- "double"; v <- v[is.finite(v)]
  if (!length(v)) return("empty")
  if (min(v) < 0) return("already_log")
  if (median(v) > 0 && max(v) / median(v) <= 50) return("already_log")
  if (median(v) <= 0 && max(v) <= 25) return("already_log")
  "linear_needs_log2"
}

# matrix (probes x samples) -> data/expr/<ACC>.rds ; applies log2 only if linear
write_expr <- function(acc, m, note = "") {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  if (anyDuplicated(rownames(m))) m <- m[!duplicated(rownames(m)), , drop = FALSE]
  m <- m[!is.na(rownames(m)) & nzchar(rownames(m)), , drop = FALSE]
  v0 <- scale_verdict(m)
  if (v0 == "linear_needs_log2") {
    if (min(m, na.rm = TRUE) < 0) stop(acc, ": verdict linear but negative values - refusing")
    m <- log2(m + 1)
    log_msg(acc, ": scale=linear -> applied log2(x+1)")
  } else log_msg(acc, ": scale=", v0, " -> no transform")
  out <- data.frame(ID_REF = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  saveRDS(out, file.path(ROOT, "data/expr", paste0(acc, ".rds")))
  v <- m[is.finite(m)]
  log_msg(acc, ": expr ", nrow(m), " probes x ", ncol(m), " samples | range ",
          sprintf("%.3f-%.3f", min(v), max(v)), " median ", sprintf("%.3f", median(v)))
  invisible(out)
}

# patient-level surv table; rownames = patient id (must match expr column names)
write_surv <- function(acc, df, token, idcol = NULL) {
  stopifnot(token %in% c("OS", "DSS", "DFS", "RFS", "PFS", "MFS", "EFS"))
  st <- paste0(token, "_status"); tm <- paste0(token, "_time")
  if (!st %in% names(df) || !tm %in% names(df)) stop(acc, ": need ", st, " and ", tm)
  ids <- if (is.null(idcol)) rownames(df) else as.character(df[[idcol]])
  keep_c <- setdiff(names(df), c(st, tm))
  out <- data.frame(df[, keep_c, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  out[[st]] <- as.integer(df[[st]]); out[[tm]] <- as.numeric(df[[tm]])
  # endpoint columns first, matching the existing surv files
  out <- out[, c(st, tm, keep_c), drop = FALSE]
  rownames(out) <- ids
  ok <- !is.na(out[[st]]) & !is.na(out[[tm]]) & is.finite(out[[tm]]) &
        out[[st]] %in% c(0L, 1L) & out[[tm]] > 0
  out <- out[ok, , drop = FALSE]
  saveRDS(out, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))
  log_msg(acc, ": surv ", nrow(out), " patients | ", token, " events ",
          sum(out[[st]] == 1L), " | time(y) ", sprintf("%.2f-%.2f", min(out[[tm]]), max(out[[tm]])))
  invisible(out)
}

write_pheno <- function(acc, df) {
  saveRDS(df, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
}

# two-row MAGE-TAB matrix reader (row1 = sample ids, row2 = value descriptor, col1 = probe)
read_magetab_matrix <- function(path, value_col = NULL) {
  ln <- readLines(path, n = 2, warn = FALSE)
  h1 <- strsplit(ln[1], "\t")[[1]]; h2 <- strsplit(ln[2], "\t")[[1]]
  ids <- h1[-1]; desc <- h2[-1]
  want <- if (is.null(value_col)) seq_along(ids) else which(desc %in% value_col)
  if (!length(want)) stop("no columns match value_col=", paste(value_col, collapse = ","),
                          " | descriptors: ", paste(unique(desc), collapse = " / "))
  cols <- c(1L, want + 1L)
  dt <- data.table::fread(path, sep = "\t", header = FALSE, skip = 2, select = cols,
                          data.table = FALSE, check.names = FALSE, quote = "",
                          colClasses = "character")
  rn <- dt[[1]]; m <- as.data.frame(lapply(dt[, -1, drop = FALSE], num), check.names = FALSE)
  rownames(m) <- rn
  colnames(m) <- ids[want]
  # a sample can own several value columns -> collapse by rowMeans if duplicated
  if (anyDuplicated(colnames(m))) {
    m <- t(vapply(split(seq_len(ncol(m)), colnames(m)), function(j)
      rowMeans(as.matrix(m[, j, drop = FALSE]), na.rm = TRUE), numeric(nrow(m))))
    m <- as.data.frame(m, check.names = FALSE)
  }
  m
}

# processed matrices whose header row lists the sample names but has NO placeholder
# for the leading probe-id column (E-MTAB-4032 training/validation, E-MTAB-6134)
read_unlabelled_tsv <- function(path, tag = basename(path)) {
  h <- strsplit(readLines(path, n = 1, warn = FALSE), "\t")[[1]]
  h <- gsub('^"|"$', "", h)
  x <- data.table::fread(path, sep = "\t", header = FALSE, skip = 1, data.table = FALSE,
                         check.names = FALSE, quote = "\"", colClasses = "character")
  ids <- as.character(x[[1]])
  m <- as.data.frame(lapply(x[, -1, drop = FALSE], num), check.names = FALSE)
  if (length(h) == ncol(m) + 1L) h <- h[-1]      # header carries an id placeholder
  if (length(h) != ncol(m)) stop(tag, ": header fields ", length(h),
                                 " != value columns ", ncol(m))
  colnames(m) <- h
  if (anyDuplicated(ids)) {
    log_msg(tag, ": ", sum(duplicated(ids)), " duplicated probe rows -> de-duplicated")
    keep <- !duplicated(ids); m <- m[keep, , drop = FALSE]; ids <- ids[keep]
  }
  rownames(m) <- ids
  m
}

# ============================== cohort builders ==============================

build_E_MEXP_2780 <- function() {
  acc <- "E-MEXP-2780"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "final-gene-expression-data-file.txt.magetab")
  stopifnot(file.exists(f))
  m <- read_magetab_matrix(f, value_col = "RMA")
  sid <- pick_col(s, "^Source Name$"); tt <- pick_col(s, "TimeOfSurvival")
  ss <- pick_col(s, "OverallSurvival")
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex   <- s[[pick_col(s, "Sex", FALSE)]]
  surv$disease <- s[[pick_col(s, "DiseaseState", FALSE)]]
  surv$OS_status <- ifelse(tolower(s[[ss]]) %in% c("dead", "1", "yes"), 1L,
                    ifelse(tolower(s[[ss]]) %in% c("alive", "0", "no"), 0L, NA_integer_))
  surv$OS_time <- num(s[[tt]]) / 365.25          # SDRF unit = day
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

build_E_MTAB_1719 <- function() {
  acc <- "E-MTAB-1719"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "MESO-Jaurand-EXP-RMA-normalized-data.txt")
  stopifnot(file.exists(f))
  m <- read_magetab_matrix(f, value_col = "log2 RMA")
  sid <- pick_col(s, "^Source Name$"); tt <- pick_col(s, "OS delay"); ss <- pick_col(s, "OS event")
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "Sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "Age", FALSE)]])
  surv$histology <- s[[pick_col(s, "histology", FALSE)]]
  surv$OS_status <- ifelse(num(s[[ss]]) == 1, 1L, ifelse(num(s[[ss]]) == 0, 0L, NA_integer_))
  surv$OS_time <- num(s[[tt]]) / 12              # SDRF unit = month
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

build_E_MTAB_1205 <- function() {
  acc <- "E-MTAB-1205"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "cog.exp.meanadj2.txt")
  stopifnot(file.exists(f))
  m <- read_magetab_matrix(f)
  sid <- pick_col(s, "^Source Name$")
  tt <- pick_col(s, "\\[Relapse-Free Survival\\]$"); ss <- pick_col(s, "Relapse-Free Survival Status")
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "Sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "Age at Diagnosis", FALSE)]])
  surv$disease_state <- s[[pick_col(s, "DiseaseState", FALSE)]]
  surv$wbc_at_diagnosis <- num(s[[pick_col(s, "White Blood Cell Count", FALSE)]])
  surv$RFS_status <- ifelse(num(s[[ss]]) == 1, 1L, ifelse(num(s[[ss]]) == 0, 0L, NA_integer_))
  surv$RFS_time <- num(s[[tt]]) / 365.25         # SDRF unit = day
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "RFS"); write_pheno(acc, surv)
}

build_E_TABM_346 <- function() {
  acc <- "E-TABM-346"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "E-TABM-346-processed-data-2496486478.txt")
  stopifnot(file.exists(f))
  m <- read_magetab_matrix(f, value_col = "Affymetrix:CHPSignal")   # 53 samples x 6 fields
  sid <- pick_col(s, "^Source Name$")
  # the matrix is keyed by the SDRF `Scan Name` (KL_B8_U133A), not by Source Name
  sc <- pick_col(s, "Scan Name", required = FALSE)
  if (!is.na(sc)) {
    map <- setNames(s[[sid]], s[[sc]])
    hit <- intersect(colnames(m), names(map))
    if (length(hit) >= 0.9 * nrow(s)) { colnames(m) <- unname(map[colnames(m)]) }
  }
  gv <- function(pat) {
    c <- pick_col(s, pat, required = FALSE); if (is.na(c)) return(rep(NA_character_, nrow(s)))
    sub("^[^:]*:\\s*", "", s[[c]])
  }
  os  <- gv("Overall survival"); ds <- gv("Death status")
  efs <- gv("Event free survival"); ps <- gv("Progression status")
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "Sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "Age", FALSE)]])
  surv$disease_state <- s[[pick_col(s, "DiseaseState", FALSE)]]
  surv$clinical_treatment <- s[[pick_col(s, "ClinicalTreatment", FALSE)]]
  surv$tumor_grading <- num(s[[pick_col(s, "TumorGrading", FALSE)]])
  surv$OS_status  <- as.integer(num(ds))
  surv$OS_time    <- num(os)            # "2.349 years"
  surv$PFS_status <- as.integer(num(ps))
  surv$PFS_time   <- as.numeric(NA)
  surv$EFS_time   <- num(efs)
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m)
  write_surv(acc, surv, "OS")
  write_pheno(acc, surv)
}

build_E_MTAB_3580 <- function() {
  acc <- "E-MTAB-3580"; s <- read_sdrf(acc)
  d <- file.path(SUP, acc)
  rcc <- list.files(d, pattern = "\\.RCC$", full.names = TRUE)
  stopifnot(length(rcc) > 0)
  read_rcc <- function(p) {
    ln <- readLines(p, warn = FALSE)
    i <- grep("^<Code_Summary>", ln); j <- grep("^</Code_Summary>", ln)
    if (!length(i) || !length(j)) stop("no Code_Summary in ", p)
    body <- ln[(i[1] + 1):(j[1] - 1)]
    body <- body[nzchar(body)]
    hdr <- strsplit(body[1], ",")[[1]]
    dt <- read.csv(text = paste(body[-1], collapse = "\n"), header = FALSE,
                   col.names = hdr, stringsAsFactors = FALSE, check.names = FALSE)
    dt <- dt[tolower(dt$CodeClass) == "endogenous", c("Name", "Count")]
    v <- as.numeric(dt$Count); names(v) <- dt$Name
    v
  }
  lst <- lapply(rcc, read_rcc)
  genes <- Reduce(union, lapply(lst, names))
  m <- matrix(NA_real_, nrow = length(genes), ncol = length(lst),
              dimnames = list(genes, basename(rcc)))
  for (k in seq_along(lst)) m[names(lst[[k]]), k] <- lst[[k]]
  # map each RCC to its SDRF row through `Array Data File`, then to the PATIENT id
  adf <- pick_col(s, "Array Data File", required = FALSE)
  sid <- pick_col(s, "^Source Name$")
  # The SDRF `Source Name` here is the generic row label `SampleN`; the patient /
  # RNA identifier is `Characteristics[sampleid]` (68 distinct for 84 rows - the
  # remaining 16 rows are repeat nCounter runs of the same RNA with identical
  # time+event). The audit's "84 patients" is a row-count artefact, so this
  # builder de-duplicates to the true patient level.
  pid_col <- pick_col(s, "sampleid", required = FALSE)
  if (is.na(pid_col)) pid_col <- sid
  idx <- if (!is.na(adf)) match(basename(rcc), s[[adf]]) else NA_integer_
  if (any(is.na(idx))) {
    stem <- sub("\\.RCC$", "", basename(rcc))
    idx2 <- match(stem, sub("\\.CEL$", "", s[[adf]]))
    idx[is.na(idx)] <- idx2[is.na(idx)]
  }
  if (any(is.na(idx))) stop(acc, ": ", sum(is.na(idx)), " RCC files not matched to the SDRF")
  colnames(m) <- s[[pid_col]][idx]
  if (anyDuplicated(colnames(m))) {                     # collapse repeat runs
    mm <- vapply(split(seq_len(ncol(m)), colnames(m)), function(j)
      rowMeans(as.matrix(m[, j, drop = FALSE]), na.rm = TRUE), numeric(nrow(m)))
    m <- as.data.frame(mm, check.names = FALSE)         # genes x patients
  }
  tt <- pick_col(s, "overall survival time"); ss <- pick_col(s, "overall survival event")
  surv <- data.frame(pid = s[[pid_col]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "age", FALSE)]])
  surv$histology <- s[[pick_col(s, "histology", FALSE)]]
  surv$irs_stage <- s[[pick_col(s, "irs stage", FALSE)]]
  surv$OS_status <- as.integer(num(s[[ss]]))
  surv$OS_time   <- num(s[[tt]])
  surv <- surv[!duplicated(surv$pid), , drop = FALSE]        # repeat rows agree
  rownames(surv) <- surv$pid; surv$pid <- NULL
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

# CEL-name key: drop punctuation and leading zeros *per numeric field only*
# (E-MTAB-4032's processed matrices write MM-10-002 as MM-10-02, MM-11-027 as MM-11-27)
cel_key <- function(s) {
  s <- tolower(as.character(s))
  parts <- regmatches(s, gregexpr("[0-9]+|[^0-9]+", s))
  vapply(parts, function(p) {
    p <- ifelse(grepl("^[0-9]+$", p), as.character(as.numeric(p)), gsub("[^a-z]", "", p))
    paste(p, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

build_E_MTAB_4032 <- function() {
  acc <- "E-MTAB-4032"; s <- read_sdrf(acc)
  d <- file.path(SUP, acc)
  rd <- function(f) read_unlabelled_tsv(file.path(d, f), paste0(acc, "/", f))
  tr <- rd("training_set.txt"); va <- rd("validation_set.txt")
  common <- intersect(rownames(tr), rownames(va))
  m <- cbind(tr[common, , drop = FALSE], va[common, , drop = FALSE])
  log_msg(acc, ": matrices merged on ", length(common), " common probes (",
          ncol(tr), " + ", ncol(va), " samples)")
  adf <- pick_col(s, "Array Data File")
  sid <- pick_col(s, "^Source Name$")
  key2row <- setNames(seq_len(nrow(s)), cel_key(s[[adf]]))
  idx <- unname(key2row[cel_key(colnames(m))])
  if (any(is.na(idx))) stop(acc, ": ", sum(is.na(idx)), " matrix columns unmatched to the SDRF")
  colnames(m) <- s[[sid]][idx]
  if (anyDuplicated(colnames(m))) stop(acc, ": duplicate patient column after mapping")
  # OS column is labelled |[years]| but is really MONTHS - proven against E-MTAB-1038
  # (73 shared patients, ratio 4032/1038 == 1/(365.25/12) == 1/30.4375, i.e. months vs days)
  co <- pick_col(s, "from diagnosis"); ca <- pick_col(s, "patient alive")
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "age", FALSE)]])
  surv$iss_stage <- s[[pick_col(s, "iss stage", FALSE)]]
  surv$treatment_line <- s[[pick_col(s, "treatment line", FALSE)]]
  surv$cohort_set <- s[[pick_col(s, "sample set", FALSE)]]
  surv$ttp_months <- num(s[[pick_col(s, "ttp", FALSE)]])
  surv$OS_status <- ifelse(tolower(s[[ca]]) == "no", 1L,
                    ifelse(tolower(s[[ca]]) == "yes", 0L, NA_integer_))
  surv$OS_time <- num(s[[co]]) / 12          # months -> years
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

# ---- map a set of sample-level ids (CEL / matrix columns) onto SDRF Source Names
map_samples <- function(keys, s, adf, sid) {
  idx <- match(keys, s[[adf]])
  if (anyNA(idx)) {
    idx[is.na(idx)] <- match(cel_key(keys)[is.na(idx)], cel_key(s[[adf]]))
  }
  if (anyNA(idx)) {
    idx[is.na(idx)] <- match(keys, s[[sid]])[is.na(idx)]
  }
  if (anyNA(idx)) stop("unmatched sample keys: ", paste(utils::head(keys[is.na(idx)], 5), collapse = ", "))
  s[[sid]][idx]
}

# Affymetrix CEL cohort: RMA with the platform CDF (affy + hgu133plus2cdf for GPL570)
build_rma_affy <- function(acc) {
  s <- read_sdrf(acc)
  d <- file.path(SUP, acc)
  cels <- list.files(d, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)
  if (!length(cels)) stop(acc, ": no CEL files under ", d)
  log_msg(acc, ": RMA over ", length(cels), " CEL files (affy)")
  suppressMessages(library(affy))
  ab <- affy::ReadAffy(filenames = cels)
  m <- Biobase::exprs(affy::rma(ab))
  adf <- pick_col(s, "Array Data File", required = FALSE)
  sid <- pick_col(s, "^Source Name$")
  colnames(m) <- map_samples(basename(cels), s, adf, sid)
  list(expr = m, sdrf = s, sid = sid)
}

# Affymetrix ST cohort: RMA with oligo + the pd annotation package
build_rma_oligo <- function(acc, pd_pkg) {
  s <- read_sdrf(acc)
  d <- file.path(SUP, acc)
  cels <- list.files(d, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)
  if (!length(cels)) stop(acc, ": no CEL files under ", d)
  log_msg(acc, ": RMA over ", length(cels), " CEL files (oligo / ", pd_pkg, ")")
  suppressMessages({ library(oligo); library(Biobase) })
  ab <- oligo::read.celfiles(cels, pkgname = pd_pkg, verbose = FALSE)
  m <- Biobase::exprs(oligo::rma(ab, target = "core"))
  adf <- pick_col(s, "Array Data File", required = FALSE)
  sid <- pick_col(s, "^Source Name$")
  colnames(m) <- map_samples(basename(cels), s, adf, sid)
  list(expr = m, sdrf = s, sid = sid)
}

# the Unit column that immediately follows `value_col` (SDRFs repeat `Unit[time unit]`,
# and the first one belongs to a different characteristic - E-MTAB-6134 col 13 is the
# unit of `ffpeblock age` = "year", while col 17 is the real unit of `os.delay` = "month")
unit_for <- function(s, value_col) {
  i <- match(value_col, names(s))
  if (is.na(i) || i >= ncol(s)) return(rep(NA_character_, nrow(s)))
  nm <- names(s)[i + 1L]
  if (!is.na(nm) && grepl("unit", nm, ignore.case = TRUE)) return(s[[i + 1L]])
  rep(NA_character_, nrow(s))
}

# time column -> years, driven by the SDRF's own Unit column (default_unit otherwise)
to_years <- function(v, unit, default_unit = "day") {
  u <- if (is.null(unit)) rep(default_unit, length(v)) else tolower(unit)
  u[is.na(u) | !nzchar(trimws(u))] <- default_unit
  out <- v
  out[grepl("year", u)] <- out[grepl("year", u)]
  out[grepl("month", u)] <- out[grepl("month", u)] / 12
  out[grepl("day", u)]   <- out[grepl("day", u)] / 365.25
  out[grepl("week", u)]  <- out[grepl("week", u)] / 52.1786
  out
}

build_E_MTAB_6134 <- function() {
  acc <- "E-MTAB-6134"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "ProcessedExpression.tsv")
  stopifnot(file.exists(f))
  m <- read_unlabelled_tsv(f, acc)
  sid <- pick_col(s, "^Source Name$")
  hit <- intersect(colnames(m), s[[sid]])
  if (length(hit) < 0.9 * nrow(s)) stop(acc, ": only ", length(hit), " matrix columns match Source Name")
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$disease <- s[[pick_col(s, "disease", FALSE)]]
  surv$tumor_grading <- s[[pick_col(s, "tumor grading", FALSE)]]
  surv$tnm <- s[[pick_col(s, "TNM", FALSE)]]
  surv$resection_margin <- s[[pick_col(s, "resection margin", FALSE)]]
  surv$clinical_center <- s[[pick_col(s, "clinical center", FALSE)]]
  dc <- pick_col(s, "dfs.delay", FALSE)
  surv$dfs_time <- to_years(num(s[[dc]]), unit_for(s, dc), "month")
  surv$dfs_event <- suppressWarnings(as.integer(num(s[[pick_col(s, "dfs.event", FALSE)]])))
  surv$OS_status <- suppressWarnings(as.integer(num(s[[pick_col(s, "os.event")]])))
  oc <- pick_col(s, "os.delay")
  surv$OS_time <- to_years(num(s[[oc]]), unit_for(s, oc), "month")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

build_E_MTAB_6389 <- function() {
  acc <- "E-MTAB-6389"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "data_exp_icc.txt")
  stopifnot(file.exists(f))
  m <- read_unlabelled_tsv(f, acc)
  sid <- pick_col(s, "^Source Name$")
  ind <- pick_col(s, "individual", required = FALSE)
  if (is.na(ind)) ind <- sid
  if (!length(intersect(colnames(m), s[[sid]]))) stop(acc, ": matrix columns match no Source Name")
  pid <- s[[ind]][match(colnames(m), s[[sid]])]
  keep <- !is.na(pid); m <- m[, keep, drop = FALSE]; pid <- pid[keep]
  if (anyDuplicated(pid)) {
    log_msg(acc, ": ", sum(duplicated(pid)), " duplicate patient columns (tumour/normal) -> collapsed")
    mm <- vapply(split(seq_len(ncol(m)), pid), function(j)
      rowMeans(as.matrix(m[, j, drop = FALSE]), na.rm = TRUE), numeric(nrow(m)))
    m <- as.data.frame(mm, check.names = FALSE); pid <- colnames(mm)
  }
  colnames(m) <- pid
  surv <- data.frame(row.names = s[[ind]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$disease <- s[[pick_col(s, "disease", FALSE)]]
  surv$organism_part <- s[[pick_col(s, "organism part", FALSE)]]
  surv$sampling_site <- s[[pick_col(s, "sampling site", FALSE)]]
  surv$cirrhosis <- s[[pick_col(s, "cirrhosis", FALSE)]]
  surv$vascular_invasion <- s[[pick_col(s, "vascular invasion", FALSE)]]
  surv$OS_status <- as.integer(num(s[[pick_col(s, "event death")]]))
  oc <- pick_col(s, "overall survival")
  surv$OS_time <- to_years(num(s[[oc]]), unit_for(s, oc), "month")
  surv <- surv[!duplicated(rownames(surv)), , drop = FALSE]
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

build_E_MTAB_4321 <- function() {
  acc <- "E-MTAB-4321"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "UROMOL_gene_fpkm_gtf.txt")
  stopifnot(file.exists(f))
  # 4 leading annotation columns (tracking_id / gene.type / gene.status / gene.name)
  # followed by one column per sample
  hdr <- strsplit(readLines(f, n = 1, warn = FALSE), "\t")[[1]]
  nann <- max(grep("^(tracking_id|gene\\.type|gene\\.status|gene\\.name)$", hdr))
  x <- data.table::fread(f, sep = "\t", header = TRUE, data.table = FALSE,
                         check.names = FALSE, quote = "",
                         select = c(1L, (nann + 1L):length(hdr)))
  m <- as.data.frame(lapply(x[, -1, drop = FALSE], num), check.names = FALSE)
  rownames(m) <- as.character(x[[1]])
  colnames(m) <- hdr[(nann + 1L):length(hdr)]
  log_msg(acc, ": FPKM matrix ", nrow(m), " genes x ", ncol(m), " samples")
  sid <- pick_col(s, "^Source Name$")
  if (!length(intersect(colnames(m), s[[sid]]))) {
    j <- match(cel_key(colnames(m)), cel_key(s[[sid]]))
    if (anyNA(j)) stop(acc, ": ", sum(is.na(j)), " matrix columns unmatched")
    colnames(m) <- s[[sid]][j]
  }
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "age", FALSE)]])
  surv$disease <- s[[pick_col(s, "disease", FALSE)]]
  surv$tumor_grading <- s[[pick_col(s, "tumor grading", FALSE)]]
  surv$disease_staging <- s[[pick_col(s, "disease staging", FALSE)]]
  surv$eortc_risk_score <- s[[pick_col(s, "EORTC risk score", FALSE)]]
  surv$bcg_treatment <- s[[pick_col(s, "BCG treatment", FALSE)]]
  surv$cystectomy <- s[[pick_col(s, "cystectomy", FALSE)]]
  ev <- tolower(s[[pick_col(s, "progression to T2")]])
  surv$PFS_status <- ifelse(grepl("^yes", ev), 1L, ifelse(grepl("^no", ev), 0L,
                     suppressWarnings(as.integer(num(ev)))))
  pc <- pick_col(s, "progression free survival")
  surv$PFS_time <- to_years(num(s[[pc]]), unit_for(s, pc), "month")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "PFS"); write_pheno(acc, surv)
}

build_E_TABM_1202 <- function() {
  acc <- "E-TABM-1202"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "E-TABM-1202.eSet.r")
  stopifnot(file.exists(f))
  suppressMessages({ library(Biobase); library(affy) })
  e <- new.env(); load(f, envir = e)
  obj <- get(ls(e)[1], envir = e)
  if (inherits(obj, "AffyBatch")) {
    # the deposit ships the *raw* AffyBatch (A-AFFY-44 = GPL570, cdfName HG-U133_Plus_2)
    log_msg(acc, ": eSet.r holds an AffyBatch (", length(sampleNames(obj)),
            " CELs, cdf ", obj@cdfName, ") -> RMA with hgu133plus2cdf")
    es <- affy::rma(obj)
  } else if (inherits(obj, "ExpressionSet")) {
    es <- obj
  } else stop(acc, ": unsupported object class in ", f, ": ", paste(class(obj), collapse = ","))
  m <- Biobase::exprs(es)
  log_msg(acc, ": ", nrow(m), " probes x ", ncol(m), " samples")
  sid <- pick_col(s, "^Source Name$"); adf <- pick_col(s, "Array Data File", required = FALSE)
  if (length(intersect(colnames(m), s[[sid]])) < 0.9 * nrow(s)) {
    # eSet sampleNames == SDRF `Scan Name`; `Array Data File` is the same stem + .CEL
    sc <- pick_col(s, "Scan Name", required = FALSE)
    mp <- if (!is.na(sc)) setNames(s[[sid]], s[[sc]]) else NULL
    if (!is.null(mp) && mean(colnames(m) %in% names(mp)) > 0.9) colnames(m) <- unname(mp[colnames(m)])
    else if (!is.na(adf)) colnames(m) <- map_samples(colnames(m), s, adf, sid)
    else { j <- match(cel_key(colnames(m)), cel_key(s[[sid]])); colnames(m) <- s[[sid]][j] }
  }
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "Sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "Age", FALSE)]])
  surv$histology <- s[[pick_col(s, "Histology", FALSE)]]
  surv$clinical_group <- s[[pick_col(s, "ClinicalGroup", FALSE)]]
  surv$irs_tnm_stage <- s[[pick_col(s, "IRS_TNM_Stage", FALSE)]]
  surv$metastasis <- s[[pick_col(s, "Metastasis", FALSE)]]
  surv$lymph_node <- s[[pick_col(s, "LymphNode", FALSE)]]
  surv$chromosomal_translocation <- s[[pick_col(s, "chromosomal_translocation", FALSE)]]
  st <- tolower(s[[pick_col(s, "\\[Status\\]")]])
  surv$OS_status <- ifelse(grepl("dead", st), 1L, ifelse(grepl("alive", st), 0L,
                     suppressWarnings(as.integer(num(st)))))
  tc <- pick_col(s, "survival_time")
  surv$OS_time <- to_years(num(s[[tc]]), unit_for(s, tc), "year")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

build_E_MTAB_3892 <- function() {
  acc <- "E-MTAB-3892"; s <- read_sdrf(acc)
  r <- build_rma_affy(acc); m <- r$expr; s <- r$sdrf; sid <- r$sid
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "age", FALSE)]])
  surv$histology <- s[[pick_col(s, "histology", FALSE)]]
  surv$histology_grade <- s[[pick_col(s, "histology_grade", FALSE)]]
  surv$disease_staging <- s[[pick_col(s, "disease staging", FALSE)]]
  surv$codeletion_1p19q <- s[[pick_col(s, "1p/19q", FALSE)]]
  surv$idh1_mutation <- s[[pick_col(s, "idh1", FALSE)]]
  surv$OS_status <- suppressWarnings(as.integer(num(s[[pick_col(s, "os.event")]])))
  oc <- pick_col(s, "os.delay")
  surv$OS_time <- to_years(num(s[[oc]]), unit_for(s, oc), "month")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

build_E_MTAB_1727 <- function() {
  acc <- "E-MTAB-1727"; s <- read_sdrf(acc)
  r <- build_rma_affy(acc); m <- r$expr; s <- r$sdrf; sid <- r$sid
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "Sex", FALSE)]]
  surv$disease <- s[[pick_col(s, "Disease", FALSE)]]
  surv$tnm <- s[[pick_col(s, "TNM", FALSE)]]
  surv$disease_staging <- s[[pick_col(s, "disease staging", FALSE)]]
  surv$scc_histology <- s[[pick_col(s, "SCC Histology", FALSE)]]
  surv$tobacco_smoking <- s[[pick_col(s, "Tobacco Smoking", FALSE)]]
  surv$age_at_diagnosis <- num(s[[pick_col(s, "Age At Diagnosis", FALSE)]])
  rc <- pick_col(s, "Relaspe-Free Survival In Month", FALSE)
  surv$rfs_time <- to_years(num(s[[rc]]), unit_for(s, rc), "month")
  surv$rfs_event <- suppressWarnings(as.integer(num(s[[pick_col(s, "Relaspe-Free Survival Event", FALSE)]])))
  surv$OS_status <- suppressWarnings(as.integer(num(s[[pick_col(s, "Overall Survival Event")]])))
  oc <- pick_col(s, "\\[Overall Survival\\]")
  surv$OS_time <- to_years(num(s[[oc]]), unit_for(s, oc), "month")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

build_E_MTAB_3267 <- function() {
  acc <- "E-MTAB-3267"; s <- read_sdrf(acc)
  r <- build_rma_oligo(acc, "pd.hugene.1.0.st.v1"); m <- r$expr; s <- r$sdrf; sid <- r$sid
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "age", FALSE)]])
  surv$histology_type <- s[[pick_col(s, "histology type", FALSE)]]
  surv$clinical_treatment <- s[[pick_col(s, "clinical treatment", FALSE)]]
  surv$sunitinib_response <- s[[pick_col(s, "sunitinib response", FALSE)]]
  surv$PFS_status <- suppressWarnings(as.integer(num(s[[pick_col(s, "\\[progression\\]")]])))
  pc <- pick_col(s, "progression free survival")
  surv$PFS_time <- to_years(num(s[[pc]]), unit_for(s, pc), "month")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "PFS"); write_pheno(acc, surv)
}

build_E_MTAB_6877 <- function() {
  acc <- "E-MTAB-6877"; s <- read_sdrf(acc)
  r <- build_rma_oligo(acc, "pd.hugene.2.0.st"); m <- r$expr; s <- r$sdrf; sid <- r$sid
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex <- s[[pick_col(s, "sex", FALSE)]]
  surv$age <- num(s[[pick_col(s, "age", FALSE)]])
  surv$disease <- s[[pick_col(s, "disease", FALSE)]]
  surv$histology <- s[[pick_col(s, "histology", FALSE)]]
  surv$tnm_stage <- s[[pick_col(s, "tnm.stage", FALSE)]]
  surv$asbestos_exposure <- s[[pick_col(s, "asbestos.exposure", FALSE)]]
  surv$OS_status <- as.integer(num(s[[pick_col(s, "event death")]]))
  oc <- pick_col(s, "overall survival")
  surv$OS_time <- to_years(num(s[[oc]]), unit_for(s, oc), "month")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

BUILDERS <- list(
  "E-MEXP-2780" = build_E_MEXP_2780,
  "E-MTAB-1719" = build_E_MTAB_1719,
  "E-MTAB-1205" = build_E_MTAB_1205,
  "E-TABM-346" = build_E_TABM_346,
  "E-MTAB-3580" = build_E_MTAB_3580,
  "E-MTAB-4032" = build_E_MTAB_4032,
  "E-MTAB-6134" = build_E_MTAB_6134,
  "E-MTAB-6389" = build_E_MTAB_6389,
  "E-MTAB-4321" = build_E_MTAB_4321,
  "E-TABM-1202" = build_E_TABM_1202,
  "E-MTAB-3892" = build_E_MTAB_3892,
  "E-MTAB-1727" = build_E_MTAB_1727,
  "E-MTAB-3267" = build_E_MTAB_3267,
  "E-MTAB-6877" = build_E_MTAB_6877
)

todo <- if (WANT == "all") names(BUILDERS) else strsplit(WANT, ",")[[1]]
for (a in todo) {
  if (!a %in% names(BUILDERS)) { log_msg("no builder for ", a, " (skipped)"); next }
  log_msg("=== building ", a, " ===")
  res <- tryCatch({ BUILDERS[[a]](); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
  log_msg(a, " -> ", res)
}
log_msg("done")
}

# ---------------------------------------------------------------------------
# run_106_build_embl_step2()  <-  verbatim pipeline/R/106_build_embl_step2.R
# ---------------------------------------------------------------------------
run_106_build_embl_step2 <- function() {
# 106_build_embl_step2.R ------------------------------------------------------
# Build the LAST THREE EMBL-EBI (ArrayExpress / BioStudies) cohorts:
#   E-MTAB-1727  Lung Cancer  OS  GPL570    (CEL-only -> RMA with affy)
#   E-MTAB-6389  Liver Cancer OS  GPL17585  (processed HTA-2.0 matrix)
#   E-MTAB-6877  Mesothelioma OS  GPL16686  (CEL-only -> RMA with oligo)
#
# Conventions are IDENTICAL to the precedents 92/93/100/104 (helpers copied verbatim
# from 104_build_embl_cohorts.R so the products are interchangeable):
#   data/expr/<ACC>.rds                 ID_REF + one column per sample (probe level,
#                                       de-duplicated probe rows, NO collapsing)
#   data/processed/surv/<ACC>_surv.rds  <TOK>_status(0/1) + <TOK>_time(YEARS) + clinical
#   data/pheno/<ACC>.rds                patient-level clinical
# log2 is applied ONLY when the matrix is genuinely linear (same scale_verdict rule as
# 18_log2_transform.R); `18_log2_transform.R --apply` is never invoked.
#
# Usage: Rscript pipeline/R/106_build_embl_step2.R <ACC|all> [ROOT]
# ----------------------------------------------------------------------------
suppressMessages({ library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
WANT <- if (length(args) >= 1) args[1] else "all"
ROOT <- if (length(args) >= 2) args[2] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
SUP  <- file.path(ROOT, "data/suppl/embl_step2")
dir.create(file.path(ROOT, "data/expr"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/processed/surv"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "data/pheno"), showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ---------- shared helpers (verbatim from 104_build_embl_cohorts.R) ----------
read_sdrf <- function(acc) {
  for (p in c(file.path(ROOT, "data/suppl/embi_scan/downloads", paste0(acc, ".sdrf.txt")),
              file.path(ROOT, "data/suppl/embi_scan/AUDIT/sdrf", paste0(acc, ".sdrf.txt")),
              file.path(SUP, acc, paste0(acc, ".sdrf.txt")))) {
    if (file.exists(p)) {
      x <- read.delim(p, sep = "\t", header = TRUE, check.names = FALSE, quote = "",
                      stringsAsFactors = FALSE, colClasses = "character")
      names(x) <- trimws(names(x))
      names(x) <- make.unique(names(x), sep = "__")
      return(x)
    }
  }
  stop(acc, ": SDRF not found")
}

pick_col <- function(df, pats, required = TRUE) {
  nm <- names(df)
  for (p in pats) {
    i <- grep(p, nm, ignore.case = TRUE)
    if (length(i)) return(nm[i[1]])
  }
  if (required) stop("no column matching: ", paste(pats, collapse = " | "),
                     " | have: ", paste(nm, collapse = ", "))
  NA_character_
}

num <- function(x) {
  x <- as.character(x); x <- trimws(x)
  eu <- grepl("^[+-]?[0-9]+,[0-9]{1,2}$", x)
  x[eu] <- sub(",", ".", x[eu])
  x <- gsub("(?<=[0-9])[ \u00a0](?=[0-9]{3}(\\D|$))", "", x, perl = TRUE)
  out <- suppressWarnings(as.numeric(x))
  pat <- "[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?"
  bad <- is.na(out) & !is.na(x) & nzchar(x) & x != "NA"
  if (any(bad)) {
    xx <- x[bad]; pos <- regexpr(pat, xx)
    mm <- rep(NA_character_, length(xx)); okm <- pos > 0
    if (any(okm)) mm[okm] <- regmatches(xx, pos)
    out[bad] <- suppressWarnings(as.numeric(mm))
  }
  out
}

scale_verdict <- function(m) {
  v <- as.matrix(m); storage.mode(v) <- "double"; v <- v[is.finite(v)]
  if (!length(v)) return("empty")
  if (min(v) < 0) return("already_log")
  if (median(v) > 0 && max(v) / median(v) <= 50) return("already_log")
  if (median(v) <= 0 && max(v) <= 25) return("already_log")
  "linear_needs_log2"
}

write_expr <- function(acc, m, note = "") {
  m <- as.matrix(m); storage.mode(m) <- "double"
  if (anyDuplicated(rownames(m))) m <- m[!duplicated(rownames(m)), , drop = FALSE]
  m <- m[!is.na(rownames(m)) & nzchar(rownames(m)), , drop = FALSE]
  v0 <- scale_verdict(m)
  if (v0 == "linear_needs_log2") {
    if (min(m, na.rm = TRUE) < 0) stop(acc, ": verdict linear but negative values - refusing")
    m <- log2(m + 1)
    log_msg(acc, ": scale=linear -> applied log2(x+1)")
  } else log_msg(acc, ": scale=", v0, " -> no transform")
  out <- data.frame(ID_REF = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  saveRDS(out, file.path(ROOT, "data/expr", paste0(acc, ".rds")))
  v <- m[is.finite(m)]
  log_msg(acc, ": expr ", nrow(m), " probes x ", ncol(m), " samples | range ",
          sprintf("%.3f-%.3f", min(v), max(v)), " median ", sprintf("%.3f", median(v)))
  invisible(out)
}

write_surv <- function(acc, df, token, idcol = NULL) {
  stopifnot(token %in% c("OS", "DSS", "DFS", "RFS", "PFS", "MFS", "EFS"))
  st <- paste0(token, "_status"); tm <- paste0(token, "_time")
  if (!st %in% names(df) || !tm %in% names(df)) stop(acc, ": need ", st, " and ", tm)
  ids <- if (is.null(idcol)) rownames(df) else as.character(df[[idcol]])
  keep_c <- setdiff(names(df), c(st, tm))
  out <- data.frame(df[, keep_c, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  out[[st]] <- as.integer(df[[st]]); out[[tm]] <- as.numeric(df[[tm]])
  out <- out[, c(st, tm, keep_c), drop = FALSE]
  rownames(out) <- ids
  ok <- !is.na(out[[st]]) & !is.na(out[[tm]]) & is.finite(out[[tm]]) &
        out[[st]] %in% c(0L, 1L) & out[[tm]] > 0
  out <- out[ok, , drop = FALSE]
  saveRDS(out, file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds")))
  log_msg(acc, ": surv ", nrow(out), " patients | ", token, " events ",
          sum(out[[st]] == 1L), " | time(y) ", sprintf("%.2f-%.2f", min(out[[tm]]), max(out[[tm]])))
  invisible(out)
}

write_pheno <- function(acc, df) saveRDS(df, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))

# time column -> years, driven by the SDRF's own Unit column (verbatim from 104)
unit_for <- function(s, value_col) {
  i <- match(value_col, names(s))
  if (is.na(i) || i >= ncol(s)) return(rep(NA_character_, nrow(s)))
  nm <- names(s)[i + 1L]
  if (!is.na(nm) && grepl("unit", nm, ignore.case = TRUE)) return(s[[i + 1L]])
  rep(NA_character_, nrow(s))
}

to_years <- function(v, unit, default_unit = "day") {
  u <- if (is.null(unit)) rep(default_unit, length(v)) else tolower(unit)
  u[is.na(u) | !nzchar(trimws(u))] <- default_unit
  out <- v
  out[grepl("year", u)]  <- out[grepl("year", u)]
  out[grepl("month", u)] <- out[grepl("month", u)] / 12
  out[grepl("day", u)]   <- out[grepl("day", u)] / 365.25
  out[grepl("week", u)]  <- out[grepl("week", u)] / 52.1786
  out
}

# ============================== cohort builders ==============================

# ---- E-MTAB-1727 : Lung Cancer, OS, GPL570, CEL-only ------------------------
# SDRF: 173 assay rows / 109 Source Names. 64 patients carry a SECOND assay =
# `*_DNA.txt` (genomic DNA, SNP/CGH), whose survival fields are byte-identical to
# the patient's CEL row. -> keep the expression (CEL) row per patient = 109 samples,
# de-duplicated to distinct patients; 90 usable OS pairs / 71 deaths.
build_E_MTAB_1727 <- function() {
  acc <- "E-MTAB-1727"; s <- read_sdrf(acc)
  d <- file.path(SUP, acc)
  cels <- list.files(d, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)
  if (length(cels) != 109) stop(acc, ": expected 109 CEL, found ", length(cels))
  log_msg(acc, ": RMA over ", length(cels), " CEL files (affy + hgu133plus2cdf)")
  suppressMessages(library(affy))
  ab <- affy::ReadAffy(filenames = cels)
  m <- Biobase::exprs(affy::rma(ab))
  sid <- pick_col(s, "^Source Name$")
  adf <- pick_col(s, "Array Data File", required = FALSE)
  is_cel_row <- grepl("\\.CEL$", s[[adf]], ignore.case = TRUE)
  log_msg(acc, ": SDRF rows ", nrow(s), " | CEL rows ", sum(is_cel_row),
          " | distinct Source Name ", length(unique(s[[sid]])),
          " | non-CEL (genomic DNA) rows ", sum(!is_cel_row))
  s_cel <- s[is_cel_row, , drop = FALSE]
  colnames(m) <- map_samples_1727(basename(cels), s_cel, adf, sid)
  # patient identity = Source Name; one CEL row per patient by construction
  if (anyDuplicated(colnames(m))) {
    log_msg(acc, ": collapsing ", sum(duplicated(colnames(m))),
            " repeated patient CEL columns by rowMeans")
    mm <- vapply(split(seq_len(ncol(m)), colnames(m)), function(j)
      rowMeans(as.matrix(m[, j, drop = FALSE]), na.rm = TRUE), numeric(nrow(m)))
    m <- mm; colnames(m) <- colnames(mm)
  }
  surv <- data.frame(row.names = s_cel[[sid]], stringsAsFactors = FALSE)
  surv$sex              <- s_cel[[pick_col(s_cel, "Sex", FALSE)]]
  surv$disease          <- s_cel[[pick_col(s_cel, "Disease", FALSE)]]
  surv$organ_part       <- s_cel[[pick_col(s_cel, "Organism Part", FALSE)]]
  surv$tnm              <- s_cel[[pick_col(s_cel, "TNM", FALSE)]]
  surv$disease_staging  <- s_cel[[pick_col(s_cel, "disease staging", FALSE)]]
  surv$scc_histology    <- s_cel[[pick_col(s_cel, "SCC Histology", FALSE)]]
  surv$tobacco_smoking  <- s_cel[[pick_col(s_cel, "Tobacco Smoking", FALSE)]]
  surv$age_at_diagnosis <- num(s_cel[[pick_col(s_cel, "Age At Diagnosis", FALSE)]])
  rc <- pick_col(s_cel, "Relaspe-Free Survival In Month", FALSE)
  surv$rfs_time  <- to_years(num(s_cel[[rc]]), unit_for(s_cel, rc), "month")
  surv$rfs_event <- suppressWarnings(as.integer(num(s_cel[[pick_col(s_cel, "Relaspe-Free Survival Event", FALSE)]])))
  surv$OS_status <- suppressWarnings(as.integer(num(s_cel[[pick_col(s_cel, "Overall Survival Event")]])))
  oc <- pick_col(s_cel, "\\[Overall Survival\\]")
  surv$OS_time <- to_years(num(s_cel[[oc]]), unit_for(s_cel, oc), "month")
  surv <- surv[!duplicated(rownames(surv)), , drop = FALSE]
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

# map CEL basenames onto SDRF Source Names; 1727's CEL filenames drop the cohort
# prefix that the Source Name keeps ("CBA_114_U133_2.CEL" -> "LUTN_CBR_HS_114"),
# so match on the SDRF `Array Data File` first and fall back to Assay/Extract Name.
map_samples_1727 <- function(keys, s, adf, sid) {
  idx <- match(keys, s[[adf]])
  if (anyNA(idx)) {
    an <- pick_col(s, "Assay Name", required = FALSE)
    if (!is.na(an)) idx[is.na(idx)] <- match(cel_key(keys)[is.na(idx)], cel_key(s[[an]]))
  }
  if (anyNA(idx)) stop("unmatched CEL keys: ", paste(utils::head(keys[is.na(idx)], 5), collapse = ", "))
  s[[sid]][idx]
}

cel_key <- function(x) tolower(gsub("\\.cel$", "", basename(x), ignore.case = TRUE))

# ---- E-MTAB-6389 : Liver Cancer (intrahepatic cholangiocarcinoma), OS, GPL17585
# SDRF: 109 rows (ICC###G tumour / ICC###NT normal), 82 distinct individuals.
# Expression = data_exp_icc.txt (HTA-2.0 PSR probesets). Patient-level.
build_E_MTAB_6389 <- function() {
  acc <- "E-MTAB-6389"; s <- read_sdrf(acc)
  f <- file.path(SUP, acc, "data_exp_icc.txt")
  stopifnot(file.exists(f))
  log_msg(acc, ": reading processed matrix ", f)
  m <- read_unlabelled_tsv(f, acc)
  log_msg(acc, ": matrix ", nrow(m), " probes x ", ncol(m), " samples")
  sid <- pick_col(s, "^Source Name$")
  ind <- pick_col(s, "individual", required = FALSE)
  if (is.na(ind)) ind <- sid
  j <- match(colnames(m), s[[sid]])
  if (all(is.na(j))) stop(acc, ": matrix columns match no Source Name")
  if (anyNA(j)) { log_msg(acc, ": dropping ", sum(is.na(j)), " matrix columns absent from the SDRF")
                  m <- m[, !is.na(j), drop = FALSE]; j <- j[!is.na(j)] }
  pid <- s[[ind]][j]
  colnames(m) <- pid
  if (anyDuplicated(colnames(m))) {
    nd <- sum(duplicated(colnames(m)))
    log_msg(acc, ": ", nd, " duplicate patient columns (tumour+normal of the same individual) -> rowMeans per patient")
    mm <- vapply(split(seq_len(ncol(m)), colnames(m)), function(k)
      rowMeans(as.matrix(m[, k, drop = FALSE]), na.rm = TRUE), numeric(nrow(m)))
    m <- mm; colnames(m) <- colnames(mm)
  }
  # one survival record per patient (the individual's fields are identical across its rows)
  s2 <- s[!duplicated(s[[ind]]), , drop = FALSE]
  surv <- data.frame(row.names = s2[[ind]], stringsAsFactors = FALSE)
  surv$sex               <- s2[[pick_col(s2, "sex", FALSE)]]
  surv$disease           <- s2[[pick_col(s2, "disease", FALSE)]]
  surv$organism_part     <- s2[[pick_col(s2, "organism part", FALSE)]]
  surv$sampling_site     <- s2[[pick_col(s2, "sampling site", FALSE)]]
  surv$cirrhosis         <- s2[[pick_col(s2, "cirrhosis", FALSE)]]
  surv$vascular_invasion <- s2[[pick_col(s2, "vascular invasion", FALSE)]]
  surv$OS_status <- suppressWarnings(as.integer(num(s2[[pick_col(s2, "event death")]])))
  oc <- pick_col(s2, "overall survival")
  surv$OS_time <- to_years(num(s2[[oc]]), unit_for(s2, oc), "month")
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

# processed matrix whose header row has NO placeholder for the leading probe column
# (verbatim from 104)
read_unlabelled_tsv <- function(path, tag = basename(path)) {
  h <- strsplit(readLines(path, n = 1, warn = FALSE), "\t")[[1]]
  h <- gsub('^"|"$', "", h)
  x <- data.table::fread(path, sep = "\t", header = FALSE, skip = 1, data.table = FALSE,
                         check.names = FALSE, quote = "\"", colClasses = "character")
  ids <- as.character(x[[1]])
  m <- as.data.frame(lapply(x[, -1, drop = FALSE], num), check.names = FALSE)
  if (length(h) == ncol(m) + 1L) h <- h[-1]
  if (length(h) != ncol(m)) stop(tag, ": header fields ", length(h),
                                 " != value columns ", ncol(m))
  colnames(m) <- h
  if (anyDuplicated(ids)) {
    log_msg(tag, ": ", sum(duplicated(ids)), " duplicated probe rows -> de-duplicated")
    keep <- !duplicated(ids); m <- m[keep, , drop = FALSE]; ids <- ids[keep]
  }
  rownames(m) <- ids
  m
}

# ---- E-MTAB-6877 : Mesothelioma, OS, GPL16686, CEL-only ---------------------
# SDRF: 67 rows (63 malignant pleural mesothelioma + 4 normal pleura). Survival
# fields: `Characteristics[overall survival]` + `[event death]` (dead/alive/NA).
# No Unit column -> the 1..164 range is MONTHS (median 16). 4 normals + 1 OS=0 row
# have no usable pair -> 59 patients / 50 deaths. RMA with oligo + pd.hugene.2.0.st.
build_E_MTAB_6877 <- function() {
  acc <- "E-MTAB-6877"; s <- read_sdrf(acc)
  d <- file.path(SUP, acc)
  cels <- list.files(d, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)
  if (length(cels) != 67) stop(acc, ": expected 67 CEL, found ", length(cels))
  log_msg(acc, ": RMA over ", length(cels), " CEL files (oligo + pd.hugene.2.0.st)")
  suppressMessages({ library(oligo); library(Biobase) })
  ab <- oligo::read.celfiles(cels, pkgname = "pd.hugene.2.0.st", verbose = FALSE)
  m <- Biobase::exprs(oligo::rma(ab, target = "core"))
  sid <- pick_col(s, "^Source Name$")
  adf <- pick_col(s, "Array Data File", required = FALSE)
  idx <- match(basename(cels), s[[adf]])
  if (anyNA(idx)) idx[is.na(idx)] <- match(cel_key(cels)[is.na(idx)], cel_key(s[[adf]]))
  if (anyNA(idx)) stop(acc, ": unmatched CEL keys: ",
                       paste(utils::head(basename(cels)[is.na(idx)], 5), collapse = ", "))
  colnames(m) <- s[[sid]][idx]
  surv <- data.frame(row.names = s[[sid]], stringsAsFactors = FALSE)
  surv$sex              <- s[[pick_col(s, "sex", FALSE)]]
  surv$age              <- num(s[[pick_col(s, "age", FALSE)]])
  surv$disease          <- s[[pick_col(s, "disease", FALSE)]]
  surv$histology        <- s[[pick_col(s, "histology", FALSE)]]
  surv$tnm_stage        <- s[[pick_col(s, "tnm.stage", FALSE)]]
  surv$asbestos_exposure<- s[[pick_col(s, "asbestos.exposure", FALSE)]]
  ev <- tolower(trimws(s[[pick_col(s, "event death")]]))
  surv$OS_status <- ifelse(ev == "dead", 1L, ifelse(ev == "alive", 0L, NA_integer_))
  oc <- pick_col(s, "overall survival")
  surv$OS_time <- to_years(num(s[[oc]]), unit_for(s, oc), "month")   # unit col absent -> default month
  log_msg(acc, ": SDRF event counts dead=", sum(ev == "dead"), " alive=", sum(ev == "alive"),
          " unparsed/NA=", sum(!ev %in% c("dead", "alive")))
  m <- m[, intersect(colnames(m), rownames(surv)), drop = FALSE]
  surv <- surv[colnames(m), , drop = FALSE]
  write_expr(acc, m); write_surv(acc, surv, "OS"); write_pheno(acc, surv)
}

BUILDERS <- list(
  "E-MTAB-1727" = build_E_MTAB_1727,
  "E-MTAB-6389" = build_E_MTAB_6389,
  "E-MTAB-6877" = build_E_MTAB_6877
)

todo <- if (WANT == "all") names(BUILDERS) else strsplit(WANT, ",")[[1]]
for (a in todo) {
  if (!a %in% names(BUILDERS)) { log_msg("no builder for ", a, " (skipped)"); next }
  log_msg("=== building ", a, " ===")
  res <- tryCatch({ BUILDERS[[a]](); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
  log_msg(a, " -> ", res)
}
log_msg("done")
}
