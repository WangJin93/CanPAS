# ===========================================================================
# CanPAS curation pipeline -- 08_qc_screen_and_verify
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 04_qc_report.R
#   2. 16_verify_catalog_mirror.R
#   3. 36_screen_geo_candidates.R
#   4. 37_geo_expansion_screen.R
#   5. 113_sweep_gpl_mirror_divergence.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "08_qc_screen_and_verify.R", package = "CanPAS"))
#   run_04_qc_report()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================


# ---------------------------------------------------------------------------
# run_04_qc_report()  <-  verbatim pipeline/R/04_qc_report.R
# ---------------------------------------------------------------------------
run_04_qc_report <- function() {
# 04_qc_report.R
# -----------------------------------------------------------------------------
# CanPAS 数据处理 第 4 步（QC）：汇总一个/多个数据集的处理状态，输出
#   pipeline/out/QC_<date>.csv  每数据集：样本数、事件数、时间范围、
#     表达谱探针数、平台注释覆盖探针数、产物文件是否存在
#   并在控制台打印可读摘要。
# 用法: Rscript 04_qc_report.R [<ProjectRoot>] [<ACC1> <ACC2> ...]
# -----------------------------------------------------------------------------
suppressMessages({library(dplyr)})
args <- commandArgs(trailingOnly=TRUE)
ROOT <- if (length(args) >= 1) args[1] else "~/data/Project/CanPAS"
ROOT <- path.expand(ROOT)
accs <- if (length(args) >= 2) args[-1] else c("GSE52903", "GSE44001")

summary_row <- function(acc) {
  expr <- tryCatch(readRDS(file.path(ROOT,"data/expr", paste0(acc,".rds"))),
                   error=function(e) NULL)
  pheno <- tryCatch(readRDS(file.path(ROOT,"data/pheno", paste0(acc,".rds"))),
                    error=function(e) NULL)
  surv <- NULL
  if (file.exists(file.path(ROOT,"data/processed/surv", paste0(acc,"_surv.rds"))))
    surv <- tryCatch(readRDS(file.path(ROOT,"data/processed/surv", paste0(acc,"_surv.rds"))),
                     error=function(e) NULL)
  info <- tryCatch(read.csv(file.path(ROOT,"data/dataset_info.csv")), error=function(e) NULL)
  gpl <- if (!is.null(info)) {
    g <- info$GPL[info$Accession == acc]
    if (!length(g)) NA_character_ else as.character(g[1])
  } else NA_character_
  gmap <- if (!is.na(gpl))
    tryCatch(readRDS(file.path(ROOT,"data/processed/gpl", paste0(gpl,".rds"))),
             error=function(e) NULL) else NULL
  n_probe_mapped <- if (!is.null(expr) && !is.null(gmap))
    sum(expr$ID_REF %in% rownames(gmap)) else NA_integer_
  n_probe_gene <- if (!is.null(expr) && !is.null(gmap))
    sum(expr$ID_REF %in% rownames(gmap)[!is.na(gmap$gene_id) & gmap$gene_id != ""]) else NA_integer_

  surv_endpoints <- if (!is.null(surv)) {
    sc <- grep("_(status|time)$", colnames(surv), value=TRUE)
    paste(unique(sub("_(status|time)$","",sc)), collapse="|")
  } else "none"
  ev <- if (!is.null(surv)) {
    sc <- grep("_status$", colnames(surv), value=TRUE)
    if (length(sc)) sum(surv[[sc[1]]] == 1, na.rm=TRUE) else NA_integer_
  } else NA_integer_
  tm <- if (!is.null(surv)) {
    tc <- grep("_time$", colnames(surv), value=TRUE)
    if (length(tc)) paste(round(range(surv[[tc[1]]], na.rm=TRUE),3), collapse=" - ") else NA_character_
  } else NA_character_

  data.frame(
    Accession=acc, GPL=gpl, N_pheno=if(!is.null(pheno)) nrow(pheno) else NA_integer_,
    N_expr_samples=if(!is.null(expr)) ncol(expr)-1 else NA_integer_,
    N_probes=if(!is.null(expr)) nrow(expr) else NA_integer_,
    probes_in_gpl_map=n_probe_mapped,
    probes_with_gene=n_probe_gene,
    surv_file=file.exists(file.path(ROOT,"data/processed/surv", paste0(acc,"_surv.rds"))),
    surv_endpoint=surv_endpoints, events=ev, time_range_yr=tm,
    raw_file=file.exists(file.path(ROOT,"data/raw", paste0(acc,".gz"))),
    stringsAsFactors=FALSE)
}

res <- do.call(rbind, lapply(accs, summary_row))
outf <- file.path(ROOT, "pipeline/out", paste0("QC_", format(Sys.Date(), "%Y%m%d"), ".csv"))
write.csv(res, outf, row.names=FALSE)
print(res, row.names=FALSE)
message("QC table written to ", outf)
}

# ---------------------------------------------------------------------------
# run_16_verify_catalog_mirror()  <-  verbatim pipeline/R/16_verify_catalog_mirror.R
# ---------------------------------------------------------------------------
run_16_verify_catalog_mirror <- function() {
# 16_verify_catalog_mirror.R -------------------------------------------------
# 自动校验 catalog 与镜像/本地数据的一致性,可在 CI 中运行(不一致时退出码非 0)。
#
# 检查项:
#   1. Accession 命名规范:不得出现 `-GPL` 破折号写法(API 按 catalog 字面查库)
#   2. Accession 重复
#   3. 有终点标注但镜像内无表达/无生存表(不可分析的行)
#   4. 无终点标注的行(仅登记,不算错误,但计数并列出)
#   5. 家族覆盖:OS/DSS/DFS/PFS/MFS 各自行数
#   6. 镜像中不属于 catalog 的表(孤儿表)
#   7. 镜像表名里出现破折号(API 把表名直接拼进 SQL,破折号表一律 500,
#      实测:`X-Y` 表存在也返回 500,`X_Y` 同内容表返回 200;catalog 里带破折号的
#      accession 由 get_data() 归一化为下划线后查询,镜像必须存下划线名)
#   8. 生存表的状态编码是否为 0/1(抽样每表,全表检查事件列取值集合)
#   9. 终点列整组缺失:SurvivalTypes 已登记但 EndpointFamilies / EndpointPrimary /
#      EP_OS..EP_MFS 全为 NA —— 即"漏跑 11_endpoint_families.R"。报 ERROR。
#  10. 表达表在镜像里、其平台(GPL)注释表却不在镜像里 —— 06_upload_db.R 旧版把 GPL
#      上传挂在 n_expr>79 上留下的缺陷类(GPL1223/GPL16686/GPL27956)。报 ERROR。
#      (只查镜像;TCGA 行不进镜像,已排除)
#
# 输出: pipeline/out/catalog_mirror_check.csv / REPORT_catalog_mirror_check.md
# 用法: CPAS_DB_PASSWORD=... Rscript pipeline/R/16_verify_catalog_mirror.R [--no-db]
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(RMySQL)
})
root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
args <- commandArgs(trailingOnly = TRUE)
use_db <- !"--no-db" %in% args

di <- get(load(file.path(root, "CanPAS/data/dataset_info.rda")))
acc <- as.character(di$Accession)
acc_db <- gsub("-", "_", acc, fixed = TRUE)

loc_expr <- sub("\\.rds$", "", list.files(file.path(root, "data/expr"), pattern = "\\.rds$"))
loc_surv <- sub("_surv\\.rds$", "", list.files(file.path(root, "data/processed/surv"),
                                              pattern = "_surv\\.rds$"))

expr_db <- surv_db <- character(0); gpl_db <- character(0); db_ok <- FALSE
if (use_db) {
  pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
  if (!nzchar(pw)) {
    message("CPAS_DB_PASSWORD is not set: running the local-file checks only.")
  } else {
    con <- tryCatch(dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                              dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                              password = pw), error = function(e) e)
    if (inherits(con, "error")) {
      message("Could not connect to the mirror: ", conditionMessage(con))
    } else {
      tb <- dbGetQuery(con, "SHOW TABLES")[[1]]
      surv_db <- tb[grepl("_surv$", tb)]
      gpl_db  <- tb[grepl("^GPL[0-9]", tb)]
      expr_db <- setdiff(setdiff(tb, surv_db), gpl_db)
      db_ok <- TRUE
    }
  }
}

chk <- data.frame(
  Accession = acc, Type = di$Type, N = di$N,
  endpoint_annotated = !is.na(di$SurvivalTypes),
  expr_local = acc %in% loc_expr | acc_db %in% loc_expr |
    gsub("_", "-", acc, fixed = TRUE) %in% loc_expr,
  surv_local = acc %in% loc_surv | acc_db %in% loc_surv |
    gsub("_", "-", acc, fixed = TRUE) %in% loc_surv,
  expr_mirror = if (db_ok) acc_db %in% expr_db else NA,
  surv_mirror = if (db_ok) paste0(acc_db, "_surv") %in% surv_db else NA,
  stringsAsFactors = FALSE)
chk$dash_spelling <- grepl("-GPL", acc)
chk$usable <- if (db_ok) (chk$expr_mirror & chk$surv_mirror) | startsWith(acc, "TCGA-") else NA

# 9. 终点列整组缺失：SurvivalTypes 已登记，但 EndpointFamilies / EndpointPrimary /
#    EP_OS..EP_MFS 全为 NA。这是"登记了生存类型、却漏跑 11_endpoint_families.R"的
#    真实缺陷（EMBL-EBI 批次曾留下 11 行），它让该队列在家族配对/整合 KM 里不可见：
#    SurvivalTypes 非空会让第 4 项把该行算作 endpoint_annotated，于是既不在
#    "without endpoint annotation" 计数里，又不贡献任何家族成员，静默漏掉。
#    判定为 ERROR；本地即可判定（db 不可用时只取 (a) 支）。
ep_cols <- grep("^EP_(OS|DSS|DFS|PFS|MFS)$", colnames(di), value = TRUE)
anno_elsewhere <- !is.na(di$SurvivalTypes) & nzchar(as.character(di$SurvivalTypes))
all_ep_na <- if (length(ep_cols))
  rowSums(!is.na(di[, ep_cols, drop = FALSE])) == 0L else rep(TRUE, nrow(di))
if ("EndpointFamilies" %in% colnames(di))
  all_ep_na <- all_ep_na & is.na(di$EndpointFamilies)
if ("EndpointPrimary" %in% colnames(di))
  all_ep_na <- all_ep_na & is.na(di$EndpointPrimary)
# (b) 支：镜像里确实有该队列的生存表（本地有 surv 文件或用 --no-db 时的本地代理）
has_surv_any <- chk$surv_local | (!is.na(chk$surv_mirror) & chk$surv_mirror)
ep_missing <- all_ep_na & (anno_elsewhere | has_surv_any)

issues <- list()
add_issue <- function(level, item, n, detail) {
  issues[[length(issues) + 1L]] <<- data.frame(level = level, check = item,
                                               n = n, detail = detail,
                                               stringsAsFactors = FALSE)
}
add_issue("error", "duplicated accessions", sum(duplicated(acc)),
          paste(unique(acc[duplicated(acc)]), collapse = ", "))
add_issue("error", "'-GPL' dash spelling", sum(chk$dash_spelling),
          paste(acc[chk$dash_spelling], collapse = ", "))
# 第 3 项（endpoint_annotated 的定义只看 SurvivalTypes）看不见「终点列整组 NA」，
# 所以这条必须独立报错。detail 里带上 SurvivalTypes，便于直接看出漏跑的是哪一步。
add_issue("error",
          "endpoint columns all NA while the row is endpoint-annotated or has a survival table",
          sum(ep_missing),
          if (any(ep_missing))
            paste0(acc[ep_missing], " (SurvivalTypes=", as.character(di$SurvivalTypes[ep_missing]), ")",
                   collapse = ", ")
          else "")
if (db_ok) {
  miss <- chk[chk$endpoint_annotated & !chk$expr_mirror & !startsWith(chk$Accession, "TCGA-"), ]
  add_issue("error", "annotated but no expression in the mirror", nrow(miss),
            paste(miss$Accession, collapse = ", "))
  miss2 <- chk[chk$endpoint_annotated & !chk$surv_mirror & !startsWith(chk$Accession, "TCGA-"), ]
  add_issue("error", "annotated but no survival table in the mirror", nrow(miss2),
            paste(miss2$Accession, collapse = ", "))
  # 10. 表达表在镜像里、平台(GPL)注释表却不在镜像里。这是 06_upload_db.R 旧版把 GPL
  #     上传挂在 n_expr > 79 上的真实遗留缺陷：n_expr<=79 的小队列表达/生存表进了镜像、
  #     注释表被静默跳过，于是 ID_map -> GPL -> expr 的基因查询路径走不通，只能靠手工
  #     脚本补救（GPL1223/GSE1379、GPL16686/E-MTAB-6877、GPL27956/GSE205209 各一次）。
  #     GPL 修好之后 06_upload_db.R 的规则是"镜像缺表就上传、已有表永不覆盖"，所以本项
  #     应恒为 0；非 0 就意味着又出现了注释表没进镜像的队列。判定为 ERROR。
  #     只用镜像表清单，故与其它镜像检查一样在 --no-db 下不参与判定（--no-db 仍退出 0）。
  #     平台表名的归一化与 06_upload_db.R 一致（GPL 里的 '-' 写成 '_'），并且是在**全部**
  #     镜像表名里查，而不是只看 '^GPL[0-9]' —— 否则 CGGA_*_PLAT / E-MTAB-4321_PLAT
  #     这类真实存在的平台表会被误判为缺失。TCGA 行不进镜像（fetch-on-demand，共用
  #     TCGA_HiSeqV2），显式排除。
  gpl_tab <- gsub("-", "_", as.character(di$GPL), fixed = TRUE)
  gpl_known <- !is.na(gpl_tab) & nzchar(gpl_tab)
  gpl_missing <- !startsWith(acc, "TCGA-") & chk$expr_mirror &
    !(gpl_known & gpl_tab %in% tb)
  add_issue("error",
            "expression table in the mirror but platform/GPL table is missing",
            sum(gpl_missing),
            if (any(gpl_missing))
              paste0(acc[gpl_missing], " (GPL=", as.character(di$GPL[gpl_missing]), ")",
                     collapse = ", ")
            else "")
  orph <- setdiff(expr_db, acc_db)
  orph <- orph[!grepl("_PLAT$", orph)]
  add_issue("warn", "mirror tables not in the catalog", length(orph),
            paste(orph, collapse = ", "))
  # 破折号表名：API 把表名直接拼进 SQL，破折号一律 500（实测），因此镜像里
  # 出现破折号表名就是"上传路径没有归一化"的信号，必须报错而不是留到运行时。
  dash_tb <- grep("-", tb, value = TRUE, fixed = TRUE)
  add_issue("error", "mirror table name contains a dash (unreachable through the API)",
            length(dash_tb), paste(dash_tb, collapse = ", "))
  if (length(surv_db)) {
    bad <- character(0)
    for (t in surv_db) {
      cols <- tryCatch(dbGetQuery(con, sprintf("SHOW COLUMNS FROM `%s`", t))$Field,
                       error = function(e) character(0))
      # only real endpoint pairs (<TOKEN>_time + <TOKEN>_status) are survival
      # status columns; biomarker columns such as Kras_status (Mut/Wt) are not
      times <- sub("_time$", "", grep("_time$", cols, value = TRUE))
      sc <- intersect(paste0(times, "_status"), grep("_status$", cols, value = TRUE))
      for (cc in sc) {
        v <- tryCatch(dbGetQuery(con, sprintf("SELECT DISTINCT `%s` v FROM `%s` WHERE `%s` IS NOT NULL",
                                             cc, t, cc))$v,
                      error = function(e) numeric(0))
        if (length(v) && !all(v %in% c(0, 1))) bad <- c(bad, paste0(t, ".", cc))
      }
    }
    add_issue("error", "survival status columns not coded 0/1", length(bad),
              paste(utils::head(bad, 10), collapse = ", "))
  }
  # Sample-size columns vs the mirror. N used to hold a legacy planned size
  # (e.g. GSE31210: 133 instead of 226), so it is re-checked here against the
  # actual table dimensions to stop it drifting again.
  if (all(c("N", "n_expr", "n_surv", "n_events") %in% colnames(di))) {
    bad_n <- ev0 <- na_n <- character(0)
    for (i in seq_len(nrow(di))) {
      a <- as.character(di$Accession[i])
      if (startsWith(a, "TCGA-")) next
      adb <- gsub("-", "_", a, fixed = TRUE)
      ne <- if (adb %in% expr_db) length(dbListFields(con, adb)) - 1L else NA_integer_
      ns <- if (paste0(adb, "_surv") %in% surv_db)
        dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s_surv`", adb))$n else NA_integer_
      if (!is.na(ne) && !is.na(di$n_expr[i])   && ne != di$n_expr[i]) bad_n <- c(bad_n, sprintf("%s expr %s!=%s", a, di$n_expr[i], ne))
      if (!is.na(ns) && !is.na(di$n_surv[i])   && ns != di$n_surv[i]) bad_n <- c(bad_n, sprintf("%s surv %s!=%s", a, di$n_surv[i], ns))
      if (!is.na(di$N[i]) && !is.na(ne) && di$N[i] > ne) bad_n <- c(bad_n, sprintf("%s N>n_expr", a))
      if (!is.na(di$N[i]) && !is.na(ns) && di$N[i] > ns) bad_n <- c(bad_n, sprintf("%s N>n_surv", a))
      if (!is.na(di$N[i]) && !is.na(di$n_events[i]) && di$n_events[i] == 0L) ev0 <- c(ev0, a)
      if (is.na(di$N[i]) && !is.na(ne)) na_n <- c(na_n, a)
    }
    add_issue("error", "catalog N/n_expr/n_surv disagree with the mirror",
              length(bad_n), paste(utils::head(bad_n, 10), collapse = ", "))
    # 新增两类曾真实发生过的缺陷：
    #   (a) 生存表样本 ID 与表达表无法连接（如 GSE4573/GSE3494 的 row_names 是行号）
    #   (b) 终点状态列整列为空（如 GSE48075 的 OS_status，源里是 "os censor"）
    nojoin <- empty_st <- character(0)
    for (i in seq_len(nrow(di))) {
      a <- as.character(di$Accession[i]); if (startsWith(a, "TCGA-")) next
      adb <- gsub("-", "_", a, fixed = TRUE); st_t <- paste0(adb, "_surv")
      if (!st_t %in% surv_db) next
      cols <- dbListFields(con, st_t)
      sv <- tryCatch(dbGetQuery(con, sprintf("SELECT * FROM `%s`", st_t)), error = function(e) NULL)
      if (is.null(sv) || !nrow(sv)) next
      if (adb %in% expr_db) {
        ids <- as.character(sv[[1]]); ex <- setdiff(dbListFields(con, adb), "row_names")
        if (!length(intersect(ids, ex))) nojoin <- c(nojoin, a)
      }
      for (cc in grep("_status$", cols, value = TRUE)) {
        v <- sv[[cc]]
        if (!any(!is.na(v)) && !any(is.na(v))) next          # 全空表跳过
        if (all(is.na(v))) empty_st <- c(empty_st, paste0(a, ".", cc))
      }
    }
    add_issue("error", "survival table shares no sample ID with the expression table",
              length(nojoin), paste(nojoin, collapse = ", "))
    add_issue("warn", "endpoint status column is empty (source field not parsed?)",
              length(empty_st), paste(empty_st, collapse = ", "))
    add_issue("warn", "primary endpoint has 0 events", length(ev0),
              paste(ev0, collapse = ", "))
    add_issue("info", "has expression but N is NA (no primary endpoint data)",
              length(na_n), paste(na_n, collapse = ", "))
  }
}
add_issue("info", "rows without endpoint annotation", sum(!chk$endpoint_annotated),
          paste(acc[!chk$endpoint_annotated], collapse = ", "))

fams <- c("OS", "DSS", "DFS", "PFS", "MFS")
fam_counts <- vapply(fams, function(f) {
  col <- paste0("EP_", f)
  if (!col %in% colnames(di)) return(0L)
  v <- as.character(di[[col]]); sum(!is.na(v) & nzchar(v))
}, integer(1))

res <- do.call(rbind, issues)
out_csv <- file.path(root, "pipeline/out/catalog_mirror_check.csv")
utils::write.csv(res, out_csv, row.names = FALSE)
utils::write.csv(chk, file.path(root, "pipeline/out/catalog_mirror_check_detail.csv"), row.names = FALSE)

rep_lines <- c(
  "# catalog <-> mirror consistency check", "",
  sprintf("Run: %s | database checks: %s", format(Sys.time()), if (db_ok) "yes" else "no (no DB access)"),
  sprintf("catalog rows: %d | with endpoints: %d | fully usable: %s", nrow(chk),
          sum(chk$endpoint_annotated),
          if (db_ok) as.character(sum(chk$usable, na.rm = TRUE)) else "n/a"),
  "", "## Family coverage", "",
  paste0("| family | cohorts |", collapse = ""),
  "|---|---|",
  sprintf("| %s | %d |", fams, fam_counts),
  "", "## Checks", "",
  "| level | check | n | detail |", "|---|---|---|---|",
  sprintf("| %s | %s | %d | %s |", res$level, res$check, res$n,
          substr(gsub("\\|", "/", res$detail), 1, 400)))
out_md <- file.path(root, "pipeline/out/REPORT_catalog_mirror_check.md")
writeLines(rep_lines, out_md)
cat(paste(rep_lines, collapse = "\n"), "\n")

if (db_ok) dbDisconnect(con)
errs <- sum(res$n[res$level == "error"])
cat(sprintf("\nERRORS: %d | warnings: %d | info: %d\n", errs,
            sum(res$n[res$level == "warn"]), sum(res$n[res$level == "info"])))
if (errs > 0) quit(status = 1) else quit(status = 0)
}

# ---------------------------------------------------------------------------
# run_36_screen_geo_candidates()  <-  verbatim pipeline/R/36_screen_geo_candidates.R
# ---------------------------------------------------------------------------
run_36_screen_geo_candidates <- function() {
# 36_screen_geo_candidates.R -------------------------------------------------
# 目的（离线，不下载任何东西）
#   扫描 data/raw/*.gz 里的**每一份 GEO series matrix**，把生存信息尽量完整地
#   挖出来，产出"有效样本数 > 50 且尚未编目"的候选清单，供后续按 6 步流程入库。
#
# 为什么
#   前面的候选脚本只看 `!Sample_characteristics_ch1` 的少数行，漏掉两类关键来源：
#     (a) `!Sample_description` 里的自由文本（例如
#         "Age: 72, ... RFS event: no relapse, RFS months: 109,
#          Overall survival event: alive, Overall survival months: 109"），
#     (b) 特征行里的 key:value / key=value 混写（例如
#         "Clinical information: SampleID=25;...;DFS=75;Status=recur"、
#         "[SURTIM=69.24 (Follow-up time in months ...)]"、
#         "Did_Patient_Die(0=No,1=Yes) = 1"、
#         "status: Deceased after 27 months after diagnosis"）。
#   本脚本同时登记 `!Series_supplementary_file` 的文件名（**不下载**）。
#
# 实现要点（踩过的坑）
#   * GEO series matrix 是 CRLF 结尾：一律用 trimws() 去掉 \r（不要自造 whitespace=
#     正则），再去掉单元格两端成对的引号，否则 `$` 锚点会静默失效。
#   * 只读 header（到 `!series_matrix_table_begin` 为止），不解析表达谱矩阵，
#     因此 160+ 份文件可在一两分钟内扫完。
#   * 解析顺序：先按**顶层** `;`/`,` 切段（括号内的分隔符不切，避免
#     "Did_Patient_Die(0=No,1=Yes)"、"(Follow-up ...; integer)" 被切坏），
#     段内取**第一个顶层** `=`（无 `=` 再取 `:`）作为 key/value 分隔符。
#   * 端点 token 由 key 判定（小写）：DSS/CSS->DSS，DFS/DFI/EFS->DFS，
#     RFS/复发->RFS，PFS/PFI/进展->PFS，MFS/转移->MFS，其余 survival/随访/
#     生死/vital->OS（与 11_endpoint_families.R 的家族口径一致）。
#   * 角色（时间/事件）**由值决定**：值能解析出数字就进"时间"候选，值能解码成
#     0/1 就进"事件"候选；同一个 key 可以两者都是（如 "status: 0:alive" 会被
#     时间侧的全 0/1 惩罚规则剔除）。
#   * 没有明确端点的随访时间（如 SURTIM）视为**通用时间**，可配任意 token 的
#     事件列；端点专属时间只配同 token 的事件列。
#   * 事件值解码：NA 词 -> NA；纯 0/1 -> 0/1；否定词（no/alive/...）先于肯定词
#     （dead/recur/...）；其它 -> NA。
#   * 有效样本数 effective_n = 同一 token 下"时间可解析 + 事件可解码"的样本数；
#     events = 其中事件=1 的样本数。
#
# 输出
#   pipeline/out/geo_candidates_gt50.csv   主清单（9 列，见下）
#   pipeline/out/geo_candidates_gt50.md    按癌种分组的中文报告
#   pipeline/out/geo_scan_all.csv          全部文件的扫描明细（含未达阈值者，便于复核）
# CSV 列: local_file, series, platform, n_samples, effective_n, events, token,
#         source_of_survival_info, in_catalog
#
# 用法: Rscript pipeline/R/36_screen_geo_candidates.R [ROOT]
# ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

ROOT <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(ROOT)) ROOT[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
RAW <- file.path(ROOT, "data/raw")
OUT <- file.path(ROOT, "pipeline/out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
MIN_EFF <- 50L

di <- read.csv(file.path(ROOT, "data/dataset_info.csv"), check.names = FALSE, stringsAsFactors = FALSE)

# ---------------------------------------------------------------- 基础工具
strip_quotes <- function(x) {
  x <- trimws(x)                                  # 去掉 CRLF 的 \r（关键）
  x <- sub('^"', '', x); x <- sub('"$', '', x)
  x <- gsub('""', '"', x, fixed = TRUE)
  trimws(x)
}

split_top <- function(s, seps = c(";", ",")) {     # 顶层（括号外）切分
  ch <- strsplit(s, "", fixed = TRUE)[[1]]
  if (!length(ch)) return(character(0))
  depth <- 0L; cur <- character(0); out <- character(0)
  for (cc in ch) {
    if (cc == "(") depth <- depth + 1L
    else if (cc == ")") depth <- max(0L, depth - 1L)
    if (depth == 0L && cc %in% seps) { out <- c(out, paste(cur, collapse = "")); cur <- character(0) }
    else cur <- c(cur, cc)
  }
  out <- trimws(c(out, paste(cur, collapse = "")))
  out[nzchar(out)]
}

# 审计结论（pipeline/out/geo_expansion_AUDIT.md §4）：
#   * SURV_HINT 曾当**硬门槛**（!has_hint(cell) -> 整个单元格丢弃），会静默丢掉真实
#     队列：GSE304059 的 `futime`/`fustat`、GSE22138 的 `metastasis` —— key 名里不含
#     任何 hint 词，于是这些行从未进入解析器。现在 SURV_HINT 只作**同分排序提示**
#     （见 infer_pairs()），不再当门槛。
#   * CRLF：审计实测本批 743 份 header 已由抓取步骤归一为 LF（只有 2 份在
#     !Sample_data_processing 的**值**里含字面 \r，不是行结尾），故此处继续用普通
#     trimws() 去两端空白，不另造 whitespace 正则 —— 该缺陷在本批未被数据触发。
SURV_HINT <- "(surv|follow|vital|death|dead|die|died|recur|relaps|progress|event|status|indicator|month|year|day|dfs|rfs|pfs|dss|mfs|\\bos\\b|surtim|surind|outcome|metasta|futime|fu[ _-]?stat)"
has_hint <- function(x) !is.na(x) & grepl(SURV_HINT, x, ignore.case = TRUE, perl = TRUE)

parse_pairs <- function(cell) {                    # 单元格 -> data.frame(key, value)
  cell <- strip_quotes(cell)
  if (!nzchar(cell)) return(NULL)
  k <- character(0); v <- character(0)
  for (s in split_top(cell)) {
    ch <- strsplit(s, "", fixed = TRUE)[[1]]
    depth <- 0L; eqs <- integer(0); cols <- integer(0)
    for (i in seq_along(ch)) {
      cc <- ch[i]
      if (cc == "(") depth <- depth + 1L
      else if (cc == ")") depth <- max(0L, depth - 1L)
      else if (depth == 0L && cc == "=") eqs <- c(eqs, i)
      else if (depth == 0L && cc == ":") cols <- c(cols, i)
    }
    pos <- if (length(eqs)) eqs[1] else if (length(cols)) cols[1] else NA_integer_
    if (is.na(pos) || pos < 2L || pos >= length(ch)) next
    key <- trimws(paste(ch[seq_len(pos - 1)], collapse = ""))
    val <- trimws(paste(ch[(pos + 1):length(ch)], collapse = ""))
    key <- gsub("^[\\[{(]+", "", key); key <- gsub("[\\]}]+$", "", key)
    val <- gsub("^[\\[{(]+", "", val); val <- gsub("[\\]}]+$", "", val)
    # "[SURTIM=69.24 (Follow-up time in months ...)]" 这类：保留前导数字，丢掉括号里的说明文档
    val <- sub("^([-+]?[0-9]+(?:\\.[0-9]+)?)\\s*\\(.*$", "\\1", val)
    if (nzchar(key)) { k <- c(k, key); v <- c(v, trimws(val)) }
  }
  if (!length(k)) NULL else data.frame(key = k, value = v, stringsAsFactors = FALSE)
}

TOKEN_RE <- list(                                  # 顺序即优先级：先具体后泛化
  DSS = "(^|[^a-z])(dss|css|bcss)([^a-z]|$)|disease[ _-]?specific|disease[ _-]?related|cause[ _-]?specific",
  DFS = "(^|[^a-z])(dfs|dfi|efs)([^a-z]|$)|disease[ _-]?free",
  RFS = "(^|[^a-z])(rfs)([^a-z]|$)|recurrence[ _-]?free|relapse[ _-]?free|recurrence|relaps|recur",
  PFS = "(^|[^a-z])(pfs|pfi)([^a-z]|$)|progression[ _-]?free|progress",
  MFS = "(^|[^a-z])(mfs|drfs)([^a-z]|$)|metasta|distant[ _-]?recurrence",
  OS  = "(^|[^a-z])(os)([^a-z]|$)|overall[ _-]?survival|survival|surv|vital|death|dead|died|die|surtim|surind")
token_of_key <- function(key) {
  k <- tolower(key)
  for (tk in names(TOKEN_RE)) if (grepl(TOKEN_RE[[tk]], k, perl = TRUE)) return(tk)
  NA_character_
}

EXCLUDE_RE <- paste0("age|braf|margin|residual|menopause|race|ethnic|smok|hpv|histolog|",
                     "diagnos|\\bdx\\b|pathologic|birth|surgery|operation|treat|therap|drug|dose|",
                     "batch|plate|scan|center|institut|country|gender|^sex|height|weight|stage|grade|",
                     "site|tumor[ _-]?status|tissue|cell|culture|passage|protocol|keyword|date|",
                     "group|subgrp|cluster|class|subtype|response|responder|censor[ _-]?reason|",
                     "sample[ _-]?id|record[ _-]?id|patient[ _-]?id|^id$|number|prior[ _-]?lines")
excluded_key <- function(key) grepl(EXCLUDE_RE, tolower(key), perl = TRUE)

TIME_RE  <- "time|month|year|day|follow"
EVENT_RE <- "event|status|vital|death|dead|died|die|recur|relaps|progress|indicator|surv|outcome|censor|metasta"
BARE_RE  <- "^(dfs|rfs|os|pfs|pfi|dss|mfs|dfi|efs|drfs|surv|survival|surtim|surind|status|event)$"
# 通用（非 token 专属）时间/事件列的准入：必须带生存语义，避免把
# "doubling time (days)"、"ER status"、"braf status" 之类当成随访时间/事件
# 通用（非 token 专属）时间/事件列的准入。审计修复（§4）：把「随访时间」缩写补进词表
# （futime / ftime / months / years / days / last fu / follow-up），否则 GSE304059 的
# `futime`+`fustat`、GSE22138 的 `months to endpoint`+`metastasis` 这类真实随访对在
# 通用袋里拿不到角色。判定前先剥掉括号内的说明文档（如 "doubling time (days)"），
# 以免把非随访的 "... (days)" 误当随访时间。
GEN_TIME_OK  <- "surv|follow|vital|surtim|overall|disease|recur|relaps|progress|metasta|dfs|rfs|pfs|dss|(^|[^a-z])os([^a-z]|$)|^(time|months|years|days|follow[ _-]?up)([[:space:]]*\\(.*\\))?$|(^|[^a-z])(futime|ftime|months?|years?|days?|fu|last[ _-]?fu|follow[ _-]?up)([^a-z]|$)"
GEN_EVENT_OK <- "surv|vital|death|dead|die|recur|relaps|progress|event|censor|follow|^(status|event|indicator)$|(^|[^a-z])(fustat|stat|fu[ _-]?stat|dead|died)([^a-z]|$)"
EVENT_VAL_MAX <- 40L
# 时间/事件袋打分：全 0/1 的列不算时间（否则 "status: 0:alive" 会被当时间）
score_time  <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return(0L); if (all(v %in% c(0, 1))) return(0L); length(v) }
score_event <- function(v) sum(!is.na(v))     # 事件值超过这个长度视为说明性文字，不当 0/1 解析

NA_VAL_RE <- "^(na|n/?a|unknown|not available|not applicable|no more data|none|null|nd|--?|-|\\?|not determined|undetermined|not done|missing|no data)$"
NEG_RE <- "alive|^no\\b|^n$|^false|negative|free of|no evidence|not recur|not progress|did not|without|^0$|no relapse|no recurrence|no progression|no death|ned|tumor free|^non[ _-]?(recur|relaps|progress|event)"
POS_RE <- "dead|death|deceased|died|^yes|^y$|^true|positive|recur|relaps|progress|event|^1$|recurrence|metasta|decease"
parse_event <- function(v) {
  v <- tolower(strip_quotes(v))
  if (!nzchar(v) || grepl(NA_VAL_RE, v)) return(NA_real_)
  if (grepl("^-?[0-9.]+$", v)) { n <- suppressWarnings(as.numeric(v)); return(if (!is.na(n) && n %in% c(0, 1)) n else NA_real_) }
  if (grepl(NEG_RE, v)) return(0)
  if (grepl(POS_RE, v)) return(1)
  NA_real_
}
parse_time <- function(v) {
  v <- strip_quotes(v)
  m <- regmatches(v, regexpr("[-+]?[0-9]*\\.?[0-9]+", v))
  if (!length(m) || !nzchar(m)) return(NA_real_)
  suppressWarnings(as.numeric(m))
}
COMBINED_RE <- "(alive|dead|deceased|died|no evidence of disease|ned|recurrence|relapse|recurred|relapsed|progression|progressed|metastasis|metastatic)[^0-9]{0,30}?([0-9]+(\\.[0-9]+)?)\\s*(months?|years?|days?|weeks?)"
parse_combined <- function(v) {
  s <- tolower(strip_quotes(v))
  r <- regmatches(s, regexec(COMBINED_RE, s, perl = TRUE))[[1]]
  if (!length(r)) return(NULL)
  word <- r[2]; num <- suppressWarnings(as.numeric(r[3])); unit <- r[6]
  ev <- if (grepl("alive|no evidence|ned", word)) 0 else 1
  tk <- if (grepl("recur|relaps", word)) "RFS" else if (grepl("progress", word)) "PFS" else
        if (grepl("metasta", word)) "MFS" else "OS"
  list(token = tk, time = num, unit = unit, event = ev)
}

# ---------------------------------------------------------------- header 读取
read_header <- function(path, max_lines = 300000L) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  out <- character(0)
  repeat {
    chunk <- readLines(con, n = 5000L, warn = FALSE)
    if (!length(chunk)) break
    hit <- grep("^!series_matrix_table_begin", chunk, ignore.case = TRUE)
    if (length(hit)) { out <- c(out, chunk[seq_len(hit[1] - 1L)]); break }
    out <- c(out, chunk)
    if (length(out) > max_lines) break
  }
  out
}

scan_file <- function(path) {
  h <- read_header(path)
  tag <- sub("\t.*$", "", h)
  body <- sub("^[^\t]*\t?", "", h)
  split_vals <- function(i) {
    v <- strsplit(body[i], "\t", fixed = TRUE)[[1]]
    if (!length(v)) character(0) else vapply(v, strip_quotes, character(1), USE.NAMES = FALSE)
  }
  idx <- function(pattern) grep(pattern, tag, ignore.case = TRUE)
  first_val <- function(pattern) { i <- idx(pattern); if (!length(i)) NA_character_ else (split_vals(i[1])[1] %||% NA_character_) }

  series <- first_val("^!Series_geo_accession$")
  if (is.na(series) || !nzchar(series)) series <- sub("\\.gz$", "", basename(path))
  platform <- sub("^.*-(GPL[0-9]+)$", "\\1", sub("\\.gz$", "", basename(path)))
  if (!grepl("^GPL[0-9]+$", platform)) {
    p <- unique(split_vals(idx("^!Sample_platform_id$")[1]))
    platform <- if (length(p)) p[1] else NA_character_
  }
  title <- first_val("^!Series_title$")
  si <- idx("^!Series_summary$")
  summ <- if (!length(si)) "" else paste(vapply(si, function(j) (split_vals(j)[1] %||% ""), character(1)), collapse = " ")
  fi <- idx("^!Series_supplementary_file")
  supp <- if (!length(fi)) character(0) else unique(unlist(lapply(fi, split_vals)))
  gi <- idx("^!Sample_geo_accession$")
  n_samples <- if (length(gi)) length(split_vals(gi[1])) else NA_integer_
  ti <- idx("^!Sample_title$")
  titles <- if (length(ti)) split_vals(ti[1]) else character(0)
  sni <- idx("^!Sample_source_name_ch1$")
  src_names <- if (length(sni)) split_vals(sni[1]) else character(0)
  if (is.na(n_samples) || !n_samples) n_samples <- length(titles)
  if (!n_samples) return(NULL)

  field_idx <- sort(c(idx("^!Sample_characteristics_ch[0-9]*$"), idx("^!Sample_description$")))
  token_time <- list(); token_event <- list(); gen_time <- list(); gen_event <- list()
  src <- character(0); ref_like <- rep(FALSE, n_samples)
  if (length(src_names) == n_samples)
    ref_like <- ref_like | grepl("reference|universal|stratagene|pool", src_names, ignore.case = TRUE)
  if (length(titles) == n_samples)
    ref_like <- ref_like | grepl("reference|universal|stratagene", titles, ignore.case = TRUE)

  put <- function(bag, key, s, val) {
    stopifnot(is.character(key), length(key) == 1L)
    val <- if (is.numeric(val) && length(val) == 1L) val else NA_real_
    if (is.null(bag[[key]]) || !is.numeric(bag[[key]]) || length(bag[[key]]) != n_samples)
      bag[[key]] <- rep(NA_real_, n_samples)
    bag[[key]][s] <- val
    bag
  }
  for (i in field_idx) {
    vals <- split_vals(i)
    if (length(vals) != n_samples) next
    # 审计修复：不再用 has_hint(vals)/has_hint(cell) 做硬门槛（会丢掉 futime/fustat 行）
    is_desc <- grepl("^!Sample_description$", tag[i], ignore.case = TRUE)
    where <- if (is_desc) "description" else "characteristics"
    for (s in seq_len(n_samples)) {
      cell <- vals[s]
      if (!nzchar(cell)) next
      if (grepl("reference|universal|stratagene|pool", cell, ignore.case = TRUE)) ref_like[s] <- TRUE
      pr <- parse_pairs(cell)
      if (!is.null(pr)) for (r in seq_len(nrow(pr))) {
        key <- pr$key[r]; val <- pr$value[r]
        if (excluded_key(key)) next
        k <- tolower(key)
        tk <- token_of_key(key)
        # ---- 角色由「解出的 token/家族 + 值」决定，不再只看 key 文本（审计 §4 第 2 条）
        # 旧实现 is_time_k = TIME_RE(key) || BARE_RE(key) 之类，使 code_os / drfs /
        # tt.drfs / pfs_index / fustat 这类 key 虽已被 token_of_key() 识别出端点，
        # 却一个角色都拿不到。新规则：
        #   纯时间词 key（time/month/year/day/...） -> 只做时间；
        #   纯事件词 key（status/vital/death/...）   -> 只做事件；
        #   其余                                      -> 交给值决定（能解出数字进时间候选，
        #                                               能解出 0/1 进事件候选；全 0/1 的列
        #                                               由 score_time() 的惩罚挡在时间袋外）。
        k_np <- gsub("\\([^)]*\\)", " ", k)     # 括号内是说明文档，不参与通用词表判定
        time_hint  <- grepl(TIME_RE, k)
        event_hint <- grepl(EVENT_RE, k)
        tv <- if (time_hint  || !event_hint) parse_time(val) else NA_real_
        ev <- if ((event_hint || !time_hint) && nchar(val) <= EVENT_VAL_MAX) parse_event(val) else NA_real_
        if (!is.na(tk)) {
          if (!is.na(tv)) { token_time[[tk]] <- put(token_time[[tk]], key, s, tv); src <- c(src, sprintf("%s[%s]", key, where)) }
          if (!is.na(ev)) { token_event[[tk]] <- put(token_event[[tk]], key, s, ev); src <- c(src, sprintf("%s[%s]", key, where)) }
        } else {
          if (!is.na(tv) && grepl(GEN_TIME_OK, k_np, perl = TRUE)) {
            gen_time <- put(gen_time, key, s, tv); src <- c(src, sprintf("%s[%s]", key, where))
          }
          if (!is.na(ev) && grepl(GEN_EVENT_OK, k_np, perl = TRUE)) {
            gen_event <- put(gen_event, key, s, ev); src <- c(src, sprintf("%s[%s]", key, where))
          }
        }
      }
      cb <- parse_combined(cell)                     # "Deceased after 27 months"
      if (!is.null(cb)) {
        k <- paste0("value_text[", cb$token, "]")
        token_time[[cb$token]]  <- put(token_time[[cb$token]], k, s, cb$time)
        token_event[[cb$token]] <- put(token_event[[cb$token]], k, s, cb$event)
        src <- c(src, sprintf("value_text[%s]", where))
      }
    }
  }

  drop_empty <- function(bag) bag[vapply(bag, function(v) any(!is.na(v)), logical(1))]
  token_time <- lapply(token_time, drop_empty); token_event <- lapply(token_event, drop_empty)
  gen_time <- drop_empty(gen_time); gen_event <- drop_empty(gen_event)
  # 「最大一致时间+状态对」选择器（审计 §4 第 3 条）：
  #   旧实现各自挑一个时间袋/事件袋，且只在 token 袋为空时才回退到通用袋 —— 只要
  #   token 事件袋里有 1 个可解码的 0/1 值（GSE183088 的 survival (months)=1），
  #   就不再回退，结果报成 1/1（实为 86/75）。
  #   现在把「该 token 的专属袋 ∪ 通用袋」全部两两组合，取**重叠样本数最大**的一对；
  #   平手时优先 token 专属 > 通用，其次 has_hint()（SURV_HINT 仅作排序提示），最后按键名。
  infer_pairs <- function(tbags, ebags) {
    Tb <- c(tbags, gen_time); Eb <- c(ebags, gen_event)
    if (!length(Tb) || !length(Eb)) return(NULL)
    tpref <- c(rep(0L, length(tbags)), rep(1L, length(gen_time)))
    epref <- c(rep(0L, length(ebags)), rep(1L, length(gen_event)))
    tsc <- vapply(Tb, score_time, integer(1)); esc <- vapply(Eb, score_event, integer(1))
    ta <- which(tsc > 0L); ea <- which(esc > 0L)
    if (!length(ta) || !length(ea)) return(NULL)
    Tm <- vapply(Tb[ta], function(v) !is.na(v), logical(n_samples)) + 0
    Em <- vapply(Eb[ea], function(v) !is.na(v), logical(n_samples)) + 0
    if (is.null(dim(Tm))) Tm <- matrix(Tm, ncol = 1L)
    if (is.null(dim(Em))) Em <- matrix(Em, ncol = 1L)
    ov <- crossprod(Tm, Em)
    th <- as.integer(has_hint(names(Tb)[ta])); eh <- as.integer(has_hint(names(Eb)[ea]))
    grid <- expand.grid(ii = seq_along(ta), jj = seq_along(ea))
    grid$n <- ov[cbind(grid$ii, grid$jj)]
    grid <- grid[grid$n > 0L, , drop = FALSE]
    if (!nrow(grid)) return(NULL)
    grid$pref <- tpref[ta][grid$ii] + epref[ea][grid$jj]
    grid$hint <- th[grid$ii] + eh[grid$jj]
    grid <- grid[order(-grid$n, grid$pref, -grid$hint, ta[grid$ii], ea[grid$jj]), , drop = FALSE]
    a <- ta[grid$ii[1]]; b <- ea[grid$jj[1]]
    ok <- !is.na(Tb[[a]]) & !is.na(Eb[[b]])
    list(effective_n = sum(ok), events = sum(Eb[[b]][ok] == 1),
         time_key = names(Tb)[a], event_key = names(Eb)[b])
  }

  per_token <- list()
  for (tk in c("OS", "DSS", "DFS", "RFS", "PFS", "MFS")) {
    pt <- infer_pairs(if (is.null(token_time[[tk]])) list() else token_time[[tk]],
                      if (is.null(token_event[[tk]])) list() else token_event[[tk]])
    if (is.null(pt)) next
    per_token[[tk]] <- pt
  }
  list(series = series, platform = platform, n_samples = n_samples, title = title %||% "",
       summary = summ, supp = supp, per_token = per_token, sources = unique(src),
       ref_like = sum(ref_like))
}

# ---------------------------------------------------------------- 主循环
files <- sort(list.files(RAW, pattern = "\\.gz$", full.names = FALSE))
cat("scanning", length(files), "files under", RAW, "\n")
res <- vector("list", length(files)); errs <- character(0)
for (i in seq_along(files)) {
  f <- files[i]
  r <- tryCatch(scan_file(file.path(RAW, f)),
                error = function(e) { errs <<- c(errs, sprintf("%s: %s", f, conditionMessage(e))); NULL })
  if (is.null(r)) { errs <- c(errs, sprintf("%s: not a parsable series matrix", f)); next }
  best <- NULL
  for (tk in names(r$per_token)) {
    pt <- r$per_token[[tk]]
    if (is.null(best) || pt$effective_n > best$effective_n ||
        (pt$effective_n == best$effective_n && pt$events > best$events)) best <- c(pt, list(token = tk))
  }
  res[[i]] <- data.frame(
    local_file = f, series = r$series, platform = r$platform, n_samples = r$n_samples,
    effective_n = if (is.null(best)) 0L else best$effective_n,
    events = if (is.null(best)) 0L else best$events,
    token = if (is.null(best)) NA_character_ else best$token,
    tokens_all = if (length(r$per_token)) paste(sprintf("%s(%d/%d)", names(r$per_token),
                                                       vapply(r$per_token, function(z) z$effective_n, integer(1)),
                                                       vapply(r$per_token, function(z) z$events, integer(1))), collapse = ";") else NA_character_,
    time_key = if (is.null(best)) NA_character_ else best$time_key,
    event_key = if (is.null(best)) NA_character_ else best$event_key,
    source_of_survival_info = if (is.null(best)) NA_character_ else
      sprintf("time=%s; event=%s", best$time_key, best$event_key),
    sources_all = paste(unique(r$sources), collapse = ";"),
    title = r$title, summary = substr(r$summary, 1, 600),
    supp_files = paste(r$supp, collapse = ";"),
    ref_like_samples = r$ref_like, stringsAsFactors = FALSE)
}
ok <- !vapply(res, is.null, logical(1))
all_res <- do.call(rbind, res[ok])

gpl_tokens <- function(x) unique(trimws(unlist(strsplit(gsub("[+/,]", " ", as.character(x)), " "))))
# in_catalog = 该 GSE 已在 catalog 里（<series>_<GPL> 精确命中，或 <series> 本身命中）。
# 后者常见于"本地文件带 -GPL 后缀、catalog 只写 series"的情况；若 catalog 的 GPL 与本地
# 平台不一致，记到 catalog_gpl 里，供报告标注（重复 deposit 或 catalog 标签待核）。
all_res$in_catalog <- vapply(seq_len(nrow(all_res)), function(i) {
  s <- all_res$series[i]; p <- all_res$platform[i]
  if (!is.na(p) && paste0(s, "_", p) %in% di$Accession) return(TRUE)
  s %in% di$Accession
}, logical(1))
all_res$catalog_gpl <- vapply(seq_len(nrow(all_res)), function(i) {
  s <- all_res$series[i]; p <- all_res$platform[i]
  hit <- di$Accession == s | (!is.na(p) & di$Accession == paste0(s, "_", p))
  paste(setdiff(unique(di$GPL[hit]), NA), collapse = ",")
}, character(1))
all_res$gpl_mismatch <- all_res$in_catalog & vapply(seq_len(nrow(all_res)), function(i) {
  p <- all_res$platform[i]
  is.na(p) || !p %in% gpl_tokens(all_res$catalog_gpl[i])
}, logical(1))

# 既往排除清单（pipeline/out/excluded.csv）：候选若在里面，标注原因
excl <- NULL
if (file.exists(file.path(OUT, "excluded.csv"))) {
  ex <- read.csv(file.path(OUT, "excluded.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("Accession", "Reason") %in% colnames(ex))) excl <- ex
}
excl_reason <- function(s) if (is.null(excl)) NA_character_ else {
  r <- excl$Reason[excl$Accession == s]; if (length(r)) r[1] else NA_character_ }

# ---------------------------------------------------------------- 交叉核对既有候选清单
extra <- list()
for (cf in c("candidates_final_20260907.csv", "candidates_surv_scan.csv")) {
  fp <- file.path(OUT, cf)
  if (!file.exists(fp)) next
  cc <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  gse_col <- if ("GSE" %in% colnames(cc)) "GSE" else if ("series" %in% colnames(cc)) "series" else NA
  if (is.na(gse_col)) next
  gp <- if ("GPL" %in% colnames(cc)) cc$GPL else if ("platform" %in% colnames(cc)) cc$platform else rep(NA, nrow(cc))
  for (i in seq_len(nrow(cc))) {
    gse <- trimws(as.character(cc[[gse_col]][i])); gpls <- gpl_tokens(gp[i])
    gpls <- gpls[grepl("^GPL[0-9]+$", gpls)]
    for (g in if (length(gpls)) gpls else NA_character_) {
      hit <- which(all_res$series == gse & (is.na(g) | all_res$platform == g))
      if (!length(hit)) hit <- which(all_res$series == gse)
      extra[[length(extra) + 1L]] <- data.frame(
        list_name = cf, series = gse, platform = g, 
        local_file = if (length(hit)) paste(unique(all_res$local_file[hit]), collapse = ";") else NA_character_,
        scanned_effective_n = if (length(hit)) max(all_res$effective_n[hit]) else NA_integer_,
        in_catalog = (paste0(gse, "_", g) %in% di$Accession || gse %in% di$Accession),
        stringsAsFactors = FALSE)
    }
  }
}
extra_df <- if (length(extra)) do.call(rbind, extra) else NULL

# ---------------------------------------------------------------- 输出主清单
sel <- all_res[all_res$effective_n > MIN_EFF & !all_res$in_catalog, , drop = FALSE]
sel <- sel[order(-sel$effective_n, sel$series), , drop = FALSE]
out_csv <- data.frame(
  local_file = sel$local_file, series = sel$series, platform = sel$platform,
  n_samples = sel$n_samples, effective_n = sel$effective_n, events = sel$events,
  token = sel$token, source_of_survival_info = sel$source_of_survival_info,
  in_catalog = sel$in_catalog, stringsAsFactors = FALSE)
write.csv(out_csv, file.path(OUT, "geo_candidates_gt50.csv"), row.names = FALSE)
write.csv(all_res, file.path(OUT, "geo_scan_all.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 癌种推断
cancer_of <- function(txt) {
  t <- tolower(txt)
  pat <- c(
    "Glioma/Brain" = "gliom|astrocyt|oligodendro|gbm|brain|medulloblast|meningioma|ependymom",
    "Breast" = "breast|mammary|brca",
    "Lung" = "lung|pulmonary|nsclc|luad|lusc",
    "Colorectal" = "colorect|colon|rectal|rectum",
    "Gastric/Esophageal" = "gastric|stomach|esophag|oesophag|barrett",
    "Liver" = "hepatocell|liver|hcc|cholangiocarcin",
    "Pancreatic" = "pancrea",
    "Leukemia" = "leukemi|leukem|aml|cml|myelodysplas|myeloid",
    "Lymphoma/Myeloma" = "lymphoma|myeloma|hodgkin|dlbcl|follicular lymph",
    "Melanoma/Skin" = "melanoma|cutaneous",
    "Renal/Kidney" = "renal|kidney|nephro|wilms",
    "Bladder/Urothelial" = "bladder|urothelial|transitional cell",
    "Ovarian" = "ovarian|ovary",
    "Endometrial/Cervical/Uterine" = "endometri|uterine|uterus|cervix|cervical",
    "Prostate" = "prostate",
    "Head and Neck" = "head and neck|oral|tongue|larynx|pharyn|hnscc",
    "Sarcoma/Bone" = "sarcoma|osteosarc|bone|ewing|liposarc",
    "Thyroid" = "thyroid",
    "Neuroendocrine/Other" = "neuroblast|retinoblast|pheochrom|paragangliom|adrenocortical|mesothelio|testicular|germ cell|thymoma|gastroenteropancreatic|neuroendocrine")
  for (nm in names(pat)) if (grepl(pat[[nm]], t)) return(nm)
  "unknown"
}
sel$cancer <- vapply(seq_len(nrow(sel)), function(i) cancer_of(paste(sel$series[i], sel$title[i])), character(1))

# ---------------------------------------------------------------- markdown
md <- c("# data/raw 离线生存信息扫描：有效样本 > 50 且未编目的 GEO 候选",
        "", sprintf("扫描时间: %s | 扫描 data/raw/*.gz: %d 份 | 解析失败: %d 份",
                    format(Sys.time()), length(files), length(errs)),
        sprintf("阈值: effective_n > %d（同一 token 下\"时间可解析 + 事件可解码\"的样本数）且 Accession 不在 data/dataset_info.csv", MIN_EFF),
        "", "## 0. 方法与局限", "",
        "- 只看 series matrix 的 header（`!Sample_characteristics_ch1` 各行 + `!Sample_description` 各行 + 系列级信息），不解析表达谱矩阵。",
        "- key 判定端点 token，值决定角色：能解析出数字 -> 时间候选；能解码成 0/1 -> 事件候选。",
        "- **未做单位归一**（月/年/天混用）：本表只报有效样本数与事件数，入库前需在 `03_surv_table.R` 的 specs 里写清单位换算。",
        "- `token` 列是多 token 中有效样本数最大者；完整分布见 `geo_scan_all.csv` 的 `tokens_all` 列（格式 `TOKEN(有效n/事件数)`）。",
        "- `!Series_supplementary_file` 只登记文件名，**未下载**；若临床表在补充文件里（本批命中多数不是），需另行人工下载。",
        "- 本表只覆盖 data/raw 已有的 series matrix；既有候选清单里若有本地没有文件的 GSE，无法离线给有效样本数（见第 4 节）。",
        "", "## 1. 命中清单（按有效样本数降序）", "",
        "| local_file | series | platform | n_samples | effective_n | events | token | 生存信息来源 |",
        "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(sel)))
  md <- c(md, sprintf("| %s | %s | %s | %d | %d | %d | %s | %s |",
                      sel$local_file[i], sel$series[i], sel$platform[i], sel$n_samples[i],
                      sel$effective_n[i], sel$events[i], sel$token[i], sel$source_of_survival_info[i]))
md <- c(md, "", "## 2. 按癌种分组（由 series title/summary 关键词推断，仅用于排序）", "")
for (cn in sort(unique(sel$cancer))) {
  sub <- sel[sel$cancer == cn, ]
  md <- c(md, sprintf("### %s（%d 份）", cn, nrow(sub)), "",
          "| series | platform | effective_n | events | token | title |", "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(sub)))
    md <- c(md, sprintf("| %s | %s | %d | %d | %s | %s |", sub$series[i], sub$platform[i],
                        sub$effective_n[i], sub$events[i], sub$token[i],
                        substr(gsub("\\|", "/", sub$title[i]), 1, 90)))
  md <- c(md, "")
}

md <- c(md, "## 3. 建议排除 / 需人工确认", "",
        "### 3.1 同一 GSE 的多平台重复 deposit（同一批病人，只需保留注释最全的一个平台）", "")
dup_series <- unique(sel$series[sel$series %in% names(which(table(all_res$series) > 1))])
if (length(dup_series)) for (s_ in dup_series) {
  sub <- all_res[all_res$series == s_, ]
  md <- c(md, sprintf("- **%s**: %s", s_, paste(sprintf("%s（%s）=有效%d", sub$platform, sub$local_file, sub$effective_n), collapse = ", ")))
} else md <- c(md, "- 无未编目的多文件 series")
md <- c(md, "- 人工已知同源对：GSE1378（microdissected）与 GSE1379（whole tissue）是同一批 60 例他莫昔芬患者的两种组织制备，只应取其一。")

md <- c(md, "", "### 3.2 两色参考设计（样本 source_name/title 含 reference/universal/pool）——表达值是 log ratio，建议排除", "")
refhits <- all_res[all_res$effective_n > MIN_EFF & all_res$ref_like_samples > 0.5 * all_res$n_samples &
                      (!all_res$in_catalog | all_res$gpl_mismatch), ]
if (nrow(refhits)) for (i in seq_len(nrow(refhits)))
  md <- c(md, sprintf("- **%s**（%s）：n=%d，其中 %d 个样本标注为参考 RNA/Universal Reference",
                      refhits$local_file[i], refhits$series[i], refhits$n_samples[i], refhits$ref_like_samples[i])) else
  md <- c(md, "- 无")

md <- c(md, "", "### 3.3 已在 catalog（series 命中，但本地平台与 catalog 的 GPL 标签不一致）——视为重复，不建议再登记", "")
mism <- all_res[all_res$effective_n > MIN_EFF & all_res$gpl_mismatch, ]
if (nrow(mism)) for (i in seq_len(nrow(mism)))
  md <- c(md, sprintf("- %s：本地 %s（%s）=有效%d，catalog 记为 GPL=%s",
                      mism$series[i], mism$platform[i], mism$local_file[i], mism$effective_n[i], mism$catalog_gpl[i])) else
  md <- c(md, "- 无")

md <- c(md, "", "### 3.4 本次扫到生存信息、但 `pipeline/out/excluded.csv` 里曾判定为排除的 series —— 需复核是否翻案", "")
rev_ <- sel[!is.na(vapply(sel$series, function(x) excl_reason(x), character(1))), ]
if (nrow(rev_)) for (i in seq_len(nrow(rev_)))
  md <- c(md, sprintf("- **%s**（%s, %s）：本次 effective_n=%d / events=%d / token=%s；既往排除理由：%s",
                      rev_$series[i], rev_$platform[i], rev_$local_file[i], rev_$effective_n[i], rev_$events[i],
                      rev_$token[i], excl_reason(rev_$series[i]))) else
  md <- c(md, "- 无")

md <- c(md, "", "### 3.5 events == effective_n（事件列可能是稀疏编码：只给事件病例写了值）——数值需人工核对", "")
sparse <- sel[sel$events == sel$effective_n & sel$effective_n > 0, ]
if (nrow(sparse)) for (i in seq_len(nrow(sparse)))
  md <- c(md, sprintf("- %s（%s, %s）：effective_n=events=%d（token=%s）",
                      sparse$local_file[i], sparse$series[i], sparse$token[i], sparse$effective_n[i], sparse$token[i])) else
  md <- c(md, "- 无")

md <- c(md, "", "### 3.6 effective_n 达标但 events < 10（可入库，但 KM/多因素分析功效有限）", "")
fewev <- sel[sel$events < 10, ]
if (nrow(fewev)) for (i in seq_len(nrow(fewev)))
  md <- c(md, sprintf("- %s（%s, %s）：effective_n=%d，events=%d", fewev$local_file[i], fewev$series[i], fewev$token[i], fewev$effective_n[i], fewev$events[i])) else
  md <- c(md, "- 无")

md <- c(md, "", "### 3.7 命中清单涉及的 series 级补充文件（未下载，按需人工取）", "")
wit <- sel[nzchar(sel$supp_files), ]
if (nrow(wit)) for (i in seq_len(nrow(wit)))
  md <- c(md, sprintf("- %s: %s", wit$series[i], gsub(";", " ; ", wit$supp_files[i]))) else
  md <- c(md, "- 命中清单里的文件都没有列 supplementary file（生存信息只在 series matrix 内）")

if (!is.null(extra_df)) {
  md <- c(md, "", "## 4. 与既有候选清单交叉核对", "",
          "`candidates_final_20260907.csv` / `candidates_surv_scan.csv` 里的每个 (GSE, GPL)。",
          "第 5 列含义：本地无文件 / 本次扫描的有效样本数。", "",
          "| 清单 | series | platform | 本地文件 | 本次扫描 effective_n | 已在 catalog |", "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(extra_df)))
    md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |", extra_df$list_name[i], extra_df$series[i],
                        extra_df$platform[i], extra_df$local_file[i],
                        ifelse(is.na(extra_df$local_file[i]), "本地无文件",
                               ifelse(is.na(extra_df$scanned_effective_n[i]), "未扫描到", extra_df$scanned_effective_n[i])),
                        extra_df$in_catalog[i]))
  inl <- extra_df[!extra_df$in_catalog & !is.na(extra_df$scanned_effective_n) & extra_df$scanned_effective_n > MIN_EFF, ]
  md <- c(md, "", sprintf("其中满足 >%d 且未编目的: %d 条记录 / %d 个 series（%s）",
                          MIN_EFF, nrow(inl), length(unique(inl$series)),
                          if (nrow(inl)) paste(unique(inl$series), collapse = ", ") else "无"),
          "若某条目未出现在第 1 节：本地没有该 series matrix（无法离线评估），或本地有文件但 effective_n 未达阈值。")
}
if (length(errs)) md <- c(md, "", "## 5. 解析失败/跳过的文件", "", paste0("- ", errs))
writeLines(md, file.path(OUT, "geo_candidates_gt50.md"))

cat(sprintf("scanned=%d  errors=%d  candidates(>%d, not catalogued)=%d\n", nrow(all_res), length(errs), MIN_EFF, nrow(sel)))
cat("wrote pipeline/out/geo_candidates_gt50.csv | geo_candidates_gt50.md | geo_scan_all.csv\n")
if (nrow(sel)) print(sel[, c("series", "platform", "n_samples", "effective_n", "events", "token")])
}

# ---------------------------------------------------------------------------
# run_37_geo_expansion_screen()  <-  verbatim pipeline/R/37_geo_expansion_screen.R
# ---------------------------------------------------------------------------
run_37_geo_expansion_screen <- function() {
# 37_geo_expansion_screen.R ---------------------------------------------------
# 目的
#   为 data/dataset_info.csv 里"有 TCGA 队列但没有 GEO 队列"的癌种（15 个）筛选
#   GEO 候选队列，产出可审计的候选清单：
#     pipeline/out/geo_expansion_candidates.csv   一行 = 一个 (series, platform) 候选
#     pipeline/out/geo_expansion_candidates.md     按癌种分组的报告
#
# 输入（全部为**已获取的元数据**，脚本本身不联网、不下载任何矩阵）
#   pipeline/out/geo_expansion_headers/*.header.txt
#       由 /tmp/geo/fetch_all.py 用 HTTP Range + 部分 gzip 解压取得的 series matrix
#       **头部**（到 `!series_matrix_table_begin` 为止）。表达谱矩阵体从未下载。
#   pipeline/out/geo_expansion_series_meta.csv
#       NCBI E-utilities esummary(db=gds) 的系列级元数据 + 该系列命中的癌种/是否命中
#       survival 检索式。
#
# 解析器
#   直接照搬 pipeline/R/36_screen_geo_candidates.R 的 header 解析逻辑
#   （strip_quotes / split_top / parse_pairs / token_of_key / parse_event /
#     parse_time / parse_combined），以保证与既有扫描口径一致。
#   关键坑：GEO series matrix 为 CRLF 结尾，一律 trimws() 去 \r，再剥成对引号，
#   不要自造 whitespace 正则，否则 $ 锚点静默失效。此处输入已是普通文本文件。
#
# 用法: Rscript pipeline/R/37_geo_expansion_screen.R [ROOT]
# ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

ROOT <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(ROOT)) ROOT[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
OUT  <- file.path(ROOT, "pipeline/out")
HDR  <- file.path(OUT, "geo_expansion_headers")
MIN_EFF <- 50L                       # 门槛：有效样本数必须 **大于** 50

di <- read.csv(file.path(ROOT, "data/dataset_info.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cat_names <- unique(di$Type)
# 15 个缺口癌种（catalog 里有 TCGA、没有 GEO）
GAP_TYPES <- c("Kidney Cancer", "Adrenocortical Cancer", "Endometrial Cancer",
               "Head and Neck Cancer", "Melanoma", "Mesothelioma", "Pheochromocytoma",
               "Sarcoma", "Testicular Cancer", "Thymoma", "Thyroid Cancer",
               "Uterine Carcinosarcoma", "Uveal Melanoma",
               "Cholangiocarcinoma", "Lymphoma")
smeta <- read.csv(file.path(OUT, "geo_expansion_series_meta.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
rownames(smeta) <- smeta$gse

# ---------------------------------------------------------------- 基础工具（照搬 36）
strip_quotes <- function(x) {
  x <- trimws(x)                                  # 去掉 CRLF 的 \r（关键）
  x <- sub('^"', '', x); x <- sub('"$', '', x)
  x <- gsub('""', '"', x, fixed = TRUE)
  trimws(x)
}

split_top <- function(s, seps = c(";", ",")) {
  ch <- strsplit(s, "", fixed = TRUE)[[1]]
  if (!length(ch)) return(character(0))
  depth <- 0L; cur <- character(0); out <- character(0)
  for (cc in ch) {
    if (cc == "(") depth <- depth + 1L
    else if (cc == ")") depth <- max(0L, depth - 1L)
    if (depth == 0L && cc %in% seps) { out <- c(out, paste(cur, collapse = "")); cur <- character(0) }
    else cur <- c(cur, cc)
  }
  out <- trimws(c(out, paste(cur, collapse = "")))
  out[nzchar(out)]
}

# 审计结论（pipeline/out/geo_expansion_AUDIT.md §4）：
#   * SURV_HINT 曾当**硬门槛**（!has_hint(cell) -> 整个单元格丢弃），会静默丢掉真实
#     队列：GSE304059 的 `futime`/`fustat`、GSE22138 的 `metastasis` —— key 名里不含
#     任何 hint 词，于是这些行从未进入解析器。现在 SURV_HINT 只作**同分排序提示**
#     （见 infer_pairs()），不再当门槛。
#   * CRLF：审计实测本批 743 份 header 已由抓取步骤归一为 LF（只有 2 份在
#     !Sample_data_processing 的**值**里含字面 \r，不是行结尾），故此处继续用普通
#     trimws() 去两端空白，不另造 whitespace 正则 —— 该缺陷在本批未被数据触发。
SURV_HINT <- "(surv|follow|vital|death|dead|die|died|recur|relaps|progress|event|status|indicator|month|year|day|dfs|rfs|pfs|dss|mfs|\\bos\\b|surtim|surind|outcome|metasta|futime|fu[ _-]?stat)"
has_hint <- function(x) !is.na(x) & grepl(SURV_HINT, x, ignore.case = TRUE, perl = TRUE)

parse_pairs <- function(cell) {
  cell <- strip_quotes(cell)
  if (!nzchar(cell)) return(NULL)
  k <- character(0); v <- character(0)
  for (s in split_top(cell)) {
    ch <- strsplit(s, "", fixed = TRUE)[[1]]
    depth <- 0L; eqs <- integer(0); cols <- integer(0)
    for (i in seq_along(ch)) {
      cc <- ch[i]
      if (cc == "(") depth <- depth + 1L
      else if (cc == ")") depth <- max(0L, depth - 1L)
      else if (depth == 0L && cc == "=") eqs <- c(eqs, i)
      else if (depth == 0L && cc == ":") cols <- c(cols, i)
    }
    pos <- if (length(eqs)) eqs[1] else if (length(cols)) cols[1] else NA_integer_
    if (is.na(pos) || pos < 2L || pos >= length(ch)) next
    key <- trimws(paste(ch[seq_len(pos - 1)], collapse = ""))
    val <- trimws(paste(ch[(pos + 1):length(ch)], collapse = ""))
    key <- gsub("^[\\[{(]+", "", key); key <- gsub("[\\]}]+$", "", key)
    val <- gsub("^[\\[{(]+", "", val); val <- gsub("[\\]}]+$", "", val)
    val <- sub("^([-+]?[0-9]+(?:\\.[0-9]+)?)\\s*\\(.*$", "\\1", val)
    if (nzchar(key)) { k <- c(k, key); v <- c(v, trimws(val)) }
  }
  if (!length(k)) NULL else data.frame(key = k, value = v, stringsAsFactors = FALSE)
}

TOKEN_RE <- list(
  DSS = "(^|[^a-z])(dss|css|bcss)([^a-z]|$)|disease[ _-]?specific|disease[ _-]?related|cause[ _-]?specific",
  DFS = "(^|[^a-z])(dfs|dfi|efs)([^a-z]|$)|disease[ _-]?free",
  RFS = "(^|[^a-z])(rfs)([^a-z]|$)|recurrence[ _-]?free|relapse[ _-]?free|recurrence|relaps|recur",
  PFS = "(^|[^a-z])(pfs|pfi)([^a-z]|$)|progression[ _-]?free|progress",
  MFS = "(^|[^a-z])(mfs|drfs)([^a-z]|$)|metasta|distant[ _-]?recurrence",
  OS  = "(^|[^a-z])(os)([^a-z]|$)|overall[ _-]?survival|survival|surv|vital|death|dead|died|die|surtim|surind")
token_of_key <- function(key) {
  k <- tolower(key)
  for (tk in names(TOKEN_RE)) if (grepl(TOKEN_RE[[tk]], k, perl = TRUE)) return(tk)
  NA_character_
}

EXCLUDE_RE <- paste0("age|braf|margin|residual|menopause|race|ethnic|smok|hpv|histolog|",
                     "diagnos|\\bdx\\b|pathologic|birth|surgery|operation|treat|therap|drug|dose|",
                     "batch|plate|scan|center|institut|country|gender|^sex|height|weight|stage|grade|",
                     "site|tumor[ _-]?status|tissue|cell|culture|passage|protocol|keyword|date|",
                     "group|subgrp|cluster|class|subtype|response|responder|censor[ _-]?reason|",
                     "sample[ _-]?id|record[ _-]?id|patient[ _-]?id|^id$|number|prior[ _-]?lines")
excluded_key <- function(key) grepl(EXCLUDE_RE, tolower(key), perl = TRUE)

TIME_RE  <- "time|month|year|day|follow"
EVENT_RE <- "event|status|vital|death|dead|died|die|recur|relaps|progress|indicator|surv|outcome|censor|metasta"
BARE_RE  <- "^(dfs|rfs|os|pfs|pfi|dss|mfs|dfi|efs|drfs|surv|survival|surtim|surind|status|event)$"
# 通用（非 token 专属）时间/事件列的准入。审计修复（§4）：把「随访时间」缩写补进词表
# （futime / ftime / months / years / days / last fu / follow-up），否则 GSE304059 的
# `futime`+`fustat`、GSE22138 的 `months to endpoint`+`metastasis` 这类真实随访对在
# 通用袋里拿不到角色。判定前先剥掉括号内的说明文档（如 "doubling time (days)"），
# 以免把非随访的 "... (days)" 误当随访时间。
GEN_TIME_OK  <- "surv|follow|vital|surtim|overall|disease|recur|relaps|progress|metasta|dfs|rfs|pfs|dss|(^|[^a-z])os([^a-z]|$)|^(time|months|years|days|follow[ _-]?up)([[:space:]]*\\(.*\\))?$|(^|[^a-z])(futime|ftime|months?|years?|days?|fu|last[ _-]?fu|follow[ _-]?up)([^a-z]|$)"
GEN_EVENT_OK <- "surv|vital|death|dead|die|recur|relaps|progress|event|censor|follow|^(status|event|indicator)$|(^|[^a-z])(fustat|stat|fu[ _-]?stat|dead|died)([^a-z]|$)"
EVENT_VAL_MAX <- 40L
# 时间/事件袋打分：全 0/1 的列不算时间（否则 "status: 0:alive" 会被当时间）
score_time  <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return(0L); if (all(v %in% c(0, 1))) return(0L); length(v) }
score_event <- function(v) sum(!is.na(v))

NA_VAL_RE <- "^(na|n/?a|unknown|not available|not applicable|no more data|none|null|nd|--?|-|\\?|not determined|undetermined|not done|missing|no data)$"
NEG_RE <- "alive|^no\\b|^n$|^false|negative|free of|no evidence|not recur|not progress|did not|without|^0$|no relapse|no recurrence|no progression|no death|ned|tumor free|^non[ _-]?(recur|relaps|progress|event)"
POS_RE <- "dead|death|deceased|died|^yes|^y$|^true|positive|recur|relaps|progress|event|^1$|recurrence|metasta|decease"
parse_event <- function(v) {
  v <- tolower(strip_quotes(v))
  if (!nzchar(v) || grepl(NA_VAL_RE, v)) return(NA_real_)
  if (grepl("^-?[0-9.]+$", v)) { n <- suppressWarnings(as.numeric(v)); return(if (!is.na(n) && n %in% c(0, 1)) n else NA_real_) }
  if (grepl(NEG_RE, v)) return(0)
  if (grepl(POS_RE, v)) return(1)
  NA_real_
}
parse_time <- function(v) {
  v <- strip_quotes(v)
  m <- regmatches(v, regexpr("[-+]?[0-9]*\\.?[0-9]+", v))
  if (!length(m) || !nzchar(m)) return(NA_real_)
  suppressWarnings(as.numeric(m))
}
COMBINED_RE <- "(alive|dead|deceased|died|no evidence of disease|ned|recurrence|relapse|recurred|relapsed|progression|progressed|metastasis|metastatic)[^0-9]{0,30}?([0-9]+(\\.[0-9]+)?)\\s*(months?|years?|days?|weeks?)"
parse_combined <- function(v) {
  s <- tolower(strip_quotes(v))
  r <- regmatches(s, regexec(COMBINED_RE, s, perl = TRUE))[[1]]
  if (!length(r)) return(NULL)
  word <- r[2]; num <- suppressWarnings(as.numeric(r[3])); unit <- r[6]
  ev <- if (grepl("alive|no evidence|ned", word)) 0 else 1
  tk <- if (grepl("recur|relaps", word)) "RFS" else if (grepl("progress", word)) "PFS" else
        if (grepl("metasta", word)) "MFS" else "OS"
  list(token = tk, time = num, unit = unit, event = ev)
}

# ---------------------------------------------------------------- header 解析（改自 36 scan_file）
read_header <- function(path) {
  if (grepl("\\.gz$", path)) { con <- gzfile(path, "rt"); on.exit(close(con)); return(readLines(con, warn = FALSE)) }
  readLines(path, warn = FALSE)
}

scan_header <- function(path) {
  h <- read_header(path)
  tag <- sub("\t.*$", "", h)
  body <- sub("^[^\t]*\t?", "", h)
  split_vals <- function(i) {
    v <- strsplit(body[i], "\t", fixed = TRUE)[[1]]
    if (!length(v)) character(0) else vapply(v, strip_quotes, character(1), USE.NAMES = FALSE)
  }
  idx <- function(pattern) grep(pattern, tag, ignore.case = TRUE)
  first_val <- function(pattern) { i <- idx(pattern); if (!length(i)) NA_character_ else (split_vals(i[1])[1] %||% NA_character_) }

  series <- first_val("^!Series_geo_accession$")
  if (is.na(series) || !nzchar(series)) series <- sub("\\.header\\.txt(\\.gz)?$", "", basename(path))
  title <- first_val("^!Series_title$")
  si <- idx("^!Series_summary$")
  summ <- if (!length(si)) "" else paste(vapply(si, function(j) (split_vals(j)[1] %||% ""), character(1)), collapse = " ")
  fi <- idx("^!Series_supplementary_file")
  supp <- if (!length(fi)) character(0) else unique(unlist(lapply(fi, split_vals)))
  ri <- idx("^!Series_relation")
  rel <- if (!length(ri)) character(0) else unique(unlist(lapply(ri, split_vals)))
  oi <- idx("^!Series_overall_design$")
  design <- if (!length(oi)) "" else paste(vapply(oi, function(j) (split_vals(j)[1] %||% ""), character(1)), collapse = " ")
  sti <- idx("^!Series_type$")
  stype <- if (!length(sti)) "" else paste(vapply(sti, function(j) (split_vals(j)[1] %||% ""), character(1)), collapse = " ")
  gi <- idx("^!Sample_geo_accession$")
  n_samples <- if (length(gi)) length(split_vals(gi[1])) else NA_integer_
  ti <- idx("^!Sample_title$")
  titles <- if (length(ti)) split_vals(ti[1]) else character(0)
  sni <- idx("^!Sample_source_name_ch1$")
  src_names <- if (length(sni)) split_vals(sni[1]) else character(0)
  if (is.na(n_samples) || !n_samples) n_samples <- length(titles)
  if (!n_samples) return(NULL)

  field_idx <- sort(c(idx("^!Sample_characteristics_ch[0-9]*$"), idx("^!Sample_description$")))
  token_time <- list(); token_event <- list(); gen_time <- list(); gen_event <- list()
  src <- character(0); ref_like <- rep(FALSE, n_samples)
  ref_words <- "reference|universal|stratagene|pool|common reference|control rna|reference rna"
  if (length(src_names) == n_samples)
    ref_like <- ref_like | grepl(ref_words, src_names, ignore.case = TRUE)
  if (length(titles) == n_samples)
    ref_like <- ref_like | grepl(ref_words, titles, ignore.case = TRUE)

  put <- function(bag, key, s, val) {
    stopifnot(is.character(key), length(key) == 1L)
    val <- if (is.numeric(val) && length(val) == 1L) val else NA_real_
    if (is.null(bag[[key]]) || !is.numeric(bag[[key]]) || length(bag[[key]]) != n_samples)
      bag[[key]] <- rep(NA_real_, n_samples)
    bag[[key]][s] <- val
    bag
  }
  for (i in field_idx) {
    vals <- split_vals(i)
    if (length(vals) != n_samples) next
    # 审计修复：不再用 has_hint(vals)/has_hint(cell) 做硬门槛（会丢掉 futime/fustat 行）
    is_desc <- grepl("^!Sample_description$", tag[i], ignore.case = TRUE)
    where <- if (is_desc) "description" else "characteristics"
    for (s in seq_len(n_samples)) {
      cell <- vals[s]
      if (!nzchar(cell)) next
      if (grepl(ref_words, cell, ignore.case = TRUE)) ref_like[s] <- TRUE
      pr <- parse_pairs(cell)
      if (!is.null(pr)) for (r in seq_len(nrow(pr))) {
        key <- pr$key[r]; val <- pr$value[r]
        if (excluded_key(key)) next
        k <- tolower(key)
        tk <- token_of_key(key)
        # ---- 角色由「解出的 token/家族 + 值」决定，不再只看 key 文本（审计 §4 第 2 条）
        # 旧实现 is_time_k = TIME_RE(key) || BARE_RE(key) 之类，使 code_os / drfs /
        # tt.drfs / pfs_index / fustat 这类 key 虽已被 token_of_key() 识别出端点，
        # 却一个角色都拿不到。新规则：
        #   纯时间词 key（time/month/year/day/...） -> 只做时间；
        #   纯事件词 key（status/vital/death/...）   -> 只做事件；
        #   其余                                      -> 交给值决定（能解出数字进时间候选，
        #                                               能解出 0/1 进事件候选；全 0/1 的列
        #                                               由 score_time() 的惩罚挡在时间袋外）。
        k_np <- gsub("\\([^)]*\\)", " ", k)     # 括号内是说明文档，不参与通用词表判定
        time_hint  <- grepl(TIME_RE, k)
        event_hint <- grepl(EVENT_RE, k)
        tv <- if (time_hint  || !event_hint) parse_time(val) else NA_real_
        ev <- if ((event_hint || !time_hint) && nchar(val) <= EVENT_VAL_MAX) parse_event(val) else NA_real_
        if (!is.na(tk)) {
          if (!is.na(tv)) { token_time[[tk]] <- put(token_time[[tk]], key, s, tv); src <- c(src, sprintf("%s[%s]", key, where)) }
          if (!is.na(ev)) { token_event[[tk]] <- put(token_event[[tk]], key, s, ev); src <- c(src, sprintf("%s[%s]", key, where)) }
        } else {
          if (!is.na(tv) && grepl(GEN_TIME_OK, k_np, perl = TRUE)) {
            gen_time <- put(gen_time, key, s, tv); src <- c(src, sprintf("%s[%s]", key, where))
          }
          if (!is.na(ev) && grepl(GEN_EVENT_OK, k_np, perl = TRUE)) {
            gen_event <- put(gen_event, key, s, ev); src <- c(src, sprintf("%s[%s]", key, where))
          }
        }
      }
      cb <- parse_combined(cell)
      if (!is.null(cb)) {
        k <- paste0("value_text[", cb$token, "]")
        token_time[[cb$token]]  <- put(token_time[[cb$token]], k, s, cb$time)
        token_event[[cb$token]] <- put(token_event[[cb$token]], k, s, cb$event)
        src <- c(src, sprintf("value_text[%s]", where))
      }
    }
  }

  drop_empty <- function(bag) bag[vapply(bag, function(v) any(!is.na(v)), logical(1))]
  token_time <- lapply(token_time, drop_empty); token_event <- lapply(token_event, drop_empty)
  gen_time <- drop_empty(gen_time); gen_event <- drop_empty(gen_event)
  # 「最大一致时间+状态对」选择器（审计 §4 第 3 条）：
  #   旧实现各自挑一个时间袋/事件袋，且只在 token 袋为空时才回退到通用袋 —— 只要
  #   token 事件袋里有 1 个可解码的 0/1 值（GSE183088 的 survival (months)=1），
  #   就不再回退，结果报成 1/1（实为 86/75）。
  #   现在把「该 token 的专属袋 ∪ 通用袋」全部两两组合，取**重叠样本数最大**的一对；
  #   平手时优先 token 专属 > 通用，其次 has_hint()（SURV_HINT 仅作排序提示），最后按键名。
  infer_pairs <- function(tbags, ebags) {
    Tb <- c(tbags, gen_time); Eb <- c(ebags, gen_event)
    if (!length(Tb) || !length(Eb)) return(NULL)
    tpref <- c(rep(0L, length(tbags)), rep(1L, length(gen_time)))
    epref <- c(rep(0L, length(ebags)), rep(1L, length(gen_event)))
    tsc <- vapply(Tb, score_time, integer(1)); esc <- vapply(Eb, score_event, integer(1))
    ta <- which(tsc > 0L); ea <- which(esc > 0L)
    if (!length(ta) || !length(ea)) return(NULL)
    Tm <- vapply(Tb[ta], function(v) !is.na(v), logical(n_samples)) + 0
    Em <- vapply(Eb[ea], function(v) !is.na(v), logical(n_samples)) + 0
    if (is.null(dim(Tm))) Tm <- matrix(Tm, ncol = 1L)
    if (is.null(dim(Em))) Em <- matrix(Em, ncol = 1L)
    ov <- crossprod(Tm, Em)
    th <- as.integer(has_hint(names(Tb)[ta])); eh <- as.integer(has_hint(names(Eb)[ea]))
    grid <- expand.grid(ii = seq_along(ta), jj = seq_along(ea))
    grid$n <- ov[cbind(grid$ii, grid$jj)]
    grid <- grid[grid$n > 0L, , drop = FALSE]
    if (!nrow(grid)) return(NULL)
    grid$pref <- tpref[ta][grid$ii] + epref[ea][grid$jj]
    grid$hint <- th[grid$ii] + eh[grid$jj]
    grid <- grid[order(-grid$n, grid$pref, -grid$hint, ta[grid$ii], ea[grid$jj]), , drop = FALSE]
    a <- ta[grid$ii[1]]; b <- ea[grid$jj[1]]
    ok <- !is.na(Tb[[a]]) & !is.na(Eb[[b]])
    list(effective_n = sum(ok), events = sum(Eb[[b]][ok] == 1),
         time_key = names(Tb)[a], event_key = names(Eb)[b])
  }

  per_token <- list()
  for (tk in c("OS", "DSS", "DFS", "RFS", "PFS", "MFS")) {
    pt <- infer_pairs(if (is.null(token_time[[tk]])) list() else token_time[[tk]],
                      if (is.null(token_event[[tk]])) list() else token_event[[tk]])
    if (is.null(pt)) next
    per_token[[tk]] <- pt
  }
  any_time <- length(gen_time) > 0 || any(vapply(token_time, length, integer(1)) > 0)
  any_event <- length(gen_event) > 0 || any(vapply(token_event, length, integer(1)) > 0)
  list(series = series, n_samples = n_samples, title = title %||% "",
       summary = summ, design = design, stype = stype, supp = supp, rel = rel,
       any_time = any_time, any_event = any_event,
       per_token = per_token, sources = unique(src), ref_like = sum(ref_like),
       src_tab = if (length(src_names)) {
         tb <- sort(table(src_names), decreasing = TRUE)
         paste(sprintf("%s(%d)", names(tb), as.integer(tb)), collapse = " | ")
       } else "")
}

# ---------------------------------------------------------------- 主循环
files <- sort(list.files(HDR, pattern = "\\.header\\.txt(\\.gz)?$", full.names = FALSE))
cat("scanning", length(files), "series-matrix headers under", HDR, "\n")
res <- vector("list", length(files)); errs <- character(0)
for (i in seq_along(files)) {
  f <- files[i]
  r <- tryCatch(scan_header(file.path(HDR, f)),
                error = function(e) { errs <<- c(errs, sprintf("%s: %s", f, conditionMessage(e))); NULL })
  if (is.null(r)) { errs <- c(errs, sprintf("%s: not a parsable header", f)); next }
  gse <- r$series
  gpls <- character(0)
  if (gse %in% rownames(smeta)) {
    g <- as.character(smeta[gse, "gpl_list"])
    if (!is.na(g) && nzchar(g)) gpls <- unique(trimws(strsplit(g, "[+;, ]")[[1]]))
  }
  gpls <- paste0("GPL", gsub("^GPL", "", gpls)); gpls <- gpls[grepl("^GPL[0-9]+$", gpls)]
  if (!length(gpls)) gpls <- NA_character_

  best <- NULL
  for (tk in names(r$per_token)) {
    pt <- r$per_token[[tk]]
    if (is.null(best) || pt$effective_n > best$effective_n ||
        (pt$effective_n == best$effective_n && pt$events > best$events)) best <- c(pt, list(token = tk))
  }
  for (g in gpls) res[[length(res) + 1L]] <- data.frame(
    gse = gse, platform = g, header_file = f,
    n_samples = r$n_samples,
    effective_n = if (is.null(best)) 0L else best$effective_n,
    events = if (is.null(best)) 0L else best$events,
    token = if (is.null(best)) NA_character_ else best$token,
    tokens_all = if (length(r$per_token)) paste(sprintf("%s(%d/%d)", names(r$per_token),
                                                       vapply(r$per_token, function(z) z$effective_n, integer(1)),
                                                       vapply(r$per_token, function(z) z$events, integer(1))), collapse = ";") else NA_character_,
    time_key = if (is.null(best)) NA_character_ else best$time_key,
    event_key = if (is.null(best)) NA_character_ else best$event_key,
    source_of_survival_info = if (is.null(best)) NA_character_ else
      sprintf("time=%s; event=%s", best$time_key, best$event_key),
    sources_all = paste(unique(r$sources), collapse = ";"),
    any_time_column = r$any_time, any_event_column = r$any_event,
    title = r$title, summary = substr(r$summary, 1, 400),
    overall_design = substr(r$design, 1, 300), series_type = substr(r$stype, 1, 120),
    source_names = substr(r$src_tab, 1, 500),
    supp_files = paste(r$supp, collapse = ";"),
    series_relations = paste(head(r$rel, 20), collapse = ";"),
    ref_like_samples = r$ref_like, stringsAsFactors = FALSE)
}
if (!length(res)) stop("no headers parsed")
all_res <- do.call(rbind, res)
cat("series-platform rows:", nrow(all_res), " errors:", length(errs), "\n")

# ---------------------------------------------------------------- 附加判定
lookup <- function(gse, col) if (gse %in% rownames(smeta)) as.character(smeta[gse, col]) else NA_character_

all_res$gap_types   <- vapply(all_res$gse, lookup, character(1), col = "gap_types")
all_res$in_surv_query <- vapply(all_res$gse, lookup, character(1), col = "in_surv_query")
all_res$gdstype     <- vapply(all_res$gse, lookup, character(1), col = "gdstype")
all_res$pdat        <- vapply(all_res$gse, lookup, character(1), col = "pdat")

# 已在 catalog？
all_res$in_catalog <- vapply(seq_len(nrow(all_res)), function(i) {
  s <- all_res$gse[i]; p <- all_res$platform[i]
  (paste0(s, "_", p) %in% di$Accession) || (s %in% di$Accession)
}, logical(1))
all_res$catalog_type <- vapply(all_res$gse, function(s) {
  hit <- di$Type[di$Accession == s | sub("_.*$", "", di$Accession) == s]
  if (length(hit)) paste(unique(hit), collapse = ",") else NA_character_
}, character(1))

# （参考设计 / QC / 平台注释 / 排除标记 / 多平台去重 的判定统一放在下面的"平台技术分类"一节，
#   因为那些判定需要平台表和 taxon 字段。）

# ---------------------------------------------------------------- 平台技术分类
plat <- NULL
pf <- file.path(OUT, "geo_expansion_platforms.csv")
if (file.exists(pf)) plat <- read.csv(pf, check.names = FALSE, stringsAsFactors = FALSE)
all_res$platform_title <- ""
all_res$platform_technology <- ""
all_res$platform_organism <- ""
if (!is.null(plat)) for (i in seq_len(nrow(all_res))) {
  h <- which(plat$gpl == all_res$platform[i])
  if (!length(h)) next
  h <- h[1]
  pick <- function(col) { v <- plat[[col]][h]; if (length(v) == 1L && !is.na(v)) as.character(v) else "" }
  all_res$platform_title[i]      <- pick("platform_title")
  all_res$platform_technology[i] <- pick("technology")
  all_res$platform_organism[i]   <- pick("organism")
}
all_res$taxon <- vapply(all_res$gse, lookup, character(1), col = "taxon")
all_res$human <- grepl("Homo sapiens", all_res$taxon)
all_res$series_gdstype <- all_res$gdstype

# ---- 平台分类：先看平台标题/技术，再看 series 的 gdstype 兜底
pt <- tolower(paste(all_res$platform_title, all_res$platform_technology))
cls <- rep(NA_character_, nrow(all_res))
set_cls <- function(re, v) cls[is.na(cls) & grepl(re, pt, perl = TRUE)] <<- v
set_cls("methylat|450k|450 k|27k|epic|cpg|infinium|humanmethylation|help_promoter|promoter array", "methylation")
set_cls("\\bsnp\\b|cytosnp|omni|a?cgh|comparative genomic|copy number|genome tiling|genotyping|bac array", "snp_cgh")
set_cls("mirna|micro ?rna|mirbase|lna|mircury|mir-", "mirna")
set_cls("nanostring|ncounter", "ncounter")
set_cls("taqman|openarray|fluidigm|stem ?loop|realtime pcr|real-time pcr|rt-?pcr|qpcr|high.?throughput real", "qpcr")
set_cls("hiseq|nextseq|novaseq|genome analyzer|miseq|ion torrent|\\bbgi\\b|dnbseq|helicos|\\b454\\b|solid|pacbio|nanopore|sequencing", "rna_seq")
set_cls("hgu|hg-?u|u133|u95|primeview|hugene|clariom|humanht|humanwg|sentrix|beadchip|sureprint|agilent-0|whole genome|cdna|oligonucleotide|spotted|affymetrix|illumina", "expression_array")
# gdstype 兜底
gdt <- tolower(all_res$series_gdstype)
cls[is.na(cls) & grepl("methylation", gdt)] <- "methylation"
cls[is.na(cls) & grepl("genome variation|snp", gdt)] <- "snp_cgh"
cls[is.na(cls) & grepl("non-coding rna", gdt)] <- "mirna"
cls[is.na(cls) & grepl("high throughput sequencing", gdt)] <- "rna_seq"
cls[is.na(cls) & grepl("expression profiling by array", gdt)] <- "expression_array"
# 纠正：series 级 gdstype 若完全不含 "Expression profiling"，以 gdstype 为准
# （平台 technology 字段常写 "spotted DNA/cDNA"/"oligonucleotide"，会把 aCGH/HELP 甲基化平台误判成表达谱）
gdt2cls <- rep(NA_character_, nrow(all_res))
gdt2cls[grepl("methylation", gdt)] <- "methylation"
gdt2cls[grepl("genome variation|snp|comparative genomic", gdt)] <- "snp_cgh"
gdt2cls[grepl("non-coding rna", gdt)] <- "mirna"
gdt2cls[grepl("protein|antibody|genome binding|occupancy", gdt)] <- "other"
no_expr <- !grepl("expression profiling", gdt)
cls[no_expr & !is.na(gdt2cls)] <- gdt2cls[no_expr & !is.na(gdt2cls)]
cls[no_expr & is.na(gdt2cls)] <- "other"
cls[is.na(cls)] <- "unknown"
all_res$platform_class <- cls
EXPR_CLASSES <- c("expression_array", "rna_seq", "ncounter")
all_res$expression_assay <- ifelse(all_res$platform_class %in% EXPR_CLASSES,
  sprintf("YES(%s)", all_res$platform_class),
  sprintf("NO(%s: %s)", all_res$platform_class, substr(all_res$platform_title, 1, 40)))
all_res$expression_assay_ok <- all_res$platform_class %in% EXPR_CLASSES

# 参考设计：样本级 + series 级（two-colour / common reference / ratio data）
REF_SERIES_RE <- paste0("two[- ]colou?r|two[- ]channel|common reference|universal reference|",
                        "equal (mixture|combination) of all samples|ratio data|log ?ratio|",
                        "reference (sample|rna)|stratagene|pooled reference")
all_res$series_reference_text <- grepl(REF_SERIES_RE,
  paste(all_res$overall_design, all_res$summary, all_res$series_type), ignore.case = TRUE)
all_res$reference_design <- (all_res$ref_like_samples > 0.5 * all_res$n_samples) | all_res$series_reference_text

# SuperSeries（其样本就是各 SubSeries 的并集 → 与 SubSeries 重复计一次）
all_res$is_superseries <- grepl("SuperSeries is composed of the SubSeries", all_res$summary, ignore.case = TRUE)
all_res$subseries_of <- vapply(all_res$series_relations, function(z) {
  if (is.na(z) || !nzchar(z)) return("")
  m <- regmatches(z, regexpr("SubSeries of: *GSE[0-9]+", z))
  if (length(m)) sub(".*(GSE[0-9]+)$", "\\1", m[1]) else ""
}, character(1))
all_res$superseries_of <- vapply(all_res$series_relations, function(z) {
  if (is.na(z) || !nzchar(z)) return("")
  m <- regmatches(z, gregexpr("SuperSeries of: *GSE[0-9]+", z))[[1]]
  if (length(m)) paste(unique(sub(".*(GSE[0-9]+)$", "\\1", m)), collapse = ";") else ""
}, character(1))
# SuperSeries 只有在"它的某个 SubSeries 也在本表里且有效样本>50"时才算重复 deposit
sup_covered <- vapply(seq_len(nrow(all_res)), function(i) {
  if (!all_res$is_superseries[i] || !nzchar(all_res$superseries_of[i])) return(FALSE)
  subs <- strsplit(all_res$superseries_of[i], ";")[[1]]
  any(all_res$gse %in% subs & all_res$effective_n > MIN_EFF)
}, logical(1))
all_res$superseries_covered <- sup_covered

# 同一课题的重复 deposit：按规范化标题分组，**只在同样是可用表达平台的候选之间**判重
# （同一课题的表达臂与甲基化臂标题相同，但甲基化臂本来就被排除，不应连带排除表达臂）
norm_title <- function(t) { t <- tolower(t); t <- gsub("\\[[^]]*\\]", " ", t)
  t <- gsub("\\(v[0-9]+\\)|\\bv[0-9]+\\b", " ", t); t <- gsub("[^a-z0-9 ]", " ", t)
  t <- gsub("\\s+", " ", trimws(t)); t }
all_res$study_key <- vapply(all_res$title, norm_title, character(1))
all_res$study_group <- ""
all_res$same_study_as <- ""
viable <- all_res$human & all_res$expression_assay_ok
for (k in unique(all_res$study_key[viable])) {
  ii <- which(all_res$study_key == k & viable)
  if (length(ii) < 2) next
  ord <- ii[order(-all_res$effective_n[ii])]
  all_res$study_group[ord] <- k
  for (j in ord[-1]) all_res$same_study_as[j] <- all_res$gse[ord[1]]
}

# 泛癌 panel / 细胞系资源（被疾病关键词命中的公共资源，不是某一癌种的临床队列）
#   致命：明确是公共资源库；仅"摘要提到 cell line"的只做警告，因为很多真实患者队列的摘要也会提到细胞系验证
PANEL_RE <- paste0("\\bCCLE\\b|cancer cell line encyclopedia|expression project for oncology|",
                   "\\bexpO\\b|SAGE librar|cancer genome anatomy project|NCI-?60|",
                   "sanger cell line|cell line panel|cell lines from .* cancer")
all_res$panel_resource <- grepl(PANEL_RE, paste(all_res$title, all_res$summary), ignore.case = TRUE)
all_res$cell_line_mentioned <- grepl("cell line|cell lines|xenograft",
                                     paste(all_res$title, all_res$summary, all_res$source_names), ignore.case = TRUE)

# QC / 里程碑类人工制品（MAQC 式 dummy outcome）
QC_RE <- "microarray quality control|\\bMAQC\\b|quality control (project|study|phase)|spike[- ]?in|milestone|dilution series|titration series|technical replicate|platform comparison"
all_res$qc_artefact <- grepl(QC_RE, paste(all_res$title, all_res$summary), ignore.case = TRUE)

# 平台注释可用性：本地 data/gpl/<GPL>.rds，否则由 38脚本用 AnnoProbe::idmap 检查
gpl_rds <- file.path(ROOT, "data/gpl", paste0(all_res$platform, ".rds"))
all_res$annotation_local <- file.exists(gpl_rds)
af <- file.path(OUT, "geo_expansion_annotation_check.csv")
ann <- if (file.exists(af)) read.csv(af, check.names = FALSE, stringsAsFactors = FALSE) else
  data.frame(gpl = character(0), status = character(0), detail = character(0))
all_res$annotation_status <- vapply(seq_len(nrow(all_res)), function(i) {
  if (isTRUE(all_res$annotation_local[i])) return("local_rds")
  hit <- ann$status[ann$gpl == all_res$platform[i]]
  if (length(hit)) hit[1] else "not_tested"
}, character(1))

# 生存信息位置：矩阵内 / 补充文件里可能带临床表 / 补充文件只有 RAW.tar 之类
all_res$supp_files_n <- ifelse(nzchar(all_res$supp_files),
                               vapply(strsplit(all_res$supp_files, ";"), length, integer(1)), 0L)
CLIN_SUPP_RE <- "clinical|surviv|follow|patient|outcome|pheno|annotation|sample.?info|meta.?data|characteristic|os_|_os|dfs|rfs|pfs"
supp_has_clin <- vapply(seq_len(nrow(all_res)), function(i) {
  f <- all_res$supp_files[i]
  if (!nzchar(f)) return(FALSE)
  fs <- strsplit(f, ";")[[1]]
  fs <- fs[!grepl("RAW\\.tar|_raw\\.tar|matrix|norm|fpkm|rpkm|counts|expression|CEL|IDAT", fs, ignore.case = TRUE)]
  length(fs) > 0 && any(grepl(CLIN_SUPP_RE, fs, ignore.case = TRUE))
}, logical(1))
all_res$supp_clinical_hint <- supp_has_clin
all_res$survival_location <- ifelse(all_res$effective_n > 0, "series_matrix_header",
  ifelse(supp_has_clin, "supplement_possible_clinical",
  ifelse(all_res$supp_files_n > 0, "supplement_matrix_only(no clinical-hinted file)", "absent")))

# 排除标记
all_res$exclude_flag <- ""
all_res$exclude_reason <- ""
set_flag <- function(i, flag, reason) {
  all_res$exclude_flag[i] <<- if (nzchar(all_res$exclude_flag[i])) paste(all_res$exclude_flag[i], flag, sep = ";") else flag
  all_res$exclude_reason[i] <<- if (nzchar(all_res$exclude_reason[i])) paste(all_res$exclude_reason[i], reason, sep = " | ") else reason
}
for (i in seq_len(nrow(all_res))) {
  if (all_res$in_catalog[i]) set_flag(i, "DUPLICATE_CATALOG",
      sprintf("series already in data/dataset_info.csv (Type=%s)", all_res$catalog_type[i]))
  if (all_res$panel_resource[i]) set_flag(i, "PANEL_RESOURCE",
      "pan-cancer cell-line / public panel resource (not a single-cancer clinical cohort)")
  if (all_res$qc_artefact[i]) set_flag(i, "QC_ARTEFACT", "title/summary looks like a QC / technical / MAQC-style study")
  if (all_res$reference_design[i]) set_flag(i, "REFERENCE_DESIGN",
      sprintf("%d/%d samples reference-like; series text has two-colour/common-reference wording",
              all_res$ref_like_samples[i], all_res$n_samples[i]))
  if (!all_res$human[i]) set_flag(i, "NON_HUMAN", paste0("series organism = ", all_res$taxon[i]))
  if (!all_res$expression_assay_ok[i]) set_flag(i, "NON_EXPRESSION_ASSAY",
      sprintf("assay class = %s (from series gdstype / platform); platform = %s",
              all_res$platform_class[i], substr(all_res$platform_title[i], 1, 50)))
  if (all_res$is_superseries[i]) set_flag(i, "SUPERSERIES",
      if (isTRUE(all_res$superseries_covered[i]))
        sprintf("SuperSeries whose SubSeries (%s) is also >50 in this table -> duplicate", all_res$superseries_of[i])
      else sprintf("SuperSeries aggregating SubSeries %s (samples may be counted more than once if subseries are added later)",
                   all_res$superseries_of[i]))
  if (nzchar(all_res$same_study_as[i]) && all_res$same_study_as[i] != all_res$gse[i])
    set_flag(i, "SAME_STUDY_DUPLICATE",
            sprintf("same study/title as %s (overlapping patients likely)", all_res$same_study_as[i]))
}
# 同一 GSE 多平台：仅保留一个平台行，其余标 DUPLICATE_PLATFORM
all_res$keep_platform <- TRUE
for (s in unique(all_res$gse)) {
  ii <- which(all_res$gse == s)
  if (length(ii) < 2) next
  ord <- ii[order(-all_res$effective_n[ii], all_res$annotation_local[ii] == FALSE, all_res$platform[ii])]
  for (j in ord[-1]) { all_res$keep_platform[j] <- FALSE
    set_flag(j, "DUPLICATE_PLATFORM", sprintf("another platform (%s) of the same series is the primary row", all_res$platform[ord[1]])) }
}

# ---------------------------------------------------------------- 门槛与排序
# 门槛：effective_n > 50，且必须是一个真正可用的 0/1 事件列（events >= 1）
# 其余为致命排除：catalog 重复 / 参考设计 / QC / 非人 / 非表达谱 / SuperSeries
all_res$event_column_valid <- all_res$events >= 1
all_res$passes_gate <- all_res$effective_n > MIN_EFF & all_res$event_column_valid &
  !all_res$in_catalog & !all_res$reference_design & !all_res$qc_artefact &
  all_res$human & all_res$expression_assay_ok & !all_res$panel_resource &
  !(all_res$is_superseries & all_res$superseries_covered) &
  all_res$keep_platform & !nzchar(all_res$same_study_as)
# 定向 panel（nCounter / qPCR）基因覆盖有限：不是致命排除，但要标注
all_res$gene_panel_warning <- ""
all_res$gene_panel_warning[all_res$platform_class %in% c("ncounter", "qpcr")] <-
  sprintf("targeted panel (%s) - gene coverage limited, verify the panel size", all_res$platform_title[all_res$platform_class %in% c("ncounter", "qpcr")])
# 本次抓取的 series 是表内另一 SuperSeries 的 SubSeries
all_res$is_subseries_of_listed <- vapply(seq_len(nrow(all_res)), function(i) {
  z <- all_res$subseries_of[i]
  nzchar(z) && z %in% all_res$gse
}, logical(1))
for (i in which(all_res$is_subseries_of_listed))
  set_flag(i, "SUBSERIES_OF_LISTED_SUPER", sprintf("SubSeries of %s, which is also in this table (same patients)", all_res$subseries_of[i]))
# 只有时间列没有事件列（如 GSE22541 的 DFI）——不能入组，但值得记下来
all_res$time_without_event <- all_res$any_time_column & !all_res$any_event_column & all_res$effective_n == 0
for (i in which(all_res$time_without_event))
  set_flag(i, "TIME_ONLY_NO_EVENT", "survival time column present but no decodable 0/1 status column in the matrix header")

# 平台注释是否可得（决定"能否用现有流水线建库"）
ANN_OK <- c("local_rds", "annoprobe_pipe", "geo_annot_symbol", "gene_level_no_probe_map_needed")
all_res$annotation_ok <- all_res$annotation_status %in% ANN_OK |
  grepl("^annoprobe_pipe|^geo_annot_symbol|^gene_level_no_probe_map_needed", all_res$annotation_status)

for (i in which(all_res$effective_n > MIN_EFF & !all_res$event_column_valid))
  set_flag(i, "EVENT_COLUMN_UNPARSEABLE",
           sprintf("effective_n=%d but decodable events = 0: the status column is not a plain 0/1 (multi-level code) or is fully censored",
                   all_res$effective_n[i]))
all_res$few_events  <- all_res$passes_gate & all_res$events < 10
all_res$sparse_event_coding <- all_res$passes_gate & all_res$events == all_res$effective_n & all_res$effective_n > 0

# 癌种排序键：取该 series 命中的缺口癌种（多命中时按固定顺序合并）
GAP_ORDER <- GAP_TYPES
type_key <- function(gt) {
  if (is.na(gt) || !nzchar(gt)) return("zz_other")
  ts <- unique(sub("\\(surv\\)$", "", trimws(strsplit(gt, ";")[[1]])))
  ts <- ts[ts %in% GAP_ORDER]
  if (!length(ts)) return("zz_other")
  paste(ts[order(match(ts, GAP_ORDER))], collapse = "+")
}
all_res$cancer_type <- vapply(all_res$gap_types, type_key, character(1))
all_res$cancer_type_verified <- all_res$cancer_type
all_res$type_note <- ""
all_res$manual_drop <- FALSE

# 人工复核（脚本内显式登记，便于审计）：关键词检索会误命中（摘要里顺带提到别的瘤种）
MANUAL <- list(
  "GSE17118" = list(ct = "Mesothelioma+Sarcoma", drop = TRUE, drop_reason = "mixed-histology matrix (20 peritoneal mesothelioma + 24 liposarcoma + 16 MPNST in one 60-sample series): no single cancer type reaches the >50 gate",
    note = "MIXED HISTOLOGY: 20 peritoneal mesothelioma + 24 liposarcoma + 16 MPNST in one 60-sample matrix; per-type n is below the >50 gate for every type (source_name_ch1 breakdown)"),
  "GSE54303" = list(ct = "Lymphoma", note = "aCGH (genome variation) on 332 lymphoma biopsies - NOT an expression cohort"),
  "GSE53870" = list(ct = "Cholangiocarcinoma", note = "miRNA profiling by array (non-coding RNA), not mRNA expression"),
  "GSE102166" = list(ct = "Melanoma", note = "plasma EV miRNA qPCR-array; 53 samples = 28 patients x ~2 timepoints; 'Other' platform; PFS parsed from free text with events==effective_n"),
  "GSE119810" = list(ct = "Sarcoma+Uterine Carcinosarcoma", note = "CANINE (dog) mammary tumours, not human sarcoma/UCS"),
  "GSE302459" = list(ct = "Lymphoma", note = "CANINE (dog) DLBCL PBMC NanoString cohort, not human"),
  "GSE131961" = list(ct = "Kidney Cancer", note = "miRNA RT-PCR array; status column is a 3-level code (A/C/O) that decodes to 0 events -> event column unusable"),
  "GSE131960" = list(ct = "Kidney Cancer", note = "miRNA RT-PCR array; same 3-level A/C/O status code -> event column unusable"),
  "GSE52793" = list(ct = "Head and Neck Cancer", note = "oral-rinse DNA methylation array (GPL13534), not expression"),
  "GSE68717" = list(ct = "Head and Neck Cancer", note = "SNP array (LOH/copy number), not expression"),
  "GSE277573" = list(ct = "Head and Neck Cancer", note = "pretreatment blood methylation cytometry, not a tumour expression cohort"),
  "GSE144487" = list(ct = "Melanoma", note = "DNA methylation array-based immune microenvironment study, not expression"),
  "GSE320217" = list(ct = "Sarcoma", note = "cartilage-tumour methylome profiling, not expression; only 5 events anyway"),
  "GSE237299" = list(ct = "Lymphoma", note = "DNA methylation classifier for haematolymphoid neoplasms (GPL13534/GPL23976), not expression"),
  "GSE62371" = list(ct = "Melanoma", note = "miRNA two-colour array with a common reference made of all samples -> log-ratio data, not a usable expression cohort"),
  "GSE39040" = list(ct = "Sarcoma", note = "miRNA profiling by array of osteosarcoma (non-coding RNA)"),
  "GSE244807" = list(ct = "Cholangiocarcinoma", note = "iCCA RNA-seq; 137 biopsies + 109 surgical specimens; os/death in matrix"),
  "GSE108712" = list(ct = "Head and Neck Cancer", drop = TRUE, drop_reason = "NanoString nCounter 4-gene signature panel (MMP1/COL4A1/P4HA2/THBS2): effective_n=522 looks large but the matrix has only 4 genes - cannot support a transcriptome-wide CanPAS cohort",
    note = "HNSCC surgical-margin cohort (294 tumour + 277 normal, PAIRED -> effective_n double-counts patients), BUT GPL24460 is a NanoString nCounter **4-gene** signature panel (MMP1, COL4A1, P4HA2, THBS2) - not a transcriptome-wide cohort; exclude for CanPAS"),
  "GSE17679" = list(ct = "Sarcoma", note = "SuperSeries aggregating SubSeries GSE17618 + GSE17674 (each only ~43 usable pairs). The aggregate (117 samples / 86 usable) is the only Ewing-sarcoma cohort above the gate; keep as ONE cohort, do not also take the subseries"),
  "GSE17618" = list(ct = "Sarcoma", note = "SubSeries of GSE17679 (Ewing sarcoma) - same patients as the aggregate; below the gate on its own"),
  "GSE17674" = list(ct = "Sarcoma", note = "SubSeries of GSE17679 (Ewing sarcoma) - same patients as the aggregate; below the gate on its own"),
  "GSE218006" = list(ct = "Melanoma", note = "SuperSeries of GSE218004/GSE218005 (ACT-treated melanoma, TIL/T-cell samples); the OS columns are patient-level, so effective_n counts samples per patient - check the tumour-vs-TIL composition before building"),
  "GSE22153" = list(ct = "Melanoma", note = "Test set (57) of the GSE22155 SuperSeries (test+validation, same patients); use only one of GSE22153/GSE22155"),
  "GSE22155" = list(ct = "Melanoma", note = "SuperSeries = GSE22153 (test) + GSE22154 (validation); overlaps GSE22153"),
  "GSE69051" = list(ct = "Lymphoma", note = "DASL v4 arm (156 new patients + 71 replicates from the v3 arm GSE69049/GSE32918 + 12 replicates of GSE32918) - PATIENT OVERLAP with GSE32918; GSE69053 is the SuperSeries of GSE69049+GSE69051"),
  "GSE32918" = list(ct = "Lymphoma", note = "DASL lymphoma cohort (249 usable / 137 events); 12 of its samples are replicated inside GSE69051 - do not treat the two as independent validation sets"),
  "GSE23501" = list(ct = "Lymphoma", note = "Expression arm (GPL570) of the 69-sample DLBCL study whose methylation arm is GSE23967 - same patients, different assay; only one of the two can be used"),
  "GSE23967" = list(ct = "Lymphoma", note = "HELP promoter methylation arm of the same 69 DLBCL patients as GSE23501 - not an expression cohort, and same patients"),
  "GSE12945" = list(ct = "zz_other", note = "OUT OF SCOPE bonus: 62 colorectal cancers (Charite CRC cohort) with OS in GEO; it only matched the lymphoma query because the summary mentions 'leukemia and lymphoma' (Wiskott-Aldrich syndrome). Not already in the catalog (catalog Colorectal Cancer rows are different GSEs) - flag for the Colorectal bucket, not for the 15 gap types"),
  "GSE22138" = list(ct = "Melanoma+Uveal Melanoma", note = "Uveal melanoma primary tumours (n=63, GPL570): GEO characteristics carry tissue/age/gender/eye/tumour location/size but NO survival or follow-up column, and the only supplementary file is GSE22138_RAW.tar (CEL files) - metastasis-free survival is not obtainable from GEO"),
  "GSE22541" = list(ct = "Kidney Cancer", note = "ccRCC: GEO has 'Disease free interval (DFI) in months' but NO status/event column -> cannot form a usable pair; survival table would have to come from the publication"),
  "GSE19949" = list(ct = "Kidney Cancer", note = "RCC subgroups (n=133): GEO characteristics are grade/stage/ICD-O diagnosis only, no survival"),
  "GSE11151" = list(ct = "Kidney Cancer+Sarcoma", note = "renal tumours + normal kidneys: only 'renal tumor' source annotation, no clinical/survival fields"),
  "GSE167093" = list(ct = "Kidney Cancer", note = "Largest kidney series found (n=656, GPL10558): GEO characteristics are age/grade/histo/Sex/Stage/tissue/type only - NO survival or follow-up column; supplements are normalised expression matrices, no clinical table"),
  "GSE163722" = list(ct = "Mesothelioma", note = "Malignant pleural mesothelioma cohort (n=131, GPL11532) used for a RERG female-survival-advantage paper: GEO characteristics are compartment/organ/Sex only and the sole supplementary file is GSE163722_RAW.tar -> the survival table is in the publication, not in GEO"),
  "GSE162739" = list(ct = "Melanoma+Mesothelioma+Uveal Melanoma", note = "NOT a clinical cohort: BAP1-knockout cell-line RNA-seq/ChIP study; it matched the queries only because the summary lists melanoma/mesothelioma/uveal melanoma cell lines"),
  "GSE2109" = list(ct = "zz_other", note = "Expression Project for Oncology (expO) pan-cancer resource (n=2158): matched almost every disease query; excluded as PANEL_RESOURCE with no survival"),
  "GSE3218" = list(ct = "Testicular Cancer", note = "101 adult male germ cell tumours + 5 normal testis: GEO characteristics are histology only, no survival/follow-up"),
  "GSE205209" = list(ct = "Endometrial Cancer", note = "uterine serous cancer, OS in matrix (NanoString PanCancer IO360 targeted panel)")
)
for (g in names(MANUAL)) {
  ii <- which(all_res$gse == g)
  if (!length(ii)) next
  all_res$cancer_type_verified[ii] <- MANUAL[[g]]$ct
  all_res$type_note[ii] <- MANUAL[[g]]$note
  if (isTRUE(MANUAL[[g]]$drop)) {
    all_res$manual_drop[ii] <- TRUE
    for (j in ii) set_flag(j, "MANUAL_EXCLUDE", MANUAL[[g]]$drop_reason)
  }
}
# 人工排除生效（必须在 MANUAL 循环之后）
all_res$passes_gate <- all_res$passes_gate & !all_res$manual_drop
all_res$few_events  <- all_res$passes_gate & all_res$events < 10
all_res$sparse_event_coding <- all_res$passes_gate & all_res$events == all_res$effective_n & all_res$effective_n > 0
all_res <- all_res[order(all_res$cancer_type_verified, -all_res$effective_n, all_res$gse, all_res$platform), ]

out_cols <- c("cancer_type_verified", "cancer_type", "gse", "platform", "n_samples", "effective_n", "events", "token",
              "tokens_all", "time_key", "event_key", "source_of_survival_info", "sources_all",
              "survival_location", "supp_clinical_hint", "supp_files", "supp_files_n", "annotation_status", "annotation_local", "annotation_ok",
              "in_catalog", "catalog_type", "reference_design", "ref_like_samples",
              "qc_artefact", "panel_resource", "cell_line_mentioned", "human", "taxon", "series_gdstype", "platform_class", "expression_assay",
              "manual_drop",
              "platform_title", "platform_technology", "platform_organism",
              "is_superseries", "same_study_as", "study_group",
              "exclude_flag", "exclude_reason", "keep_platform",
              "passes_gate", "event_column_valid", "few_events", "sparse_event_coding",
              "any_time_column", "any_event_column", "time_without_event", "gene_panel_warning",
              "is_subseries_of_listed", "superseries_of", "subseries_of", "superseries_covered",
              "in_surv_query", "pdat", "source_names", "header_file",
              "series_relations", "type_note", "title", "overall_design", "summary")
write.csv(all_res[, out_cols], file.path(OUT, "geo_expansion_candidates.csv"), row.names = FALSE)
saveRDS(all_res, file.path(OUT, "geo_expansion_candidates.rds"))
cat("wrote pipeline/out/geo_expansion_candidates.csv  rows=", nrow(all_res),
    " pass>", MIN_EFF, "=", sum(all_res$passes_gate), "\n", sep = "")
cat("parse errors:", length(errs), "\n")
if (length(errs)) writeLines(errs, file.path(OUT, "geo_expansion_parse_errors.txt"))

# ---------------------------------------------------------------- markdown 报告
qlog <- NULL
qf <- file.path(OUT, "geo_expansion_queries.csv")
if (file.exists(qf)) qlog <- read.csv(qf, check.names = FALSE, stringsAsFactors = FALSE)
nfiles <- length(files)
nseries <- length(unique(all_res$gse))
md <- c(
 "# GEO expansion candidates — filling the 15 gap cancer types",
 "",
 sprintf("Generated: %s by `pipeline/R/37_geo_expansion_screen.R`", format(Sys.time())),
 "",
 "## 0. What this is / how it was produced",
 "",
 "`data/dataset_info.csv` (155 rows) contains 15 cancer types that have TCGA cohorts but **zero GEO cohorts**:",
 paste0("`", paste(GAP_TYPES, collapse = "`, `"), "`."),
 "",
 "This file lists GEO candidate cohorts for those 15 types. It is a **search + triage** artifact only:",
 "no cohort was built, no expression matrix and no supplementary file was downloaded.",
 "",
 "### Pipeline",
 "",
 "1. **GEO search** — NCBI E-utilities on `db=gds` (GEO DataSets), `gse[Entry Type]` filter,",
 "   sample-count filter `51:100000[Number of Samples]`. For every cancer type two queries were run:",
 "   a *disease-only* query and the same disease query ANDed with a survival/prognosis vocabulary",
 "   (`survival OR prognosis OR prognostic OR outcome OR \"follow-up\" OR recurrence OR relapse OR",
 "   mortality OR \"clinical characteristics\" OR \"disease-free\" OR \"progression-free\"`).",
 "   The two result sets were unioned per type (a series can belong to more than one type).",
 "   Exact query strings: `pipeline/out/geo_expansion_queries.csv` (also reproduced in §5 below).",
 sprintf("   Universe: **%d GSE series** (>50 samples) over the 15 types.", length(unique(smeta$gse))),
 "2. **Header-only retrieval** — for each series in the *survival* result set, plus every series in the",
 "   disease-only set whose title/summary matches a survival-endpoint regex, the **series-matrix header**",
 "   was fetched with an HTTP `Range` request on the `.gz` and a *partial* gzip inflate, stopping at",
 "   `!series_matrix_table_begin`. The expression matrix body was never transferred.",
 sprintf("   Headers retrieved: **%d series** (%d series-platform rows after expansion over `gpl_list`).",
         nseries, nrow(all_res)),
 "   Files kept under `pipeline/out/geo_expansion_headers/` (gzipped after the run).",
 "3. **Survival extraction** — the header parser of `pipeline/R/36_screen_geo_candidates.R` was reused",
 "   verbatim (CRLF-safe `trimws()` + quote stripping; top-level `;`/`,` splitting that respects",
 "   parentheses; first top-level `=` else `:` as key/value separator; endpoint token from the key;",
 "   role (time vs event) from the value; `parse_combined()` for free-text `!Sample_description`",
 "   sentences such as `\"Age: 72, ... RFS months: 109\"`; `value_text[TOKEN]` rows come from that path).",
 "   Both `!Sample_characteristics_ch*` and `!Sample_description` are parsed.",
 "4. **Gate** — a candidate *passes* iff",
 "   `effective_n > 50` AND `events >= 1` AND the series is not fatally excluded.",
 "   `effective_n` counts samples in the same endpoint token whose **time parses to a number AND status",
 "   decodes to 0/1**; `events` is the number of `status == 1` among them.",
 "   Fatal exclusions (each recorded in `exclude_flag` / `exclude_reason`):",
 "   already in `data/dataset_info.csv`; two-colour / common-reference design; QC or MAQC-style artefact;",
 "   non-human organism; platform is not an expression assay (methylation / SNP / aCGH / miRNA / qPCR);",
 "   SuperSeries whose SubSeries is itself in the table; a non-primary platform row of the same series;",
 "   same-study duplicate title; and a status column that is not a plain 0/1 (`events == 0` after decoding,",
 "   e.g. the A/C/O code in GSE131960/GSE131961).",
 "   Two manual exclusions are registered explicitly in the script (`MANUAL`): GSE108712 (NanoString",
 "   **4-gene** panel) and GSE17118 (mixed-histology matrix, no single type reaches 50).",
 "   The `>50` gate is applied to **all** samples of the matrix; the tumour/normal composition is reported",
 "   per series in `source_names` so the auditor can re-derive a per-type n.",
 "",
 "### Column dictionary (`geo_expansion_candidates.csv`)",
 "",
 "| column | meaning |",
 "|---|---|",
 "| `cancer_type_verified` | cancer type(s) used for grouping: the query-based `cancer_type` with the manual corrections applied |",
 "| `cancer_type` | raw query-based type(s) this series matched (`+` = matched several) |",
 "| `type_note` | per-series manual note (mixed histology, overlap with another deposit, why survival is missing...) |",
 "| `gse` / `platform` | GEO series / GPL (from E-utilities esummary `gpl`) |",
 "| `n_samples` | samples in the series matrix header |",
 "| `effective_n` | samples with time **and** status parsed, same endpoint token (the gate metric) |",
 "| `events` | `status == 1` count among those samples |",
 "| `token` / `tokens_all` | best endpoint token; all tokens as `TOKEN(effective_n/events)` |",
 "| `time_key` / `event_key` / `source_of_survival_info` | exact characteristic keys (or `value_text[...]`) used |",
 "| `sources_all` | every survival-bearing key seen, with `[characteristics]` or `[description]` origin |",
 "| `survival_location` | `series_matrix_header` (usable pair inside the matrix) / `supplement_possible_clinical` (a supplementary file name hints at clinical data) / `supplement_matrix_only(...)` (only RAW.tar / matrix files) / `absent` |",
 "| `supp_files` | `!Series_supplementary_file` names (**recorded only, never downloaded**) |",
 "| `annotation_status` | `local_rds` (a `data/gpl/<GPL>.rds` exists) or AnnoProbe::idmap outcome |",
 "| `in_catalog` / `catalog_type` | already represented in `data/dataset_info.csv` |",
 "| `reference_design` / `ref_like_samples` | two-colour Universal-Reference design flag |",
 "| `qc_artefact` | MAQC-style / QC / spike-in / milestone series flag |",
 "| `exclude_flag` / `exclude_reason` | `DUPLICATE_CATALOG`, `DUPLICATE_PLATFORM`, `REFERENCE_DESIGN`, `QC_ARTEFACT`, `NON_HUMAN`, `NON_EXPRESSION_ASSAY`, `SUPERSERIES`, `SUBSERIES_OF_LISTED_SUPER`, `SAME_STUDY_DUPLICATE`, `EVENT_COLUMN_UNPARSEABLE`, `TIME_ONLY_NO_EVENT`, `MANUAL_EXCLUDE` |",
 "| `platform_class` | `expression_array` / `rna_seq` / `ncounter` (usable) vs `methylation` / `snp_cgh` / `mirna` / `qpcr` / `other` / `unknown` (fatal) |",
 "| `any_time_column` / `any_event_column` | whether the parser found **any** survival time / any decodable status column at all |",
 "| `time_without_event` | time present but no status column anywhere (e.g. GSE22541 DFI) |",
 "| `gene_panel_warning` | targeted panel (nCounter/qPCR): gene coverage limited |",
 "| `source_names` | `!Sample_source_name_ch1` breakdown (tumour vs normal vs cell line) |",
 "| `keep_platform` | FALSE when another platform row of the same series is the primary one |",
 "| `passes_gate` | `effective_n > 50` and no fatal exclusion |",
 "| `few_events` | passes the >50 gate but `events < 10` (power warning) |",
 "| `sparse_event_coding` | `events == effective_n` — event column likely sparse/one-sided, verify by hand |",
 "| `in_surv_query` | series was returned by the *disease AND survival* query |",
 "| `series_relations` | `!Series_relation` lines (SuperSeries/SubSeries/BioProject — duplicate-deposit risk) |",
 "")

# ---- §1 per-type summary
md <- c(md, "## 1. Summary by cancer type", "",
        "| cancer type | series screened (>50 samples) | candidates in table | pass >50 gate | buildable (pass + annotation + >=10 events) | blocked by | recommendation |",
        "|---|---|---|---|---|---|---|")
blocked_txt <- function(gse_ok, all_gse) {
  # 只统计"未通过"的原因（针对该癌种）
  sel <- all_res[all_res$gse %in% all_gse & !all_res$gse %in% gse_ok, , drop = FALSE]
  if (!nrow(sel)) return("—")
  tags <- character(0)
  add <- function(mask, label) if (any(mask)) {
    gs <- unique(sel$gse[mask]); shown <- head(gs, 5)
    tags <<- c(tags, sprintf("%s: %s%s", label, paste(shown, collapse = ","),
                             if (length(gs) > length(shown)) sprintf(" (+%d more)", length(gs) - length(shown)) else ""))
  }
  add(!sel$human, "non-human")
  add(!sel$expression_assay_ok, "non-expression assay")
  add(sel$panel_resource, "pan-cancer cell-line/panel resource")
  add(sel$in_catalog, "already in catalog")
  add(sel$reference_design, "reference design")
  add(sel$qc_artefact, "QC artefact")
  add(sel$is_superseries, "SuperSeries")
  add(nzchar(sel$same_study_as) & sel$same_study_as != sel$gse, "same-study duplicate")
  add(sel$effective_n > MIN_EFF & !sel$event_column_valid, "status column not a usable 0/1")
  add(sel$passes_gate %in% c(TRUE, "TRUE") & !sel$annotation_ok, "passes the gate but the probe->gene map is not obtainable")
  low <- sel[sel$effective_n > 0 & sel$effective_n <= MIN_EFF, ]
  if (nrow(low)) tags <- c(tags, sprintf("%d series with 1-50 usable pairs", length(unique(low$gse))))
  zero <- sel[sel$effective_n == 0 & !sel$panel_resource, ]
  if (nrow(zero)) tags <- c(tags, sprintf("%d series with no time+status pair in the matrix header (survival may be in a supplement, or none)", length(unique(zero$gse))))
  if (!length(tags)) tags <- "—"
  paste(tags, collapse = "; ")
}
for (ct in GAP_TYPES) {
  inct <- all_res[grepl(paste0("(^|\\+)", gsub("([()])", "\\\\\\1", ct), "(\\+|$)"), all_res$cancer_type_verified), , drop = FALSE]
  ok <- inct[inct$passes_gate, , drop = FALSE]
  ok10 <- ok[ok$events >= 10 & ok$annotation_ok, , drop = FALSE]
  allg <- unique(inct$gse)
  scr <- length(unique(smeta$gse[vapply(smeta$gap_types, function(z) !is.na(z) && grepl(paste0("(^|;)", gsub("([()])", "\\\\\\1", ct), "(\\(surv\\))?(;|$)"), z), logical(1))]))
  suppc <- inct[inct$effective_n == 0 & inct$supp_clinical_hint & inct$keep_platform &
                 inct$human & inct$expression_assay_ok & !inct$panel_resource, , drop = FALSE]
  rec <- if (nrow(ok10) >= 3) "**usable** — several cohorts >50 pairs with >=10 events; recommend auditing the top rows"
         else if (nrow(ok10) >= 1) "**partially usable** — 1-2 cohorts pass; GEO coverage thin, TCGA remains the backbone"
         else if (nrow(ok) >= 1) "**marginal** — passes the >50 gate but events < 10"
         else if (nrow(suppc) >= 1) sprintf("**supplement-only route** — no >50 cohort inside a matrix, but %s list a clinical-looking supplementary file (the only route to a GEO cohort for this type)",
                                            paste(sprintf("%s (n=%d)", head(unique(suppc$gse), 3), head(suppc$n_samples[match(unique(suppc$gse), suppc$gse)], 3)), collapse = ", "))
         else "**no GEO cohort found** — keep TCGA-only"
  md <- c(md, sprintf("| %s | %d | %d | %d | %d | %s | %s |", ct, scr, length(unique(inct$gse)),
                      length(unique(ok$gse)), length(unique(ok10$gse)),
                      blocked_txt(unique(ok$gse), allg), rec))
}

# ---- §1b top candidates overall
top <- all_res[all_res$passes_gate & all_res$annotation_ok, , drop = FALSE]
top <- top[order(-top$effective_n), , drop = FALSE]
md <- c(md, "", "### 1b. Top candidates overall (pass the >50 gate, annotation obtainable)", "",
        "Sorted by effective n. These are the rows to audit first.", "",
        "| rank | cancer type | gse | platform | platform title | n_samples | effective_n | events | token |",
        "|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(min(nrow(top), 15))) md <- c(md, sprintf("| %d | %s | %s | %s | %s | %d | %d | %d | %s |",
    i, top$cancer_type_verified[i], top$gse[i], top$platform[i],
    gsub("\\|", "/", substr(top$platform_title[i], 1, 44)), top$n_samples[i],
    top$effective_n[i], top$events[i], top$token[i]))
md <- c(md, "")

# ---- §1c supplement-only candidates
soc <- all_res[all_res$effective_n == 0 & all_res$supp_clinical_hint & all_res$human &
               all_res$expression_assay_ok & !all_res$panel_resource & all_res$keep_platform, , drop = FALSE]
soc <- soc[order(-soc$n_samples), , drop = FALSE]
md <- c(md, "", "### 1c. Supplement-only candidates (no usable pair in the matrix, but a supplementary file name hints at clinical data)", "",
        "These cannot be built without downloading the named file; they are the second-best tier after §1b.", "",
        "| cancer type | gse | platform | n_samples | supplementary file |", "|---|---|---|---|---|")
if (nrow(soc)) for (i in seq_len(min(nrow(soc), 20)))
  md <- c(md, sprintf("| %s | %s | %s | %d | %s |", soc$cancer_type_verified[i], soc$gse[i], soc$platform[i],
                      soc$n_samples[i], gsub("\\|", "/", substr(soc$supp_files[i], 1, 130)))) else
  md <- c(md, "| - | - | - | - | none |")
md <- c(md, "")

# ---- §2 per-type candidate detail
md <- c(md, "", "## 2. Candidates by cancer type", "")
for (ct in GAP_TYPES) {
  inct <- all_res[grepl(paste0("(^|\\+)", gsub("([()])", "\\\\\\1", ct), "(\\+|$)"), all_res$cancer_type_verified), , drop = FALSE]
  ok <- inct[inct$passes_gate, , drop = FALSE]
  suppc <- inct[inct$effective_n == 0 & inct$supp_clinical_hint & inct$keep_platform &
                inct$human & inct$expression_assay_ok & !inct$panel_resource, , drop = FALSE]
  md <- c(md, sprintf("### %s", ct), "",
          sprintf("- series screened (>50 samples, both queries unioned): **%d**", length(unique(inct$gse))),
          sprintf("- series-platform rows in table: **%d**", nrow(inct)),
          sprintf("- **pass >50 gate: %d** (of which >=10 events: %d)", length(unique(ok$gse)), length(unique(ok$gse[ok$events >= 10]))),
          sprintf("- blocked: %s", blocked_txt(unique(ok$gse), unique(inct$gse))),
          if (nrow(suppc))
            sprintf("- supplement-only candidates (clinical-looking supplementary file, no pair in the matrix): %s",
                    paste(sprintf("%s (n=%d)", suppc$gse, suppc$n_samples), collapse = ", ")) else NULL, "")
  if (nrow(ok)) {
    md <- c(md, "| gse | platform | platform title | n_samples | effective_n | events | token | survival source (key names) | location | annotation | warnings |",
            "|---|---|---|---|---|---|---|---|---|---|---|")
    for (i in seq_len(nrow(ok))) md <- c(md, sprintf("| %s | %s | %s | %d | %d | %d | %s | %s | %s | %s | %s |",
        ok$gse[i], ok$platform[i], gsub("\\|", "/", substr(ok$platform_title[i], 1, 46)),
        ok$n_samples[i], ok$effective_n[i], ok$events[i], ok$token[i],
        gsub("\\|", "/", ok$source_of_survival_info[i]), ok$survival_location[i], ok$annotation_status[i],
        paste(c(if (ok$few_events[i]) "few events(<10)", if (ok$sparse_event_coding[i]) "sparse event coding",
                if (nzchar(ok$type_note[i])) "see type_note in CSV",
                if (ok$in_surv_query[i] == "FALSE") "found by disease-only query"), collapse = "; ")))
    md <- c(md, "", sprintf("Sample composition (`!Sample_source_name_ch1`): %s", 
              paste(sprintf("%s -> %s", ok$gse, substr(ok$source_names, 1, 150)), collapse = " ; ")))
  } else md <- c(md, "_No candidate passes the >50 gate._")
  # 该癌种"最值得追"的受阻 series（按样本数排序，只列头部）
  nm0 <- inct[!inct$passes_gate & inct$keep_platform, , drop = FALSE]
  nm0 <- nm0[order(-nm0$n_samples), , drop = FALSE]
  if (nrow(nm0)) {
    md <- c(md, "", "Largest blocked/excluded series for this type (by `n_samples`):", "",
            "| gse | platform | n_samples | effective_n | why blocked |", "|---|---|---|---|---|")
    for (i in seq_len(min(nrow(nm0), 6))) md <- c(md, sprintf("| %s | %s | %d | %d | %s |",
        nm0$gse[i], nm0$platform[i], nm0$n_samples[i], nm0$effective_n[i],
        gsub("\\|", "/", if (nzchar(nm0$exclude_flag[i])) nm0$exclude_flag[i]
             else if (nm0$effective_n[i] == 0) "no time+status pair in the matrix header (survival may be in a supplement)"
             else "1-50 usable pairs (below the >50 gate)")))
    md <- c(md, "")
  }
  # blocked / near-miss list
  nm <- inct[!inct$passes_gate, , drop = FALSE]
  nm <- nm[order(-nm$effective_n), , drop = FALSE]
  if (nrow(nm)) {
    md <- c(md, "", sprintf("<details><summary>Near misses / excluded for %s (%d rows)</summary>", ct, nrow(nm)), "",
            "| gse | platform | n_samples | effective_n | events | flag | reason |", "|---|---|---|---|---|---|---|")
    for (i in seq_len(min(nrow(nm), 25)))
      md <- c(md, sprintf("| %s | %s | %d | %d | %d | %s | %s |", nm$gse[i], nm$platform[i], nm$n_samples[i],
              nm$effective_n[i], nm$events[i], nm$exclude_flag[i],
              substr(gsub("\\|", "/", ifelse(nzchar(nm$exclude_reason[i]), nm$exclude_reason[i],
                   ifelse(nm$effective_n[i] == 0, "no usable time+status pair in the series-matrix header", "1-50 usable pairs (below gate)"))), 1, 170)))
    md <- c(md, "", "</details>", "")
  } else md <- c(md, "")
}

# ---- §3 explicit exclusions
md <- c(md, "## 3. Explicitly excluded (and why)", "",
        "### 3.0 Non-human or non-expression series (silently fatal for this pipeline)", "")
nh <- all_res[!all_res$human | !all_res$expression_assay_ok, , drop = FALSE]
if (nrow(nh)) {
  md <- c(md, "| gse | platform | platform title | n_samples | effective_n | why excluded |", "|---|---|---|---|---|---|")
  for (i in seq_len(min(nrow(nh), 40))) md <- c(md, sprintf("| %s | %s | %s | %d | %d | %s |",
      nh$gse[i], nh$platform[i], gsub("\\|", "/", substr(nh$platform_title[i], 1, 44)), nh$n_samples[i], nh$effective_n[i],
      gsub("\\|", "/", paste(c(if (!nh$human[i]) paste0("organism=", nh$taxon[i]),
                               if (!nh$expression_assay_ok[i]) nh$expression_assay[i]), collapse = "; "))))
  if (nrow(nh) > 40) md <- c(md, "", sprintf("_... and %d more rows (see the CSV `exclude_flag` column)._", nrow(nh) - 40))
} else md <- c(md, "- none")
md <- c(md, "", "### 3.0c Pan-cancer cell-line / public panel resources (matched by disease keywords, not single-cancer cohorts)", "")
pr <- unique(all_res$gse[all_res$panel_resource])
if (length(pr)) {
  md <- c(md, "| gse | platform | n_samples | title |", "|---|---|---|---|")
  for (g in sort(pr)) { s <- all_res[all_res$gse == g, ][1, ]
    md <- c(md, sprintf("| %s | %s | %d | %s |", g, s$platform, s$n_samples, gsub("\\|", "/", substr(s$title, 1, 80)))) }
} else md <- c(md, "- none")
md <- c(md, "", "### 3.0b Usable sample count but a status column that is not a plain 0/1 (decodes to 0 events)", "")
ez <- all_res[all_res$effective_n > MIN_EFF & !all_res$event_column_valid, , drop = FALSE]
if (nrow(ez)) for (i in seq_len(nrow(ez)))
  md <- c(md, sprintf("- **%s** (%s, n=%d, effective_n=%d): %s", ez$gse[i], ez$platform[i], ez$n_samples[i],
                      ez$effective_n[i], gsub("\\|", "/", ez$exclude_reason[ez$gse[i] == ez$gse][1]))) else md <- c(md, "- none")
md <- c(md, "", "### 3.1 Reference-array designs (two-colour, all/most samples Universal Reference)", "")
re <- unique(all_res$gse[all_res$reference_design])
if (length(re)) for (g in re) { s <- all_res[all_res$gse == g, ][1, ]
  md <- c(md, sprintf("- **%s** (%s, n=%d): %d/%d samples reference-like, series text two-colour/common-reference=%s — expression values are log-ratios, unusable as a survival cohort.",
                      g, s$platform, s$n_samples, s$ref_like_samples, s$n_samples, s$series_reference_text)) } else md <- c(md, "- none")
md <- c(md, "", "### 3.2 QC / milestone / technical artefacts (MAQC-style dummy outcomes)", "")
qc <- unique(all_res$gse[all_res$qc_artefact])
if (length(qc)) for (g in qc) { s <- all_res[all_res$gse == g, ][1, ]
  md <- c(md, sprintf("- **%s** (%s, n=%d): %s", g, s$platform, s$n_samples, substr(s$title, 1, 110))) } else md <- c(md, "- none")
md <- c(md, "", "### 3.3 Duplicate deposits of series already in the catalog", "")
du <- unique(all_res$gse[all_res$in_catalog])
if (length(du)) for (g in du) { s <- all_res[all_res$gse == g, ][1, ]
  md <- c(md, sprintf("- **%s** (%s): catalogued as `%s`", g, s$platform, s$catalog_type)) } else md <- c(md, "- none")
md <- c(md, "", "### 3.4 Same series on several platforms (only one platform row is primary)", "")
dp <- all_res[!all_res$keep_platform, , drop = FALSE]
if (nrow(dp)) for (i in seq_len(min(nrow(dp), 40)))
  md <- c(md, sprintf("- %s %s — %s", dp$gse[i], dp$platform[i], dp$exclude_reason[i]))
if (nrow(dp) > 40) md <- c(md, sprintf("- _... and %d more duplicate-platform rows._", nrow(dp) - 40))
if (!nrow(dp)) md <- c(md, "- none")
md <- c(md, "", "### 3.5 Survival only in a supplementary file (matrix header has no usable pair)", "")
so <- all_res[all_res$effective_n == 0 & all_res$supp_clinical_hint &
              !all_res$in_catalog & !all_res$qc_artefact & all_res$human & all_res$expression_assay_ok, , drop = FALSE]
so <- so[order(so$cancer_type_verified, -so$n_samples), , drop = FALSE]
if (nrow(so)) {
  md <- c(md, "These series are plausible cohorts but the survival columns are **not** in the series matrix header;",
          "they would need the named supplement (not downloaded here).", "",
          "| cancer type | gse | platform | n_samples | supplementary file(s) |", "|---|---|---|---|---|")
  for (i in seq_len(min(nrow(so), 80)))
    md <- c(md, sprintf("| %s | %s | %s | %d | %s |", so$cancer_type_verified[i], so$gse[i], so$platform[i],
                        so$n_samples[i], gsub("\\|", "/", substr(so$supp_files[i], 1, 120))))
} else md <- c(md, "- none")

md <- c(md, "", "### 3.6 SuperSeries whose samples are the union of their SubSeries", "")
ss <- unique(all_res$gse[all_res$is_superseries])
if (length(ss)) for (g in ss) { s <- all_res[all_res$gse == g, ][1, ]
  md <- c(md, sprintf("- **%s** (%s, n=%d, effective_n=%d) — %s", g, s$platform, s$n_samples, s$effective_n,
                      gsub("\\|", "/", substr(s$summary, 1, 120)))) } else md <- c(md, "- none")

md <- c(md, "", "### 3.7 Same-study deposits (identical normalised title; overlapping patients likely)", "")
sg <- unique(all_res$gse[nzchar(all_res$same_study_as) & all_res$same_study_as != all_res$gse])
if (length(sg)) for (g in sort(sg)) { s <- all_res[all_res$gse == g, ][1, ]
  md <- c(md, sprintf("- %s (%s, effective_n=%d) — same study as %s: %s", g, s$platform, s$effective_n,
                      s$same_study_as, substr(s$title, 1, 100))) } else md <- c(md, "- none")

# ---- §4 could not verify
md <- c(md, "", "## 4. Could not verify / limitations", "",
        "- **Event semantics of a few columns** — where `time_key == event_key` or `value_text[TOKEN]` was used,",
        "  the pair came from a free-text `!Sample_description` sentence; the exact key names are in the CSV.",
        "  Rows with `sparse_event_coding = TRUE` (`events == effective_n`) may be one-sided encodings and must be",
        "  checked against the raw header before building a cohort.",
        "- **Unit normalisation (months / days / years) was NOT applied.** `37_geo_expansion_screen.R` reports counts only;",
        "  the unit conversion must be written into the `03_surv_table.R` specs at build time.",
        "- **Supplement contents were not inspected** (no download). `supplement_possible_clinical` only means that a",
        "  supplementary **file name** hints at clinical data (clinical/patient/follow-up/annotation/...); whether an event",
        "  column actually exists there is unverified. Files that are clearly matrices (`RAW.tar`, `*norm*`, `*fpkm*`,",
        "  `*counts*`) are not counted as clinical.",
        "- **Sample-level composition (tumour vs adjacent normal vs cell line) was not used to filter `effective_n`.**",
        "  It is reported per series in `source_names` and repeated under each type in §2. Series with paired",
        "  tumour/normal designs (e.g. GSE108712 = 294 tumour + 277 normal) double-count patients in `effective_n`,",
        "  the same caveat the existing catalog records for `GSE53625` and `GSE102238`.",
        "- **`n_samples` is the whole matrix, not the cancer-type-specific subset.** For mixed-histology matrices",
        "  (GSE17118) or multi-disease matrices, `effective_n` overstates the usable cohort for any single type;",
        "  those rows carry a `type_note`.",
        "- **No survival data exists in the matrix header for the majority of series.** Only series whose GEO text",
        "  mentions a survival/clinical endpoint were header-fetched; a cohort whose survival lives purely in a",
        "  supplement with no endpoint word in the title/summary cannot be ruled out (see limitations of the query).",
        "- **Patient-level overlap between these GEO cohorts and the TCGA cohorts could not be checked**;",
        "  `series_relations` records SuperSeries/SubSeries/BioProject links, the only duplicate signal available",
        "  without downloading clinical tables.",
        "",
        "### 4.1 Platform annotation status (`annotation_status`)",
        "",
        "| `annotation_status` | meaning |", "|---|---|",
        "| `local_rds` | a `data/gpl/<GPL>.rds` probe→gene map already exists in this project |",
        "| `annoprobe_pipe` | `AnnoProbe::idmap(<GPL>, type=\"pipe\")` returns a `probe_id,symbol` table (one download at build time) |",
        "| `gene_level_no_probe_map_needed` | RNA-seq / NanoString: matrix rows are already gene identifiers, no probe map required |",
        "| `geo_annot_symbol` | AnnoProbe has no map, but GEO's own platform table (`<GPL>.annot.gz` or `<GPL>_family.soft.gz`) has a gene column — read from the table header only (HTTP Range, no full download) |",
        "| `geo_annot_no_symbol` | GEO's platform table has no gene column (e.g. GPL22757 Lymphochip = `ID,SPOT_ID` only) → **not buildable** without an external clone map |",
        "| `idmap_failed:<msg>` | neither route worked |",
        "| `not_tested` | the platform never reached a passing candidate row |",
        "`annotation_ok` (CSV column) is TRUE for the first four.",
        "")
ann_tab <- if (nrow(ann)) table(ann$status) else integer(0)
if (length(ann_tab)) md <- c(md, paste0("- probe outcomes: ", paste(sprintf("`%s`=%d", names(ann_tab), as.integer(ann_tab)), collapse = ", ")), "")
pass_plat <- unique(all_res[all_res$passes_gate, c("platform", "platform_title", "annotation_status")])
if (nrow(pass_plat)) {
  md <- c(md, "Platforms used by the passing candidates:", "",
          "| gpl | platform title | annotation |", "|---|---|---|")
  for (i in seq_len(nrow(pass_plat))) md <- c(md, sprintf("| %s | %s | %s |", pass_plat$platform[i],
        gsub("\\|", "/", substr(pass_plat$platform_title[i], 1, 60)), pass_plat$annotation_status[i]))
  md <- c(md, "")
}
if (length(errs)) md <- c(md, "### Header files that failed to parse", "", paste0("- `", errs, "`"), "")

# ---- §5 queries + artifacts
md <- c(md, "## 5. Reproducibility — exact GEO queries", "",
        "E-utilities endpoints (2026-09; NCBI rate limit respected, 4-10 concurrent workers):",
        "",
        "```",
        "esearch.fcgi?db=gds&term=<query>&retmax=10000&retmode=json",
        "esummary.fcgi?db=gds&id=<uids>&retmode=json",
        "```",
        "")
if (!is.null(qlog)) for (ct in GAP_TYPES) {
  sub <- qlog[qlog$cancer_type == ct, , drop = FALSE]
  if (!nrow(sub)) next
  md <- c(md, sprintf("**%s**", ct), "", "```")
  for (i in seq_len(nrow(sub))) md <- c(md, sprintf("# %s (hits: %s)", sub$query_kind[i], sub$n_hits[i]), sub$query_string[i], "")
  md <- c(md, "```", "")
}
md <- c(md, "### Artifacts written by this step", "",
        "- `pipeline/out/geo_expansion_candidates.csv` — the candidate table (this report's machine-readable twin)",
        "- `pipeline/out/geo_expansion_candidates.md` — this report",
        "- `pipeline/out/geo_expansion_candidates.rds` — full R object (all series-platform rows, all flags)",
        "- `pipeline/out/geo_expansion_series_meta.csv` — E-utilities esummary metadata + gap-type membership + survival-query membership",
        "- `pipeline/out/geo_expansion_platforms.csv` — `!Platform_title` / technology / organism for every GPL in the table",
        "- `pipeline/out/geo_expansion_queries.csv` — the exact esearch query strings with hit counts",
        "- `pipeline/out/geo_expansion_annotation_check.csv` — AnnoProbe::idmap probe results per GPL",
        "- `pipeline/out/geo_expansion_headers/` — gzipped series-matrix **headers only** (the raw evidence for every `effective_n`)",
        "- `pipeline/R/37_geo_expansion_screen.R` — the screening script",
        "- `pipeline/R/38_geo_expansion_annotation_check.R` — the AnnoProbe probe script",
        "")
writeLines(md, file.path(OUT, "geo_expansion_candidates.md"))
cat("wrote pipeline/out/geo_expansion_candidates.md\n")
}

# ---------------------------------------------------------------------------
# run_113_sweep_gpl_mirror_divergence()  <-  verbatim pipeline/R/113_sweep_gpl_mirror_divergence.R
# ---------------------------------------------------------------------------
run_113_sweep_gpl_mirror_divergence <- function() {
# 113_sweep_gpl_mirror_divergence.R ------------------------------------------
# READ-ONLY sweep (report only; no writes anywhere).
#
# 目的：对**每一个**被 catalog 引用、且在镜像中存在的平台表（GPL* 与 *_PLAT），
# 把镜像行数与对应本地注释表（data/processed/gpl/<GPL>.rds）的
# **已注释探针行数（annotated row_names）** 对比，列出全部真实分歧。
#
# 口径警告（本次事故的根因）：本地 rds 的 gene_id 可能是 "A /// B /// C"，
# 99_extend_gpl_db.R 会按 " /// " 展开成多行 —— 展开后的行数是「探针-基因对」数，
# 不是「已注释探针」数。GPL16686 的 53,981 探针里只有 30,766 个带 gene_id，
# 展开后 35,637 行；把 35,637 当成镜像应有多少行就是本次的假分歧。
# 本脚本同时输出两个口径，判决用 annotated 口径。
#
# 用法: Rscript pipeline/R/113_sweep_gpl_mirror_divergence.R [ROOT]
# 输出: pipeline/out/gpl_divergence_sweep.csv
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(RMySQL); library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args)) args[1] else
  Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
ROOT <- path.expand(ROOT)

dl <- read.csv(file.path(ROOT, "data/dataset_info.csv"),
               stringsAsFactors = FALSE, check.names = FALSE)
if (!"GPL" %in% names(dl)) stop("dataset_info.csv 里没有 GPL 列")
cat_plats <- sort(unique(unlist(strsplit(as.character(dl$GPL[!is.na(dl$GPL) & dl$GPL != ""]),
                                          "[,;[:space:]]+"))))
cat_plats <- cat_plats[cat_plats != ""]
cat("catalog rows:", nrow(dl), "| distinct platforms referenced:", length(cat_plats), "\n")

con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = Sys.getenv("CPAS_DB_PASSWORD"), client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)
tabs <- dbListTables(con)
mirror_tabs <- sort(unique(c(cat_plats, tabs[grepl("^GPL", tabs) | grepl("_PLAT$", tabs)])))

gpl_dir <- file.path(ROOT, "data/processed/gpl")
res <- vector("list", 0)
for (tab in mirror_tabs) {
  f <- file.path(gpl_dir, paste0(tab, ".rds"))
  has_local <- file.exists(f)
  in_cat <- tab %in% cat_plats
  in_db <- tab %in% tabs
  mirror_n <- if (in_db) dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", tab))$n else NA_integer_
  mirror_d <- if (in_db) dbGetQuery(con, sprintf("SELECT COUNT(DISTINCT row_names) n FROM `%s`", tab))$n else NA_integer_

  n_probe <- annotated <- expanded <- NA_integer_
  dbonly <- localonly <- n_extra_pairs <- NA_integer_
  if (has_local) {
    g <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(g)) {
      rn <- rownames(g)
      n_probe <- length(rn)
      gid <- if ("gene_id" %in% names(g)) as.character(g$gene_id) else
        if ("ENTREZ_GENE_ID" %in% names(g)) as.character(g$ENTREZ_GENE_ID) else NA_character_
      keep <- !is.na(gid) & gid != "" & gid != "---"
      annotated <- sum(keep)
      rnk <- rn[keep]; gidk <- gid[keep]
      expanded <- sum(vapply(strsplit(gidk, " /// ", fixed = TRUE), length, 1L))
      if (in_db) {
        db2 <- dbGetQuery(con, sprintf("SELECT row_names, gene_id FROM `%s`", tab))
        dbset <- unique(db2$row_names)
        lset  <- unique(rnk)
        dbonly    <- length(setdiff(dbset, lset))
        localonly <- length(setdiff(lset, dbset))
        n_extra_pairs <- sum(grepl(" /// ", gidk, fixed = TRUE))  # probes whose split would add rows
      }
    }
  }
  # 判决：只用 annotated 口径
  delta <- if (!is.na(annotated) && !is.na(mirror_d)) annotated - mirror_d else NA_integer_
  verdict <- if (is.na(delta)) "no local map / not in mirror" else if (delta == 0 && dbonly == 0 && localonly == 0)
    "in sync (annotated)" else if (delta > 0 && dbonly == 0) "local newer: additive top-up possible" else
    if (delta < 0) "mirror has more rows than the local map (mirror newer)" else
    "divergent (needs manual review)"
  res[[length(res) + 1L]] <- data.frame(
    table = tab, in_catalog = in_cat, in_mirror = in_db, local_map = has_local,
    local_probes = n_probe, local_annotated_probes = annotated,
    local_expanded_pairs = expanded, mirror_rows = mirror_n, mirror_distinct_row_names = mirror_d,
    delta_annotated_minus_mirror = delta, mirror_only_row_names = dbonly,
    local_only_row_names = localonly, probes_with_multi_entrez = n_extra_pairs,
    additive_topup_possible = ifelse(!is.na(delta) & delta > 0 & dbonly == 0, "yes",
                                     ifelse(!is.na(delta) & delta == 0, "n/a (already complete)", "no")),
    verdict = verdict, stringsAsFactors = FALSE)
}
out <- do.call(rbind, res)
dir.create(file.path(ROOT, "pipeline/out"), showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(ROOT, "pipeline/out/gpl_divergence_sweep.csv"), row.names = FALSE)

cat("\n=== platforms referenced by the catalog ===\n")
print(out[out$in_catalog, c("table","local_probes","local_annotated_probes","local_expanded_pairs",
                            "mirror_rows","delta_annotated_minus_mirror",
                            "mirror_only_row_names","local_only_row_names","verdict")], row.names = FALSE)
cat("\n=== divergences among catalog platforms ===\n")
dv <- out[out$in_catalog & (is.na(out$delta_annotated_minus_mirror) | out$delta_annotated_minus_mirror != 0 |
                            out$mirror_only_row_names != 0 | out$local_only_row_names != 0), ]
print(dv[, c("table","local_annotated_probes","mirror_rows","delta_annotated_minus_mirror",
             "mirror_only_row_names","local_only_row_names","verdict")], row.names = FALSE)
cat("\n=== mirror tables that are NOT referenced by the catalog ===\n")
print(out[!out$in_catalog, c("table","mirror_rows","local_annotated_probes",
                             "delta_annotated_minus_mirror","verdict")], row.names = FALSE)
}
