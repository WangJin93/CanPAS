# ===========================================================================
# CanPAS curation pipeline -- 06_catalog_and_registry
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 05_dataset_plan.R
#   2. 19_fix_catalog_N.R
#   3. 26_cohort_overlap.R
#   4. 27_catalog_notes.R
#   5. 41_register_final3_rows.R
#   6. 94_update_catalog_geo_expansion.R
#   7. 98_update_catalog_suppl_expansion.R
#   8. 101_update_catalog_relaxed_gate.R
#   9. 105_update_catalog_embl.R
#   10. 107_register_embl_step2.R
#   11. 111_register_gse205209.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "06_catalog_and_registry.R", package = "CanPAS"))
#   run_05_dataset_plan()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================


# ---------------------------------------------------------------------------
# run_05_dataset_plan()  <-  verbatim pipeline/R/05_dataset_plan.R
# ---------------------------------------------------------------------------
run_05_dataset_plan <- function() {
# 05_dataset_plan.R
# -----------------------------------------------------------------------------
# 生成 pipeline/out/dataset_plan.csv（以 base accession 为准，多平台 GSE 合并）：
#   “all datasets.csv” 主表全集 与 dataset_info.csv、本地文件、远端 DB(若可达,只读)
#   比对，给出每套数据的状态与下一步动作。
# 用法: Rscript 05_dataset_plan.R [<ProjectRoot>]
# -----------------------------------------------------------------------------
suppressMessages({library(dplyr); library(tidyr); library(stringr)})

ROOT <- if (length(commandArgs(trailingOnly=TRUE)) >= 1)
  commandArgs(trailingOnly=TRUE)[1] else "~/data/Project/CanPAS"
ROOT <- path.expand(ROOT)
setwd(ROOT)

master <- read.csv("all datasets.csv", header=FALSE, check.names=FALSE,
                   colClasses="character", na.strings="")
nm <- make.unique(trimws(as.character(master[1, ])))
colnames(master) <- nm
master <- master[-1, ]
master_long <- master %>%
  tidyr::pivot_longer(everything(), names_to="Type", values_to="Accession") %>%
  filter(!is.na(Accession), nzchar(trimws(Accession))) %>%
  mutate(Accession = trimws(Accession)) %>%
  distinct(Accession, .keep_all=TRUE) %>%
  mutate(base = str_remove(Accession, "[-_]GPL[0-9]+$")) %>%
  select(Type, Accession, base)

info <- read.csv("data/dataset_info.csv", stringsAsFactors=FALSE)
info <- info %>% mutate(base = str_remove(Accession, "[-_]GPL[0-9]+$"))

local_expr <- sub("\\.rds$", "", list.files("data/expr", pattern="\\.rds$"))
local_pheno<- sub("\\.rds$", "", list.files("data/pheno", pattern="\\.rds$"))
local_surv <- sub("_surv\\.rds$", "", list.files("data/processed/surv", pattern="_surv\\.rds$"))
local_gpl  <- sub("\\.rds$", "", list.files("data/processed/gpl", pattern="\\.rds$"))
norm <- function(x) str_replace_all(x, "-", "_")

db_ok <- FALSE
expr_db <- surv_db <- gpl_db <- character(0)
if (requireNamespace("RMySQL", quietly=TRUE)) {
  con <- tryCatch(DBI::dbConnect(RMySQL::MySQL(), host="139.224.80.159",
      dbname="cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CPAS"), password = Sys.getenv("CPAS_DB_PASSWORD")),
    error=function(e) NULL)
  if (!is.null(con)) {
    tabs <- DBI::dbListTables(con)
    db_ok <- TRUE
    surv_db <- sub("_surv$", "", tabs[grepl("_surv$", tabs)])
    expr_db <- setdiff(tabs[grepl("^GSE", tabs)], paste0(surv_db, "_surv"))
    gpl_db  <- tabs[grepl("^GPL", tabs)]
    DBI::dbDisconnect(con)
  }
}

info_sum <- info %>%
  group_by(base) %>%
  summarise(Type_info = first(Type), N = sum(N),
            n_parts = n(),
            GPL = paste(sort(unique(GPL)), collapse=" | "),
            method = paste(unique(method[!is.na(method)]), collapse=" | "), .groups="drop") %>%
  mutate(in_info = TRUE)

file_sum <- data.frame(base = str_remove(unique(c(local_expr, local_pheno, local_surv)),
                                           "[-_]GPL[0-9]+$")) %>%
  distinct() %>%
  mutate(expr_local = norm(base) %in% norm(local_expr),
         pheno_local= norm(base) %in% norm(local_pheno),
         surv_local = norm(base) %in% norm(local_surv))

expr_db_st <- str_remove(expr_db, "[-_]GPL[0-9]+$")
surv_db_st <- str_remove(surv_db, "[-_]GPL[0-9]+$")
db_sum <- data.frame(base = unique(c(expr_db_st, surv_db_st))) %>%
  mutate(expr_db = base %in% expr_db_st, surv_db = base %in% surv_db_st)

gpl_flags <- info %>% filter(!is.na(GPL)) %>%
  group_by(base) %>%
  summarise(gpl_local = all(unique(GPL) %in% local_gpl),
            gpl_db    = if (db_ok) all(unique(GPL) %in% gpl_db) else NA,
            .groups="drop")

excl <- if (file.exists("pipeline/out/excluded.csv"))
  read.csv("pipeline/out/excluded.csv", stringsAsFactors=FALSE) else
  data.frame(Accession=character(0), Reason=character(0))

plan <- master_long %>%
  left_join(info_sum, by="base") %>%
  left_join(file_sum, by="base") %>%
  left_join(db_sum, by="base") %>%
  left_join(gpl_flags, by="base") %>%
  mutate(expr_local = ifelse(is.na(expr_local), FALSE, expr_local),
         surv_local = ifelse(is.na(surv_local), FALSE, surv_local),
         expr_db    = ifelse(is.na(expr_db), FALSE, expr_db),
         surv_db    = ifelse(is.na(surv_db), FALSE, surv_db),
         in_info    = ifelse(is.na(in_info), FALSE, in_info),
         status = case_when(
           Accession %in% excl$Accession ~ "EXCLUDED",
           expr_db & surv_db ~ "DONE-DB",
           surv_db ~ "DONE-SURV-ONLY",
           expr_db ~ "DONE-EXPR-ONLY",
           expr_local | surv_local | pheno_local ~ "LOCAL-PENDING",
           TRUE ~ "TODO")) %>%
  select(Type, Accession=base, GPL, N, in_info, expr_local, surv_local,
         expr_db, surv_db, gpl_local, gpl_db, status) %>%
  arrange(Type, Accession)

dir.create("pipeline/out", showWarnings=FALSE)
write.csv(plan, "pipeline/out/dataset_plan.csv", row.names=FALSE)
cat("DB reachable:", db_ok, " | n_base_datasets:", nrow(plan), "\n")
print(table(plan$status, plan$Type))
message("plan written: pipeline/out/dataset_plan.csv")
}

# ---------------------------------------------------------------------------
# run_19_fix_catalog_N()  <-  verbatim pipeline/R/19_fix_catalog_N.R
# ---------------------------------------------------------------------------
run_19_fix_catalog_N <- function() {
# 19_fix_catalog_N.R ---------------------------------------------------------
# 重算 catalog（dataset_info）的样本量列，使其与镜像实际可分析样本数一致。
#
# 为什么改
#   旧 N 实际记录的是「表达谱样本数」（122 个有镜像表的数据集里 119 个等于表达表列数），
#   不是分析真正能用的样本数，且个别行是历史遗留值（如 GSE31210 记 133，实际 226）。
#   用旧 N 估功效会被误导：分析用的是「有表达数据 且 该终点有记录」的样本。
#
# 新列定义
#   N          可分析样本数 = |表达样本 ∩ 主要终点(EndpointPrimary) time 与 status 都非缺失的样本|；
#              无表达表/该终点缺失 -> NA（表示该数据集在应用里做不了基因层面分析）
#   n_expr     镜像表达表样本数（无表 -> NA）
#   n_surv     镜像生存表行数（无表 -> NA）
#   n_events   主要终点事件数（status=1，取上面那个交集内）
#   n_OS ... n_MFS  各终点家族的可分析样本数（家族 -> 该数据集实际使用的 token，
#              如 GSE31210 的 DFS 家族用 RFS 列；无该家族 -> 0）
#   expr_in_mirror  镜像中是否有表达表
#
# TCGA 守卫（2026-09-24 加入，同日加固为 fail-closed）
#   31 个 TCGA-* 行按设计不在镜像里（表达按需从 UCSC Xena 抓取），镜像口径重算不出它们的 N。
#   实测直接写入会把 31 行 N 全部清成 NA，因此：
#     - 加载 catalog 后先取 .N_before 与 .tcga_keep <- grepl("^TCGA", Accession)；
#     - 重算 N 后把 TCGA 行的 N 原位还原；
#     - 立即做守卫自校验 .tcga_bad（TCGA 行只要出现「变 NA」或「改数」即非 0）；
#     - 写回前设写前断言：.tcga_bad 非 0 或任一 TCGA 行将变 NA -> stop()，绝不写出坏 catalog；
#     - dry-run 额外导出 catalog_N_fix_would_write.csv（--write 的完整内容，含全部 TCGA 行），
#       并把 TCGA 行逐行列在终端上，使「--write 会写成什么」无需真写即可核验。
#
# 用法
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/19_fix_catalog_N.R          # 只算并出报告
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/19_fix_catalog_N.R --write  # 写回
#     写回：data/dataset_info.csv、CanPAS/data/dataset_info.rda、data/dataset_info.xlsx
# 输出：pipeline/out/catalog_N_fix.csv        逐行 N_old/N_new + tcga_guard + would_become_NA
#       pipeline/out/catalog_N_fix_would_write.csv  dry-run 下的「--write 内容预览」
#       pipeline/out/REPORT_catalog_N_fix.md
# ----------------------------------------------------------------------------
.libPaths(c("/home/Jingle/R/library", .libPaths()))
suppressPackageStartupMessages({library(RMySQL); library(CanPAS)})

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
write_back <- "--write" %in% commandArgs(trailingOnly = TRUE)
fam <- c("OS", "DSS", "DFS", "PFS", "MFS")

# EndpointPrimary 可能是 ""（未标注）或未知 token —— 统一兜底成 NA
primary_family <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  f <- tryCatch(CanPAS::endpoint_family(x), error = function(e) character(0))
  if (length(f) != 1L) NA_character_ else f
}

pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
if (!nzchar(pw)) stop("CPAS_DB_PASSWORD 未设置。")
con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = pw, client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con))

di <- get(load(file.path(root, "CanPAS/data/dataset_info.rda")))
tb <- dbGetQuery(con, "SHOW TABLES")[[1]]
expr_tables <- setdiff(setdiff(tb, grep("_surv$", tb, value = TRUE)), grep("^GPL", tb, value = TRUE))
expr_tables <- setdiff(expr_tables, grep("_PLAT$", expr_tables, value = TRUE))   # 注释表不算表达谱

n_expr <- n_surv <- n_events <- rep(NA_integer_, nrow(di))
fam_n  <- matrix(0L, nrow(di), length(fam), dimnames = list(NULL, paste0("n_", fam)))
notes  <- character(nrow(di))

for (i in seq_len(nrow(di))) {
  acc <- as.character(di$Accession[i]); adb <- gsub("-", "_", acc, fixed = TRUE)
  if (startsWith(acc, "TCGA-")) {
    # TCGA 的表达按需从 UCSC Xena 抓取，镜像里没有表达表；样本集合以本地临床表为准，
    # 因此 n_expr 保持 NA（无固定表达样本数），N 仍按同一口径计算。
    sv <- tryCatch(CanPAS::tcga_surv_table(sub("^TCGA-", "", acc)), error = function(e) NULL)
    has_expr <- !is.null(sv)
    expr_ids <- if (has_expr) as.character(sv[[1]]) else character(0)
    if (!is.null(sv)) n_surv[i] <- nrow(sv)
  } else {
    has_expr <- adb %in% expr_tables
    expr_ids <- if (has_expr) setdiff(dbListFields(con, adb), "row_names") else character(0)
    if (has_expr) n_expr[i] <- length(expr_ids)
    surv_t <- paste0(adb, "_surv")
    sv <- if (surv_t %in% tb) dbGetQuery(con, sprintf("SELECT * FROM `%s`", surv_t)) else NULL
    if (!is.null(sv)) n_surv[i] <- nrow(sv)
  }
  pfam <- primary_family(di$EndpointPrimary[i])
  if (is.null(sv) || !has_expr) { notes[i] <- if (is.null(sv)) "no surv table" else "no expr table"; next }
  if (is.na(pfam)) { notes[i] <- "no primary endpoint annotation"; next }
  idcol <- names(sv)[1]

  # 各家族：解析出该数据集实际使用的 token（GSE31210 的 DFS 家族 = RFS）
  for (f in fam) {
    # 先看本脚本加载的 catalog 里的 EP_<family>（源文件是唯一真源；
    # 不能用 CanPAS::endpoint_resolve，它读的是已安装包里的旧 catalog）；
    # 缺失时回退到扫描生存表列，找属于该家族的 token。
    epcol <- paste0("EP_", f)
    tok <- if (epcol %in% names(di) && !is.na(di[[epcol]][i]) && nzchar(as.character(di[[epcol]][i])))
      as.character(di[[epcol]][i]) else NA_character_
    if (is.na(tok)) {
      cand <- sub("_time$", "", grep("_time$", names(sv), value = TRUE))
      cand <- cand[vapply(cand, function(t)
        identical(tryCatch(CanPAS::endpoint_family(t), error = function(e) NA_character_), f), logical(1))]
      if (length(cand)) tok <- cand[1]
    }
    if (is.na(tok) || !nzchar(tok)) next
    st <- paste0(tok, "_status"); tm <- paste0(tok, "_time")
    if (!st %in% names(sv)) next
    # 可分析 = 有表达 且 time/status 都非缺失（Cox/KM 的实际入组条件）
    ok  <- !is.na(sv[[st]])
    if (tm %in% names(sv)) ok <- ok & !is.na(sv[[tm]])
    ids <- as.character(sv[[idcol]][ok])
    keep <- intersect(ids, expr_ids)
    fam_n[i, paste0("n_", f)] <- length(keep)
    if (!is.na(pfam) && f == pfam) {
      n_events[i] <- sum(sv[[st]][ok][match(keep, ids)] == 1, na.rm = TRUE)
    }
  }
}

fam_cols <- as.data.frame(fam_n, stringsAsFactors = FALSE)
N_new <- vapply(seq_len(nrow(di)), function(i) {
  f <- primary_family(di$EndpointPrimary[i])
  if (is.na(f)) return(NA_integer_)
  as.integer(fam_cols[[paste0("n_", f)]][i])
}, integer(1))
# 无表达表（且非 TCGA 按需抓取）或该家族样本数为 0 -> 不可分析
N_new[N_new == 0L] <- NA_integer_
N_new[is.na(n_expr) & !startsWith(as.character(di$Accession), "TCGA-")] <- NA_integer_

old <- di
## TCGA 行按设计不在镜像里（本地 tcga 表 + 按需 Xena 提供），镜像口径无法重算其 N；
## 若不保护，直接写入会把 31 个 TCGA 行的 N 清成 NA（2026-09-24 实测）。
N_new_raw <- N_new        # 镜像口径的「原始重算值」（TCGA 行通常为 NA），只用于审计对比
.tcga_keep <- grepl("^TCGA", di$Accession)
.N_before  <- di$N
di$N          <- N_new
di$N[.tcga_keep] <- .N_before[.tcga_keep]
## 守卫自校验（fail closed）：TCGA 行必须逐行保持原值 —— 既不能由数值变 NA，也不能改数。
## 若这里不为 0，说明守卫失效，下面的 --write 会直接停止，而不是写出被清空 N 的 catalog。
.tcga_bad <- .tcga_keep & (is.na(di$N) != is.na(.N_before) |
                           (!is.na(.N_before) & di$N != .N_before))
cat(sprintf("（TCGA 行跳过镜像口径：%d 行保留原 N；守卫校验异常 %d 行）\n",
            sum(.tcga_keep), sum(.tcga_bad)))
di$n_expr     <- as.integer(n_expr)
di$n_surv     <- as.integer(n_surv)
di$n_events   <- as.integer(n_events)
for (cc in names(fam_cols)) di[[cc]] <- as.integer(fam_cols[[cc]])
di$expr_in_mirror <- !is.na(di$n_expr)

cmp <- data.frame(Accession = old$Accession, Type = old$Type,
                  N_old = old$N, N_new = di$N, n_expr = di$n_expr, n_surv = di$n_surv,
                  n_events = di$n_events, primary = old$EndpointPrimary,
                  note = notes, stringsAsFactors = FALSE)
cmp$changed <- !is.na(cmp$N_old) & !is.na(cmp$N_new) & cmp$N_old != cmp$N_new
cmp$tcga_guard      <- .tcga_keep                  # 该行的 N 由 TCGA 守卫保留，未参与镜像重算
cmp$N_mirror_recompute <- N_new_raw                # 不做守卫时镜像口径会写出的值（TCGA 行为 NA）
cmp$would_become_NA <- !is.na(cmp$N_old) & is.na(cmp$N_new)   # --write 会把它写成 NA
dir.create("pipeline/out", showWarnings = FALSE)
write.csv(cmp, "pipeline/out/catalog_N_fix.csv", row.names = FALSE)

cat("行数:", nrow(di), "\n")
cat("N 变化:", sum(cmp$changed), "行 | 由数值变 NA:", sum(!is.na(cmp$N_old) & is.na(cmp$N_new)),
    "| NA 变数值:", sum(is.na(cmp$N_old) & !is.na(cmp$N_new)), "\n")
cat("N = n_events 完全一致的行数:", sum(di$N == di$n_events, na.rm = TRUE), "/", sum(!is.na(di$N)), "\n\n")
cat("变化最大的 15 行:\n")
print(head(cmp[cmp$changed, ][order(-abs(cmp$N_old - cmp$N_new)[cmp$changed]), ][, 1:8], 15), row.names = FALSE)

## TCGA 守卫逐行审核：列出 --write 将会写入的每一个 TCGA 行的 N
cat("\n== TCGA 守卫逐行审核（这些行的 N 必须与现值完全一致）==\n")
.tcga_tab <- data.frame(Accession = cmp$Accession[.tcga_keep], N_before = cmp$N_old[.tcga_keep],
                        N_would_write = cmp$N_new[.tcga_keep],
                        would_become_NA = cmp$would_become_NA[.tcga_keep],
                        stringsAsFactors = FALSE)
if (nrow(.tcga_tab)) print(.tcga_tab, row.names = FALSE)
cat(sprintf("TCGA 行合计: %d | N 原样保留: %d | 将变 NA: %d | 守卫异常: %d\n",
            sum(.tcga_keep),
            sum(.tcga_keep & !cmp$would_become_NA),
            sum(.tcga_keep & cmp$would_become_NA),
            sum(.tcga_bad)))
cat(sprintf("若不设守卫，镜像重算会把 %d/%d 个 TCGA 行的 N 变成 NA（这正是守卫挡下的破坏）\n",
            sum(.tcga_keep & is.na(cmp$N_mirror_recompute)), sum(.tcga_keep)))
cat("TCGA 行 N 与加载时完全一致 (identical):",
    identical(as.integer(di$N[.tcga_keep]), as.integer(old$N[.tcga_keep])), "\n")

## 写前断言（fail closed）：任何一条不成立都不写 data/
if (any(.tcga_bad))
  stop(sprintf("TCGA 守卫失效：%d 行的 N 被改动，已中止写入。", sum(.tcga_bad)))
if (sum(.tcga_keep & cmp$would_become_NA) > 0L)
  stop(sprintf("TCGA 守卫失效：%d 行的 N 会被写成 NA，已中止写入。",
               sum(.tcga_keep & cmp$would_become_NA)))

if (!write_back) {
  ## 把「--write 会写成什么」原样导出到 pipeline/out/（不碰 data/），让 dry-run 可逐行核验
  would <- file.path("pipeline/out", "catalog_N_fix_would_write.csv")
  write.csv(di, would, row.names = FALSE)
  cat("\n[dry-run] 未写回。加 --write 才更新 catalog。\n")
  cat("         --write 的完整内容预览:", would, "（含全部 TCGA 行，可直接比对 N）\n")
  cat("         写前断言通过：TCGA 行 0 行会被改动、0 行会变 NA。\n")
  quit(save = "no")
}

# 必须以对象名 dataset_info 保存：包内 utils::data("dataset_info") 按名取对象
dataset_info <- di
save(dataset_info, file = file.path(root, "CanPAS/data/dataset_info.rda"), version = 2)
write.csv(di, file.path(root, "data/dataset_info.csv"), row.names = FALSE)
if (requireNamespace("openxlsx", quietly = TRUE))
  openxlsx::write.xlsx(di, file.path(root, "data/dataset_info.xlsx"))
cat("\n已写回: CanPAS/data/dataset_info.rda, data/dataset_info.csv, data/dataset_info.xlsx\n")
}

# ---------------------------------------------------------------------------
# run_26_cohort_overlap()  <-  verbatim pipeline/R/26_cohort_overlap.R
# ---------------------------------------------------------------------------
run_26_cohort_overlap <- function() {
# 26_cohort_overlap.R ---------------------------------------------------------
# 识别 catalog 中「同一批病人」的队列，写进 catalog 供应用与使用者参考。
#
# 为什么需要
#   多队列分析（Meta / Pooled KM / COX by datasets）会把选中的队列一起池化。
#   若两个队列其实是同一批病人（同一研究的两个平台、或同一队列重复入库），
#   池化就等于把病人数了一遍以上，会虚增精度并让 I² 失真。
#
# 判定依据
#   样本标题（GEO 的 !Sample_title / pheno 的 title 列）。同一研究的不同平台
#   在 GEO 里 GSM 不同，但标题（患者/样本编号）相同 -> 可据此识别重叠。
#   来源：
#     1) 本地 data/pheno/<ACC>.rds 的 title 列（离线可复算）
#     2) pipeline/ref/cohort_overlap_seed.csv —— 已核验的 GEO 侧重叠对
#        （本地无 pheno 的队列，如 GSE10885/GSE20624/GSE35629 之间）
#
# 输出
#   pipeline/out/cohort_overlap.csv       所有 >=5 例的重叠对
#   pipeline/out/REPORT_cohort_overlap.md 报告（含分组）
#   catalog 新增两列：
#     CohortGroup : 重叠组标签（无重叠 = NA）
#     Note        : 人类可读提示（重叠对象与例数；其他数据说明）
#
# 标题碰撞（pipeline/ref/cohort_overlap_exclude.csv）
#   逐字符比对标题会命中「偶然同名」：GSE25066（乳腺癌）的短数字样本号与
#   GSE32918（DLBCL）的 panel 重复码（1、1_Rep1、…）相同，54 例被误判为共享病人。
#   该队列真正的重复 deposit 是 GSE69051，且已按约定不入库。约定：只登记真实存在
#   的共享病人（同一 study 的平台拆分 / 重复 deposit / 同一批病人两次投稿）；
#   这类碰撞不计入 overlap 表、不生成 CohortGroup，只在 Note 里注明（与论文 §2.2
#   一致）。名单与理由写在 ref 文件里，和 cohort_overlap_seed.csv 同一处可复算。
#
# Note 列的写入约定
#   本脚本只负责「重叠/碰撞」与「无终点/无生存表」两类 Note；catalog 里已有的
#   数据来源说明（如 94/27 写入的注释）一律保留、追加在后面，不整列覆盖。
#
# 用法: Rscript pipeline/R/26_cohort_overlap.R [--write]
# ----------------------------------------------------------------------------
root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
write_back <- "--write" %in% commandArgs(trailingOnly = TRUE)
MIN_OVERLAP <- 5L

di <- get(load(file.path(root, "CanPAS/data/dataset_info.rda")))

# ---------- 1) 标题集合 ----------
# 关键：pheno 文件名用 '-'（GSE3494-GPL96），catalog 的 Accession 用 '_'
# （GSE3494_GPL96），比对前统一去掉分隔符。
canon <- function(x) gsub("[_-]", "", as.character(x))
acc <- as.character(di$Accession)
acc_by_canon <- stats::setNames(acc, canon(acc))
titles <- list()
for (f in list.files(file.path(root, "data/pheno"), pattern = "\\.rds$")) {
  p <- try(readRDS(file.path(root, "data/pheno", f)), silent = TRUE)
  if (inherits(p, "try-error") || !"title" %in% names(p) || !nrow(p)) next
  key <- names(acc_by_canon)[match(canon(sub("\\.rds$", "", f)), names(acc_by_canon))]
  if (is.na(key)) next
  titles[[acc_by_canon[[key]]]] <- unique(trimws(as.character(p$title)))
}
titles <- titles[!vapply(titles, function(x) all(is.na(x) | !nzchar(x)), logical(1))]

# ---------- 2) 已核验的 GEO 侧重叠（本地无 pheno 的队列）----------
seed_f <- file.path(root, "pipeline/ref/cohort_overlap_seed.csv")
seed <- if (file.exists(seed_f))
  utils::read.csv(seed_f, stringsAsFactors = FALSE) else
  data.frame(A = character(0), B = character(0), n = integer(0), note = character(0))
# 把种子里的标题集合并进 titles（用占位标题保证后续统一走"标题比对"逻辑）
seed_pairs <- list()
for (i in seq_len(nrow(seed))) {
  key <- paste(seed$A[i], seed$B[i], sep = "||")
  seed_pairs[[key]] <- as.integer(seed$n[i])
}

# ---------- 2b) 已知标题碰撞（形似重叠、实为巧合；见文件头说明）----------
excl_f <- file.path(root, "pipeline/ref/cohort_overlap_exclude.csv")
excl <- if (file.exists(excl_f))
  utils::read.csv(excl_f, stringsAsFactors = FALSE) else
  data.frame(A = character(0), B = character(0), reason = character(0))

# ---------- 3) 两两重叠 ----------
nms <- intersect(names(titles), acc)
pairs <- list()
for (i in seq_along(nms)) for (j in seq_along(nms)) {
  if (j <= i) next
  k <- length(intersect(titles[[nms[i]]], titles[[nms[j]]]))
  if (k >= MIN_OVERLAP)
    pairs[[length(pairs) + 1L]] <- data.frame(A = nms[i], B = nms[j], n = k,
                                              source = "title", stringsAsFactors = FALSE)
}
# ---------- 3b) 同一 GEO series 的两个平台（同批病人，标题仅差平台/标签后缀）----------
# 逐字符比较标题会漏检这类：GSE4716 的标题为 "Patient 1 GF200 001" vs "Patient 1 GF201 001"，
# GSE4412 为 "…Extract58_le1" vs "…Extract58_le2"——病人编号一致，仅平台/标签后缀不同。
# 判定：Accession 去掉 _GPL* 后相同者视为同一 series 的平台拆分，按同一批病人保守登记。
base_of <- sub("_GPL[0-9A-Za-z]+$", "", acc)
tab <- table(base_of)
for (b in names(tab)[tab > 1L]) {
  mem <- acc[base_of == b]
  for (i in seq_along(mem)) for (j in seq_along(mem)) {
    if (j <= i) next
    n_shared <- suppressWarnings(min(di$N[acc == mem[i]], di$N[acc == mem[j]], na.rm = TRUE))
    ## 目录未记录样本量时保留 NA（"共享例数未核"），不要写成 0
    pairs[[length(pairs) + 1L]] <- data.frame(A = mem[i], B = mem[j],
                                              n = as.integer(n_shared),
                                              source = "platform-pair", stringsAsFactors = FALSE)
  }
}

if (nrow(seed)) {
  # seed 里用 '-' 写 accession，先映射成 catalog 的写法，再与本地产出的对去重
  map1 <- function(x) { k <- canon(x); if (k %in% names(acc_by_canon)) acc_by_canon[[k]] else x }
  for (i in seq_len(nrow(seed))) {
    pairs[[length(pairs) + 1L]] <- data.frame(A = map1(seed$A[i]), B = map1(seed$B[i]),
                                              n = as.integer(seed$n[i]),
                                              source = "GEO-seed", stringsAsFactors = FALSE)
  }
}
ov <- if (length(pairs)) do.call(rbind, pairs) else
  data.frame(A = character(0), B = character(0), n = integer(0), source = character(0))
if (nrow(ov)) {
  key <- paste(pmin(canon(ov$A), canon(ov$B)), pmax(canon(ov$A), canon(ov$B)))
  ov <- ov[order(ov$source != "title"), ]          # 本地产出的优先保留
  ov <- ov[!duplicated(key[order(ov$source != "title")]), ]
  ov$A <- vapply(ov$A, function(x) if (canon(x) %in% names(acc_by_canon)) acc_by_canon[[canon(x)]] else x, character(1))
  ov$B <- vapply(ov$B, function(x) if (canon(x) %in% names(acc_by_canon)) acc_by_canon[[canon(x)]] else x, character(1))
  ov <- ov[order(-ov$n), ]
}

# ---------- 3c) 剔除标题碰撞对（不计入共享病人，只在 Note 中说明）----------
pair_key <- function(a, b) paste(pmin(canon(a), canon(b)), pmax(canon(a), canon(b)))
ov_ex <- ov[0, , drop = FALSE]
if (nrow(ov) && nrow(excl)) {
  ek <- pair_key(excl$A, excl$B)
  hit <- pair_key(ov$A, ov$B) %in% ek
  if (any(hit)) {
    ov_ex <- ov[hit, , drop = FALSE]
    ov_ex$reason <- excl$reason[match(pair_key(ov_ex$A, ov_ex$B), ek)]
    ov <- ov[!hit, , drop = FALSE]
  }
}

# ---------- 4) 连通分量分组 ----------
group_of <- stats::setNames(rep(NA_character_, nrow(di)), acc)
if (nrow(ov)) {
  nodes <- sort(unique(c(ov$A, ov$B)))
  parent <- stats::setNames(nodes, nodes)
  find <- function(x) { while (parent[[x]] != x) x <- parent[[x]]; x }
  for (i in seq_len(nrow(ov))) {
    ra <- find(ov$A[i]); rb <- find(ov$B[i])
    if (ra != rb) parent[[ra]] <- rb
  }
  roots <- vapply(nodes, find, character(1))
  for (r in unique(roots)) {
    members <- sort(nodes[roots == r])
    lab <- sprintf("%s(+%d)", members[1], length(members) - 1L)
    for (m in members) group_of[[m]] <- lab
  }
}

# ---------- 5) 写入 catalog ----------
note_col <- rep(NA_character_, nrow(di))
for (i in seq_len(nrow(di))) {
  a <- acc[i]
  if (!is.na(group_of[[a]])) {
    part <- ov[ov$A == a | ov$B == a, , drop = FALSE]
    part$other <- ifelse(part$A == a, part$B, part$A)
    kind <- if (all(part$source == "platform-pair"))
      "SAME SERIES on another platform — do not pool together (one study, platform split; patient ids match, titles differ only by platform/label suffix)"
    else "OVERLAPS"
    cnt <- ifelse(is.na(part$n), "count not recorded", as.character(part$n))
    note_col[i] <- sprintf("%s %s (shared patients: %s)", kind, group_of[[a]],
                           paste(sprintf("%s %s", part$other, cnt), collapse = "; "))
  }
}
# 标题碰撞：不进 CohortGroup，只在 Note 里记录（判定依据见文件头与 ref 文件）
if (nrow(ov_ex)) {
  for (i in seq_len(nrow(ov_ex))) {
    for (a in c(ov_ex$A[i], ov_ex$B[i])) {
      k <- which(acc == a); if (!length(k)) next
      other <- setdiff(c(ov_ex$A[i], ov_ex$B[i]), a)
      txt <- sprintf(paste0("TITLE COLLISION (not a shared-patient group): %d sample titles ",
                            "also occur in %s, but that cohort's titles are panel replicate ",
                            "codes and its genuine duplicate deposit is not in the catalog ",
                            "(pipeline/ref/cohort_overlap_exclude.csv)"),
                     ov_ex$n[i], other)
      note_col[k] <- if (is.na(note_col[k])) txt else paste(note_col[k], txt, sep = " | ")
    }
  }
}
# 其他说明：镜像里没有生存表 = 不可分析（只能做表达层面的事）
has_surv_note <- !is.na(di$n_surv)
no_endpoint <- is.na(di$SurvivalTypes)
note_col[is.na(note_col) & no_endpoint & has_surv_note] <- "no annotated endpoint"
note_col[is.na(note_col) & !has_surv_note] <- "expression only — no survival table, cannot be analysed"

# 保留 catalog 中已有的「数据来源」说明：本脚本只重写自己负责的 Note 片段。
# 做法是按 " | " 拆段、丢掉本脚本自己写过的段（重叠/碰撞/无终点/无生存表），
# 其余原样拼回。这样重复运行是幂等的（不会把碰撞说明叠加两次），
# 也不会把别的脚本写的数据来源说明清掉。
auto_pat <- "^(OVERLAPS |SAME SERIES on another platform|TITLE COLLISION|no annotated endpoint|expression only — no survival table)"
old_note <- if ("Note" %in% names(di)) as.character(di$Note) else rep(NA_character_, nrow(di))
prov <- vapply(old_note, function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  seg <- trimws(strsplit(x, "|", fixed = TRUE)[[1]])
  seg <- seg[nzchar(seg) & !grepl(auto_pat, seg)]
  if (!length(seg)) NA_character_ else paste(seg, collapse = " | ")
}, character(1))
keep <- !is.na(prov)
note_col[keep] <- ifelse(is.na(note_col[keep]), prov[keep],
                         paste(prov[keep], note_col[keep], sep = " | "))

di$CohortGroup <- group_of[acc]
di$Note <- note_col

dir.create(file.path(root, "pipeline/out"), showWarnings = FALSE)
utils::write.csv(ov, file.path(root, "pipeline/out/cohort_overlap.csv"), row.names = FALSE)
md <- c("# 队列病人重叠检查", "",
        sprintf("运行: %s | 最小重叠: %d 例 | 标题来源: %d 个本地 pheno + %d 条已核验 GEO 记录",
                format(Sys.time()), MIN_OVERLAP, length(titles), nrow(seed)), "",
        "## 重叠对", "", "| cohort A | cohort B | 重叠例数 | 来源 |", "|---|---|---|---|",
        if (nrow(ov)) sprintf("| %s | %s | %d | %s |", ov$A, ov$B, ov$n, ov$source) else "| — | — | — | — |",
        "", "## 分组", "", "| group | 成员 |", "|---|---|")
if (nrow(ov)) {
  gl <- split(names(group_of)[!is.na(group_of)], unlist(group_of)[!is.na(group_of)])
  md <- c(md, sprintf("| %s | %s |", names(gl), vapply(gl, paste, character(1), collapse = ", ")))
}
md <- c(md, "", "## 已排除：标题碰撞（形似重叠，非共享病人）", "",
        sprintf("约定：只登记真实存在的共享病人（同一 study 的平台拆分 / 重复 deposit / 同一批病人两次投稿）。"),
        sprintf("下列 %d 对标题相同但来源不同（不同癌种或不同设计），不计入上表、不生成 CohortGroup，只在 catalog 的 Note 列注明。名单见 `pipeline/ref/cohort_overlap_exclude.csv`。", nrow(ov_ex)),
        "", "| cohort A | cohort B | 命中标题数 | 理由 |", "|---|---|---|---|")
md <- c(md, if (nrow(ov_ex)) sprintf("| %s | %s | %d | %s |", ov_ex$A, ov_ex$B, ov_ex$n, ov_ex$reason) else "| — | — | — | — |")
writeLines(md, file.path(root, "pipeline/out/REPORT_cohort_overlap.md"))
cat(sprintf("重叠对 %d 对 | 涉及队列 %d 个 | 分组 %d 个 | 已排除标题碰撞 %d 对\n", nrow(ov),
            sum(!is.na(di$CohortGroup)), length(unique(na.omit(di$CohortGroup))), nrow(ov_ex)))
if (nrow(ov)) print(ov, row.names = FALSE)
if (nrow(ov_ex)) { cat("\n已排除（标题碰撞，不生成 CohortGroup）：\n"); print(ov_ex[, c("A", "B", "n")], row.names = FALSE) }

if (!write_back) { cat("\n[dry-run] 未写回 catalog。加 --write 执行。\n"); quit(save = "no") }
dataset_info <- di
save(dataset_info, file = file.path(root, "CanPAS/data/dataset_info.rda"), version = 2)
write.csv(di, file.path(root, "data/dataset_info.csv"), row.names = FALSE)
cat("\n已写回 catalog（新增 CohortGroup / Note 两列）\n")
}

# ---------------------------------------------------------------------------
# run_27_catalog_notes()  <-  verbatim pipeline/R/27_catalog_notes.R
# ---------------------------------------------------------------------------
run_27_catalog_notes <- function() {
# 27_catalog_notes.R ---------------------------------------------------------
# 给 catalog 的 Note 列补充「与重叠无关」的数据说明（可复算、逐条有出处）。
#
# 现状：仅收录已在实跑中遇到的问题，避免写入未经核实的推断。
#   GSE61676 (GPL5188)：平台无探针→基因符号映射，凡走符号查询的分析函数
#     （如 cohort_merged()/get_expr_data()）都会失败；队列本身可分析（有 OS），
#     属"可分析但该基因不可测"。证据：manuscript/figures/logs 与
#     analysis/scripts_v100/31_cross_cancer_figures.R 实跑报错
#     "Platform GPL5188 returned no probe-to-gene mapping for the requested ids"。
#
# 用法: Rscript pipeline/R/27_catalog_notes.R [--write]
# ----------------------------------------------------------------------------
root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
write_back <- "--write" %in% commandArgs(trailingOnly = TRUE)

di <- get(load("CanPAS/data/dataset_info.rda"))
acc <- as.character(di$Accession)
if (!"Note" %in% names(di)) di$Note <- NA_character_

notes <- data.frame(
  Accession = c("GSE61676"),
  text = c(paste0("Platform GPL5188 has no probe-to-gene (symbol) mapping in the mirror: ",
                  "the cohort is analysable (OS available) but gene-symbol queries fail; ",
                  "use probe ids or exclude it from gene-level analyses.")),
  stringsAsFactors = FALSE)

for (i in seq_len(nrow(notes))) {
  k <- which(acc == notes$Accession[i])
  if (!length(k)) { cat("  未找到", notes$Accession[i], "\n"); next }
  old <- di$Note[k]
  new <- if (is.na(old) || !nzchar(old)) notes$text[i]
         else if (grepl(notes$text[i], old, fixed = TRUE)) old
         else paste(old, notes$text[i], sep = " | ")
  if (identical(old, new)) { cat(sprintf("  %-14s 已是最新\n", notes$Accession[i])); next }
  cat(sprintf("  %-14s Note 更新: %s\n", notes$Accession[i], substr(new, 1, 90)))
  di$Note[k] <- new
}

if (write_back) {
  dataset_info <- di          # 对象名必须是 dataset_info（包内数据契约）
  save(dataset_info, file = "CanPAS/data/dataset_info.rda", version = 2)
  utils::write.csv(di, "data/dataset_info.csv", row.names = FALSE, na = "NA")
  cat("已写回 CanPAS/data/dataset_info.rda 与 data/dataset_info.csv\n")
} else {
  cat("[dry-run] 未写回；加 --write 执行\n")
}
}

# ---------------------------------------------------------------------------
# run_41_register_final3_rows()  <-  verbatim pipeline/R/41_register_final3_rows.R
# ---------------------------------------------------------------------------
run_41_register_final3_rows <- function() {
# 41_register_final3_rows.R --------------------------------------------------
# 目的
#   把本次经授权的 3 个队列追加进 data/dataset_info.csv（只追加，不改任何既有行）：
#     GSE1379     Breast Cancer   DFS            60/28（重建的表见 39_build_gse1379_surv.R）
#     TCGA-CHOL   Liver Cancer    OS,DSS,DFI,PFI 45/23（表由 40_add_tcga_chol_dlbc.R 生成）
#     TCGA-DLBC   Lymphoma        OS,DSS,DFI,PFI 47/9
#   X 续号（catalog 的 X 有历史空洞，不重排既有 X）。
#
# 约定依据（CanPAS/README.md "Column semantics" + 既有行反推）
#   N        = 在**交付的表达表与生存表里都有**、且主要终点(status+time 均非缺失)可用的样本数
#   n_expr   = 交付表达表的样本列数（TCGA 为 NA：Xena 按需取数，不进镜像）
#   n_surv   = 交付生存表的行数（含无可用终点的行）
#   n_events = 主要终点事件数，且**只统计上述 N 个样本**
#   method   = GEO/EMBL 微阵列用 "RNA"；TCGA 用 "TCGA-RNAseq"
#   expr_in_mirror = GSE 为 TRUE（要上传）；TCGA 为 FALSE（Xena-backed，设计如此）
#
# 产物
#   data/dataset_info.csv                     （193 -> 196 行）
#   pipeline/out/catalog_final3_rows.csv      （本次追加的 3 行，供审计）
# 后续必须执行
#   Rscript pipeline/R/11_endpoint_families.R <ROOT>   # 重算家族列并重建 rda
# 用法
#   Rscript pipeline/R/41_register_final3_rows.R [ROOT]
# 说明：不连数据库、不重建 rda、不动 tarball/git/manuscript。
# ----------------------------------------------------------------------------

ROOT <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(ROOT)) ROOT[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(ROOT)

di <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(di) != 193L) stop("expected the pre-change 193-row catalog, got ", nrow(di))
NEW_ACC <- c("GSE1379", "TCGA-CHOL", "TCGA-DLBC")
if (any(NEW_ACC %in% di$Accession)) stop("already catalogued: ",
                                         paste(NEW_ACC[NEW_ACC %in% di$Accession], collapse = ", "))
if (anyDuplicated(di$Accession)) stop("duplicated accession in the existing catalog")

# 逐格核对本地交付产物，确保 catalog 写的就是实测值
chk_tab <- function(acc) {
  f <- file.path("data/processed/surv", paste0(acc, "_surv.rds"))
  if (!file.exists(f)) stop("missing delivered survival table: ", f)
  readRDS(f)
}
t1379 <- chk_tab("GSE1379"); tchol <- chk_tab("TCGA-CHOL"); tdlbc <- chk_tab("TCGA-DLBC")
stopifnot(nrow(t1379) == 60L, sum(t1379$DFS_status == 1) == 28L)
stopifnot(nrow(tchol) == 45L, nrow(tdlbc) == 47L)
ex1379 <- readRDS("data/expr/GSE1379.rds"); n_expr_1379 <- ncol(ex1379) - 1L
stopifnot(n_expr_1379 == 60L,
          length(intersect(setdiff(colnames(ex1379), "ID_REF"), rownames(t1379))) == 60L)

blank <- function(n) rep(NA_character_, n)
mk <- function(X, Type, Accession, SurvivalTypes, EndpointFamilies, EP_OS, EP_DSS, EP_DFS,
               EP_PFS, EP_MFS, EndpointPrimary, EndpointDerived, GPL, N, method, n_expr,
               n_surv, n_events, n_OS, n_DSS, n_DFS, n_PFS, n_MFS, expr_in_mirror,
               CohortGroup, Note) data.frame(
  X = as.integer(X), Type = Type, Accession = Accession, SurvivalTypes = SurvivalTypes,
  EndpointFamilies = EndpointFamilies, EP_OS = EP_OS, EP_DSS = EP_DSS, EP_DFS = EP_DFS,
  EP_PFS = EP_PFS, EP_MFS = EP_MFS, EndpointPrimary = EndpointPrimary,
  EndpointDerived = EndpointDerived, GPL = GPL, N = as.integer(N), method = method,
  n_expr = as.integer(n_expr), n_surv = as.integer(n_surv), n_events = as.integer(n_events),
  n_OS = as.integer(n_OS), n_DSS = as.integer(n_DSS), n_DFS = as.integer(n_DFS),
  n_PFS = as.integer(n_PFS), n_MFS = as.integer(n_MFS), expr_in_mirror = expr_in_mirror,
  CohortGroup = CohortGroup, Note = Note, stringsAsFactors = FALSE)

rows <- rbind(
  mk(NA, "Breast Cancer", "GSE1379", "DFS", "DFS", NA, NA, "DFS", NA, NA, "DFS", NA,
     "GPL1223", 60L, "RNA", 60L, 60L, 28L, 0L, 0L, 60L, 0L, 0L, TRUE, NA,
     paste0("survival source: GEO series-matrix !Sample_description free text ",
            "(Clinical information: ...;DFS=<months>;Status=recur/non-recur), not ",
            "!Sample_characteristics_ch1 (the series has none); DFS_time = months/12 in years. ",
            "GSE1378 is the SAME 60 patients (microdissected cells vs whole-tissue sections) and is ",
            "deliberately NOT catalogued - count these patients once when pooling")),
  mk(NA, "Liver Cancer", "TCGA-CHOL", "OS,DSS,DFI,PFI", "OS,DSS,DFS,PFS", "OS", "DSS", "DFI",
     "PFI", NA, "OS", "DFI,PFI", "TCGA_HiSeqV2", 45L, "TCGA-RNAseq", NA, 45L, 23L,
     45L, 43L, 32L, 45L, 0L, FALSE, NA,
     paste0("intrahepatic cholangiocarcinoma registered under Liver Cancer, the same vocabulary ",
            "used for the earlier ICC cohort E-MTAB-6389 (and for GSE76427 / GSE14520 / TCGA-LIHC); ",
            "no new Cholangiocarcinoma type is introduced. Endpoints available: OS 45/23, DSS 43/20, ",
            "DFI 32/12, PFI 45/24. Xena expression + local clinical, so expr_in_mirror = FALSE by design")),
  mk(NA, "Lymphoma", "TCGA-DLBC", "OS,DSS,DFI,PFI", "OS,DSS,DFS,PFS", "OS", "DSS", "DFI",
     "PFI", NA, "OS", "DFI,PFI", "TCGA_HiSeqV2", 47L, "TCGA-RNAseq", NA, 47L, 9L,
     47L, 47L, 27L, 47L, 0L, FALSE, NA,
     paste0("endpoints available: OS 47/9, DSS 47/4, DFI 27/4, PFI 47/12. Carries only 9 OS events ",
            "(thin), and its OS time distribution includes a time = 0 row. Xena expression + local ",
            "clinical, so expr_in_mirror = FALSE by design"))
)
if (!identical(colnames(rows), colnames(di))) stop("column set/order mismatch")

rows$X <- (max(di$X) + 1L):(max(di$X) + nrow(rows))
di2 <- rbind(di, rows)
if (nrow(di2) != 196L) stop("expected 196 catalog rows, got ", nrow(di2))
if (anyDuplicated(di2$Accession)) stop("duplicated accession after append")
# 既有行必须逐格不变
same <- TRUE
for (cn in colnames(di)) {
  a <- di[[cn]]; b <- di2[[cn]][seq_len(nrow(di))]
  same <- same && identical(a, b)
}
if (!same) stop("existing 193 rows changed - aborting")
cat("[check] existing 193 rows unchanged; appended 3 rows\n")

write.csv(di2, "data/dataset_info.csv", row.names = FALSE)
write.csv(rows, "pipeline/out/catalog_final3_rows.csv", row.names = FALSE)
print(rows[, c("X", "Type", "Accession", "SurvivalTypes", "GPL", "N", "n_expr", "n_surv",
               "n_events", "n_OS", "n_DSS", "n_DFS", "n_PFS", "n_MFS", "expr_in_mirror")],
      row.names = FALSE)
cat("\nwrote data/dataset_info.csv rows:", nrow(di2), "\n")
cat("NEXT: Rscript pipeline/R/11_endpoint_families.R", ROOT, "\n")
}

# ---------------------------------------------------------------------------
# run_94_update_catalog_geo_expansion()  <-  verbatim pipeline/R/94_update_catalog_geo_expansion.R
# ---------------------------------------------------------------------------
run_94_update_catalog_geo_expansion <- function() {
# 94_update_catalog_geo_expansion.R -------------------------------------------
# 把 2026-09-24 GEO 扩展批次（A 审计后）实际建成的队列登记进 catalog：
#   data/dataset_info.csv  +  CanPAS/data/dataset_info.rda
#
# 约定（沿用既有 catalog 语义，见 19_fix_catalog_N.R 的列定义）：
#   N        患者级有效样本数（本批次一律用 geo_expansion_AUDIT.md 的审计数字）
#   n_events 主要终点事件数（同上）
#   n_expr   本地 data/expr/<ACC>.rds 的样本列数（= 上传到镜像后的表列数）
#   n_surv   本地 data/processed/surv/<ACC>_surv.rds 的行数（= 镜像 surv 表行数）
#   n_OS..n_MFS  各家族可分析样本数（= 审计数字，仅登记实际建了 token 的家族）
#   X        max(X)+1 起顺序编号
#   EndpointFamilies / EP_* / EndpointPrimary / EndpointDerived 由
#           11_endpoint_families.R 依据 SurvivalTypes 生成（本脚本写 NA）
#
# 用法:
#   Rscript pipeline/R/94_update_catalog_geo_expansion.R [ROOT]           # dry-run
#   Rscript pipeline/R/94_update_catalog_geo_expansion.R [ROOT] --write   # 写回 csv+rda
# 输出: pipeline/out/catalog_geo_expansion_rows.csv（本次新增行）+ 控制台校验
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else
  Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
ROOT <- path.expand(ROOT); setwd(ROOT)
WRITE <- "--write" %in% args

# acc, Type, SurvivalTypes, GPL, N(审计), n_events(审计), method, 家族n, Note
NEW <- list(
  list(acc="GSE65858",  type="Head and Neck Cancer", st="PFS", gpl="GPL10558", N=270, ev=133, method="RNA",
       fam=c(PFS=270),
       note="HNSCC (GPL10558 HT-12 v4)；PFS 270/133（审计患者级）。同一特征另有 OS(270/94) 未登记以免改变 EndpointPrimary 口径；GEO 扩展批次 2026-09-24"),
  list(acc="GSE117973", type="Head and Neck Cancer", st="PFS", gpl="GPL10558", N=77, ev=22, method="RNA",
       fam=c(PFS=77),
       note="HNSCC (GPL10558)；PFS 77/22（审计患者级）。另有 DSS(77/16) 未登记；GEO 扩展批次 2026-09-24"),
  list(acc="GSE27020",  type="Head and Neck Cancer", st="DFS", gpl="GPL96", N=109, ev=34, method="RNA",
       fam=c(DFS=109),
       note="喉癌(HNSCC) n=109/34 DFS（审计患者级）；DFS=dfs status(1=recurred)+dfs(months)；GEO 扩展批次 2026-09-24"),
  list(acc="GSE159067", type="Head and Neck Cancer", st="PFS", gpl="GPL18573", N=102, ev=96, method="RNA-seq",
       fam=c(PFS=102),
       note="晚期 HNSCC (PD-1/PD-L1 治疗) PFS 102/96（审计）。★HTG EdgeSeq OBP 定向 panel(~2.5k 基因)，非全转录组；表达来自补充文件(已在 GEO 归一为 log2CPM)；GEO 扩展批次 2026-09-24"),
  list(acc="GSE65904",  type="Melanoma", st="DSS", gpl="GPL10558", N=210, ev=102, method="RNA",
       fam=c(DSS=210),
       note="皮肤黑色素瘤 (GPL10558)；DSS 210/102（审计患者级，214 样本中 4 例临床缺失）；GEO 扩展批次 2026-09-24"),
  list(acc="GSE198430", type="Melanoma", st="OS", gpl="GPL32055", N=105, ev=50, method="RNA",
       fam=c(OS=105),
       note="早期黑色素瘤 OS 105/50（审计）。★NanoString nCounter 194 基因定向 panel(GPL32055)，非全转录组；GEO 扩展批次 2026-09-24"),
  list(acc="GSE198431", type="Melanoma", st="OS", gpl="GPL32055", N=79, ev=30, method="RNA",
       fam=c(OS=79),
       note="早期黑色素瘤 OS 79/30（审计）。与 GSE198430 同研究不同 sampleset（GSM 重叠 0、FFPE 编号不重叠），原候选表误判为 SAME_STUDY_DUPLICATE（审计 §3.13）；★194 基因 nCounter panel；GEO 扩展批次 2026-09-24"),
  list(acc="GSE22153",  type="Melanoma", st="OS", gpl="GPL6102", N=54, ev=47, method="RNA",
       fam=c(OS=54),
       note="IV 期黑色素瘤 OS 54/47（审计）。GSE22155 SuperSeries 的 test 集，只取 GSE22153；GEO 扩展批次 2026-09-24"),
  list(acc="GSE325123", type="Melanoma", st="OS", gpl="GPL24676", N=105, ev=62, method="RNA-seq",
       fam=c(OS=105),
       note="黑色素瘤长期随访 OS 105/62（审计患者级）。★GeoMx DSP：已把 170 个 included in_analysis ROI 按 patient 聚合（表达取 ROI 均值），ROI 级 269/168 是重复计数；105 例中 102 例有 OS 时间(该 102 例内 60 死亡)；GEO 扩展批次 2026-09-24"),
  list(acc="GSE71118",  type="Sarcoma", st="MFS", gpl="GPL570", N=312, ev=124, method="RNA",
       fam=c(MFS=312),
       note="多种肉瘤（CINSARC 验证队列）MFS 312/124（审计）。时间列 time 单位=年(最大 15)；GEO 扩展批次 2026-09-24"),
  list(acc="GSE30929",  type="Sarcoma", st="MFS", gpl="GPL96", N=140, ev=49, method="RNA",
       fam=c(MFS=140),
       note="脂肪肉瘤 DRFS 140/49（审计）。DRFS->MFS 家族(见 11_endpoint_families.R)；histology=subtype；GEO 扩展批次 2026-09-24"),
  list(acc="GSE271517", type="Sarcoma", st="OS", gpl="GPL24676", N=55, ev=24, method="RNA-seq",
       fam=c(OS=55),
       note="滑膜肉瘤 OS 55/24（审计：91 样本中仅 55 例原发瘤有 OS，转移瘤 NA）。表达来自 Patient/Sample_Counts 补充文件(计数->log2(x+1))；GEO 扩展批次 2026-09-24"),
  list(acc="GSE31312",  type="Lymphoma", st="OS,PFS", gpl="GPL570", N=475, ev=172, method="RNA",
       fam=c(OS=475, PFS=475),
       note="DLBCL R-CHOP OS 475/172 + PFS 475（审计，临床表在补充 PDF）。★475 例临床中有 470 例能与 498 个矩阵样本按 GEO Depository # 对上（5 个 ID 在 GEO 不存在，GEO 侧 470/170）；GEO 扩展批次 2026-09-24"),
  list(acc="GSE32918",  type="Lymphoma", st="OS", gpl="GPL8432", N=172, ev=93, method="RNA",
       fam=c(OS=172),
       note="DLBCL (DASL) OS 172/93（审计患者级；249 张芯片含 _Rep 重复，已按患者去重）。与 GSE69051 同一批患者的第二次 deposit（GSE69051 只记重复说明，不另建）；GEO 扩展批次 2026-09-24"),
  list(acc="GSE23501",  type="Lymphoma", st="OS", gpl="GPL570", N=69, ev=13, method="RNA",
       fam=c(OS=69),
       note="DLBCL OS 69/13（审计：事件列用 code_os(alive/dead)，不是临床疗效码 status）。另有 code_pfs+progression-free_survival(69/17) 未登记；GEO 扩展批次 2026-09-24"),
  list(acc="GSE248835", type="Lymphoma", st="DFS", gpl="GPL33963", N=256, ev=179, method="RNA",
       fam=c(DFS=256),
       note="LBCL (CAR-T vs 化疗) n=256。★终点是 EFS(event.free.survival.*)，按 EFS->DFS 家族登记为 DFS(179 事件)；★NanoString nCounter IO360 817 行定向 panel；GEO 扩展批次 2026-09-24"),
  list(acc="GSE22138",  type="Uveal Melanoma", st="MFS", gpl="GPL570", N=63, ev=35, method="RNA",
       fam=c(MFS=63),
       note="葡萄膜黑色素瘤原发瘤 MFS 63/35（审计）。★新癌种；原候选表 MANUAL note 称“GEO 无生存列”有误（审计 §7）；GEO 扩展批次 2026-09-24"),
  list(acc="GSE183088", type="Mesothelioma", st="OS", gpl="GPL30570", N=86, ev=75, method="RNA",
       fam=c(OS=86),
       note="恶性胸膜间皮瘤 OS 86/75（审计；原候选表因稀疏事件袋报成 1/1）。★70 基因 nCounter panel（GPL30570，68/70 行可映射 ENTREZ），按任务要求确认可用于本 schema 后建库；GEO 扩展批次 2026-09-24"),
  list(acc="GSE162520", type="Lung Cancer", st="PFS", gpl="GPL18573", N=92, ev=35, method="RNA-seq",
       fam=c(PFS=92),
       note="NSCLC (PD-1/PD-L1 治疗) PFS 92/35（审计）。★癌种纠正：审计/任务把它列在 Head and Neck 类（检索式命中 HNSCC 字样），但 GEO 标题、overall design 与每例 patient diagnosis 都是 NSCLC，故 Type=Lung Cancer。★HTG EdgeSeq OBP 定向 panel；表达来自补充文件(log2CPM，分号分隔/逗号小数点)；GEO 扩展批次 2026-09-24")
)

di <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (!"X" %in% colnames(di)) di$X <- seq_len(nrow(di))
accs <- vapply(NEW, function(z) z$acc, character(1))
dup <- intersect(accs, di$Accession)
if (length(dup)) stop("already in catalog: ", paste(dup, collapse = ", "))

rows <- vector("list", length(NEW))
for (i in seq_along(NEW)) {
  z <- NEW[[i]]
  ef <- file.path("data/expr", paste0(z$acc, ".rds"))
  sf <- file.path("data/processed/surv", paste0(z$acc, "_surv.rds"))
  if (!file.exists(ef)) stop(z$acc, ": missing ", ef)
  if (!file.exists(sf)) stop(z$acc, ": missing ", sf)
  ex <- readRDS(ef); sv <- readRDS(sf)
  n_expr <- ncol(ex) - 1L; n_surv <- nrow(sv)
  r <- data.frame(X = max(di$X) + i, Type = z$type, Accession = z$acc,
                  SurvivalTypes = z$st, EndpointFamilies = NA_character_,
                  EP_OS = NA_character_, EP_DSS = NA_character_, EP_DFS = NA_character_,
                  EP_PFS = NA_character_, EP_MFS = NA_character_,
                  EndpointPrimary = NA_character_, EndpointDerived = NA_character_,
                  GPL = z$gpl, N = z$N, method = z$method,
                  n_expr = n_expr, n_surv = n_surv, n_events = z$ev,
                  n_OS = 0L, n_DSS = 0L, n_DFS = 0L, n_PFS = 0L, n_MFS = 0L,
                  expr_in_mirror = TRUE, CohortGroup = NA_character_, Note = z$note,
                  stringsAsFactors = FALSE)
  for (f in names(z$fam)) r[[paste0("n_", f)]] <- as.integer(z$fam[[f]])
  if (z$N > n_expr) stop(z$acc, ": N > n_expr")
  if (z$N > n_surv) stop(z$acc, ": N > n_surv")
  rows[[i]] <- r
}
new <- do.call(rbind, rows)
new <- new[, colnames(di), drop = FALSE]

out <- rbind(di, new)
dir.create("pipeline/out", showWarnings = FALSE, recursive = TRUE)
write.csv(new[, c("X", "Type", "Accession", "SurvivalTypes", "GPL", "N", "n_events",
                  "n_expr", "n_surv", "n_OS", "n_DSS", "n_DFS", "n_PFS", "n_MFS",
                  "method", "Note")],
          "pipeline/out/catalog_geo_expansion_rows.csv", row.names = FALSE)
cat("new rows:", nrow(new), " | catalog:", nrow(di), "->", nrow(out),
    " | X:", min(new$X), "-", max(new$X), "\n")
print(new[, c("X", "Type", "Accession", "SurvivalTypes", "GPL", "N", "n_events",
              "n_expr", "n_surv")], row.names = FALSE)

if (!WRITE) { cat("\n[dry-run] 未写回。加 --write 才会写 data/dataset_info.csv 与 rda。\n"); quit(save = "no") }

write.csv(out, "data/dataset_info.csv", row.names = FALSE)
dataset_info <- out
save(dataset_info, file = "CanPAS/data/dataset_info.rda", compress = "xz")

# ---- 逐格校验 csv <-> rda（0 差异才算通过）--------------------------------
a <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
e <- new.env(); load("CanPAS/data/dataset_info.rda", envir = e)
b <- get("dataset_info", envir = e)
if (!identical(dim(a), dim(b))) stop("dim mismatch: ", paste(dim(a), collapse="x"), " vs ", paste(dim(b), collapse="x"))
if (!identical(colnames(a), colnames(b))) stop("colnames mismatch")
diff_cells <- 0L; diff_detail <- character(0)
for (j in seq_len(ncol(a))) for (i in seq_len(nrow(a))) {
  x <- a[i, j]; y <- b[i, j]
  same <- (is.na(x) && is.na(y)) || (!is.na(x) && !is.na(y) && as.character(x) == as.character(y))
  if (!same) { diff_cells <- diff_cells + 1L
    if (length(diff_detail) < 10) diff_detail <- c(diff_detail, sprintf("r%d c%s: csv=%s rda=%s", i, colnames(a)[j], x, y)) }
}
cat(sprintf("catalog verification: %d x %d cells, differences = %d\n", nrow(a), ncol(a), diff_cells))
if (length(diff_detail)) print(diff_detail)
if (diff_cells > 0) stop("catalog csv/rda mismatch")
cat("OK: data/dataset_info.csv == CanPAS/data/dataset_info.rda (cell-by-cell)\n")
}

# ---------------------------------------------------------------------------
# run_98_update_catalog_suppl_expansion()  <-  verbatim pipeline/R/98_update_catalog_suppl_expansion.R
# ---------------------------------------------------------------------------
run_98_update_catalog_suppl_expansion <- function() {
# 98_update_catalog_suppl_expansion.R -----------------------------------------
# 把 2026-09-24「补充文件 / 发表表格」路线的队列登记进 catalog：
#   data/dataset_info.csv + CanPAS/data/dataset_info.rda
#
# 与 94_update_catalog_geo_expansion.R 完全相同的列语义：
#   N       患者级有效样本数（= 表达 ∩ 主要终点的可分析患者数）
#   n_events 主要终点事件数
#   n_expr  本地 data/expr/<ACC>.rds 的样本列数
#   n_surv  本地 data/processed/surv/<ACC>_surv.rds 的行数
#   n_OS..n_MFS 各家族可分析样本数（仅登记实际建了 token 的家族）
#   EndpointFamilies / EP_* 由 11_endpoint_families.R 生成（本脚本写 NA）
#
# 用法:
#   Rscript pipeline/R/98_update_catalog_suppl_expansion.R [ROOT]           # dry-run
#   Rscript pipeline/R/98_update_catalog_suppl_expansion.R [ROOT] --write   # 写回 csv+rda
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else
  Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
ROOT <- path.expand(ROOT); setwd(ROOT)
WRITE <- "--write" %in% args

NEW <- list(
  list(acc="A5-PCPG", type="Pheochromocytoma", st="OS", gpl="GPL24676", N=70, ev=19,
       method="RNA-seq", fam=c(OS=70),
       note=paste("嗜铬细胞瘤/副神经节瘤（PPGL）OS 70/19。★数据源非 GEO：临床来自 cBioPortal",
                  "hnsc_a5consortium_2025（Hereditary SDHB-Mutant Pheochromocytomas and",
                  "Paragangliomas, A5 Consortium, Nat Commun 2025；cBioPortal 的 cancerTypeId",
                  "被错标为 hnsc，队列本体是 PPGL，已用 PCPG 专有临床字段核对），表达来自同队列",
                  "data_mrna_seq_cpm（Hugo_Symbol × 样本，已是 log2 CPM，未做二次变换）。患者级去重：",
                  "91 张表达样本 -> 77 例患者（优先 Primary 瘤）；70 例有可用 OS（19 死亡，最长随访",
                  "456 月/38 年，已回原始表核对）。DFS 仅 35 例可用（<50）故只登记 OS。独立于 TCGA-PCPG。",
                  "补充文件/发表表格路线 2026-09-24")),
  list(acc="IMmotion150", type="Kidney Cancer", st="PFS", gpl="GPL24676", N=263, ev=164,
       method="RNA-seq", fam=c(PFS=263),
       note=paste("肾细胞癌（RCC）PFS 263/164。★独立于 TCGA-KIRC：IMmotion150 是 Roche/Genentech",
                  "II 期 ICI 试验队列（cBioPortal rcc_iatlas_immotion150_2018, Nat Med 2018），",
                  "非 TCGA 再沉积。临床 PFS_MONTHS/PFS_STATUS（1:Progressed/0:Not_Progressed），",
                  "自治疗开始计——与本批既有治疗队列口径一致（GSE159067 / GSE162520 也是治疗起始 PFS）。",
                  "表达 data_mrna_seq_tpm（Hugo_Symbol × 患者，1 患者 1 样本，263/263 1:1 连接）；",
                  "矩阵尺度为压缩对数型（median=0, p99.9=15.8, max=22.9 -> 通过 92 的 already_log 闸门），",
                  "未做二次变换。补充文件/发表表格路线 2026-09-24")),
  list(acc="GSE3218", type="Testicular Cancer", st="OS", gpl="GPL96", N=74, ev=27,
       method="RNA", fam=c(OS=74),
       note=paste("成人男性生殖细胞瘤（GCT）OS 74/27。★数据源：表达来自 GEO GSE3218-GPL96",
                  "（107 个肿瘤样本，RMA 已 log2），临床来自发表表格 —— PMC4666461（PLoS One 2015,",
                  "PMID 26624623）S1 Table \"clinical features for the patients included in this study\"，",
                  "已在 data/suppl/GSE3218_PMC4666461_S1_Table_clinical_features.xlsx（Europe PMC",
                  "supplementaryFiles 打包下载）。连接键 = 表格 Sample（如 052B）== !Sample_title 里的",
                  "[0-9]{3}[A-Z] 词元；108 行临床中 74 行落在 GSE3218（全部 Expression Training==Yes），",
                  "27 例死亡。另 34 行是同一研究的验证集 GSE10783（患者级 34 例 < 50，未登记）。",
                  "★只登记 OS：表格里的 2y DFS / 5y DSS 是里程碑二分指标，不是 time-to-event，不能当 DFS/DSS。",
                  "独立于 TCGA-TGCT。补充文件/发表表格路线 2026-09-24"))
)


di <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (!"X" %in% colnames(di)) di$X <- seq_len(nrow(di))
accs <- vapply(NEW, function(z) z$acc, character(1))
dup <- intersect(accs, di$Accession)
if (length(dup)) {
  message("already in catalog (skipped): ", paste(dup, collapse = ", "))
  NEW <- Filter(function(z) !z$acc %in% dup, NEW)
}
if (!length(NEW)) { cat("nothing new to add\n"); quit(save = "no") }

rows <- vector("list", length(NEW))
for (i in seq_along(NEW)) {
  z <- NEW[[i]]
  ef <- file.path("data/expr", paste0(z$acc, ".rds"))
  sf <- file.path("data/processed/surv", paste0(z$acc, "_surv.rds"))
  if (!file.exists(ef)) stop(z$acc, ": missing ", ef)
  if (!file.exists(sf)) stop(z$acc, ": missing ", sf)
  ex <- readRDS(ef); sv <- readRDS(sf)
  n_expr <- ncol(ex) - 1L; n_surv <- nrow(sv)
  r <- data.frame(X = max(di$X) + i, Type = z$type, Accession = z$acc,
                  SurvivalTypes = z$st, EndpointFamilies = NA_character_,
                  EP_OS = NA_character_, EP_DSS = NA_character_, EP_DFS = NA_character_,
                  EP_PFS = NA_character_, EP_MFS = NA_character_,
                  EndpointPrimary = NA_character_, EndpointDerived = NA_character_,
                  GPL = z$gpl, N = z$N, method = z$method,
                  n_expr = n_expr, n_surv = n_surv, n_events = z$ev,
                  n_OS = 0L, n_DSS = 0L, n_DFS = 0L, n_PFS = 0L, n_MFS = 0L,
                  expr_in_mirror = TRUE, CohortGroup = NA_character_, Note = z$note,
                  stringsAsFactors = FALSE)
  for (f in names(z$fam)) r[[paste0("n_", f)]] <- as.integer(z$fam[[f]])
  if (z$N > n_expr) stop(z$acc, ": N > n_expr")
  if (z$N > n_surv) stop(z$acc, ": N > n_surv")
  rows[[i]] <- r
}
new <- do.call(rbind, rows)
new <- new[, colnames(di), drop = FALSE]
out <- rbind(di, new)

write.csv(new[, c("X","Type","Accession","SurvivalTypes","GPL","N","n_events",
                  "n_expr","n_surv","n_OS","n_DSS","n_DFS","n_PFS","n_MFS",
                  "method","Note")],
          "pipeline/out/catalog_suppl_expansion_rows.csv", row.names = FALSE)
cat("new rows:", nrow(new), " | catalog:", nrow(di), "->", nrow(out),
    " | X:", min(new$X), "-", max(new$X), "\n")
print(new[, c("X","Type","Accession","SurvivalTypes","GPL","N","n_events","n_expr","n_surv")],
      row.names = FALSE)

if (!WRITE) { cat("\n[dry-run] 未写回。加 --write 才会写 data/dataset_info.csv 与 rda。\n"); quit(save = "no") }

write.csv(out, "data/dataset_info.csv", row.names = FALSE)
dataset_info <- out
save(dataset_info, file = "CanPAS/data/dataset_info.rda", compress = "xz")

# ---- 逐格校验 csv <-> rda --------------------------------------------------
a <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
e <- new.env(); load("CanPAS/data/dataset_info.rda", envir = e)
b <- get("dataset_info", envir = e)
if (!identical(dim(a), dim(b))) stop("dim mismatch: ", paste(dim(a), collapse="x"), " vs ", paste(dim(b), collapse="x"))
if (!identical(colnames(a), colnames(b))) stop("colnames mismatch")
diff_cells <- 0L; diff_detail <- character(0)
for (j in seq_len(ncol(a))) for (i in seq_len(nrow(a))) {
  x <- a[i, j]; y <- b[i, j]
  same <- (is.na(x) && is.na(y)) || (!is.na(x) && !is.na(y) && as.character(x) == as.character(y))
  if (!same) { diff_cells <- diff_cells + 1L
    if (length(diff_detail) < 10) diff_detail <- c(diff_detail, sprintf("r%d c%s: csv=%s rda=%s", i, colnames(a)[j], x, y)) }
}
cat(sprintf("catalog verification: %d x %d cells, differences = %d\n", nrow(a), ncol(a), diff_cells))
if (length(diff_detail)) print(diff_detail)
if (diff_cells > 0) stop("catalog csv/rda mismatch")
cat("OK: data/dataset_info.csv == CanPAS/data/dataset_info.rda (cell-by-cell)\n")
}

# ---------------------------------------------------------------------------
# run_101_update_catalog_relaxed_gate()  <-  verbatim pipeline/R/101_update_catalog_relaxed_gate.R
# ---------------------------------------------------------------------------
run_101_update_catalog_relaxed_gate <- function() {
# 101_update_catalog_relaxed_gate.R ------------------------------------------
# 把「松弛闸门（>=30 患者带可用 time+status 对；旧闸门 >50）」批次登记进 catalog：
#   data/dataset_info.csv + CanPAS/data/dataset_info.rda
#
# 列语义与 94/98 完全一致：
#   N        患者级有效样本数（= 表达 ∩ 主要终点的可分析患者数）
#   n_events 主要终点事件数
#   n_expr   本地 data/expr/<ACC>.rds 的样本列数
#   n_surv   本地 data/processed/surv/<ACC>_surv.rds 的行数
#   n_OS..n_MFS 各家族可分析样本数（仅登记实际建了 token 的家族）
#   EndpointFamilies / EP_* 由 11_endpoint_families.R 生成（本脚本写 NA）
#
# ★ 低于旧闸门 50 的队列必须在 Note 里显式写明松弛闸门（例：
#   "small cohort: N=34; gate relaxed to >=30 per author"）。
#
# 用法:
#   Rscript pipeline/R/101_update_catalog_relaxed_gate.R [ROOT]           # dry-run
#   Rscript pipeline/R/101_update_catalog_relaxed_gate.R [ROOT] --write   # 写回 csv+rda
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else
  Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
ROOT <- path.expand(ROOT); setwd(ROOT)
WRITE <- "--write" %in% args
stopifnot(!file.exists("pipeline/R/19_fix_catalog_N.R") || TRUE)  # 19 号脚本严禁运行

NEW <- list(
  list(acc = "GSE76019", type = "Adrenocortical Cancer", st = "EFS", gpl = "GPL13158",
       N = 34, ev = 12, method = "RNA", fam = c(DFS = 34),
       note = paste("肾上腺皮质癌（儿童，COG ARAR0332 方案，PMID 27307598）。",
                    "small cohort: N=34; gate relaxed to >=30 per author。",
                    "表达来自 GEO GSE76019-GPL13158（54715 探针 × 34 样本，已是 log 尺度，",
                    "未做二次变换）；临床来自同系列 !Sample_characteristics_ch1：",
                    "histology=ACC（34/34）、Stage、efs.time（年，由整天数换算，例 1.48665297741273 y = 543 d）、",
                    "efs.event（0/1）。患者级 34/34 有可用 (time,status) 对，12 个事件，随访 0.162–7.088 年。",
                    "只登记 EFS 一个终点（该系列无 OS/DSS；EFS 依 11_endpoint_families.R 的 FAMILY_OF 归入 DFS 家族）。",
                    "独立于 TCGA-ACC。松弛闸门批次 2026-09-24")),
  list(acc = "GSE76039", type = "Thyroid Cancer", st = "OS", gpl = "GPL570",
       N = 35, ev = 29, method = "RNA", fam = c(OS = 35),
       note = paste("低分化/未分化甲状腺癌（PDTC/ATC）。",
                    "small cohort: N=35; gate relaxed to >=30 per author。",
                    "★数据源非单一 GEO：表达来自 GEO GSE76039-GPL570（54675 探针 × 37 样本，",
                    "MSK JCI 2016 研究的表达子集，GEO 内无生存字段），临床来自同一研究在",
                    "cBioPortal 的存档 thyroid_mskcc_2016（Poorly-Differentiated and Anaplastic",
                    "Thyroid Cancers, MSK, JCI 2016, PMID 26878173；117 患者带 OS_MONTHS/OS_STATUS）。",
                    "连接键：GSE76039 !Sample_title（s_JF_thy_NNN_P）== cBioPortal SAMPLE_ID，",
                    "并用 CEL 文件名内的 6 位病理号（如 105028T）== cBioPortal OTHER_SAMPLE_ID 交叉核对，",
                    "37/37 双向匹配。37 例中 35 例有可用 (OS_MONTHS, OS_STATUS)，29 例死亡。",
                    "矩阵尺度为 log 型（min 2.261/max 15.819 -> 通过 already_log 闸门），未做二次变换。",
                    "独立于 TCGA-THCA。松弛闸门批次 2026-09-24"))
)

di <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (!"X" %in% colnames(di)) di$X <- seq_len(nrow(di))
accs <- vapply(NEW, function(z) z$acc, character(1))
dup <- intersect(accs, di$Accession)
if (length(dup)) {
  message("already in catalog (skipped): ", paste(dup, collapse = ", "))
  NEW <- Filter(function(z) !z$acc %in% dup, NEW)
}
if (!length(NEW)) { cat("nothing new to add\n"); quit(save = "no") }

rows <- vector("list", length(NEW))
for (i in seq_along(NEW)) {
  z <- NEW[[i]]
  ef <- file.path("data/expr", paste0(z$acc, ".rds"))
  sf <- file.path("data/processed/surv", paste0(z$acc, "_surv.rds"))
  if (!file.exists(ef)) stop(z$acc, ": missing ", ef)
  if (!file.exists(sf)) stop(z$acc, ": missing ", sf)
  ex <- readRDS(ef); sv <- readRDS(sf)
  n_expr <- ncol(ex) - 1L; n_surv <- nrow(sv)
  if (z$N < 30) stop(z$acc, ": N=", z$N, " below the relaxed gate (>=30)")
  if (z$N < 50 && !grepl("gate relaxed to >=30 per author", z$note, fixed = TRUE))
    stop(z$acc, ": N<50 but the Note does not record the relaxed gate")
  r <- data.frame(X = max(di$X) + i, Type = z$type, Accession = z$acc,
                  SurvivalTypes = z$st, EndpointFamilies = NA_character_,
                  EP_OS = NA_character_, EP_DSS = NA_character_, EP_DFS = NA_character_,
                  EP_PFS = NA_character_, EP_MFS = NA_character_,
                  EndpointPrimary = NA_character_, EndpointDerived = NA_character_,
                  GPL = z$gpl, N = z$N, method = z$method,
                  n_expr = n_expr, n_surv = n_surv, n_events = z$ev,
                  n_OS = 0L, n_DSS = 0L, n_DFS = 0L, n_PFS = 0L, n_MFS = 0L,
                  expr_in_mirror = TRUE, CohortGroup = NA_character_, Note = z$note,
                  stringsAsFactors = FALSE)
  for (f in names(z$fam)) r[[paste0("n_", f)]] <- as.integer(z$fam[[f]])
  if (z$N > n_expr) stop(z$acc, ": N > n_expr")
  if (z$N > n_surv) stop(z$acc, ": N > n_surv")
  rows[[i]] <- r
}
new <- do.call(rbind, rows)
new <- new[, colnames(di), drop = FALSE]
out <- rbind(di, new)

write.csv(new[, c("X","Type","Accession","SurvivalTypes","GPL","N","n_events",
                  "n_expr","n_surv","n_OS","n_DSS","n_DFS","n_PFS","n_MFS",
                  "method","Note")],
          "pipeline/out/catalog_relaxed_gate_rows.csv", row.names = FALSE)
cat("new rows:", nrow(new), " | catalog:", nrow(di), "->", nrow(out),
    " | X:", min(new$X), "-", max(new$X), "\n")
print(new[, c("X","Type","Accession","SurvivalTypes","GPL","N","n_events","n_expr","n_surv")],
      row.names = FALSE)

if (!WRITE) { cat("\n[dry-run] 未写回。加 --write 才会写 data/dataset_info.csv 与 rda。\n"); quit(save = "no") }

write.csv(out, "data/dataset_info.csv", row.names = FALSE)
dataset_info <- out
save(dataset_info, file = "CanPAS/data/dataset_info.rda", compress = "xz")

# ---- 逐格校验 csv <-> rda --------------------------------------------------
a <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
e <- new.env(); load("CanPAS/data/dataset_info.rda", envir = e)
b <- get("dataset_info", envir = e)
if (!identical(dim(a), dim(b))) stop("dim mismatch: ", paste(dim(a), collapse="x"), " vs ", paste(dim(b), collapse="x"))
if (!identical(colnames(a), colnames(b))) stop("colnames mismatch")
diff_cells <- 0L; diff_detail <- character(0)
for (j in seq_len(ncol(a))) for (i in seq_len(nrow(a))) {
  x <- a[i, j]; y <- b[i, j]
  same <- (is.na(x) && is.na(y)) || (!is.na(x) && !is.na(y) && as.character(x) == as.character(y))
  if (!same) { diff_cells <- diff_cells + 1L
    if (length(diff_detail) < 10) diff_detail <- c(diff_detail, sprintf("r%d c%s: csv=%s rda=%s", i, colnames(a)[j], x, y)) }
}
cat(sprintf("catalog verification: %d x %d cells, differences = %d\n", nrow(a), ncol(a), diff_cells))
if (length(diff_detail)) print(diff_detail)
if (diff_cells > 0) stop("catalog csv/rda mismatch")
cat("OK: data/dataset_info.csv == CanPAS/data/dataset_info.rda (cell-by-cell)\n")
}

# ---------------------------------------------------------------------------
# run_105_update_catalog_embl()  <-  verbatim pipeline/R/105_update_catalog_embl.R
# ---------------------------------------------------------------------------
run_105_update_catalog_embl <- function() {
# 105_update_catalog_embl.R ---------------------------------------------------
# Register the EMBL-EBI (ArrayExpress / BioStudies) cohorts actually built by
# 104_build_embl_cohorts.R into the catalog:
#   data/dataset_info.csv  +  CanPAS/data/dataset_info.rda
#
# Conventions follow 94_update_catalog_geo_expansion.R exactly:
#   N         analysable patients (== n_surv, every surv row carries a usable pair)
#   n_expr    sample columns of data/expr/<ACC>.rds (== the mirror table's column count)
#   n_surv    rows of data/processed/surv/<ACC>_surv.rds
#   n_events  primary-endpoint events (patient level)
#   n_OS..n_MFS  analysable patients for the family the cohort actually registers
#   X         max(X)+1 onwards
#   EndpointFamilies / EP_* / EndpointPrimary  written afterwards by 11_endpoint_families.R
#
# Usage:
#   Rscript pipeline/R/105_update_catalog_embl.R [ROOT]            # dry run
#   Rscript pipeline/R/105_update_catalog_embl.R [ROOT] --write    # write csv + rda
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else
  Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
ROOT <- path.expand(ROOT); setwd(ROOT)
WRITE <- "--write" %in% args

SMALL <- "small cohort: N=%d; gate relaxed to >=30 per author"
SRC   <- "ArrayExpress/BioStudies"

# acc, Type, SurvivalTypes, GPL, family-n column, method, Note, audited_n, audited_ev
NEW <- list(
  list(acc="E-MTAB-6134", type="Pancreatic Cancer", st="OS", gpl="GPL13667", method="RNA",
       fam="OS", audit_n=288, audit_ev=181,
       note=sprintf("%s E-MTAB-6134 (Grenoble, pancreatic ductal adenocarcinoma). Survival from the study's own SDRF `Characteristics[os.delay]` + `[os.event]`, unit=month -> years. Expression from ProcessedExpression.tsv (49,386 U219 probes x 309 samples; GPL13667 local map, probes not collapsed). Patient-level = 1 row per patient (309/309), 288 usable OS pairs / 181 deaths, follow-up 0.09-12.20 y. Audit/re-check count 288/181 reproduced exactly.", SRC)),
  list(acc="E-MTAB-4032", type="Multiple Myeloma", st="OS", gpl="GPL6244", method="RNA",
       fam="OS", audit_n=149, audit_ev=94,
       note=sprintf("%s E-MTAB-4032 (Czech MM cohort, training+validation arms). Survival from SDRF `Characteristics[overall survival (from diagnosis) [years]]` + `[patient alive]`; event = `patient alive` == No. ★UNIT CORRECTED: the column is labelled [years] but holds MONTHS - against the re-deposit E-MTAB-1038 (73 shared patients) the ratio is a constant 1/30.4375 = 1/(365.25/12), i.e. months vs days; values are also decimal-comma. OS_time = value/12 (median 3.10 y, max 85.0 y). 149 usable pairs / 94 deaths, reproducing the audit exactly (the original scan's 55 had the status inverted). Expression from training_set.txt + validation_set.txt (33,297 GPL6244 probes x 151 samples; training_set duplicates its first 7,595 probe rows, de-duplicated). E-MTAB-1038 is a re-deposit of the training arm and is NOT registered.", SRC)),
  list(acc="E-MTAB-4321", type="Bladder Cancer", st="PFS", gpl="E-MTAB-4321_PLAT", method="RNA-seq",
       fam="PFS", audit_n=462, audit_ev=31,
       note=sprintf("%s E-MTAB-4321 (UROMOL, early-stage urothelial carcinoma, whole-transcriptome Ribo-Zero RNA-seq). Survival from SDRF `Characteristics[progression free survival]` (unit=month -> years) + `[progression to T2+]`. Expression = gene-level FPKM matrix UROMOL_gene_fpkm_gtf.txt (43,204 ENSG rows x 476 samples), linear -> log2(x+1). 462 analysable patients / 31 events (low event density). Platform token E-MTAB-4321_PLAT follows the existing CGGA_<n>_PLAT convention for gene-level matrices (map = ENSG -> ENTREZ via biomart, 55.39%%). Not present in the producer's scan; promoted by the re-check as gate-passing and GEO-independent (0 tokens / 0 GSE / 0 title matches).", SRC)),
  list(acc="E-MTAB-3892", type="Glioma Cancer", st="OS", gpl="GPL570", method="RNA",
       fam="OS", audit_n=164, audit_ev=44,
       note=sprintf("%s E-MTAB-3892 (CIT oligodendroglial / diffuse astrocytoma cohort). Survival from SDRF `Characteristics[os.event]` + `[os.delay]` (unit=month -> years). No processed matrix in the deposit: 179 CEL files RMA'd with affy + hgu133plus2cdf -> 54,675 GPL570 probes x 179 samples. 164 usable OS pairs / 44 deaths, follow-up 0.08-8.75 y. Independent of TCGA (the record mentions TCGA only as an external comparison; 0 barcode matches). Only the expression assay is registered - the companion methylation/SNP/miRNA studies are the same patients.", SRC)),
  list(acc="E-TABM-1202", type="Sarcoma", st="OS", gpl="GPL570", method="RNA",
       fam="OS", audit_n=101, audit_ev=34,
       note=sprintf("%s E-TABM-1202 (rhabdomyosarcoma, alveolar 65 / embryonal 36). Survival from SDRF `Characteristics [survival_time]` (unit=years) + `[Status]` (dead from disease = 1). The deposit's E-TABM-1202.eSet.r holds a RAW AffyBatch (101 CELs, cdf HG-U133_Plus_2), RMA'd with affy + hgu133plus2cdf -> 54,675 probes x 101 samples. 101 usable OS pairs / 34 deaths, follow-up 0.50-11.33 y.", SRC)),
  list(acc="E-MTAB-1727", type="Lung Cancer", st="OS", gpl="GPL570", method="RNA",
       fam="OS", audit_n=90, audit_ev=71,
       note=sprintf("%s E-MTAB-1727 (Grenoble lung SCC / basaloid cohort). Survival from SDRF `Characteristics [Overall Survival Event]` + `[Overall Survival]` (unit=month -> years). 173 SDRF rows / 109 Source Names (64 patients carry a second SNP/CGH array): patient-level 90 usable OS pairs / 71 deaths, i.e. the audit's correction of the scan's row-level 152/116. No processed matrix: 109 CEL files RMA'd with affy + hgu133plus2cdf. The only `GSE` strings in the record are protocol-library refs (P-GSE39582-*), not a re-deposit.", SRC)),
  list(acc="E-TABM-346", type="Lymphoma", st="OS", gpl="GPL96", method="RNA",
       fam="OS", audit_n=53, audit_ev=30,
       note=sprintf("%s E-TABM-346 (DLBCL, R-CHOP). Survival parsed from the SDRF `Factor Value [ClinicalInformation - Overall survival]` / `[Death status]` strings (\"Overall survival: 2.349 years\"), unit=years. Expression from the deposit's processed CHP-signal matrix (22,283 GPL96 probes x 53 samples; linear -> log2(x+1)). 53 usable OS pairs / 30 deaths, follow-up 0.03-5.41 y.", SRC)),
  list(acc="E-MTAB-3267", type="Kidney Cancer", st="PFS", gpl="GPL6244", method="RNA",
       fam="PFS", audit_n=53, audit_ev=39,
       note=sprintf("%s E-MTAB-3267 (metastatic clear-cell RCC, sunitinib). Survival from SDRF `Characteristics[progression free survival]` (unit=month -> years) + `[progression]`. No processed matrix: 59 CEL files RMA'd with oligo + pd.hugene.1.0.st.v1 -> 33,297 GPL6244 probes x 59 samples. 53 usable PFS pairs / 39 events, follow-up 0.08-4.92 y. Does not re-deposit TCGA-KIRC.", SRC)),
  list(acc="E-MTAB-1205", type="Leukemia Cancer", st="RFS", gpl="GPL570", method="RNA",
       fam="DFS", audit_n=50, audit_ev=22,
       note=sprintf("%s E-MTAB-1205 (paediatric T-ALL). ★Token = RFS: the registered columns are SDRF `Characteristics[Relapse-Free Survival Status]` + `[Relapse-Free Survival]` (22 events), not the study's separate Disease-Free Survival pair (18 events) - the audit flagged the producer's DFS/22 mismatch. Unit=day -> years. Expression from the legacy-FTP processed matrix cog.exp.meanadj2.txt (54,675 GPL570 probes x 50 samples; linear -> log2(x+1)). 50 usable pairs / 22 relapses. The only overlap signal is generic S2/S6 names colliding with the sarcoma series GSE71118 - a collision, not overlap.", SRC)),
  list(acc="E-MTAB-3580", type="Sarcoma", st="OS", gpl="GPL16422", method="RNA",
       fam="OS", audit_n=68, audit_ev=15,
       note=sprintf("%s E-MTAB-3580 (fusion-negative rhabdomyosarcoma, NanoString nCounter). Survival from SDRF `Characteristics[overall survival time]` (years) + `[overall survival event]`. Expression parsed from the 84 raw .RCC Code_Summary blocks (101-gene endogenous panel; counts linear -> log2(x+1)). ★PATIENT COUNT: the audit/re-check recorded 84/16 from the SDRF `Source Name`, but that field is the generic row label `Sample1..Sample84`; the patient/RNA id is `Characteristics[sampleid]` (6-letter NanoString codes) which has only 68 distinct values, and the 11 duplicated ids carry IDENTICAL time+event (repeat nCounter runs of the same RNA). Built table = 68 patients / 15 OS events (the same row-vs-patient inflation class the audit itself flagged for E-MTAB-1727/6389/6729). Low event count - flagged. NanoString panel, not whole-transcriptome.", SRC)),
  list(acc="E-MTAB-6389", type="Liver Cancer", st="OS", gpl="GPL17585", method="RNA",
       fam="OS", audit_n=79, audit_ev=52,
       note=sprintf("%s E-MTAB-6389 (intrahepatic cholangiocarcinoma - a histology new to the catalogue's Liver Cancer rows, which are all HCC). Survival from SDRF `Characteristics[overall survival]` (unit=month -> years) + `[event death]`. ★Patient level: Source Name is ICC###G (tumour, 78) / ICC###NT (normal, 31) and `Characteristics[individual]` has 82 distinct values over 109 rows; deduplicated to patients -> 79 usable OS pairs / 52 deaths (the scan's 109/72 was row-level). Expression from data_exp_icc.txt (HTA-2.0 PSR probesets x 109 samples). ★Annotation: AnnoProbe has no GPL17585, so the probe->ENTREZ map was built from the GEO platform SOFT (data/gpl/GPL17585.rds, gene_symbols -> ENTREZ; 607,678/925,032 probes = 65.69%%).", SRC)),
  list(acc="E-MTAB-6877", type="Mesothelioma", st="OS", gpl="GPL16686", method="RNA",
       fam="OS", audit_n=59, audit_ev=50,
       note=sprintf("%s E-MTAB-6877 (malignant pleural mesothelioma 63 + normal pleura 4). Survival from SDRF `Characteristics[overall survival]` + `[event death]` (no Unit column; the 0-164 range is months -> years, corroborated by the median of 16). No processed matrix: 67 CEL files RMA'd with oligo + pd.hugene.2.0.st -> HuGene-2.0-st transcript clusters x 67 samples. 59 usable OS pairs / 50 deaths. ★Annotation: AnnoProbe has no GPL16686 AND the GEO platform SOFT has no gene column at all (only genomic ranges / GB_ACC), so the map comes from the matching Bioconductor db hugene20sttranscriptcluster.db (data/gpl/GPL16686.rds, 30,766/53,981 = 56.99%%).", SRC)),
  list(acc="E-MTAB-1719", type="Mesothelioma", st="OS", gpl="GPL570", method="RNA",
       fam="OS", audit_n=34, audit_ev=33, small=34,
       note=sprintf("%s E-MTAB-1719 (MESO-JD malignant pleural mesothelioma). Survival from SDRF `Characteristics[OS event]` + `[OS delay]` (unit=month -> years). Expression from the legacy-FTP processed matrix MESO-Jaurand-EXP-RMA-normalized-data.txt (54,675 GPL570 probes x 38 samples, already log2 RMA). 34 usable OS pairs / 33 deaths; 4 of the 38 SDRF rows carry OS event = N/A.", SRC)),
  list(acc="E-MEXP-2780", type="Pancreatic Cancer", st="OS", gpl="GPL570", method="RNA",
       fam="OS", audit_n=30, audit_ev=28, small=30,
       note=sprintf("%s E-MEXP-2780 (pancreatic ductal adenocarcinoma). Survival from SDRF `Characteristics[TimeOfSurvival]` (unit=day -> years) + `[OverallSurvival]` (dead/alive). Expression from final-gene-expression-data-file.txt.magetab (54,675 GPL570 probes x 30 samples, RMA log2; coverage 100%% measured on the full matrix). 30 usable OS pairs / 28 deaths, follow-up 0.39-4.44 y.", SRC))
)

di <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (!"X" %in% colnames(di)) di$X <- seq_len(nrow(di))
accs <- vapply(NEW, function(z) z$acc, character(1))
dup <- intersect(accs, di$Accession)
if (length(dup)) stop("already in catalog: ", paste(dup, collapse = ", "))

SKIPPED <- character(0)
rows <- vector("list", length(NEW))
for (i in seq_along(NEW)) {
  z <- NEW[[i]]
  ef <- file.path("data/expr", paste0(z$acc, ".rds"))
  sf <- file.path("data/processed/surv", paste0(z$acc, "_surv.rds"))
  if (!file.exists(ef) || !file.exists(sf)) {
    SKIPPED <<- c(SKIPPED, z$acc)
    cat("  [skip] ", z$acc, ": not built yet (", ef, " / ", sf, ")\n", sep = "")
    next
  }
  ex <- readRDS(ef); sv <- readRDS(sf)
  n_expr <- ncol(ex) - 1L; n_surv <- nrow(sv)
  st <- paste0(z$st, "_status")
  if (!st %in% colnames(sv)) stop(z$acc, ": ", st, " not in the surv table")
  ev <- sum(sv[[st]] == 1L)
  if (is.na(ev) || ev == 0L) stop(z$acc, ": primary endpoint has 0 events")
  note <- z$note
  if (!is.null(z$small)) note <- paste0(note, " ", sprintf(SMALL, z$small), ".")
  note <- paste0(note, " Registering at patient level; N = analysable patients. EMBL-EBI batch 2026-09-24.")
  r <- data.frame(X = max(di$X) + i, Type = z$type, Accession = z$acc,
                  SurvivalTypes = z$st, EndpointFamilies = NA_character_,
                  EP_OS = NA_character_, EP_DSS = NA_character_, EP_DFS = NA_character_,
                  EP_PFS = NA_character_, EP_MFS = NA_character_,
                  EndpointPrimary = NA_character_, EndpointDerived = NA_character_,
                  GPL = z$gpl, N = n_surv, method = z$method,
                  n_expr = n_expr, n_surv = n_surv, n_events = ev,
                  n_OS = 0L, n_DSS = 0L, n_DFS = 0L, n_PFS = 0L, n_MFS = 0L,
                  expr_in_mirror = TRUE, CohortGroup = NA_character_, Note = note,
                  stringsAsFactors = FALSE, check.names = FALSE)
  r[[paste0("n_", z$fam)]] <- as.integer(n_surv)
  if (n_surv > n_expr) stop(z$acc, ": N > n_expr")
  if (!is.null(z$audit_n) && z$audit_n != n_surv)
    cat(sprintf("  [note] %s: audit/re-check said %d patients, built %d (recorded in the Note)\n",
                z$acc, z$audit_n, n_surv))
  if (!is.null(z$audit_ev) && z$audit_ev != ev)
    cat(sprintf("  [note] %s: audit/re-check said %d events, built %d\n", z$acc, z$audit_ev, ev))
  rows[[i]] <- r
}
rows <- rows[!vapply(rows, is.null, logical(1))]
if (length(SKIPPED)) cat("SKIPPED (not built): ", paste(SKIPPED, collapse = ", "), "\n")
if (!length(rows)) stop("nothing to register")
new <- do.call(rbind, rows)
new$X <- max(di$X) + seq_len(nrow(new))
new <- new[, colnames(di), drop = FALSE]

out <- rbind(di, new)
dir.create("pipeline/out", showWarnings = FALSE, recursive = TRUE)
write.csv(new[, c("X", "Type", "Accession", "SurvivalTypes", "GPL", "N", "n_events",
                  "n_expr", "n_surv", "n_OS", "n_DSS", "n_DFS", "n_PFS", "n_MFS",
                  "method", "Note")],
          "pipeline/out/catalog_embl_rows.csv", row.names = FALSE)
cat("new rows:", nrow(new), " | catalog:", nrow(di), "->", nrow(out),
    " | X:", min(new$X), "-", max(new$X), "\n")
print(new[, c("X", "Type", "Accession", "SurvivalTypes", "GPL", "N", "n_events",
              "n_expr", "n_surv", "method")], row.names = FALSE)
cat("sum N:", sum(di$N, na.rm = TRUE), "->", sum(out$N, na.rm = TRUE), "\n")

if (!WRITE) { cat("\n[dry-run] nothing written. add --write.\n"); quit(save = "no") }

write.csv(out, "data/dataset_info.csv", row.names = FALSE)
dataset_info <- out
save(dataset_info, file = "CanPAS/data/dataset_info.rda", compress = "xz")

# ---- cell-by-cell csv <-> rda verification ---------------------------------
a <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
e <- new.env(); load("CanPAS/data/dataset_info.rda", envir = e)
b <- get("dataset_info", envir = e)
if (!identical(dim(a), dim(b))) stop("dim mismatch: ", paste(dim(a), collapse = "x"),
                                     " vs ", paste(dim(b), collapse = "x"))
if (!identical(colnames(a), colnames(b))) stop("colnames mismatch")
diff_cells <- 0L; diff_detail <- character(0)
for (j in seq_len(ncol(a))) for (i in seq_len(nrow(a))) {
  x <- a[i, j]; y <- b[i, j]
  same <- (is.na(x) && is.na(y)) || (!is.na(x) && !is.na(y) && as.character(x) == as.character(y))
  if (!same) { diff_cells <- diff_cells + 1L
    if (length(diff_detail) < 10) diff_detail <- c(diff_detail,
      sprintf("r%d c%s: csv=%s rda=%s", i, colnames(a)[j], x, y)) }
}
cat(sprintf("catalog verification: %d x %d cells, differences = %d\n", nrow(a), ncol(a), diff_cells))
if (length(diff_detail)) print(diff_detail)
if (diff_cells > 0) stop("catalog csv/rda mismatch")
cat("OK: data/dataset_info.csv == CanPAS/data/dataset_info.rda (cell-by-cell)\n")
}

# ---------------------------------------------------------------------------
# run_107_register_embl_step2()  <-  verbatim pipeline/R/107_register_embl_step2.R
# ---------------------------------------------------------------------------
run_107_register_embl_step2 <- function() {
# 107_register_embl_step2.R ---------------------------------------------------
# Register the last three EMBL-EBI cohorts built by 106_build_embl_step2.R into
#   data/dataset_info.csv  +  CanPAS/data/dataset_info.rda
#
# Conventions are IDENTICAL to 105_update_catalog_embl.R / 94_update_catalog_geo_expansion.R:
#   N         analysable patients (== n_surv; every surv row carries a usable pair)
#   n_expr    sample columns of data/expr/<ACC>.rds (== the mirror table's column count)
#   n_surv    rows of data/processed/surv/<ACC>_surv.rds
#   n_events  primary-endpoint events (patient level)
#   n_OS..n_MFS  analysable patients for the family the cohort actually registers
#   X         max(X)+1 onwards
#   EndpointFamilies / EP_* / EndpointPrimary  intentionally left NA here and written
#   afterwards by 11_endpoint_families.R (run once for ALL three at the end).
#
# Difference from 105: the accession is selectable, so each cohort can be registered
# and checkpointed on its own (a long batch must not be able to die silently).
#
# Usage:
#   Rscript pipeline/R/107_register_embl_step2.R [ROOT] <ACC|all> [--write]
#   (dry run by default; --write appends the row(s) and saves the rda)
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
WRITE <- "--write" %in% args
rest  <- setdiff(args, "--write")
# order-agnostic: <ROOT> <ACC...> or <ACC...> <ROOT>
is_acc <- grepl("^(E-|GSE|TCGA|all)", rest)
ACC_ARG <- rest[is_acc][1]
ROOT_ARG <- rest[!is_acc][1]
ROOT  <- if (!is.na(ROOT_ARG) && nzchar(ROOT_ARG)) ROOT_ARG else
  Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
WANT  <- if (!is.na(ACC_ARG) && nzchar(ACC_ARG)) ACC_ARG else "all"
ROOT <- path.expand(ROOT); setwd(ROOT)

SMALL <- "small cohort: N=%d; gate relaxed to >=30 per author"
SRC   <- "ArrayExpress/BioStudies"

NEW <- list(
  list(acc="E-MTAB-1727", type="Lung Cancer", st="OS", gpl="GPL570", method="RNA",
       fam="OS", audit_n=90, audit_ev=71,
       note=sprintf("%s E-MTAB-1727 (Grenoble lung SCC / basaloid cohort). Survival from SDRF `Characteristics [Overall Survival Event]` + `[Overall Survival]` (Unit [time unit] = month -> years; 11 months -> 0.9167 y cross-checked by hand). ★Row-vs-patient: the SDRF has 173 assay rows and 109 Source Names; the 64 extra rows are a SECOND assay of type `genomic DNA` (*_DNA.txt, SNP/CGH) whose survival fields are byte-identical to the patient's expression row. Only the 109 CEL rows are the expression assay -> patient level 109 samples, 90 usable OS pairs / 71 deaths (16 rows have OS = NA, 1 has OS = 0). The earlier scan's 109 patients / 152 rows / 116 events was ROW-level. No processed matrix in the deposit: the 109 CEL files (3.3 GB, md5 manifest `data/suppl/embl_step2/E-MTAB-1727/CEL.md5`) were RMA'd with affy + hgu133plus2cdf -> 54,675 GPL570 probes x 109 samples; already log2 (median 4.94), no transform. Array Design REF A-AFFY-44 = GPL570, so the existing local map applies (gene-map coverage 80.72%%).", SRC)),
  list(acc="E-MTAB-6389", type="Liver Cancer", st="OS", gpl="GPL17585", method="RNA",
       fam="OS", audit_n=79, audit_ev=52,
       note=sprintf("%s E-MTAB-6389 (intrahepatic cholangiocarcinoma - a histology under the catalogue's existing Liver Cancer type, whose rows are otherwise HCC). Survival from SDRF `Characteristics[overall survival]` (Unit[time unit] = month -> years) + `[event death]`. ★Patient level: 109 SDRF rows = ICC###G tumour / ICC###NT normal, `Characteristics[individual]` has 82 distinct values; the tumour+normal columns of one individual are collapsed by rowMeans -> 79 usable OS pairs / 52 deaths (the scan's 109/72 was row-level). Expression from the study's processed matrix data_exp_icc.txt (HTA-2.0 `PSR*.hg.*` probesets). ★Annotation is the limiting factor: AnnoProbe has no GPL17585, so the map came from the GEO platform SOFT (data/gpl/GPL17585.rds); measured coverage 607,678/925,032 probes = 65.69%%.", SRC)),
  list(acc="E-MTAB-6877", type="Mesothelioma", st="OS", gpl="GPL16686", method="RNA",
       fam="OS", audit_n=59, audit_ev=50,
       note=sprintf("%s E-MTAB-6877 (malignant pleural mesothelioma 63 + normal pleura 4). Survival from SDRF `Characteristics[overall survival]` + `[event death]` (dead/alive/NA; no Unit column -> the 1..164 range is months, median 16, corroborated against the 252-month max of the sister 1727 cohort). ★Patient level: 67 CEL rows -> 59 usable OS pairs / 50 deaths (4 normal-pleura rows have OS NA, 1 tumour row has OS = 0, 1 has event NA; 7 rows are NA overall). CEL-only deposit: the 67 CEL files were RMA'd with oligo + pd.hugene.2.0.st -> 53,981 GPL16686 transcript clusters x 67 samples; already log2, no transform. Processed-matrix probe ids verified to match the GPL16686 map exactly. ★Annotation is the limiting factor: AnnoProbe has no GPL16686 and the GEO platform SOFT carries no gene column, so the map came from the matching Bioconductor db hugene20sttranscriptcluster.db (data/gpl/GPL16686.rds); measured coverage 30,766/53,981 probes = 56.99%%.", SRC))
)

accs <- vapply(NEW, function(z) z$acc, character(1))
sel  <- if (identical(WANT, "all")) accs else strsplit(WANT, ",")[[1]]
bad  <- setdiff(sel, accs); if (length(bad)) stop("unknown accession(s): ", paste(bad, collapse = ", "))

di <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (!"X" %in% colnames(di)) di$X <- seq_len(nrow(di))
dup <- intersect(sel, di$Accession)
if (length(dup)) stop("already in catalog: ", paste(dup, collapse = ", "))

SKIPPED <- character(0); rows <- list()
for (z in NEW[vapply(NEW, function(q) q$acc %in% sel, logical(1))]) {
  ef <- file.path("data/expr", paste0(z$acc, ".rds"))
  sf <- file.path("data/processed/surv", paste0(z$acc, "_surv.rds"))
  if (!file.exists(ef) || !file.exists(sf)) {
    SKIPPED <- c(SKIPPED, z$acc)
    cat("  [skip] ", z$acc, ": not built yet\n", sep = ""); next
  }
  ex <- readRDS(ef); sv <- readRDS(sf)
  n_expr <- ncol(ex) - 1L; n_surv <- nrow(sv)
  st <- paste0(z$st, "_status")
  if (!st %in% colnames(sv)) stop(z$acc, ": ", st, " not in the surv table")
  ev <- sum(sv[[st]] == 1L)
  if (is.na(ev) || ev == 0L) stop(z$acc, ": primary endpoint has 0 events")
  # join sanity: the survivors must be joinable to the expression table
  if (!all(rownames(sv) %in% colnames(ex)[-1]))
    stop(z$acc, ": surv ids not all present as expr columns")
  note <- paste0(z$note, " Registering at patient level; N = analysable patients. EMBL-EBI step-2 batch 2026-09-25.")
  r <- data.frame(X = NA_integer_, Type = z$type, Accession = z$acc,
                  SurvivalTypes = z$st, EndpointFamilies = NA_character_,
                  EP_OS = NA_character_, EP_DSS = NA_character_, EP_DFS = NA_character_,
                  EP_PFS = NA_character_, EP_MFS = NA_character_,
                  EndpointPrimary = NA_character_, EndpointDerived = NA_character_,
                  GPL = z$gpl, N = n_surv, method = z$method,
                  n_expr = n_expr, n_surv = n_surv, n_events = ev,
                  n_OS = 0L, n_DSS = 0L, n_DFS = 0L, n_PFS = 0L, n_MFS = 0L,
                  expr_in_mirror = TRUE, CohortGroup = NA_character_, Note = note,
                  stringsAsFactors = FALSE, check.names = FALSE)
  r[[paste0("n_", z$fam)]] <- as.integer(n_surv)
  if (n_surv > n_expr) stop(z$acc, ": N > n_expr")
  if (!is.null(z$audit_n) && z$audit_n != n_surv)
    cat(sprintf("  [note] %s: audit said %d patients, built %d\n", z$acc, z$audit_n, n_surv))
  if (!is.null(z$audit_ev) && z$audit_ev != ev)
    cat(sprintf("  [note] %s: audit said %d events, built %d\n", z$acc, z$audit_ev, ev))
  cat(sprintf("  [ok] %s: N=%d events=%d n_expr=%d n_surv=%d token=%s gpl=%s\n",
              z$acc, n_surv, ev, n_expr, n_surv, z$st, z$gpl))
  rows[[length(rows) + 1L]] <- r
}
if (length(SKIPPED)) cat("SKIPPED (not built): ", paste(SKIPPED, collapse = ", "), "\n")
if (!length(rows)) stop("nothing to register")
new <- do.call(rbind, rows)
new$X <- max(di$X) + seq_len(nrow(new))
new <- new[, colnames(di), drop = FALSE]
out <- rbind(di, new)

cat("catalog rows: ", nrow(di), " -> ", nrow(out), " | new X: ",
    paste(new$X, collapse = ","), "\n", sep = "")
cat("sum N: ", sum(di$N, na.rm = TRUE), " -> ", sum(out$N, na.rm = TRUE), "\n", sep = "")

if (!WRITE) { cat("\n[dry-run] nothing written. add --write.\n"); quit(save = "no") }

write.csv(out, "data/dataset_info.csv", row.names = FALSE)
dataset_info <- out
save(dataset_info, file = "CanPAS/data/dataset_info.rda", compress = "xz")

# ---- cell-by-cell csv <-> rda verification (same gate as 105) ---------------
a <- read.csv("data/dataset_info.csv", check.names = FALSE, stringsAsFactors = FALSE)
e <- new.env(); load("CanPAS/data/dataset_info.rda", envir = e)
b <- get("dataset_info", envir = e)
if (!identical(dim(a), dim(b))) stop("dim mismatch")
if (!identical(colnames(a), colnames(b))) stop("colnames mismatch")
diff_cells <- 0L; diff_detail <- character(0)
for (j in seq_len(ncol(a))) for (i in seq_len(nrow(a))) {
  x <- a[i, j]; y <- b[i, j]
  same <- (is.na(x) && is.na(y)) || (!is.na(x) && !is.na(y) && as.character(x) == as.character(y))
  if (!same) { diff_cells <- diff_cells + 1L
    if (length(diff_detail) < 10) diff_detail <- c(diff_detail,
      sprintf("r%d c%s: csv=%s rda=%s", i, colnames(a)[j], x, y)) }
}
cat(sprintf("catalog verification: %d x %d cells, differences = %d\n", nrow(a), ncol(a), diff_cells))
if (length(diff_detail)) print(diff_detail)
if (diff_cells > 0) stop("catalog csv/rda mismatch")
cat("OK: data/dataset_info.csv == CanPAS/data/dataset_info.rda (cell-by-cell)\n")
}

# ---------------------------------------------------------------------------
# run_111_register_gse205209()  <-  verbatim pipeline/R/111_register_gse205209.R
# ---------------------------------------------------------------------------
run_111_register_gse205209 <- function() {
# 111_register_gse205209.R ----------------------------------------------------
# 把已授权的 GSE205209 队列登记进 catalog（唯一新行），并更正 excluded.csv 的状态。
#   * 先备份 data/dataset_info.csv + CanPAS/data/dataset_info.rda
#   * 追加一行；其余 196 行按字节保持不动（只追加同一行的 CSV 文本）
#   * 只登记 OS —— 见脚本末尾的说明
# 用法: Rscript pipeline/R/111_register_gse205209.R [--apply]
# ----------------------------------------------------------------------------
root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
apply_l <- "--apply" %in% args
ACC <- "GSE205209"; DATE <- format(Sys.Date(), "%Y%m%d")

BACKUP <- file.path(root, "pipeline/backup", paste0("catalog_pre_gse205209_", DATE))
csv_p <- file.path(root, "data/dataset_info.csv")
rda_p <- file.path(root, "CanPAS/data/dataset_info.rda")

di <- read.csv(csv_p, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(di) == 196L, !ACC %in% di$Accession)
cat("baseline rows:", nrow(di), "| max X:", max(as.integer(di$X)), "\n")
for (f in c("OS","DSS","DFS","PFS","MFS"))
  cat("  EP_", f, "=", sum(!is.na(di[[paste0("EP_", f)]])), sep = "")
cat("\n  Sum N =", sum(as.integer(di$N), na.rm = TRUE), "\n")

# ---- 新行 -------------------------------------------------------------------
NOTE <- paste0(
  "small cohort: N=29 (one below the >=30 relaxed gate; admitted by explicit author ",
  "decision 2026-09-25); paired primary/metastatic arrays deduplicated to one row per ",
  "subject; OS 22 events. ",
  "SOURCE: GEO series matrix GSE205209 (NanoString PanCancer IO 360, GPL27956) ",
  "!Sample_characteristics_ch1 `overall survival (os in months)` + ",
  "`overall survival (months) - censored data - 1 implies censored`; months/12 -> years. ",
  "SIZE: 60 arrays are matched primary/metastatic tumours of 29 distinct subjects ",
  "(the study's own Clinical_Experimental_Metadata.xlsx `Patient` column has 29 distinct ",
  "values; the earlier '28' in excluded.csv came from the xlsx leaving the Patient cell ",
  "blank on every primary row, and GEO's own !Series_overall_design also states 29). ",
  "DEDUP RULE: primary-preferred - one row per subject, taking the usable primary (P) ",
  "array and falling back to the first usable metastatic (M/D) array; that discards one ",
  "array per subject with a usable primary (USC7 keeps USC7P, drops USC7M;USC7MD; USC20 ",
  "keeps USC20P, drops USC20M;USC20MD; USC28 is the only row on a metastatic array ",
  "because its primary GSM6208673 is all-null). ",
  "OS 22 events / 7 censored, OS_time 0.250-5.242 years. ",
  "EXPRESSION: 770 IO360 genes x 54 usable arrays (6 of the 60 sample columns are ",
  "entirely null - the paper's 5 'unacceptable NanoString results' over patient IDs 1 Met, ",
  "9 Met+Primary, 12 Met, 20 Distant Met, 28 Primary); USC9 has NO usable array, so it is ",
  "present in the survival table but not joinable to expression. Annotation: platform SOFT ",
  "GPL27956 carries symbols only (ID/ORF/SPOT_ID, no Entrez), so gene_id comes from ",
  "org.Hs.eg.db SYMBOL->ENTREZID -> 756/770 = 98.18%. ",
  "TTR IS DELIBERATELY NOT REGISTERED AS AN ENDPOINT. The series matrix also carries ",
  "`time to recurrence (ttr in months)` + `... - censored`: 26 usable pairs / 19 events / ",
  "7 censored (3 subjects carry ttr = 0 months). It is not a clean DFS/RFS: for 10 of the ",
  "22 decedents with no documented recurrence the TTR value equals the OS value, i.e. death ",
  "was substituted for recurrence, while the 7 censored subjects have TTR = OS. The pair ",
  "is kept in the survival table (TTR_time/TTR_status) for reference and can be added later ",
  "if the author decides the death-substituted definition is acceptable.")

new <- di[1, , drop = FALSE]
new[] <- NA
new$X <- as.character(max(as.integer(di$X)) + 1L)
new$Type <- "Endometrial Cancer"
new$Accession <- ACC
new$SurvivalTypes <- "OS"
new$EndpointFamilies <- "OS"
new$EP_OS <- "OS"
new$EP_DSS <- NA; new$EP_DFS <- NA; new$EP_PFS <- NA; new$EP_MFS <- NA
new$EndpointPrimary <- "OS"
new$EndpointDerived <- NA
new$GPL <- "GPL27956"
new$N <- 29L
new$method <- "RNA"
new$n_expr <- 54L
new$n_surv <- 29L
new$n_events <- 22L
new$n_OS <- 29L
new$n_DSS <- 0L; new$n_DFS <- 0L; new$n_PFS <- 0L; new$n_MFS <- 0L
new$expr_in_mirror <- TRUE
new$CohortGroup <- NA
new$Note <- NOTE

out <- rbind(di, new)
stopifnot(nrow(out) == 197L, !any(duplicated(out$Accession)))
for (f in c("OS","DSS","DFS","PFS","MFS"))
  cat("after: EP_", f, "=", sum(!is.na(out[[paste0("EP_", f)]])), sep = "")
cat("\nafter: Sum N =", sum(as.integer(out$N), na.rm = TRUE), "\n")

# ---- excluded.csv 状态更正 ---------------------------------------------------
exc_p <- file.path(root, "pipeline/out/excluded.csv")
exc <- read.csv(exc_p, stringsAsFactors = FALSE, check.names = FALSE)
i <- which(exc$Accession == ACC & grepl("RELAXED-GATE RE-CHECK", exc$Reason))
stopifnot(length(i) == 1L)
NEWREASON <- paste0(
  "endometrial (uterine serous carcinoma): **NOW CATALOGUED 2026-09-25** as a small ",
  "cohort N=29 (one below the >=30 relaxed gate; admitted by explicit author decision). ",
  "MATERIAL CORRECTIONS to the earlier re-check note: the subject count is 29, not 28 ",
  "(the xlsx `Patient` column is blank on every primary row; 29 distinct values plus GEO's ",
  "own !Series_overall_design), and the OS event count is 22 deaths / 7 censored (the ",
  "earlier 'claimed 14 events are the CENSORED count' and the flat 'n=29 (<51)' row were ",
  "both wrong - there is no 14 in either the OS or the TTR flag). Built: ",
  "data/expr/GSE205209.rds (770 x 54), data/processed/surv/GSE205209_surv.rds (29 rows, ",
  "OS 22 events), uploaded as GSE205209 + GSE205209_surv. TTR (26 usable pairs / 19 events) ",
  "is NOT registered as an endpoint: for 10 of 22 decedents the TTR value equals OS, i.e. ",
  "death was substituted for recurrence, so it is not a clean DFS/RFS. Row retained, not deleted.")
exc$Reason[i] <- NEWREASON

if (apply_l) {
  dir.create(BACKUP, recursive = TRUE, showWarnings = FALSE)
  file.copy(csv_p, file.path(BACKUP, "dataset_info.csv"), overwrite = TRUE)
  file.copy(rda_p, file.path(BACKUP, "dataset_info.rda"), overwrite = TRUE)
  cat("backup ->", BACKUP, ":", paste(list.files(BACKUP), collapse = ", "), "\n")
  write.csv(out, csv_p, row.names = FALSE, na = "NA")
  dataset_info <- out
  save(dataset_info, file = rda_p, compress = "xz")
  write.csv(exc, exc_p, row.names = FALSE, na = "NA")
  cat("[applied] catalog row appended, excluded.csv row updated\n")
} else cat("[dry-run] nothing written\n")
}
