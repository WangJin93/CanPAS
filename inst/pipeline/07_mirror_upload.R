# ===========================================================================
# CanPAS curation pipeline -- 07_mirror_upload
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 06_upload_db.R
#   2. 09_update_db_surv.R
#   3. 13_upload_small_cohorts.R
#   4. 24_make_tables_writable.R
#   5. 42_upload_gpl1223.R
#   6. 99_extend_gpl_db.R
#   7. 112_upload_gpl27956.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "07_mirror_upload.R", package = "CanPAS"))
#   run_06_upload_db()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================


# ---------------------------------------------------------------------------
# run_06_upload_db()  <-  verbatim pipeline/R/06_upload_db.R
# ---------------------------------------------------------------------------
run_06_upload_db <- function() {
# 06_upload_db.R
# -----------------------------------------------------------------------------
# CanPAS 数据处理 第 6 步（可选，仅在你授权后执行）：把本地产物写入远端 MySQL cpas 库。
# 严格沿用既有脚本的表结构约定：
#   expr 表 : row_names=探针(VARCHAR) + GSM 列(FLOAT)     —— 同 data/processed/upload.R
#   GPL 表  : row_names=探针 + gene_id(拆 “///” 后逐行)   —— 同 data/processed/gpl_app.R
#   surv 表 : row_names=样本GSM + 生存/临床列             —— 同 data/processed/app.R
# 用法: Rscript 06_upload_db.R <Accession> [<ProjectRoot>]
#   Accession 需在 dataset_info.csv 中登记；自动附带上传其 GPL 注释表（若 DB 缺失）。
# -----------------------------------------------------------------------------
suppressMessages({library(dplyr); library(tidyr); library(stringr); library(RMySQL)})

args <- commandArgs(trailingOnly=TRUE)
ACC  <- args[1]
ROOT <- if (length(args) >= 2) args[2] else "~/data/Project/CanPAS"
ROOT <- path.expand(ROOT)
setwd(ROOT)
FORCE_ALL  <- length(args) >= 3 && args[3] == "overwrite"
FORCE_SURV <- FORCE_ALL || (length(args) >= 3 && args[3] == "overwrite-surv")
# overwrite-expr：只强制上传表达表（N<=79 的小队列要用 13_upload_small_cohorts.R 的
# 同一口径入库），但**不**强制上传 GPL 注释表 —— 避免用本地 GPL rds 覆盖镜像里
# 内容不同（例如 GPL570 镜像 48052 行 vs 本地 54675 行）的既有注释表。
FORCE_EXPR <- FORCE_ALL || (length(args) >= 3 && args[3] == "overwrite-expr")

db <- list(host="139.224.80.159", dbname="cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CPAS"),
           password = Sys.getenv("CPAS_DB_PASSWORD"))
con <- dbConnect(MySQL(), host=db$host, dbname=db$dbname, user=db$user,
                 password=db$password, client_flag=CLIENT_COMPRESS)
on.exit(dbDisconnect(con))

info <- read.csv("data/dataset_info.csv", stringsAsFactors=FALSE)
row <- info[info$Accession == ACC, ]
if (!nrow(row)) row <- info[info$Accession == gsub("_", "-", ACC), ]
if (!nrow(row)) stop("Accession not in dataset_info.csv: ", ACC)
ACC_FS <- if (file.exists(file.path(ROOT,"data/expr", paste0(ACC, ".rds")))) ACC else gsub("_","-",ACC)  # 本地文件用 '-'
rf <- function(..., ext="") { p <- file.path(...); if (file.exists(p)) p else file.path(...) }
GPL <- row$GPL[1]
# 上传闸门看「表达表样本数」，不用可分析样本数 N（N 现已改为「表达 ∩ 主要终点」
# 的可分析样本数，见 19_fix_catalog_N.R）。取值顺序：
#   1) catalog 的 n_expr（镜像实测）
#   2) 本地表达文件的样本列数（新数据集，镜像里还没有表）
#   3) 旧 catalog 的 N（兼容）
n_expr_guess <- if ("n_expr" %in% names(row)) suppressWarnings(as.integer(row$n_expr[1])) else NA_integer_
if (is.na(n_expr_guess)) {
  f_loc <- file.path(ROOT, "data/expr", paste0(ACC_FS, ".rds"))
  n_expr_guess <- if (file.exists(f_loc)) ncol(readRDS(f_loc)) - 1L else
                  suppressWarnings(as.integer(row$N[1]))
}
N_EXPR <- if (is.na(n_expr_guess)) 0L else as.integer(n_expr_guess)

# ---------- 1) GPL 注释表（镜像缺表就上传，与 n_expr 无关） ----------
# 两个闸门刻意分开，不要混用：
#   GPL 闸门   : 镜像里**没有**这张平台注释表就上传，与 n_expr 大小**无关**；镜像里
#                **已有**该表则永不上传/覆盖（大小队列一律如此，overwrite-expr 也不强制
#                GPL —— 避免用本地 rds 覆盖内容不同的既有注释表，例如 GPL570：镜像
#                48052 行 vs 本地 54675 行）。旧版把 GPL 上传挂在 n_expr>79 上，于是
#                n_expr<=79 的小队列注释表被静默跳过，GPL1223 / GPL16686 / GPL27956
#                三个平台各需一个手工补救脚本（42_ / 112_ / embl step-2 fix）。
#   expr 闸门  : n_expr <= 79 的队列跳过表达表（沿用 upload.R 惯例），唯一逃逸口是
#                overwrite / overwrite-expr。
gpl_tables <- dbListTables(con)
tab_gpl <- str_replace_all(GPL, "-", "_")
gpl_src <- file.path(ROOT, "data/processed/gpl", paste0(GPL, ".rds"))
if (!tab_gpl %in% gpl_tables) {
  if (!file.exists(gpl_src)) {
    message("GPL table absent from DB but no local map to upload from: ", tab_gpl,
            " (", gpl_src, ")")
  } else {
  g <- readRDS(gpl_src)
  g2 <- g %>%
    mutate(gene_id = ifelse(gene_id == "---", NA, gene_id)) %>%
    na.omit() %>%
    mutate(row_names = rownames(.)) %>%
    mutate(gene_id = strsplit(as.character(gene_id), " /// ")) %>%
    unnest(gene_id) %>%
    select(row_names, gene_id)
  dbWriteTable(con, name=tab_gpl, value=g2, row.names=FALSE, overwrite=TRUE)
  message("GPL table uploaded: ", tab_gpl, " (", nrow(g2), " probe-gene rows)")
  }
} else message("GPL table already in DB (never overwritten): ", tab_gpl)

# ---------- 2) 表达谱表（沿用 upload.R 惯例：仅 N > 79 的队列上传表达谱） ----------
tab_expr <- str_replace_all(ACC, "-", "_")
N <- N_EXPR                      # 表达表样本数（见文件头注释）
if (FORCE_EXPR || !tab_expr %in% dbListTables(con)) {
  if (N > 79 || FORCE_EXPR) {
    dd <- readRDS(file.path(ROOT, "data/expr", paste0(ACC_FS, ".rds")))
    g <- readRDS(file.path(ROOT, "data/processed/gpl", paste0(GPL, ".rds")))
    dd <- dd %>% filter(ID_REF %in% rownames(g)) %>%
      tibble::column_to_rownames("ID_REF")
    dd[] <- lapply(dd, function(x) as.numeric(x))
    field_types <- c("row_names" = "VARCHAR(255)",
                     setNames(rep("FLOAT", ncol(dd)), colnames(dd)))
    dbWriteTable(con, name=tab_expr, value=dd, overwrite=TRUE, row.names=TRUE,
                 field.types=field_types)
    message("expr table uploaded: ", tab_expr, " (", nrow(dd), " probes x ",
            ncol(dd), " samples)")
  } else {
    message("expr table skipped: n_expr=", N, " <= 79 (与 upload.R 惯例一致)")
  }
} else message("expr table already in DB: ", tab_expr)

# ---------- 3) 生存/临床表 ----------
tab_surv <- paste0(str_replace_all(ACC, "-", "_"), "_surv")
if (FORCE_SURV || !tab_surv %in% dbListTables(con)) {
  sv <- readRDS(file.path(ROOT, "data/processed/surv", paste0(ACC_FS, "_surv.rds")))
  sv_out <- data.frame(row_names = rownames(sv), sv, check.names=FALSE,
                       stringsAsFactors=FALSE)
  dbWriteTable(con, name=tab_surv, value=sv_out, row.names=FALSE, overwrite=TRUE)
  message("surv table uploaded: ", tab_surv, " (", nrow(sv_out), " samples)")
} else message("surv table already in DB: ", tab_surv)

# ---------- 校验 ----------
for (t in c(tab_gpl, tab_expr, tab_surv)) {
  if (t %in% dbListTables(con)) {
    n <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", t))$n
    cat(sprintf("  [ok] %-18s rows=%d\n", t, n))
  }
}
}

# ---------------------------------------------------------------------------
# run_09_update_db_surv()  <-  verbatim pipeline/R/09_update_db_surv.R
# ---------------------------------------------------------------------------
run_09_update_db_surv <- function() {
# 09_update_db_surv.R --------------------------------------------------------
# 将本地(已规范化) surv rds 同步更新到远程 MySQL cpas 数据库:
#   1) 匹配表: 本地文件 <stem>_surv.rds -> 表 <gsub('-','_',stem)>_surv
#   2) 按 DB 现有 schema 对齐列(含 row_names=样本ID; 临床列名大小写/别名兼容)
#   3) 每表事务: DELETE 旧行 -> INSERT 本地全部行 (NA -> NULL)
# TCGA-* 及本地无对应表者跳过(写入日志)。
# 输出: pipeline/out/db_surv_update_log.csv + REPORT_db_surv_update.md
# ----------------------------------------------------------------------------
suppressPackageStartupMessages(library(RMySQL))
root <- "/home/Jingle/data/Project/CPAS"
sdir <- file.path(root, "data/processed/surv")
outd <- file.path(root, "pipeline/out")
files <- list.files(sdir, pattern = "_surv\\.rds$", full.names = TRUE)

con <- dbConnect(MySQL(), host = "139.224.80.159", dbname = "cpas",
                 user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"), password = Sys.getenv("CPAS_DB_PASSWORD"))
tabs <- dbListTables(con)
surv_tabs <- tabs[grepl("_surv$", tabs)]

map_local_col <- function(local_colnames, db_col) {
  if (db_col %in% local_colnames) return(db_col)
  hit <- local_colnames[tolower(local_colnames) == tolower(db_col)]
  if (length(hit)) return(hit[1])
  if (db_col == "gender" && "sex" %in% local_colnames) return("sex")
  if (db_col == "sex" && "gender" %in% local_colnames) return("gender")
  NA_character_
}

log_rows <- list(); skipped <- list()
for (f in files) {
  nm <- basename(f)
  stem <- sub("_surv\\.rds$", "", nm)
  db_acc <- gsub("-", "_", stem)
  tbl <- paste0(db_acc, "_surv")
  if (!tbl %in% surv_tabs) { skipped[[nm]] <- "no DB table"; next }
  d <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(d)) { skipped[[nm]] <- "read fail"; next }
  sch <- dbGetQuery(con, sprintf("SHOW COLUMNS FROM `%s`", tbl))
  db_cols <- sch$Field
  idvals <- rownames(d)
  if (is.null(idvals)) idvals <- as.character(d[[1]])   # 兜底: 首列当ID
  # 构建对齐后的写入数据框 (严格 DB 列顺序)
  parts <- list()
  for (c in db_cols) {
    if (c == db_cols[1] && tolower(c) %in% c("row_names", "id")) {
      parts[[c]] <- as.character(idvals); next
    }
    src <- map_local_col(colnames(d), c)
    v <- if (!is.na(src)) d[[src]] else rep(NA, nrow(d))
    tp <- sch$Type[sch$Field == c]
    if (grepl("double|float|int", tp)) {
      parts[[c]] <- suppressWarnings(as.numeric(v))
    } else {
      parts[[c]] <- if (is.numeric(v) || is.integer(v)) as.character(v) else v
    }
  }
  out <- as.data.frame(parts, check.names = FALSE, stringsAsFactors = FALSE)
  ok <- tryCatch({
    dbBegin(con)
    dbExecute(con, sprintf("DELETE FROM `%s`", tbl))
    dbWriteTable(con, tbl, out, append = TRUE, row.names = FALSE)
    dbCommit(con); TRUE
  }, error = function(e) {
    tryCatch(dbRollback(con), error = function(e2) NULL)
    FALSE
  })
  if (ok) {
    n_db <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", tbl))$n
    log_rows[[length(log_rows) + 1L]] <- data.frame(
      table = tbl, file = nm, rows_written = nrow(out), rows_db_after = n_db,
      ok = TRUE, note = "", stringsAsFactors = FALSE)
  } else {
    skipped[[nm]] <- "write failed (transaction rolled back)"
  }
}
dbDisconnect(con)

lg <- if (length(log_rows)) do.call(rbind, log_rows) else
  data.frame(table = character(0), file = character(0), rows_written = integer(0),
             rows_db_after = integer(0), ok = logical(0), note = character(0))
write.csv(lg, file.path(outd, "db_surv_update_log.csv"), row.names = FALSE)
sk <- data.frame(file = names(skipped), reason = unlist(skipped),
                 stringsAsFactors = FALSE)
write.csv(sk, file.path(outd, "db_surv_update_skipped.csv"), row.names = FALSE)

md <- c("# DB surv 表更新报告 (标准化后本地 -> MySQL cpas)",
        sprintf("日期: %s | 成功更新: %d 张表 | 跳过: %d",
                format(Sys.time()), nrow(lg), nrow(sk)), "",
        "## 成功更新", "")
for (i in seq_len(nrow(lg)))
  md <- c(md, sprintf("- %s (%s): 写入 %d 行, DB 现有 %d 行",
                      lg$table[i], lg$file[i], lg$rows_written[i],
                      lg$rows_db_after[i]))
md <- c(md, "", "## 跳过", "")
for (i in seq_len(nrow(sk)))
  md <- c(md, sprintf("- %s: %s", sk$file[i], sk$reason[i]))
writeLines(md, file.path(outd, "REPORT_db_surv_update.md"))
cat("updated tables:", nrow(lg), "| skipped:", nrow(sk), "\n")
}

# ---------------------------------------------------------------------------
# run_13_upload_small_cohorts()  <-  verbatim pipeline/R/13_upload_small_cohorts.R
# ---------------------------------------------------------------------------
run_13_upload_small_cohorts <- function() {
# 13_upload_small_cohorts.R --------------------------------------------------
# 把「有生存表但没有表达表」的 13 个队列补入 MySQL 镜像 cpas(研究者已确认放开 N>79 门槛)。
#
#   A) 平台注释表 <GPL> (row_names, gene_id) —— 补 GPL4685 / GPL10295 / GPL3694 / GPL3696
#   B) 表达表 <ACC>(row_names + 每样本一列 float) —— 13 个队列,探针先按平台注释白名单过滤
#
# 约定与既有 upload.R / gpl_app.R 一致:row_names varchar(255),样本列 float,ENGINE=MyISAM。
# 输入: data/expr/<Accession>.rds        (第一列 ID_REF,其余列为样本)
#       data/processed/gpl/<GPL>.rds     (行名 = 探针,列 gene_id = ENTREZ)
# 输出: DB 表 + pipeline/out/db_upload_small_cohorts_log.csv
# 用法: Rscript pipeline/R/13_upload_small_cohorts.R [--dry-run]
# ----------------------------------------------------------------------------
suppressPackageStartupMessages(library(RMySQL))

root <- "/home/Jingle/data/Project/CPAS"
dry  <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

targets <- data.frame(
  Accession = c("GSE31519", "GSE7378", "GSE9195", "GSE46602", "GSE84426", "GSE17537",
                "GSE39084", "GSE57495", "GSE39055", "GSE21257", "GSE12417-GPL570",
                "GSE4716-GPL3694", "GSE4716-GPL3696"),
  GPL = c("GPL96", "GPL4685", "GPL570", "GPL570", "GPL6947", "GPL570", "GPL570",
          "GPL15048", "GPL14951", "GPL10295", "GPL570", "GPL3694", "GPL3696"),
  stringsAsFactors = FALSE)
new_gpl <- c("GPL4685", "GPL10295", "GPL3694", "GPL3696")

con <- dbConnect(MySQL(), host = "139.224.80.159", dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = Sys.getenv("CPAS_DB_PASSWORD"), client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)
existing <- dbListTables(con)

log_rows <- list()
add_log <- function(step, table, file, rows, probes, samples, status, note = "") {
  log_rows[[length(log_rows) + 1L]] <<- data.frame(
    step = step, table = table, file = file, rows_written = rows,
    probes = probes, samples = samples, status = status, note = note,
    stringsAsFactors = FALSE)
}

# ---- A) platform annotation tables ----------------------------------------
for (g in new_gpl) {
  f <- file.path(root, "data/processed/gpl", paste0(g, ".rds"))
  if (!file.exists(f)) { add_log("gpl", g, basename(f), 0, 0, NA, "missing local rds"); next }
  d <- readRDS(f)
  d <- d[!is.na(d$gene_id), , drop = FALSE]
  out <- data.frame(row_names = rownames(d), gene_id = as.character(d$gene_id),
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (dry) { add_log("gpl", g, basename(f), nrow(out), nrow(out), NA, "dry-run"); next }
  ok <- tryCatch({
    dbWriteTable(con, g, out, overwrite = TRUE, row.names = FALSE,
                 field.types = c(row_names = "TEXT", gene_id = "TEXT"))
    dbExecute(con, sprintf("ALTER TABLE `%s` ENGINE=MyISAM", g))
    TRUE
  }, error = function(e) conditionMessage(e))
  add_log("gpl", g, basename(f), if (isTRUE(ok)) nrow(out) else 0, nrow(out), NA,
          if (isTRUE(ok)) "ok" else "failed", if (isTRUE(ok)) "" else as.character(ok))
}

# ---- B) expression tables --------------------------------------------------
for (i in seq_len(nrow(targets))) {
  acc <- targets$Accession[i]; gpl <- targets$GPL[i]
  tbl <- gsub("-", "_", acc, fixed = TRUE)
  f <- file.path(root, "data/expr", paste0(acc, ".rds"))
  if (!file.exists(f)) { add_log("expr", tbl, basename(f), 0, 0, NA, "missing local rds"); next }
  gf <- file.path(root, "data/processed/gpl", paste0(gpl, ".rds"))
  if (!file.exists(gf)) { add_log("expr", tbl, basename(f), 0, 0, NA, "missing gpl rds"); next }
  dd <- readRDS(f)
  gpl_map <- readRDS(gf)
  n0 <- nrow(dd)
  keep <- as.character(dd[[1]]) %in% rownames(gpl_map)
  dd <- dd[keep, , drop = FALSE]
  ids <- as.character(dd[[1]])
  dd <- dd[, -1, drop = FALSE]
  dd[] <- lapply(dd, function(x) suppressWarnings(as.numeric(x)))
  dd <- dd[rowSums(!is.na(dd)) > 0, , drop = FALSE]
  ids <- ids[rowSums(!is.na(dd)) > 0]
  out <- as.data.frame(dd, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- ids
  ft <- rep("FLOAT", ncol(out)); names(ft) <- colnames(out)
  ft <- c(row_names = "VARCHAR(255)", ft)
  if (dry) { add_log("expr", tbl, basename(f), nrow(out), n0, ncol(out), "dry-run"); next }
  ok <- tryCatch({
    dbWriteTable(con, tbl, out, overwrite = TRUE, row.names = TRUE, field.types = ft)
    dbExecute(con, sprintf("ALTER TABLE `%s` ENGINE=MyISAM", tbl))
    TRUE
  }, error = function(e) conditionMessage(e))
  if (isTRUE(ok)) {
    n_db <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", tbl))$n
    add_log("expr", tbl, basename(f), nrow(out), n0, ncol(out), "ok",
            sprintf("db_rows=%d", n_db))
  } else {
    add_log("expr", tbl, basename(f), 0, n0, ncol(out), "failed", as.character(ok))
  }
  message(sprintf("[%d/%d] %s -> %s (%d probes x %d samples)", i, nrow(targets), acc, tbl,
                  nrow(out), ncol(out)))
}

res <- do.call(rbind, log_rows)
out_csv <- file.path(root, "pipeline/out/db_upload_small_cohorts_log.csv")
utils::write.csv(res, out_csv, row.names = FALSE)
cat("\n=== upload summary ===\n"); print(res, row.names = FALSE)
cat("\nlog:", out_csv, "\n")
}

# ---------------------------------------------------------------------------
# run_24_make_tables_writable()  <-  verbatim pipeline/R/24_make_tables_writable.R
# ---------------------------------------------------------------------------
run_24_make_tables_writable <- function() {
# 24_make_tables_writable.R --------------------------------------------------
# 把「只读」的镜像表换成可写表（不传输数据，不改变表名与内容）。
#
# 背景
#   镜像里有一部分 MyISAM 表的数据文件对 mysqld 进程不可写，任何写入都报
#   `Table 'X' is read only`（MySQL 1036）。实测 293 张表中 10 张属于这种情况（见 pipeline/out/REPORT_tables_writable.md），
#   全部是肺队列表达表（GSE11969/13213/14814/17710/30219/37745/41271/42127/
#   50081/74777/8894），其余 282 张（含全部 *_surv）可写。它们的 CREATE_TIME
#   集中在 2025-02-18 00:10–00:26 的批量导入窗口，说明是该批文件的属主/权限问题，
#   而不是账号权限（本账号对 cpas.* 有 ALL PRIVILEGES，且新建表可读可写）。
#
# 修法（客户端可完成，无需服务器 shell）
#   1) CREATE TABLE <tmp> LIKE <t>;  INSERT INTO <tmp> SELECT * FROM <t>;
#      —— 服务端内部复制，数据不出服务器，耗时与表大小无关（无网络传输）
#   2) CHECKSUM TABLE 比对原表与新表，不一致就放弃
#   3) RENAME TABLE <t> TO <old>, <tmp> TO <t>   —— 原子操作，表名不变
#   4) 再次 CHECKSUM 校验；不一致则换回
#   5) 删除 <old>，并验证 <t> 现在可写
#
# 用法
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/24_make_tables_writable.R          # dry-run
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/24_make_tables_writable.R --apply
#   （--only <表名> 可只处理一张）
# 输出：pipeline/out/REPORT_tables_writable.md
# ----------------------------------------------------------------------------
.libPaths(c("/home/Jingle/R/library", .libPaths()))
suppressPackageStartupMessages(library(RMySQL))

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
apply_db <- "--apply" %in% args
only <- if ("--only" %in% args) args[which(args == "--only") + 1L] else NULL

pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
if (!nzchar(pw)) stop("CPAS_DB_PASSWORD 未设置。")
con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = pw, client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)

log_lines <- character(0)
note <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
try_sql <- function(sql) tryCatch({ dbExecute(con, sql); TRUE },
                                 error = function(e) paste("ERR:", conditionMessage(e)))

is_readonly <- function(t) {
  c1 <- dbListFields(con, t)[1]
  r <- tryCatch({ dbExecute(con, sprintf("UPDATE `%s` SET `%s`=`%s` WHERE 1=0", t, c1, c1)); "writable" },
                error = function(e) if (grepl("read only", conditionMessage(e))) "readonly" else "error")
  r
}
checksum <- function(t) {
  r <- dbGetQuery(con, sprintf("CHECKSUM TABLE `%s`", t))
  as.character(r$Checksum[1])
}

tb <- dbGetQuery(con, "SHOW TABLES")[[1]]
targets <- character(0)
for (t in tb) if (is_readonly(t) == "readonly") targets <- c(targets, t)
if (!is.null(only)) targets <- intersect(targets, only)
note("## 只读表检测")
note(paste0("  扫描 ", length(tb), " 张表，只读 ", length(targets), " 张：",
            if (length(targets)) paste(targets, collapse = ", ") else "（无）"))

if (!apply_db || !length(targets)) {
  dir.create("pipeline/out", showWarnings = FALSE)
  writeLines(c("# 只读表检测 / 修复", "",
               paste0("- 时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               paste0("- 模式: ", if (apply_db) "--apply" else "dry-run"),
               "", "```text", log_lines, "```", ""),
             "pipeline/out/REPORT_tables_writable.md")
  if (!apply_db) note("\n[dry-run] 未做修改。加 --apply 执行换表。")
  quit(save = "no")
}

note("\n## 换表（内容不变，改为可写）")
for (t in targets) {
  tmp <- paste0("_cpas_swap_", t); old <- paste0("_cpas_old_", t)
  before <- checksum(t)
  ok <- TRUE
  for (sql in c(sprintf("DROP TABLE IF EXISTS `%s`", tmp),
                sprintf("CREATE TABLE `%s` LIKE `%s`", tmp, t),
                sprintf("INSERT INTO `%s` SELECT * FROM `%s`", tmp, t))) {
    r <- try_sql(sql)
    if (!isTRUE(r)) { note(sprintf("  %-16s 失败于 %s -> %s", t, substr(sql, 1, 40), r)); ok <- FALSE; break }
  }
  if (!ok) { try_sql(sprintf("DROP TABLE IF EXISTS `%s`", tmp)); next }
  cs_tmp <- checksum(tmp)
  if (!identical(cs_tmp, before)) {
    note(sprintf("  %-16s CHECKSUM 不一致（原 %s / 新 %s），放弃", t, before, cs_tmp))
    try_sql(sprintf("DROP TABLE IF EXISTS `%s`", tmp)); next
  }
  # 原子换名
  r <- try_sql(sprintf("RENAME TABLE `%s` TO `%s`, `%s` TO `%s`", t, old, tmp, t))
  if (!isTRUE(r)) {
    note(sprintf("  %-16s RENAME 失败：%s（已丢弃副本）", t, r))
    try_sql(sprintf("DROP TABLE IF EXISTS `%s`", tmp)); try_sql(sprintf("DROP TABLE IF EXISTS `%s`", old)); next
  }
  after <- checksum(t)
  if (!identical(after, before)) {
    note(sprintf("  %-16s 换名后 CHECKSUM 变化（%s -> %s），回滚", t, before, after))
    try_sql(sprintf("RENAME TABLE `%s` TO `%s`, `%s` TO `%s`", t, tmp, old, t)); next
  }
  writable <- is_readonly(t) == "writable"
  try_sql(sprintf("DROP TABLE IF EXISTS `%s`", old))
  note(sprintf("  %-16s 完成：CHECKSUM %s 不变 | 现在可写=%s", t, before, writable))
}

dir.create("pipeline/out", showWarnings = FALSE)
writeLines(c("# 只读表检测 / 修复", "",
             paste0("- 时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             paste0("- 模式: ", if (apply_db) "--apply" else "dry-run"),
             "", "```text", log_lines, "```", ""),
           "pipeline/out/REPORT_tables_writable.md")
note("\n已写出 pipeline/out/REPORT_tables_writable.md")
}

# ---------------------------------------------------------------------------
# run_42_upload_gpl1223()  <-  verbatim pipeline/R/42_upload_gpl1223.R
# ---------------------------------------------------------------------------
run_42_upload_gpl1223 <- function() {
# 42_upload_gpl1223.R --------------------------------------------------------
# 目的
#   把 GPL1223（Agilent-014850 Whole Human Genome 4x44K，GSE1379 的平台）的注释表
#   补进 MySQL 镜像 cpas。
#
# 为什么
#   06_upload_db.R 只在 n_expr > 79 时才补平台注释表，GSE1379 的 n_expr = 60，
#   所以按默认闸门上注释表被跳过；而镜像里原本**完全没有** GPL1223 表。
#   GSE1379 的表达表已按 `ID_REF %in% rownames(GPL1223)` 过滤后上传（11863 探针 ×
#   60 样本），没有注释表这一列就只是探针号 —— 其余每一个已编目的 GEO 队列都有平台
#   注释表。这里用与 13_upload_small_cohorts.R 完全相同的写法补上这**唯一**缺失的一张。
#
# 安全保证
#   本脚本**只在该表不存在时**创建（overwrite 仅作用于新建表）；镜像里已存在则直接跳过，
#   绝不用本地 rds 覆盖内容可能不同的既有注释表。
#
# 用法: CPAS_DB_PASSWORD=... Rscript pipeline/R/42_upload_gpl1223.R [ROOT]
# ----------------------------------------------------------------------------
suppressPackageStartupMessages(library(RMySQL))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
GPL <- "GPL1223"
f <- file.path(root, "data/processed/gpl", paste0(GPL, ".rds"))
stopifnot(file.exists(f))

con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = Sys.getenv("CPAS_DB_PASSWORD"), client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)

if (GPL %in% dbListTables(con)) {
  n <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", GPL))$n
  cat(sprintf("[skip] %s already in the mirror (%d rows) - never overwritten\n", GPL, n))
  quit(save = "no")
}

d <- readRDS(f)
d <- d[!is.na(d$gene_id) & d$gene_id != "", , drop = FALSE]
out <- data.frame(row_names = rownames(d), gene_id = as.character(d$gene_id),
                  check.names = FALSE, stringsAsFactors = FALSE)
dbWriteTable(con, GPL, out, overwrite = TRUE, row.names = FALSE,
             field.types = c(row_names = "TEXT", gene_id = "TEXT"))
dbExecute(con, sprintf("ALTER TABLE `%s` ENGINE=MyISAM", GPL))
n <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", GPL))$n
cat(sprintf("[uploaded] %s (new table): %d probe-gene rows\n", GPL, n))
}

# ---------------------------------------------------------------------------
# run_99_extend_gpl_db()  <-  verbatim pipeline/R/99_extend_gpl_db.R
# ---------------------------------------------------------------------------
run_99_extend_gpl_db <- function() {
# 99_extend_gpl_db.R ----------------------------------------------------------
# 把本地 GPL 注释表的**新增映射**追加进镜像的 <GPL> 表（只追加、不覆盖、不删行）。
#
# 背景：06_upload_db.R 只在表**不存在**时上传 GPL 表，所以 96 号脚本把新队列的
# 基因符号并入本地 data/processed/gpl/GPL24676.rds 之后，镜像里的 GPL24676 仍是
# 旧快照 —— 结果是新队列 A5-PCPG / IMmotion150 的表达矩阵在 App 的
# get_expr_data() 路径（ID_map -> <GPL> 表取 row_names -> 表达表按 row_names 取数）
# 里只有 ~21-30% 的行可查到（既有 GSE271517 也只有 40.6%）。本脚本把本地新增的
# symbol->ENTREZ 行补进镜像，使新队列可查。
#
# 安全性：只 INSERT 镜像中不存在的 row_names；既有的行一个都不动、不覆盖、不删除。
# 因此对既有队列的影响仅限于「某个 gene_id 现在多了一个探针行」——而该探针若不在
# 既有队列的表达表里，get_expr_data() 的 merge 会自动丢掉它。
#
# 用法:
#   Rscript pipeline/R/99_extend_gpl_db.R [ROOT] [GPL ...]            # dry-run
#   Rscript pipeline/R/99_extend_gpl_db.R [ROOT] --write [GPL ...]    # 追加
# 输出: pipeline/out/gpl_db_extend_log.csv
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(RMySQL); library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else
  Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
ROOT <- path.expand(ROOT); setwd(ROOT)
WRITE <- "--write" %in% args
GPLS <- setdiff(args[!startsWith(args, "--")], ROOT)
if (!length(GPLS)) GPLS <- c("GPL24676")

con <- dbConnect(MySQL(), host = "139.224.80.159", dbname = "cpas",
                 user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = Sys.getenv("CPAS_DB_PASSWORD"), client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con))

log_rows <- list()
for (gpl in GPLS) {
  tab <- str_replace_all(gpl, "-", "_")
  f <- file.path(ROOT, "data/processed/gpl", paste0(gpl, ".rds"))
  if (!file.exists(f)) { message("skip ", gpl, ": no local map"); next }
  if (!tab %in% dbListTables(con)) { message("skip ", gpl, ": table not in DB"); next }
  g <- readRDS(f)
  g2 <- data.frame(row_names = rownames(g), gene_id = as.character(g$gene_id),
                   stringsAsFactors = FALSE)
  g2 <- g2[!is.na(g2$gene_id) & g2$gene_id != "---" & g2$gene_id != "", , drop = FALSE]
  g2 <- do.call(rbind, lapply(seq_len(nrow(g2)), function(i)
    data.frame(row_names = g2$row_names[i],
               gene_id = strsplit(g2$gene_id[i], " /// ", fixed = TRUE)[[1]],
               stringsAsFactors = FALSE)))
  db_rn <- dbGetQuery(con, sprintf("SELECT row_names FROM %s", tab))$row_names
  add <- g2[!g2$row_names %in% db_rn, , drop = FALSE]
  cat(sprintf("%s: local mapped probe-gene rows=%d | DB rows=%d | to append=%d (distinct new row_names=%d)\n",
              gpl, nrow(g2), length(db_rn), nrow(add), length(unique(add$row_names))))
  if (nrow(add) && WRITE) {
    dbWriteTable(con, name = tab, value = add, row.names = FALSE, append = TRUE)
    n_after <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM %s", tab))$n
    cat(sprintf("  appended -> DB rows now %d\n", n_after))
  } else if (nrow(add)) cat("  [dry-run] 未写入\n")
  log_rows[[length(log_rows) + 1L]] <- data.frame(
    GPL = gpl, local_mapped_rows = nrow(g2), db_rows_before = length(db_rn),
    rows_appended = if (WRITE) nrow(add) else 0L,
    distinct_new_row_names = length(unique(add$row_names)),
    wrote = WRITE, stringsAsFactors = FALSE)
}
out <- do.call(rbind, log_rows)
if (!is.null(out)) {
  write.csv(out, file.path(ROOT, "pipeline/out/gpl_db_extend_log.csv"), row.names = FALSE)
  cat("\n=== 99 GPL DB extend summary ===\n"); print(out, row.names = FALSE)
}
}

# ---------------------------------------------------------------------------
# run_112_upload_gpl27956()  <-  verbatim pipeline/R/112_upload_gpl27956.R
# ---------------------------------------------------------------------------
run_112_upload_gpl27956 <- function() {
# 112_upload_gpl27956.R -------------------------------------------------------
# 目的
#   把 GPL27956（NanoString nCounter PanCancer IO 360，GSE205209 的平台）的注释表
#   补进 MySQL 镜像 cpas。
#
# 为什么
#   06_upload_db.R 只在 n_expr > 79 时才补平台注释表，GSE205209 的 n_expr = 54，
#   所以按默认闸门上注释表被跳过；而镜像里原本**完全没有** GPL27956 表
#   （本次上传前实测：API action=gpl&table=GPL27956 -> HTTP 500 "Query error"）。
#   GSE205209 的表达表已在镜像里（770 基因 x 54 样本，row_names = 基因符号），
#   没有注释表这一列，`ID_map -> GPL27956 -> GSE205209` 的基因查询路径就走不通
#   （core.R 的 .ref_probes() 先查 GPL 表拿探针，再查表达表 —— 与 E-MTAB-6877
#   上传 GPL16686 修的同一个问题）。
#
# 本地映射来源
#   data/processed/gpl/GPL27956.rds（110_build_gse205209.R 生成）：
#   平台 SOFT 只有 symbol（ID/ORF/SPOT_ID，无 Entrez），gene_id 由
#   org.Hs.eg.db SYMBOL->ENTREZID 得到；770 行 / 756 行有 Entrez = 98.18%。
#   本平台的 row_names 即基因符号（NanoString 探针就是基因），这与其它平台
#   的探针号写法不同但完全一致于镜像约定：row_names = 平台侧 ID，gene_id = Entrez。
#
# 安全保证
#   本脚本**只在该表不存在时**创建（overwrite 仅作用于新建表）；镜像里已存在则直接跳过，
#   绝不用本地 rds 覆盖内容可能不同的既有注释表。
#
# 用法: Rscript pipeline/R/112_upload_gpl27956.R [ROOT]
# ----------------------------------------------------------------------------
suppressPackageStartupMessages(library(RMySQL))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
GPL <- "GPL27956"
f <- file.path(root, "data/processed/gpl", paste0(GPL, ".rds"))
stopifnot(file.exists(f))

con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = Sys.getenv("CPAS_DB_PASSWORD"), client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)

if (GPL %in% dbListTables(con)) {
  n <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", GPL))$n
  cat(sprintf("[skip] %s already in the mirror (%d rows) - never overwritten\n", GPL, n))
  quit(save = "no")
}

d <- readRDS(f)
cat(sprintf("[local] %s: %d rows, %d with non-NA/non-empty gene_id\n",
            GPL, nrow(d), sum(!is.na(d$gene_id) & d$gene_id != "")))
d <- d[!is.na(d$gene_id) & d$gene_id != "", , drop = FALSE]
out <- data.frame(row_names = rownames(d), gene_id = as.character(d$gene_id),
                  check.names = FALSE, stringsAsFactors = FALSE)
dbWriteTable(con, GPL, out, overwrite = TRUE, row.names = FALSE,
             field.types = c(row_names = "TEXT", gene_id = "TEXT"))
dbExecute(con, sprintf("ALTER TABLE `%s` ENGINE=MyISAM", GPL))
n <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", GPL))$n
cat(sprintf("[uploaded] %s (new table): %d probe-gene rows (rows written: %d)\n",
            GPL, n, nrow(out)))
}
