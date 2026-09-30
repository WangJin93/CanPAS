# ===========================================================================
# CanPAS curation pipeline -- 05_clinical_and_scale
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 07_standardize_clinical.R
#   2. 07b_split_tnm.R
#   3. 07c_T_usage_table.R
#   4. 08_verify_normalization.R
#   5. 11_endpoint_families.R
#   6. 18_log2_transform.R
#   7. 20_log2_transform_db.R
#   8. 21_reupload_log2_tables.R
#   9. 22_report_log2_impact.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "05_clinical_and_scale.R", package = "CanPAS"))
#   run_07_standardize_clinical()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================

# --- shared helpers (glue) --------------------------------------------------
# The original 07b_split_tnm.R re-read <root>/pipeline/R/07_standardize_clinical.R
# and eval()'d its definition region to reuse the normalizers.  The consolidated
# file already carries that script verbatim, so the same definitions (MISS_TOKENS,
# is_missing_token, strip_missing, clean_text, norm_tnm) are exposed at file level,
# copied verbatim from 07_standardize_clinical.R, and its 4-line re-read preamble
# is replaced by a pointer comment.  The function bodies are unchanged.
# ----------------------------------------------------------------------------
MISS_TOKENS <- c("", "NA", "N/A", "N.A.", "#NA", "#N/A", "NULL", "NaN",
                 "-", "--", "—", "–", ".", "..", "?", "??",
                 "UNKNOWN", "UNK", "99 UNK", "NOT AVAILABLE", "NOT REPORTED",
                 "NOT SPECIFIED", "NOT APPLICABLE", "MISSING", "NONE",
                 "NON-APPLICABLE", "[UNKNOWN]", "[DISCREPANCY]",
                 "[NOT AVAILABLE]", "<NA>", "NAN")
is_missing_token <- function(x) {
  t <- trimws(x)
  tolower(t) %in% tolower(MISS_TOKENS) | grepl("^[#\\s]*n/?a[#\\s]*$", t, ignore.case = TRUE)
}
strip_missing <- function(x) {
  if (is.numeric(x)) { x[is.na(x)] <- NA; return(x) }
  x <- as.character(x)
  x[is_missing_token(x)] <- NA_character_
  x
}
clean_text <- function(x) {
  x <- strip_missing(x)
  x <- gsub("[[:space:]]+", " ", trimws(x))
  x[trimws(x) == ""] <- NA_character_
  x
}
norm_tnm <- function(x, prefix) {
  v <- clean_text(x)
  out <- rep(NA_character_, length(v))
  for (i in seq_along(v)) {
    s <- v[i]
    if (is.na(s)) next
    s <- toupper(trimws(s))
    s <- sub("^[PC](?=[TNM])", "", s, perl = TRUE)   # 去 p/c 前缀 (pT3->T3, pN0->N0)
    m <- character(0)
    if (prefix == "T") {
      m <- regmatches(s, regexpr("^(?:[PC])?T(?:IS|[0-4X])([ABC])?$", s))
      if (!length(m) && grepl("^[0-4X]$", s)) m <- paste0("T", s)
    } else if (prefix == "N") {
      m <- regmatches(s, regexpr("^N[0-3X]([ABC])?$", s))
      if (!length(m) && grepl("^[0-3X]$", s)) m <- paste0("N", s)
    } else {
      m <- regmatches(s, regexpr("^M[01X]([ABC])?$", s))
      if (!length(m) && grepl("^[01X]$", s)) m <- paste0("M", s)
    }
    if (length(m)) {
      val <- sub("([ABC])$", "", m[1])       # 去掉 a/b/c
      val <- sub("^[PC](?=T)", "", val, perl = TRUE)   # pT3 -> T3, cT1 -> T1
      val <- sub("^TIS$", "Tis", val)        # 规范 Tis
      out[i] <- val
    }
  }
  out
}


# ---------------------------------------------------------------------------
# run_07_standardize_clinical()  <-  verbatim pipeline/R/07_standardize_clinical.R
# ---------------------------------------------------------------------------
run_07_standardize_clinical <- function() {
# 07_standardize_clinical.R -------------------------------------------------
# 逐个核对并统一 data/processed/surv/ 各队列的临床信息：
#   1) 缺失值记录(NA/""/#NA/unknown/--/?/… ) -> NA(空)
#   2) stage/T/N/M 简化: 去掉 a/b/c 亚级(Ia->I, IIB->II, T1a->T1…), 统一大小写
#   3) sex -> male/female; age -> numeric; grade -> G1..G4 (well/moderate/poor并入)
#   4) 其它临床文本列: 去空白、缺失标记 -> NA
# 执行前已备份至 pipeline/backup/surv_pre_clinical_std/；结果回写 *_surv.rds
# 并在 pipeline/out 生成 变更清单 + 核对报告。
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

root <- "/home/Jingle/data/Project/CPAS"
surv_dir <- file.path(root, "data/processed/surv")
out_dir  <- cpas_out_root()
bdir <- file.path(root, "pipeline/backup", "surv_pre_clinical_std")
files <- list.files(surv_dir, pattern = "_surv\\.rds$", full.names = TRUE)

# ---------------- 缺失值 ----------------
MISS_TOKENS <- c("", "NA", "N/A", "N.A.", "#NA", "#N/A", "NULL", "NaN",
                 "-", "--", "—", "–", ".", "..", "?", "??",
                 "UNKNOWN", "UNK", "99 UNK", "NOT AVAILABLE", "NOT REPORTED",
                 "NOT SPECIFIED", "NOT APPLICABLE", "MISSING", "NONE",
                 "NON-APPLICABLE", "[UNKNOWN]", "[DISCREPANCY]",
                 "[NOT AVAILABLE]", "<NA>", "NAN")
is_missing_token <- function(x) {
  t <- trimws(x)
  tolower(t) %in% tolower(MISS_TOKENS) | grepl("^[#\\s]*n/?a[#\\s]*$", t, ignore.case = TRUE)
}
strip_missing <- function(x) {
  if (is.numeric(x)) { x[is.na(x)] <- NA; return(x) }
  x <- as.character(x)
  x[is_missing_token(x)] <- NA_character_
  x
}
clean_text <- function(x) {
  x <- strip_missing(x)
  x <- gsub("[[:space:]]+", " ", trimws(x))
  x[trimws(x) == ""] <- NA_character_
  x
}

# ---------------- stage ----------------
ROMAN <- c("0" = "0", "1" = "I", "2" = "II", "3" = "III", "4" = "IV")
norm_stage <- function(x) {
  v <- clean_text(x)
  out <- rep(NA_character_, length(v))
  for (i in seq_along(v)) {
    s <- v[i]
    if (is.na(s)) next
    s <- toupper(trimws(s))
    s <- sub("^STAGE\\s*", "", s)                   # 先大写再去 "STAGE" 前缀
    if (grepl("[TNM][0-9xX]", s) || grepl("[+ ]", s)) next      # TNM组合/描述串
    if (s %in% c("X", "TA", "TIS", "PP", "99", "UNK", "UNKNOWN")) next
    if (grepl("^[0-4][A-C]?$", s)) {                            # 数字式 1A/2/3b…
      out[i] <- unname(ROMAN[substr(s, 1, 1)]); next
    }
    m <- regmatches(s, regexpr("^IV|^III|^II|^I", s))           # 罗马式
    if (length(m)) {
      r <- m[1]; rest <- substring(s, nchar(r) + 1L)
      if (r %in% c("I", "II", "III", "IV") && rest %in% c("", "A", "B", "C"))
        out[i] <- r
    }
  }
  out
}

# ---------------- T / N / M ----------------
norm_tnm <- function(x, prefix) {
  v <- clean_text(x)
  out <- rep(NA_character_, length(v))
  for (i in seq_along(v)) {
    s <- v[i]
    if (is.na(s)) next
    s <- toupper(trimws(s))
    s <- sub("^[PC](?=[TNM])", "", s, perl = TRUE)   # 去 p/c 前缀 (pT3->T3, pN0->N0)
    m <- character(0)
    if (prefix == "T") {
      m <- regmatches(s, regexpr("^(?:[PC])?T(?:IS|[0-4X])([ABC])?$", s))
      if (!length(m) && grepl("^[0-4X]$", s)) m <- paste0("T", s)
    } else if (prefix == "N") {
      m <- regmatches(s, regexpr("^N[0-3X]([ABC])?$", s))
      if (!length(m) && grepl("^[0-3X]$", s)) m <- paste0("N", s)
    } else {
      m <- regmatches(s, regexpr("^M[01X]([ABC])?$", s))
      if (!length(m) && grepl("^[01X]$", s)) m <- paste0("M", s)
    }
    if (length(m)) {
      val <- sub("([ABC])$", "", m[1])       # 去掉 a/b/c
      val <- sub("^[PC](?=T)", "", val, perl = TRUE)   # pT3 -> T3, cT1 -> T1
      val <- sub("^TIS$", "Tis", val)        # 规范 Tis
      out[i] <- val
    }
  }
  out
}

# ---------------- sex ----------------
SEX_MAP <- c(male = "male", m = "male", man = "male", men = "male",
             female = "female", f = "female", woman = "female",
             women = "female")
norm_sex <- function(x) {
  v <- clean_text(x)
  out <- rep(NA_character_, length(v))
  for (i in seq_along(v)) {
    k <- tolower(trimws(v[i]))
    if (!is.na(k) && k %in% names(SEX_MAP)) out[i] <- unname(SEX_MAP[k])
  }
  out
}

# ---------------- grade ----------------
GRADE_WORD <- c(
  "1" = "G1", "2" = "G2", "3" = "G3", "4" = "G4",
  "I" = "G1", "II" = "G2", "III" = "G3", "IV" = "G4",
  "G1" = "G1", "G2" = "G2", "G3" = "G3", "G4" = "G4",
  "WELL" = "G1", "WELL DIFFERENTIATED" = "G1",
  "MODERATELY" = "G2", "MOD" = "G2", "MODERATE" = "G2",
  "MODERATELY DIFFERENTIATED" = "G2",
  "POORLY" = "G3", "POOR" = "G3", "POORLY DIFFERENTIATED" = "G3",
  "HIGH" = "G3", "HIGH GRADE" = "G3", "HIGH.GRADE" = "G3",
  "LOW" = "G2", "LOW GRADE" = "G2", "LOW.GRADE" = "G2",
  "UNDIFFERENTIATED" = "G4")
norm_grade <- function(x) {
  v <- clean_text(x)
  out <- rep(NA_character_, length(v))
  for (i in seq_along(v)) {
    s <- v[i]
    if (is.na(s)) next
    k <- gsub("[[:space:]]+", " ", toupper(trimws(s)))
    if (k %in% names(GRADE_WORD)) out[i] <- unname(GRADE_WORD[k])
    else if (grepl("^G[0-9]$", k)) out[i] <- k                  # G1..G9 原样
    else if (grepl("^[5-9]$|^10$", k)) out[i] <- k              # 评分型(6-10)保留
    # GX/G?/GB/… -> NA
  }
  out
}

norm_age <- function(x) suppressWarnings(as.numeric(clean_text(x)))

# ---------------- 执行 ----------------
log_rows <- list()
stage_unmapped <- list()
report <- c("# 临床信息标准化核对报告 (data/processed/surv)",
            sprintf("日期: %s | 文件数: %d", format(Sys.time()), length(files)), "",
            "## 统一规则",
            "- 缺失值记录(NA/空/#NA/unknown/--/? 等) -> NA(空)。",
            "- stage/T/N/M: 去掉 a/b/c 亚级并统一(Ia->I, IIB->II, T1a->T1, N2b->N2, M1a->M1)。",
            "  stage 中 TNM 组合串/描述串(如 pT1N0M0、cTa)无法归入临床分期 -> NA。",
            "- sex -> male/female; age -> 数值; grade -> G1..G4",
            "  (well->G1, moderate/mod->G2, poorly/poor/high grade->G3, low grade->G2);",
            "  评分型数值(5-10)保留; GX/G?/GB -> NA。",
            "- 其它临床列: 去空白、统一缺失标记 -> NA(数值列转数值)。", "")

alias <- c("Age" = "age", "Stage" = "stage", "Grade" = "grade", "gender" = "sex")

for (f in files) {
  nm <- basename(f)
  d <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(d)) { report <- c(report, sprintf("- %s: 读取失败", nm)); next }
  # 同义词列归一到规范名（规范名不存在时）
  for (a in names(alias)) {
    if (a %in% colnames(d) && !(alias[a] %in% colnames(d)))
      colnames(d)[colnames(d) == a] <- alias[a]
  }
  cn <- colnames(d)
  ep_like <- grepl("_(time|status)$", cn) &
    (sub("_(time|status)$", "", cn) %in%
       c("OS", "DSS", "DFS", "RFS", "PFS", "MFS", "DRFS", "EFS", "DFI", "PFI"))
  clin_cols <- cn[!ep_like & !cn %in% c("ID", "row_names")]
  note <- character(0)
  for (col in clin_cols) {
    x <- d[[col]]; before <- x
    y <- if (col == "age")      norm_age(x)
         else if (col == "sex") norm_sex(x)
         else if (col == "stage") norm_stage(x)
         else if (col %in% c("T", "N", "M")) norm_tnm(x, col)
         else if (col == "grade") norm_grade(x)
         else {
           y2 <- clean_text(x)
           if (is.numeric(x)) y2 <- suppressWarnings(as.numeric(y2))
           y2
         }
    n_miss_old <- sum(is.na(before)) + sum(is_missing_token(as.character(before)), na.rm = TRUE)
    n_miss_new <- sum(is.na(y))
    changed <- sum(!is.na(before) & !is.na(y) & as.character(before) != as.character(y),
                   na.rm = TRUE)
    if (n_miss_new != n_miss_old || changed > 0) {
      log_rows[[length(log_rows) + 1L]] <- data.frame(
        file = nm, column = col, n = nrow(d),
        missing_old = n_miss_old, missing_new = n_miss_new,
        changed_values = changed, stringsAsFactors = FALSE)
      note <- c(note, sprintf("%s(miss %d->%d, val %d)", col,
                              n_miss_old, n_miss_new, changed))
    }
    d[[col]] <- y
  }
  saveRDS(d, f)
  # stage 未归类记录(以备份中的原始值对照): TNM组合/描述串等无法归入分期组别
  bfile <- file.path(bdir, nm)
  if (file.exists(bfile) && "stage" %in% colnames(d)) {
    od <- readRDS(bfile)
    if ("stage" %in% colnames(od)) {
      o <- as.character(od[["stage"]]); n <- d[["stage"]]
      bad <- !is.na(o) & is.na(n)
      if (sum(bad)) stage_unmapped[[nm]] <- list(k = sum(bad),
                                                 tok = sort(unique(o[bad])))
    }
  }
  report <- c(report, sprintf("- %s: OK%s", nm,
                              if (length(note)) paste0(" [", paste(note, collapse = "; "), "]") else ""))
}

# ---------------- 汇总输出 ----------------
if (length(log_rows)) {
  logdf <- do.call(rbind, log_rows)
  write.csv(logdf, file.path(out_dir, "clinical_std_changes.csv"), row.names = FALSE)
  agg <- aggregate(cbind(missing_old, missing_new, changed_values) ~ column,
                   data = logdf, FUN = sum)
  report <- c(report, "", "## 按列汇总", "")
  for (i in seq_len(nrow(agg)))
    report <- c(report, sprintf("- %s: 缺失 %d -> %d, 值变更 %d 行",
                                agg$column[i], agg$missing_old[i],
                                agg$missing_new[i], agg$changed_values[i]))
} else report <- c(report, "", "未检测到需要变更的内容。")

# stage 未归类清单（保留供人工复核）
if (length(stage_unmapped)) {
  sdf <- do.call(rbind, lapply(names(stage_unmapped), function(nm)
    data.frame(file = nm, n_rows = stage_unmapped[[nm]]$k,
               examples = paste(head(stage_unmapped[[nm]]$tok, 6), collapse = "; "),
               stringsAsFactors = FALSE)))
  sdf <- sdf[order(-sdf$n_rows), ]
  write.csv(sdf, file.path(out_dir, "clinical_std_stage_unmapped.csv"), row.names = FALSE)
  report <- c(report, "", "## 需要人工复核：stage 列中无法归入临床分期的记录(置为 NA)",
              "这些队列的 stage 字段实际为 T/TNM 描述(如 pT2 N0Mx)或未知标记，",
              "并非分期组别；已统一置空，以下文件请按需另行处理：", "")
  for (i in seq_len(nrow(sdf)))
    report <- c(report, sprintf("- %s: %d 行，如 %s", sdf$file[i], sdf$n_rows[i],
                                sdf$examples[i]))
}

writeLines(report, file.path(out_dir, "REPORT_clinical_standardization.md"))
cat("DONE files:", length(files), "| changed-col rows:", length(log_rows), "\n")
}

# ---------------------------------------------------------------------------
# run_07b_split_tnm()  <-  verbatim pipeline/R/07b_split_tnm.R
# ---------------------------------------------------------------------------
run_07b_split_tnm <- function() {
# 07b_split_tnm.R ------------------------------------------------------------
# 对 stage 列仅含 T/TNM 描述(如 pT2 N0Mx, T1N0M0, pN0pT2)且没有独立 T/N/M 列的队列，
# 从原始值中解析出标准 T/N/M 组分并写为新列(与 norm_tnm 同一套简化规则)。
# stage 列本身仍保持统一后的临床分期(0/I/II/III/IV 或 NA)。
# 已备份; 结果回写 data/processed/surv/*_surv.rds。
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

root <- "/home/Jingle/data/Project/CPAS"
sdir <- file.path(root, "data/processed/surv")
out_dir <- cpas_out_root()
bdir <- file.path(root, "pipeline/backup", "surv_pre_clinical_std")

# reuse normalizers from 07 script
# [consolidated] the original preamble here read <root>/pipeline/R/07_standardize_clinical.R
# and eval()'d its definition region to reuse the normalizers.  That script is carried
# verbatim by this file and the same definitions are exposed at file level (see the
# shared-helper block at the top), so no re-read is needed.  Function bodies unchanged.

parse_components <- function(tok) {
  out <- list()
  if (is.na(tok) || !nzchar(trimws(tok))) return(out)
  s <- toupper(trimws(tok))
  pat <- "[PC]?([TNM])(?:IS|[0-9]+|X)([A-C])?(?=[\\s/,;_-]|[TNM]|[PC]|$)"
  m <- gregexpr(pat, s, perl = TRUE)[[1]]
  if (length(m) == 0L || m[1] == -1L) return(out)
  segs <- regmatches(s, list(m))[[1]]
  for (seg in segs) {
    mm <- regmatches(seg, regexpr("([TNM])", seg))
    letter <- substr(mm[1], 1, 1)
    if (!is.null(out[[letter]])) next
    out[[letter]] <- seg
  }
  out
}

files <- list.files(sdir, pattern = "_surv\\.rds$", full.names = TRUE)
loglist <- list()
for (f in files) {
  nm <- basename(f)
  d <- readRDS(f)
  if (!"stage" %in% colnames(d)) next
  if (any(c("T", "N", "M") %in% colnames(d))) next          # 已有独立列则跳过
  bf <- file.path(bdir, nm)
  orig <- if (file.exists(bf)) {
    od <- readRDS(bf); if ("stage" %in% colnames(od)) as.character(od[["stage"]]) else NULL
  } else NULL
  if (is.null(orig)) next
  comps <- lapply(orig, parse_components)
  have <- unique(unlist(lapply(comps, names)))
  if (!length(have)) next
  for (letter in c("T", "N", "M")) {
    if (!letter %in% have) next
    vals <- vapply(comps, function(cmp) {
      if (is.null(cmp[[letter]])) return(NA_character_)
      norm_tnm(cmp[[letter]], letter)
    }, character(1))
    d[[letter]] <- vals
  }
  saveRDS(d, f)
  nfilled <- sapply(c("T", "N", "M"), function(L)
    if (L %in% colnames(d)) sum(!is.na(d[[L]])) else 0L)
  loglist[[length(loglist) + 1L]] <- data.frame(
    file = nm, T_n = nfilled["T"], N_n = nfilled["N"], M_n = nfilled["M"],
    stringsAsFactors = FALSE)
}
if (length(loglist)) {
  lg <- do.call(rbind, loglist)
  write.csv(lg, file.path(out_dir, "clinical_std_tnm_split.csv"), row.names = FALSE)
  cat("added T/N/M columns to", nrow(lg), "files:\n")
  print(lg, row.names = FALSE)
} else cat("no files needed T/N/M split\n")
}

# ---------------------------------------------------------------------------
# run_07c_T_usage_table()  <-  verbatim pipeline/R/07c_T_usage_table.R
# ---------------------------------------------------------------------------
run_07c_T_usage_table <- function() {
# 07c_T_usage_table.R ---------------------------------------------------------
# 方案A：T 分布与使用建议表（不改动 surv 数据）
#   对 data/processed/surv/ 中所有含标准化 T 列的队列：
#   - T 分布计数 (T0..T4 / Tis / TX / NA)
#   - 该队列是否另有可用 stage 组别
#   - 使用建议: T-only 队列直接以 T 建模; 有 stage 者以 stage 为主、T 辅助
# 输出: pipeline/out/clinical_T_distribution_usage.csv
#       pipeline/out/REPORT_T_usage.md
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

root <- "/home/Jingle/data/Project/CPAS"
sdir <- file.path(root, "data/processed/surv")
outd <- cpas_out_root()
files <- list.files(sdir, pattern = "_surv\\.rds$", full.names = TRUE)
di <- read.csv(file.path(root, cpas_data("dataset_info.csv")), stringsAsFactors = FALSE)

TLEVELS <- c("T0", "T1", "T2", "T3", "T4", "Tis", "TX")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
type_of <- function(acc) {
  hit <- di$Type[di$Accession == acc]
  if (length(hit)) return(hit[1])
  pre <- di$Accession[startsWith(di$Accession, paste0(acc, "_"))]
  if (length(pre)) return(di$Type[di$Accession == pre[1]][1])
  NA_character_
}

rows <- list()
for (f in files) {
  d <- readRDS(f)
  if (!"T" %in% colnames(d)) next
  nm <- sub("_surv\\.rds$", "", basename(f))
  acc <- sub("[-_].*$", "", nm)
  tv <- d[["T"]]
  one <- data.frame(file = basename(f), accession = acc,
                    cancer_type = type_of(acc),
                    n = nrow(d),
                    stage_nonNA = if ("stage" %in% colnames(d))
                      sum(!is.na(d[["stage"]])) else 0L,
                    stringsAsFactors = FALSE)
  for (L in TLEVELS) one[[L]] <- sum(tv == L, na.rm = TRUE)
  one$NA_T <- sum(is.na(tv))
  rows[[length(rows) + 1L]] <- one
}
out <- do.call(rbind, rows)
out$suggestion <- ifelse(out$stage_nonNA > 0L,
  "有 stage 组别: 整合/分层分析以 stage(I-IV)为主, T 可作敏感性变量",
  "T-only: stage 无组别信息; 直接以 T 列作肿瘤负荷变量(COX/KM), 勿伪造分期")
out <- out[order(-out$n), ]
write.csv(out, file.path(outd, "clinical_T_distribution_usage.csv"), row.names = FALSE)

md <- c("# T 分布与使用建议表（方案 A：不改动数据）",
        sprintf("日期: %s | 覆盖含标准化 T 列的队列: %d 个", format(Sys.time()), nrow(out)), "",
        "## 统计口径",
        "- T 值已按统一规则标准化 (T0–T4 / Tis / TX；去掉 a/b/c 亚级与 p/c 前缀)。",
        "- `stage_nonNA>0` 表示该队列另有可用的临床分期组别；否则为 T-only 队列。",
        "- 本表不把 T 转换为分期组别（避免 N0M0 等未经证实的假设）。", "",
        "## 汇总表", "")
hdr <- c("数据集", "癌种", "n", "有stage非空", "T0", "T1", "T2", "T3", "T4",
         "Tis", "TX", "T缺失")
md <- c(md, paste0("| ", paste(hdr, collapse = " | "), " |"),
        paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"))
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  md <- c(md, sprintf("| %s | %s | %d | %s | %d | %d | %d | %d | %d | %d | %d | %d |",
    r$file, r$cancer_type, r$n, ifelse(r$stage_nonNA > 0, "是", "否"),
    r$T0, r$T1, r$T2, r$T3, r$T4, r$Tis, r$TX, r$NA_T))
}
md <- c(md, "",
        "## 使用建议",
        "- **T-only 队列**（stage 无组别）：直接用 T 列分组/建模，例如 T1 vs T2 vs T3+、或 T≥2 vs T1；",
        "  不可把 T 反推成 I–IV 期。",
        "- **有 stage 的队列**：跨队列整合以 stage(I–IV) 为主；T 用于敏感性或 T 特异性机制分析。",
        "- **Tis / TX**：Tis=原位癌、TX=无法评估（近缺失），建模时请显式处理（合并/剔除）。",
        "- **T 缺失较多**的队列：用非缺失子集分析，并在结果中报告缺失比例。",
        "",
        "逐队列明细见 pipeline/out/clinical_T_distribution_usage.csv")
writeLines(md, file.path(outd, "REPORT_T_usage.md"))
cat("rows:", nrow(out), "| files:", file.path(outd, "REPORT_T_usage.md"), "\n")
}

# ---------------------------------------------------------------------------
# run_08_verify_normalization()  <-  verbatim pipeline/R/08_verify_normalization.R
# ---------------------------------------------------------------------------
run_08_verify_normalization <- function() {
# 08_verify_normalization.R ---------------------------------------------------
# 复核 data/processed/surv/ 中所有 *_surv.rds 是否已完全规范化(只读, 不改文件)。
# 检查: 终点数值化; stage/T/N/M/sex/age/grade 值域; 缺失值指示词残留。
# 输出: pipeline/out/clinical_normalization_verify.csv
#       pipeline/out/REPORT_clinical_normalization_verify.md
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

root <- "/home/Jingle/data/Project/CPAS"
sdir <- file.path(root, "data/processed/surv")
outd <- cpas_out_root()
files <- list.files(sdir, pattern = "_surv\\.rds$", full.names = TRUE)

MISS <- c("", "NA", "N/A", "N.A.", "#NA", "#N/A", "NULL", "NaN", "-", "--",
          "—", "–", ".", "?", "UNKNOWN", "UNK", "99 UNK", "NOT AVAILABLE",
          "NOT REPORTED", "NOT SPECIFIED", "NOT APPLICABLE", "MISSING", "NONE",
          "[UNKNOWN]", "[DISCREPANCY]", "[NOT AVAILABLE]", "<NA>", "NAN")
is_miss <- function(x) {
  t <- trimws(x)
  tolower(t) %in% tolower(MISS) | grepl("^[#\\s]*n/?a[#\\s]*$", t, ignore.case = TRUE)
}
STD_EP <- c("OS", "DSS", "DFS", "RFS", "PFS", "MFS", "DRFS", "EFS", "DFI", "PFI")

check <- function(d) {
  issues <- character(0); infos <- character(0)
  cl <- colnames(d)
  ep_time <- cl[grepl("_(time)$", cl)]
  ep_stat <- cl[grepl("_(status)$", cl)]
  ep_stat <- ep_stat[vapply(sub("_status$", "", ep_stat),
                            function(x) x %in% STD_EP, logical(1))]
  ep <- c(ep_time, ep_stat)
  clin <- setdiff(cl, c(ep, "ID", "row_names"))

  for (cc in ep_time)
    if (!is.numeric(d[[cc]])) issues <- c(issues, paste0(cc, " 非数值"))
  for (sc in ep_stat) {
    if (!is.numeric(d[[sc]])) {
      issues <- c(issues, paste0(sc, " 非数值"))
    } else if (any(!is.na(d[[sc]]) & !d[[sc]] %in% c(0, 1))) {
      infos <- c(infos, paste0(sc, ": 状态编码为 ",
                 paste(unique(d[[sc]][!is.na(d[[sc]])]), collapse = ","),
                 " (非0/1, 上游口径, 未改动)"))
    }
  }
  if ("stage" %in% clin)
    if (any(!is.na(d$stage) & !d$stage %in% c("0", "I", "II", "III", "IV")))
      issues <- c(issues, "stage 值域异常")
  for (cc in c("T", "N", "M")) if (cc %in% clin) {
    okv <- switch(cc,
      T = c("T0", "T1", "T2", "T3", "T4", "Tis", "TX"),
      N = c("N0", "N1", "N2", "N3", "NX"),
      M = c("M0", "M1", "MX"))
    if (any(!is.na(d[[cc]]) & !d[[cc]] %in% okv))
      issues <- c(issues, paste0(cc, " 值域异常"))
  }
  if ("sex" %in% clin)
    if (any(!is.na(d$sex) & !d$sex %in% c("male", "female")))
      issues <- c(issues, "sex 值域异常")
  if ("age" %in% clin) {
    if (!is.numeric(d$age)) issues <- c(issues, "age 非数值")
    else if (any(d$age < 0 | d$age > 120, na.rm = TRUE))
      issues <- c(issues, "age 超范围")
  }
  if ("grade" %in% clin) {
    gv <- as.character(d$grade[!is.na(d$grade)])
    if (length(gv[!grepl("^G[1-4]$|^[5-9]$|^10$", gv)]))
      issues <- c(issues, "grade 值域异常")
  }
  for (cc in clin) {
    x <- d[[cc]]
    if (is.numeric(x) || is.integer(x)) next
    xc <- as.character(x)
    left <- xc[!is.na(xc) & is_miss(xc) & trimws(xc) != ""]
    if (length(left)) issues <- c(issues, paste0(cc, " 含缺失词残留(", length(left), ")"))
  }
  list(issues = unique(issues), infos = unique(infos))
}

res <- list(); info_map <- list()
for (f in files) {
  nm <- basename(f)
  d <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(d)) { res[[nm]] <- "读取失败"; next }
  ch <- check(d)
  if (length(ch$issues)) res[[nm]] <- ch$issues
  if (length(ch$infos)) info_map[[nm]] <- ch$infos
}
ok <- setdiff(basename(files), names(res))
df <- data.frame(file = c(ok, names(res)),
                 status = c(rep("OK", length(ok)), rep("ISSUE", length(res))),
                 issues = c(rep("", length(ok)),
                            vapply(res, paste, collapse = "; ", character(1))),
                 info = vapply(basename(files), function(nm)
                   if (!is.null(info_map[[nm]]))
                     paste(info_map[[nm]], collapse = "; ") else "", character(1)),
                 stringsAsFactors = FALSE)
write.csv(df, file.path(outd, "clinical_normalization_verify.csv"), row.names = FALSE)

md <- c("# 临床规范化复核报告 (data/processed/surv)",
        sprintf("日期: %s | 文件: %d | OK: %d | ISSUE: %d",
                format(Sys.time()), length(files),
                sum(df$status == "OK"), sum(df$status == "ISSUE")), "")
if (any(df$status == "ISSUE")) {
  md <- c(md, "## 存在问题的文件", "")
  bad <- df[df$status == "ISSUE", ]
  for (i in seq_len(nrow(bad)))
    md <- c(md, sprintf("- %s: %s", bad$file[i], bad$issues[i]))
} else {
  md <- c(md, "全部文件均已规范化，未发现问题。")
}
if (any(nzchar(df$info))) {
  md <- c(md, "", "## 注意（不影响规范化判定）", "")
  iv <- df[nzchar(df$info), ]
  for (i in seq_len(nrow(iv)))
    md <- c(md, sprintf("- %s: %s", iv$file[i], iv$info[i]))
}
writeLines(md, file.path(outd, "REPORT_clinical_normalization_verify.md"))
cat("verify: OK", sum(df$status == "OK"), "/", length(files),
    "| issues:", sum(df$status == "ISSUE"), "\n")
}

# ---------------------------------------------------------------------------
# run_11_endpoint_families()  <-  verbatim pipeline/R/11_endpoint_families.R
# ---------------------------------------------------------------------------
run_11_endpoint_families <- function() {
# 11_endpoint_families.R ------------------------------------------------------
# 方案B：为 dataset_info 增加生存终点"家族"映射列（不改数据库表）。
#   浏览家族(4): OS / DSS / DFS / PFS(broad, 含进展与转移类)
#   合并家族(5): OS / DSS / DFS / PFS / MFS   <- meta 与整合 KM 按此配对
# 规则:
#   单个 token 必须且只能映射到一个合并家族; 同一队列同一家族内出现多个
#   token 时按优先级择一并在报告中告警(当前数据为 0 冲突)。
# 衍生列:
#   EndpointFamilies  : 该队列可用的浏览家族(逗号分隔, OS,DSS,DFS,PFS)
#   EP_OS/EP_DSS/EP_DFS/EP_PFS/EP_MFS : 各合并家族在该队列中的具体 token(NA=无)
#   EndpointPrimary   : 默认优先家族的具体 token (OS > DSS > DFS > PFS > MFS)
#   EndpointDerived   : 其中来自 TCGA 派生的 token (DFI/PFI)
# 输出: dataset_info.csv/.rda(包内) + pipeline/out/endpoint_family_map.csv
#       + pipeline/out/REPORT_endpoint_families.md
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

root <- "/home/Jingle/data/Project/CPAS"
setwd(root)

# token -> 合并家族
FAMILY_OF <- c(
  OS = "OS",
  DSS = "DSS", CSS = "DSS", BCSS = "DSS",
  DFS = "DFS", RFS = "DFS", EFS = "DFS", DFI = "DFS",
  PFS = "PFS", PFI = "PFS",
  MFS = "MFS", DRFS = "MFS"
)
DERIVED <- c("DFI", "PFI")                       # TCGA 派生终点
BROWSE_OF <- c(OS = "OS", DSS = "DSS", DFS = "DFS",
               PFS = "PFS", MFS = "PFS")          # 合并家族 -> 浏览家族
POOL_ORDER <- c("OS", "DSS", "DFS", "PFS", "MFS") # 解析优先级

di <- read.csv(cpas_data("dataset_info.csv"), stringsAsFactors = FALSE)
if (!"SurvivalTypes" %in% colnames(di)) stop("SurvivalTypes column missing")

tokens_of <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  t <- trimws(strsplit(x, ",")[[1]]); t[nzchar(t)]
}

conflicts <- list(); unmapped <- character(0)
pool_cols <- list(OS = character(nrow(di)), DSS = character(nrow(di)),
                  DFS = character(nrow(di)), PFS = character(nrow(di)),
                  MFS = character(nrow(di)))
browse_list <- character(nrow(di)); primary <- character(nrow(di))
derived_list <- character(nrow(di))

for (i in seq_len(nrow(di))) {
  tk <- tokens_of(di$SurvivalTypes[i])
  if (!length(tk)) next
  bad <- setdiff(tk, names(FAMILY_OF))
  if (length(bad)) unmapped <- c(unmapped, bad)
  fam <- unname(FAMILY_OF[tk[tk %in% names(FAMILY_OF)]])
  # 每个家族内择一(优先级)
  for (f in POOL_ORDER) {
    hit <- tk[tk %in% names(FAMILY_OF)][fam == f]
    if (length(hit)) {
      hit <- hit[order(match(hit, names(FAMILY_OF)[FAMILY_OF == f]))]
      pool_cols[[f]][i] <- hit[1]
      if (length(hit) > 1)
        conflicts[[length(conflicts) + 1L]] <- data.frame(
          accession = di$Accession[i], family = f,
          tokens = paste(hit, collapse = ","), chosen = hit[1],
          stringsAsFactors = FALSE)
    }
  }
  avail <- POOL_ORDER[vapply(pool_cols, function(v) nzchar(v[i]), logical(1))]
  browse_list[i] <- paste(unique(unname(BROWSE_OF[avail])), collapse = ",")
  primary[i] <- if (length(avail)) pool_cols[[avail[1]]][i] else NA_character_
  dv <- tk[tk %in% DERIVED]
  derived_list[i] <- if (length(dv)) paste(dv, collapse = ",") else NA_character_
}

di$EndpointFamilies <- ifelse(nzchar(browse_list), browse_list, NA_character_)
for (f in POOL_ORDER)
  di[[paste0("EP_", f)]] <- ifelse(nzchar(pool_cols[[f]]), pool_cols[[f]], NA_character_)
di$EndpointPrimary <- primary
di$EndpointDerived <- derived_list
# 列顺序
ord <- c("X", "Type", "Accession", "SurvivalTypes", "EndpointFamilies",
         paste0("EP_", POOL_ORDER), "EndpointPrimary", "EndpointDerived",
         "GPL", "N", "method")
di <- di[, c(intersect(ord, colnames(di)), setdiff(colnames(di), ord))]
write.csv(di, cpas_data("dataset_info.csv"), row.names = FALSE)
dataset_info <- di
save(dataset_info, file = "CanPAS/data/dataset_info.rda", compress = "xz")

# ---- 审计输出 ----
mapdf <- data.frame(token = names(FAMILY_OF), pooling_family = unname(FAMILY_OF),
                    browse_family = unname(BROWSE_OF[unname(FAMILY_OF)]),
                    derived = names(FAMILY_OF) %in% DERIVED,
                    stringsAsFactors = FALSE)
write.csv(mapdf, cpas_out("endpoint_family_map.csv"), row.names = FALSE)

cov <- sapply(POOL_ORDER, function(f) sum(!is.na(di[[paste0("EP_", f)]])))
browse_cov <- sapply(c("OS", "DSS", "DFS", "PFS"), function(f)
  sum(!is.na(di$EndpointFamilies) & grepl(f, di$EndpointFamilies, fixed = TRUE)))

md <- c("# 生存终点家族映射报告 (方案 B)",
        sprintf("日期: %s | 队列: %d", format(Sys.time()), nrow(di)), "",
        "## 映射规则",
        paste0("- 合并家族(5): ",
               paste(sprintf("%s = {%s}", POOL_ORDER,
                             sapply(POOL_ORDER, function(f)
                               paste(names(FAMILY_OF)[FAMILY_OF == f], collapse = "/"))),
                     collapse = "; ")),
        "- 浏览家族(4): OS / DSS / DFS / PFS(broad = EP_PFS 或 EP_MFS)",
        paste0("- TCGA 派生终点(报告中标注): ", paste(DERIVED, collapse = ", ")), "",
        "## 覆盖队列数",
        paste0("- 浏览家族: ", paste(names(browse_cov), browse_cov, sep = "=", collapse = ", ")),
        paste0("- 合并家族: ", paste(names(cov), cov, sep = "=", collapse = ", ")),
        paste0("- 有 EndpointPrimary 默认终点的队列: ", sum(!is.na(di$EndpointPrimary))),
        paste0("- 无任何标准终点的队列: ", sum(is.na(di$EndpointFamilies))), "",
        "## 家族内多 token 冲突",
        if (length(conflicts)) "见下表(已按优先级择一)" else "无(每个队列在每个家族内至多 1 个 token)")
if (length(conflicts)) {
  cf <- do.call(rbind, conflicts)
  write.csv(cf, cpas_out("endpoint_family_conflicts.csv"), row.names = FALSE)
  for (i in seq_len(nrow(cf)))
    md <- c(md, sprintf("- %s [%s]: %s -> 采用 %s",
                        cf$accession[i], cf$family[i], cf$tokens[i], cf$chosen[i]))
}
md <- c(md, "", "## 未映射 token",
        if (length(unique(unmapped))) paste("-", paste(unique(unmapped), collapse = ", "))
        else "无", "",
        "## 迁移前后对照(可合并队列数)",
        sprintf("- DFS: 单用 DFS %d 个队列 -> DFS 家族 %d 个队列",
                sum(di$SurvivalTypes == "DFS", na.rm = TRUE), cov[["DFS"]]),
        sprintf("- 进展类: PFS %d + PFI %d -> PFS 家族 %d",
                sum(grepl("(^|,)PFS(,|$)", di$SurvivalTypes), na.rm = TRUE),
                sum(grepl("(^|,)PFI(,|$)", di$SurvivalTypes), na.rm = TRUE), cov[["PFS"]]),
        sprintf("- 转移类: MFS %d + DRFS %d -> MFS 家族 %d",
                sum(grepl("(^|,)MFS(,|$)", di$SurvivalTypes), na.rm = TRUE),
                sum(grepl("(^|,)DRFS(,|$)", di$SurvivalTypes), na.rm = TRUE), cov[["MFS"]]))
writeLines(md, cpas_out("REPORT_endpoint_families.md"))
cat("done. conflicts:", length(conflicts), "| unmapped:", length(unique(unmapped)), "\n")
print(cov); print(browse_cov)
}

# ---------------------------------------------------------------------------
# run_18_log2_transform()  <-  verbatim pipeline/R/18_log2_transform.R
# ---------------------------------------------------------------------------
run_18_log2_transform <- function() {
# 18_log2_transform.R --------------------------------------------------------
# 表达谱尺度统一：把「未经 log2 转化」的数据集做 log2(value + 1)。
#
# 背景
#   data/expr/*.rds 由 01_parse_gse.R 用 Biobase::exprs(eset) 直接保存，
#   即 GEO 提交者上传的原始尺度：部分数据集是线性荧光强度（MAS5/Agilent 等，
#   取值 10^0–10^5），部分是 log2（取值 ~0–17，常含负值，如 z-score/中心化值）。
#   两者混在同一库里会让「每 1 个表达单位」的 HR 既不可读也不可比，
#   且线性尺度下 Cox 假定 log-hazard 对原始强度线性，与 log2 数据的前提不同。
#
# ============================ 判定规则 v2（2026-09-24 加固） ============================
# v1 只用「max/median <= 50」一条比值规则，实测漏判：GSE91061 已是 log2 尺度
#   （min=0, median=0.3217, max=16.9442 -> max/median = 52.67，仅因比值略高于 50），
#   被误判为 needs_log2。若再跑一次 --apply 会被二次 log2，静默破坏数据。
# v2 改为「决定性规则 + 多证据投票」，比值只是其中一条证据。
#
# 逐条规则与理由（全部基于全矩阵精确统计，不做抽样）：
#   D1  min < 0                                  -> already_log
#       理由：线性荧光强度/RPKM 计数不可能为负；负值只可能来自 log 后的中心化/z-score。
#   D2  median > 0 且 max/median <= 50           -> already_log
#       理由：线性强度的动态范围通常跨越 3–5 个数量级（实测备份文件 213–5×10^5），
#             而 log2 尺度极少超过 50。这是 v1 唯一规则，保留但不再单独承担判定。
#   D3  median <= 0 且 max <= 25                 -> already_log
#       理由：中位数为 0 时比值退化为 Inf，比值规则失效；此时只能看绝对范围。
#   V   多证据投票 >= 3 / 4                      -> already_log
#       V1 近正的极小值: min >= 0 且 min <= 1
#           理由：log2(x+1)（x>=0）的极小值落在 0 附近；线性强度下界通常远大于 1。
#       V2 分数值占比: 非零值中「非整数」比例 >= 90%
#           理由：log2 后的连续型数据几乎处处是分数。注意必须排除 0（0 是整数且
#                 在 log 数据里占比很高，直接用全矩阵整数占比会被 0 淹没 —— GSE91061
#                 全矩阵整数占比 35.9%，但非零值里整数占比只有 0.03%）。
#       V3 99 分位有界: p99 <= 20
#           理由：log2 尺度下 99% 的探针都在 20 以内（CGGA_693 的 RSEM 极值 23.8 已接近上界）。
#       V4 绝对动态范围小: max <= 25
#           理由：同 V3 的上界依据；25 而非 20 是刻意留出 CGGA_693（max=23.8）的余量。
#   V 投票的意义：当 max/median 落在 (50, 100] 的灰区时，单看比值无法判定，
#     用 4 条彼此独立的证据投票。GSE91061 在灰区内得 4/4 票 -> already_log。
#   其余                                          -> needs_log2
#
# 误判方向的取舍（重要）：
#   本脚本的误判必须偏向「漏转」（判成 already_log 而不转），绝不能偏向「错转」
#   （把已是 log2 的数据再 log2）。漏转会被下游 QC/尺度检查发现，而双 log2 是静默
#   破坏：数值仍在合理范围内，分布形状也仍然像表达谱，几乎无法察觉。因此所有
#   边界都向 already_log 一侧放宽。
#   代价：若未来有数据集的线性值恰好全部落在 0–25（实测 167 个文件里没有这种情形），
#   会被判成 already_log 而跳过；这属于可接受的漏转。
#
# 已知局限：
#   若某数据集已是 log 尺度但 max > 25（极端离群探针）且 max/median > 100，
#   且 V1/V2 不成立，则仍会被判 needs_log2。当前 167 个文件中不存在这种情形
#   （max > 25 的 3 个文件 GSE22153 / GSE4716-GPL3696 / GSE76427 都因 min < 0 被 D1 捕获）。
#
# 来源白名单：CGGA_301/325/693 由 10_parse_cgga.R 解析，该脚本已按数据集显式声明
#   是否做 log2（RNA-seq 用 log2(RSEM+1)），不在此脚本重复处理。
#
# 幂等性：转化后 min>=0 且 max<=25，再次运行会被判为 already_log 而跳过（v2 起可靠）。
#
# ============================ --apply 安全闸门 ============================
# 自 2026-09-24 起，--apply 只有在**全部**满足下列条件时才写任何一个字节：
#   G1 必须给出显式白名单 --only=<ACC,...> / --allow=<ACC,...>（可重复、可逗号分隔）。
#      没有白名单 -> 在**读取任何 data/expr 文件之前**就拒绝执行（退出码 3），
#      只写一条拒绝记录到审计 CSV。
#      闸门目的：杜绝「无参数误跑」把整库重转一遍。
#   G2 白名单里的每个样本必须存在、不是 provenance 白名单、且判定为 needs_log2。
#      只要有一个被判为 already_log -> 整个运行拒绝执行（退出码 3）。
#      闸门目的：already_log 的数据集绝不写盘，从机制上堵死 GSE91061 这类双 log2。
#   G3 写出前先在内存里算好 log2(x+1) 并做范围自检（finite、min>=0、max<=25），
#      自检不通过就不落盘。即「先验证、后写」。
#   G4 每个数据集在写之前备份到 pipeline/backup/expr_pre_log2/<ACC>.rds（已存在则不覆盖，
#      保证回滚点是真正的「原始未转化」版本）。
#   G5 无论 dry-run 还是 --apply，都写 machine-readable 审计报告
#      pipeline/out/log2_apply_<YYYYMMDD>.csv，逐数据集记录判定、证据、动作、拒绝原因、
#      写前/写后统计与备份路径。
#   G6 无参数默认仍是 dry-run（只体检 + 出报告，不改 data/）。
#
# 用法
#   Rscript pipeline/R/18_log2_transform.R                       # 体检（默认 dry-run）
#   Rscript pipeline/R/18_log2_transform.R --only=GSE1234        # 预览白名单的判定与写后统计
#   Rscript pipeline/R/18_log2_transform.R --apply --only=GSE1234,GSE5678   # 备份 + 转化
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/18_log2_transform.R \
#       --apply --only=GSE1234 --db                              # 额外把库中对应表就地 UPDATE
# 输出
#   pipeline/out/expr_scale_survey.csv        每个数据集精确分布、判定、证据与投票
#   pipeline/out/log2_apply_<YYYYMMDD>.csv    本次运行的逐数据集审计记录（dry-run 也写）
#   pipeline/backup/expr_pre_log2/            转化前的 .rds 原件（可回滚）
#   pipeline/out/REPORT_log2_transform.md     本次 --apply 的操作记录
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

suppressPackageStartupMessages(library(RMySQL))

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
args    <- commandArgs(trailingOnly = TRUE)
apply_l <- "--apply" %in% args
use_db  <- "--db" %in% args
LIMIT   <- as.integer(Sys.getenv("CPAS_LOG2_LIMIT", unset = "0"))   # >0 时只处理前 N 个（调试用）

SKIP_PROVENANCE <- c("CGGA_301", "CGGA_325", "CGGA_693")   # 见文件头说明

# ---------- 0) 白名单解析（--only / --allow） ----------
parse_allow <- function(a) {
  vals <- character(0)
  for (i in seq_along(a)) {
    if (grepl("^--(only|allow)=", a[i])) {
      vals <- c(vals, strsplit(sub("^--(only|allow)=", "", a[i]), ",", fixed = TRUE)[[1]])
    } else if (a[i] %in% c("--only", "--allow")) {
      if (i < length(a) && !grepl("^--", a[i + 1L]))
        vals <- c(vals, strsplit(a[i + 1L], ",", fixed = TRUE)[[1]])
      else stop("--only/--allow 后需要逗号分隔的样本号，例如 --only=GSE1234,GSE5678")
    }
  }
  vals <- trimws(vals)
  unique(vals[nzchar(vals)])
}
ALLOW <- parse_allow(args)

expr_dir <- file.path(root, "data/expr")
bak_dir  <- file.path(root, "pipeline/backup/expr_pre_log2")
out_dir  <- cpas_out_root()
dir.create(bak_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 0b) 审计报告骨架（dry-run / apply / 拒绝 共用同一套列） ----------
audit_path <- file.path(out_dir, sprintf("log2_apply_%s.csv", format(Sys.Date(), "%Y%m%d")))
AUDIT_COLS <- c("mode", "timestamp", "dataset", "verdict", "evidence", "votes", "allowed",
                "action", "reason", "rows", "cols",
                "min_before", "median_before", "p99_before", "max_before", "ratio_before",
                "min_after", "median_after", "p99_after", "max_after",
                "bytes_before", "bytes_after", "backup", "db_action")
new_audit <- function(n) {
  a <- data.frame(mode = character(n), timestamp = character(n), dataset = character(n),
                  verdict = character(n), evidence = character(n), votes = integer(n),
                  allowed = logical(n), action = character(n), reason = character(n),
                  stringsAsFactors = FALSE)
  for (cc in setdiff(AUDIT_COLS, names(a))) a[[cc]] <- rep(NA, n)
  a[, AUDIT_COLS, drop = FALSE]
}

# ---------- 0c) G1 预检：--apply 必须带白名单，否则在读取任何表达文件之前就拒绝 ----------
if (apply_l && !length(ALLOW)) {
  a <- new_audit(1L)
  a$mode <- "apply-refused"; a$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  a$allowed <- FALSE; a$action <- "refused_no_allowlist"
  a$reason <- "--apply 未给出白名单（--only/--allow），已在读取任何 data/expr 文件之前拒绝"
  write.csv(a, audit_path, row.names = FALSE)
  cat("================ 拒绝执行 (REFUSED) ================\n")
  cat("--apply 需要显式白名单：--only=<ACC,...> 或 --allow=<ACC,...>。\n")
  cat("这是为了防止无参数误跑把整个表达库重转一遍。\n")
  cat("本次未读取、未写入任何 data/expr 文件。\n")
  cat("审计记录:", sub(paste0(root, "/"), "", audit_path, fixed = TRUE), "\n")
  cat("===================================================\n")
  quit(save = "no", status = 3L)
}

# ---------- 1) 全矩阵体检每个本地表达文件（精确统计，避免抽样误判） ----------
survey <- function(path) {
  x <- readRDS(path)
  if (!is.data.frame(x)) x <- as.data.frame(x)
  m <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  v <- m[is.finite(m)]
  nz <- v[v != 0]                       # V2 的证据只在非零值上算（见文件头 V2 说明）
  q <- quantile(v, c(.25, .5, .75, .99), na.rm = TRUE, names = FALSE)
  data.frame(rows = nrow(m), cols = ncol(m),
             min = min(v), q25 = q[1], median = q[2], q75 = q[3], p99 = q[4], max = max(v),
             pct_zero = round(100 * mean(v == 0), 1),
             pct_neg  = round(100 * mean(v < 0), 2),
             pct_nonint_nonzero = round(100 * (if (length(nz)) mean(abs(nz - round(nz)) >= 1e-9) else NA_real_), 3))
}

# 判定：返回 verdict / evidence / votes / reason 四元组
classify <- function(d) {
  if (is.na(d$median))
    return(list(verdict = "unreadable", evidence = "", votes = NA_integer_, reason = "文件无法读取"))
  ratio <- if (is.finite(d$median) && d$median > 0) d$max / d$median else Inf
  # 4 条独立证据
  v <- c(V1 = isTRUE(d$min >= 0 && d$min <= 1),
         V2 = isTRUE(!is.na(d$pct_nonint_nonzero) && d$pct_nonint_nonzero >= 90),
         V3 = isTRUE(d$p99 <= 20),
         V4 = isTRUE(d$max <= 25))
  nv <- sum(v)
  if (d$min < 0)
    return(list(verdict = "already_log", evidence = "D1", votes = nv,
                reason = sprintf("min=%.3f<0：线性强度不可能为负", d$min)))
  if (d$median > 0 && ratio <= 50)
    return(list(verdict = "already_log", evidence = "D2", votes = nv,
                reason = sprintf("max/median=%.2f<=50：动态范围远小于线性强度", ratio)))
  if (d$median <= 0 && d$max <= 25)
    return(list(verdict = "already_log", evidence = "D3", votes = nv,
                reason = sprintf("median=%.3f<=0 且 max=%.3f<=25", d$median, d$max)))
  if (nv >= 3)
    return(list(verdict = "already_log",
                evidence = sprintf("V%d/4[%s]", nv, paste(names(v)[v], collapse = "+")), votes = nv,
                reason = sprintf("多证据投票 %d/4（max/median=%.2f 灰区）", nv, ratio)))
  list(verdict = "needs_log2", evidence = sprintf("V%d/4", nv), votes = nv,
       reason = sprintf("线性强度：max/median=%.1f, max=%.1f, p99=%.1f, 非零分数值=%.2f%%",
                        ratio, d$max, d$p99, d$pct_nonint_nonzero))
}

files <- list.files(expr_dir, pattern = "\\.rds$", full.names = TRUE)
if (LIMIT > 0) files <- head(files, LIMIT)
message("体检 ", length(files), " 个本地表达文件 ...")
rep_rows <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- try(survey(files[i]), silent = TRUE)
  if (inherits(d, "try-error"))
    d <- data.frame(rows = NA, cols = NA, min = NA, q25 = NA, median = NA, q75 = NA,
                    p99 = NA, max = NA, pct_zero = NA, pct_neg = NA, pct_nonint_nonzero = NA)
  cl <- classify(d)
  d$dataset  <- sub("\\.rds$", "", basename(files[i]))
  d$verdict  <- cl$verdict
  d$evidence <- cl$evidence
  d$votes    <- cl$votes
  d$reason   <- cl$reason
  d$bytes    <- file.size(files[i])
  rep_rows[[i]] <- d
  if (i %% 25 == 0) message("  ...", i, "/", length(files))
}
rep <- do.call(rbind, rep_rows)
rep$ratio <- ifelse(is.finite(rep$median) & rep$median > 0, rep$max / rep$median, Inf)
rep <- rep[, c("dataset", "verdict", "evidence", "votes", "reason", "rows", "cols", "min", "q25",
               "median", "q75", "p99", "max", "ratio", "pct_zero", "pct_neg",
               "pct_nonint_nonzero", "bytes")]
write.csv(rep, file.path(out_dir, "expr_scale_survey.csv"), row.names = FALSE)
cat("\n判定结果:\n"); print(table(rep$verdict))

rep$provenance_skip <- rep$dataset %in% SKIP_PROVENANCE
todo <- rep[rep$verdict == "needs_log2" & !rep$provenance_skip, , drop = FALSE]
cat("\n需要 log2(x+1) 的数据集:", nrow(todo), "\n")
if (nrow(todo)) print(todo[, c("dataset", "rows", "cols", "min", "median", "p99", "max", "ratio", "reason")], row.names = FALSE)

# ---------- 2) 安全闸门 ----------
# 逐个白名单条目给出动作与拒绝原因（dry-run 只报告，--apply 遇拒绝即中止）
audit_action <- rep("not_selected", nrow(rep))
audit_reason <- rep("", nrow(rep))
refusals <- character(0)
if (length(ALLOW)) {
  miss <- setdiff(ALLOW, rep$dataset)
  if (length(miss)) refusals <- c(refusals, sprintf("白名单样本在 data/expr 中不存在: %s", paste(miss, collapse = ", ")))
  for (a in ALLOW) {
    i <- match(a, rep$dataset)
    if (is.na(i)) next
    if (rep$provenance_skip[i]) {
      refusals <- c(refusals, sprintf("%s 属于 provenance 白名单（由 10_parse_cgga.R 处理），不得在此转化", a))
      audit_action[i] <- "refused_provenance"; audit_reason[i] <- "provenance 白名单"
    } else if (rep$verdict[i] == "already_log") {
      # G2：这是防二次 log2 的核心闸门
      refusals <- c(refusals, sprintf("%s 判定为 already_log（证据 %s：%s），拒绝写入以免二次 log2",
                                      a, rep$evidence[i], rep$reason[i]))
      audit_action[i] <- "refused_already_log"; audit_reason[i] <- rep$reason[i]
    } else if (rep$verdict[i] == "needs_log2") {
      audit_action[i] <- "selected"; audit_reason[i] <- rep$reason[i]
    } else {
      refusals <- c(refusals, sprintf("%s 判定为 %s，无法安全转化", a, rep$verdict[i]))
      audit_action[i] <- "refused_unreadable"; audit_reason[i] <- rep$reason[i]
    }
  }
} else {
  audit_action <- ifelse(rep$verdict == "needs_log2", "skipped_no_allowlist", "skip_already_log")
  audit_reason <- ifelse(rep$verdict == "needs_log2", "未给出白名单（--only/--allow）", rep$reason)
}

TARGET <- rep$dataset[audit_action == "selected"]

# 写前/写后统计：log2(x+1) 单调，故分位数（含中位数）可解析换算，dry-run 也能给出精确的写后值
analytic_after <- function(d) log2(pmax(d + 1, 0))
audit <- new_audit(nrow(rep))
audit$mode <- if (apply_l) "apply" else "dry-run"
audit$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
audit$dataset <- rep$dataset; audit$verdict <- rep$verdict
audit$evidence <- rep$evidence; audit$votes <- as.integer(rep$votes)
audit$allowed <- rep$dataset %in% ALLOW
audit$action <- audit_action; audit$reason <- audit_reason
audit$rows <- rep$rows; audit$cols <- rep$cols
audit$min_before <- rep$min; audit$median_before <- rep$median
audit$p99_before <- rep$p99; audit$max_before <- rep$max
audit$ratio_before <- rep$ratio
audit$bytes_before <- rep$bytes
audit$db_action <- "not_requested"
sel <- audit$action == "selected"
audit$min_after[sel]    <- analytic_after(rep$min[sel])
audit$median_after[sel] <- analytic_after(rep$median[sel])
audit$p99_after[sel]    <- analytic_after(rep$p99[sel])
audit$max_after[sel]    <- analytic_after(rep$max[sel])

write.csv(audit, audit_path, row.names = FALSE)

if (nrow(todo) == 0) {
  cat("\n没有判定为 needs_log2 的数据集（全部已在 log 尺度），无需转化。\n")
}
cat("\n安全闸门状态:\n")
cat("  --apply        :", ifelse(apply_l, "是", "否（默认 dry-run）"), "\n")
cat("  白名单 (-only) :", ifelse(length(ALLOW), paste(ALLOW, collapse = ","), "<空>"), "\n")
cat("  允许写入的数据集:", ifelse(length(TARGET), paste(TARGET, collapse = ","), "<无>"), "\n")
if (length(refusals)) { cat("  拒绝原因:\n"); for (r in refusals) cat("    -", r, "\n") }
cat("  审计报告      :", sub(paste0(root, "/"), "", audit_path, fixed = TRUE), "\n")

refuse <- function(msg, code = 3L) {
  cat("\n================ 拒绝执行 (REFUSED) ================\n")
  cat(msg, "\n")
  cat("未写入任何 data/ 文件。审计记录见:", sub(paste0(root, "/"), "", audit_path, fixed = TRUE), "\n")
  cat("===================================================\n")
  quit(save = "no", status = code)
}

if (!apply_l) {
  if (length(ALLOW) && length(refusals))
    cat("\n[dry-run] 注意：按上述拒绝原因，--apply 会被闸门拦下，不会有任何写入。\n")
  cat("\n[dry-run] 未做任何修改。加 --apply 且给出 --only=<ACC,...> 才会备份并转化本地文件。\n")
  quit(save = "no", status = 0L)
}

# ---- 以下只在 --apply 生效：任何一条不满足即整轮拒绝（fail closed） ----
if (!length(ALLOW))
  refuse(paste0("--apply 需要显式白名单：--only=<ACC,...> 或 --allow=<ACC,...>。\n",
                "  这是为了防止无参数误跑把整个表达库重转一遍。\n",
                "  当前判定为 needs_log2 的数据集: ",
                ifelse(nrow(todo), paste(todo$dataset, collapse = ","), "<无>")))
if (length(refusals)) refuse(paste(refusals, collapse = "\n"))
if (!length(TARGET)) refuse("白名单未选中任何 needs_log2 数据集，无可执行动作。")
if (any(rep$min[match(TARGET, rep$dataset)] < 0))
  refuse("白名单中存在负值数据集，log2(x+1) 不适用。")

# ---------- 3) 先验证、后写：备份 + 本地转化 ----------
log_lines <- character(0)
note <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
note("## 本地文件转化（log2(value + 1)）")
note(sprintf("白名单: %s", paste(ALLOW, collapse = ", ")))
for (a in TARGET) {
  f <- file.path(expr_dir, paste0(a, ".rds"))
  b <- file.path(bak_dir,  paste0(a, ".rds"))
  x <- readRDS(f)
  ids <- x[[1]]
  m <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  m2 <- log2(m + 1)                       # 关键转化
  # G3 写前自检：不通过就不落盘
  if (any(!is.finite(m2[is.finite(m)])))
    refuse(sprintf("%s 转化后出现非有限值，已中止。", a))
  rng <- range(m2[is.finite(m2)])
  if (rng[1] < 0 || rng[2] > 25)
    refuse(sprintf("%s 转化后范围异常 (min=%.3f, max=%.3f)，超出 log2 合理区间 [0,25]，已中止。", a, rng[1], rng[2]))
  # G4 备份（已存在则不覆盖，保留真正的原始版本）
  if (!file.exists(b)) { if (!file.copy(f, b, overwrite = FALSE)) refuse(sprintf("备份失败: %s", a)) }
  out <- data.frame(ID_REF = ids, m2, check.names = FALSE)
  names(out) <- names(x)
  saveRDS(out, f)
  i <- match(a, audit$dataset)
  audit$min_after[i]    <- min(m2, na.rm = TRUE)
  audit$median_after[i] <- median(m2, na.rm = TRUE)
  audit$p99_after[i]    <- quantile(m2, .99, na.rm = TRUE, names = FALSE)
  audit$max_after[i]    <- max(m2, na.rm = TRUE)
  audit$bytes_after[i]  <- file.size(f)
  audit$backup[i]       <- sub(paste0(root, "/"), "", b, fixed = TRUE)
  note(sprintf("  %-20s rows=%d cols=%d | before min=%.3f median=%.3f max=%.3f -> after min=%.3f median=%.3f max=%.3f | backup=%s",
               a, nrow(out), ncol(out) - 1,
               rep$min[i], rep$median[i], rep$max[i],
               audit$min_after[i], audit$median_after[i], audit$max_after[i],
               audit$backup[i]))
}
note(sprintf("本地已转化 %d 个数据集，原件备份于 pipeline/backup/expr_pre_log2/", length(TARGET)))
cat("\n写后统计（逐数据集）:\n")
print(audit[audit$action == "selected",
            c("dataset", "min_before", "median_before", "max_before",
              "min_after", "median_after", "max_after")], row.names = FALSE)

# ---------- 4) 数据库就地转化（LOG2(col + 1)） ----------
if (use_db) {
  pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
  if (!nzchar(pw)) refuse("CPAS_DB_PASSWORD 未设置，无法连接镜像。")
  con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                   dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                   password = pw, client_flag = CLIENT_COMPRESS)
  on.exit(dbDisconnect(con), add = TRUE)
  tb <- dbGetQuery(con, "SHOW TABLES")[[1]]
  note("## 镜像表就地 UPDATE（LOG2(col + 1)）")
  for (a in TARGET) {                       # 只动白名单内、已通过闸门的数据集
    t <- gsub("-", "_", a)
    i <- match(a, audit$dataset)
    if (!t %in% tb) { note(sprintf("  %-20s 库中无此表，跳过（仅本地）", a)); audit$db_action[i] <- "no_table"; next }
    cols <- setdiff(dbListFields(con, t), "row_names")
    sets <- paste0("`", cols, "` = LOG2(`", cols, "` + 1)", collapse = ", ")
    nr <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", t))$n
    dbExecute(con, sprintf("UPDATE `%s` SET %s", t, sets))
    audit$db_action[i] <- "updated"
    note(sprintf("  %-20s %d 行 x %d 列 已转化", a, nr, length(cols)))
  }
}

write.csv(audit, audit_path, row.names = FALSE)   # 覆盖为含写后统计的最终版

# ---------- 5) 记录 ----------
hdr <- c("# 表达谱 log2 转化记录",
         "",
         paste0("- 时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
         paste0("- 根目录: ", root),
         paste0("- 判定规则: v2 决定性规则 D1/D2/D3 + 4 证据投票(V1 近正极小值, V2 非零分数值>=90%, ",
                "V3 p99<=20, V4 max<=25) >=3/4；见脚本文件头"),
         paste0("- 体检文件数: ", nrow(rep), "；需转化: ", nrow(todo),
                "；已是 log2: ", sum(rep$verdict == "already_log")),
         paste0("- 白名单 (-only): ", ifelse(length(ALLOW), paste(ALLOW, collapse = ","), "<空>")),
         paste0("- 实际写入: ", length(TARGET), " 个 (", ifelse(length(TARGET), paste(TARGET, collapse = ","), "-"), ")"),
         paste0("- 审计报告: pipeline/out/", basename(audit_path)),
         paste0("- 数据库就地转化: ", ifelse(use_db, "是", "否")),
         paste0("- 回滚: 用 pipeline/backup/expr_pre_log2/<ACC>.rds 覆盖回 data/expr/，",
                "库端按 06_upload_db.R 重新上传"),
         "", "```text", log_lines, "```", "")
writeLines(hdr, file.path(out_dir, "REPORT_log2_transform.md"))
cat("\n已写出: pipeline/out/REPORT_log2_transform.md, expr_scale_survey.csv 与", basename(audit_path), "\n")
}

# ---------------------------------------------------------------------------
# run_20_log2_transform_db()  <-  verbatim pipeline/R/20_log2_transform_db.R
# ---------------------------------------------------------------------------
run_20_log2_transform_db <- function() {
# 20_log2_transform_db.R -----------------------------------------------------
# 镜像端就地转化：把「线性尺度」表达表的每个样本列做 LOG2(col + 1)。
#
# 与 18_log2_transform.R 的关系
#   18 负责本地 .rds（已执行，原件在 pipeline/backup/expr_pre_log2/）；
#   本脚本把同样的变换应用到 MySQL 镜像表，不重新上传数据（列多、表大，
#   服务端 UPDATE 比整表回传快几个数量级）。变换是确定性函数，两端结果一致
#   （差异仅在 FLOAT 单精度舍入，~1e-7 相对误差）。
#
# 幂等与安全守卫
#   对每张表抽样比对：DB 值与 backup（转化前）一致 -> 执行 UPDATE；
#   DB 值与 data/expr（转化后）一致 -> 已转化，跳过；
#   两者都不一致 -> 报警并跳过，不做任何修改。
#
# 用法
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/20_log2_transform_db.R            # dry-run
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/20_log2_transform_db.R --apply
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

.libPaths(c("/home/Jingle/R/library", .libPaths()))
suppressPackageStartupMessages(library(RMySQL))

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
apply_db <- "--apply" %in% commandArgs(trailingOnly = TRUE)

bak_dir <- file.path(root, "pipeline/backup/expr_pre_log2")
accs <- sub("\\.rds$", "", list.files(bak_dir, pattern = "\\.rds$"))
if (!length(accs)) stop("备份目录为空，先跑 18_log2_transform.R --apply")
cat("备份中的数据集:", length(accs), "\n")

pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
if (!nzchar(pw)) stop("CPAS_DB_PASSWORD 未设置。")
con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = pw, client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)
tb <- dbGetQuery(con, "SHOW TABLES")[[1]]

backtick <- function(x) paste0("`", x, "`")
set.seed(42)
log_lines <- character(0)
note <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }

for (a in accs) {
  t <- gsub("-", "_", a, fixed = TRUE)
  if (!t %in% tb) { note(sprintf("  %-20s 镜像无此表（仅本地）", a)); next }
  pre  <- readRDS(file.path(bak_dir, paste0(a, ".rds")))
  post <- readRDS(file.path(root, "data/expr", paste0(a, ".rds")))
  cols <- setdiff(dbListFields(con, t), "row_names")
  # 抽样：最多 25 行 × 6 列
  ri <- if (nrow(pre) > 25) sample(nrow(pre), 25) else seq_len(nrow(pre))
  ci <- if (length(cols) > 6) sample(cols, 6) else cols
  probes <- as.character(pre[[1]][ri])
  q <- dbGetQuery(con, sprintf("SELECT row_names, %s FROM `%s` WHERE row_names IN (%s)",
                               paste(backtick(ci), collapse = ", "), t,
                               paste0("'", gsub("'", "''", probes), "'", collapse = ",")))
  if (!nrow(q)) { note(sprintf("  %-20s 抽样探针在库中不存在，跳过", a)); next }
  idx <- match(pre[[1]], q$row_names)
  got <- q[idx[!is.na(idx)], , drop = FALSE]
  ok  <- !is.na(idx)
  pre_v  <- as.matrix(pre[ok,  ci, drop = FALSE]);  storage.mode(pre_v)  <- "double"
  post_v <- as.matrix(post[ok, ci, drop = FALSE]);  storage.mode(post_v) <- "double"
  db_v   <- as.matrix(got[, ci, drop = FALSE]);     storage.mode(db_v)   <- "double"
  # FLOAT 存储：用相对误差判断（1e-5 足够区分 1 与 log2(2)=1）
  rel <- function(x, y) { d <- abs(x - y); s <- pmax(abs(y), 1e-9); max(d / s) }
  e_pre <- rel(db_v, pre_v); e_post <- rel(db_v, post_v)
  verdict <- if (is.na(e_pre) || is.na(e_post)) "unclear"
             else if (e_post <= 1e-4 && e_post <= e_pre) "already"
             else if (e_pre <= 1e-4) "linear"
             else "mismatch"
  if (verdict == "linear") {
    if (apply_db) {
      sets <- paste0(backtick(cols), " = LOG2(", backtick(cols), " + 1)", collapse = ", ")
      dbExecute(con, sprintf("UPDATE `%s` SET %s", t, sets))
      note(sprintf("  %-20s 已 UPDATE（%d 列 x %d 行；抽样相对误差 pre=%.1e post=%.1e）",
                   a, length(cols), nrow(q), e_pre, e_post))
    } else {
      note(sprintf("  %-20s [待转化] 抽样误差 pre=%.1e post=%.1e", a, e_pre, e_post))
    }
  } else if (verdict == "already") {
    note(sprintf("  %-20s 已是 log2 尺度，跳过（post 误差=%.1e）", a, e_post))
  } else {
    note(sprintf("  %-20s !! 与本地都对不上（pre=%.1e post=%.1e），未修改", a, e_pre, e_post))
  }
}

dir.create(cpas_out_root(), showWarnings = FALSE)
writeLines(c("# 镜像端 log2 转化记录", "",
             paste0("- 时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             paste0("- 模式: ", ifelse(apply_db, "--apply", "dry-run")),
             paste0("- 数据集: ", length(accs), "（来源 pipeline/backup/expr_pre_log2/）"),
             "", "```text", log_lines, "```", ""),
           cpas_out("REPORT_log2_transform_db.md"))
cat("\n已写出 pipeline/out/REPORT_log2_transform_db.md\n")
if (!apply_db) cat("[dry-run] 未修改镜像。加 --apply 执行 UPDATE。\n")
}

# ---------------------------------------------------------------------------
# run_21_reupload_log2_tables()  <-  verbatim pipeline/R/21_reupload_log2_tables.R
# ---------------------------------------------------------------------------
run_21_reupload_log2_tables <- function() {
# 21_reupload_log2_tables.R --------------------------------------------------
# 把已做 log2 转化的本地表达矩阵重新上传，覆盖镜像中对应的表。
#
# 为什么不是就地 UPDATE
#   镜像里现有表是 MyISAM，且数据文件对 mysqld 进程只读，UPDATE 直接报
#   "Table 'X' is read only"（MySQL 1036）。实测：新建/删除表、写新表都正常，
#   所以走 06_upload_db.R 一贯的 dbWriteTable(overwrite=TRUE)（DROP + CREATE +
#   LOAD DATA LOCAL INFILE）即可，且新表由本账号创建，后续可写。
#
# 几何不变式
#   行（探针）与列（样本）严格按镜像现状：若镜像表是本地文件的子集
#   （GSE21034-GPL10264 / GSE84426 / GSE84433 属此类），只上传交集，行序按镜像；
#   两者列集合必须完全一致，否则跳过该表并报警。只替换数值，不改变表结构。
#
# 用法
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/21_reupload_log2_tables.R                # dry-run
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/21_reupload_log2_tables.R --only GSE4716-GPL3694
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/21_reupload_log2_tables.R --apply
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

.libPaths(c("/home/Jingle/R/library", .libPaths()))
suppressPackageStartupMessages(library(RMySQL))

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
apply_up <- "--apply" %in% args
only <- if ("--only" %in% args) args[which(args == "--only") + 1L] else NULL

bak_dir <- file.path(root, "pipeline/backup/expr_pre_log2")
accs <- sub("\\.rds$", "", list.files(bak_dir, pattern = "\\.rds$"))
if (!length(accs)) stop("备份目录为空：先跑 18_log2_transform.R --apply")
if (!is.null(only)) accs <- intersect(accs, only)
if (!length(accs)) stop("没有匹配的数据集。")

pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
if (!nzchar(pw)) stop("CPAS_DB_PASSWORD 未设置。")
con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = pw, client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)
tb <- dbGetQuery(con, "SHOW TABLES")[[1]]

log_lines <- character(0)
note <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
set.seed(7)

for (a in accs) {
  t <- gsub("-", "_", a, fixed = TRUE)
  if (!t %in% tb) { note(sprintf("  %-20s 镜像无此表，跳过", a)); next }
  x <- readRDS(file.path(root, "data/expr", paste0(a, ".rds")))
  db_probes <- dbGetQuery(con, sprintf("SELECT row_names FROM `%s`", t))$row_names
  db_cols   <- setdiff(dbListFields(con, t), "row_names")
  loc_probes <- as.character(x[[1]]); loc_cols <- names(x)[-1]

  if (!setequal(intersect(db_probes, loc_probes), db_probes)) {
    note(sprintf("  %-20s !! 镜像探针不在本地文件中，跳过（需人工核对）", a)); next }
  if (!setequal(db_cols, loc_cols)) {
    miss <- setdiff(db_cols, loc_cols)
    note(sprintf("  %-20s !! 样本列不一致（镜像多出 %d 列，如 %s），跳过",
                 a, length(miss), paste(head(miss, 3), collapse = ","))); next }

  # 行序按镜像、列序按镜像
  m <- x[match(db_probes, loc_probes), c(1L, match(db_cols, loc_cols) + 1L)]
  dd <- as.data.frame(m, check.names = FALSE)
  names(dd) <- c("ID_REF", db_cols)
  rownames(dd) <- NULL                     # tibble 不接受继承来的行名
  dd <- tibble::column_to_rownames(dd, "ID_REF")
  dd[] <- lapply(dd, as.numeric)
  rm(m); gc(verbose = FALSE)
  v <- as.numeric(unlist(dd[seq_len(min(5, nrow(dd))), seq_len(min(5, ncol(dd)))]))
  if (all(v == 0)) { note(sprintf("  %-20s 抽样全为 0，疑似异常，跳过", a)); next }

  # 幂等：镜像已是 log2 值（与本地转化后一致）则跳过
  pr0 <- db_probes[if (length(db_probes) > 20) sample(length(db_probes), 20) else seq_along(db_probes)]
  cc0 <- db_cols[if (length(db_cols) > 6) sample(length(db_cols), 6) else seq_along(db_cols)]
  q0 <- dbGetQuery(con, sprintf("SELECT row_names, %s FROM `%s` WHERE row_names IN (%s)",
                                paste0("`", cc0, "`", collapse = ","), t,
                                paste0("'", pr0, "'", collapse = ",")))
  if (nrow(q0)) {
    i0 <- match(pr0, q0$row_names)
    a0 <- as.matrix(dd[pr0, cc0, drop = FALSE]); b0 <- as.matrix(q0[i0, cc0, drop = FALSE])
    e0 <- max(abs(a0 - b0) / pmax(abs(a0), 1e-9))
    if (is.finite(e0) && e0 < 1e-4) {
      note(sprintf("  %-20s 镜像已是 log2 尺度，跳过（相对误差 %.1e）", a, e0)); next
    }
  }

  if (!apply_up) {
    note(sprintf("  %-20s [待上传] %d 探针 x %d 样本 | 抽样值 %.3f ... %.3f",
                 a, nrow(dd), ncol(dd), min(v), max(v)))
    next
  }
  field_types <- c(row_names = "VARCHAR(255)", setNames(rep("FLOAT", ncol(dd)), colnames(dd)))
  dbWriteTable(con, name = t, value = dd, overwrite = TRUE, row.names = TRUE,
               field.types = field_types)
  n_after <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", t))$n
  # 回读抽样校验
  pr <- db_probes[if (length(db_probes) > 20) sample(length(db_probes), 20) else seq_along(db_probes)]
  cc <- db_cols[if (length(db_cols) > 6) sample(length(db_cols), 6) else seq_along(db_cols)]
  q <- dbGetQuery(con, sprintf("SELECT row_names, %s FROM `%s` WHERE row_names IN (%s)",
                               paste0("`", cc, "`", collapse = ","), t,
                               paste0("'", pr, "'", collapse = ",")))
  idx <- match(pr, q$row_names)
  lit <- as.matrix(dd[pr, cc, drop = FALSE]); lod <- as.matrix(q[idx, cc, drop = FALSE])
  rel <- max(abs(lit - lod) / pmax(abs(lit), 1e-9))
  note(sprintf("  %-20s 已上传 %d 行 x %d 列 | 回读相对误差 %.2e %s",
               a, n_after, ncol(dd), rel, ifelse(rel < 1e-4, "OK", "!! 偏大")))
  rm(dd); gc(verbose = FALSE)
}

dir.create(cpas_out_root(), showWarnings = FALSE)
writeLines(c("# 镜像表达表重新上传（log2 后）", "",
             paste0("- 时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             paste0("- 模式: ", ifelse(apply_up, "--apply", "dry-run")),
             paste0("- 数据集: ", length(accs)), "", "```text", log_lines, "```", ""),
           cpas_out("REPORT_log2_reupload.md"))
cat("\n已写出 pipeline/out/REPORT_log2_reupload.md\n")
if (!apply_up) cat("[dry-run] 未修改镜像。加 --apply 执行上传。\n")
}

# ---------------------------------------------------------------------------
# run_22_report_log2_impact()  <-  verbatim pipeline/R/22_report_log2_impact.R
# ---------------------------------------------------------------------------
run_22_report_log2_impact <- function() {
# 22_report_log2_impact.R ----------------------------------------------------
# 量化 log2 转化对分析结果的影响：同一批样本、同一个基因，比较转化前后两种尺度。
#
# 思路：log2(x+1) 可逆，x = 2^y - 1。镜像现已存 log2 值 y，反推即得转化前的
# 线性值 x（差异仅 FLOAT 舍入）。于是在完全相同的样本集上可以算出：
#   转化前：每单位(线性强度) HR 与 每 SD HR
#   转化后：每单位(log2 单位) HR 与 每 SD HR
# p 值随尺度变化（log2 是非线性变换），这正是本次改动需要公示的部分。
#
# 用法: CPAS_DB_PASSWORD 不需要（走 API）；Rscript pipeline/R/22_report_log2_impact.R
# 输出: pipeline/out/log2_impact.csv、pipeline/out/REPORT_log2_impact.md
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

.libPaths(c("/home/Jingle/R/library", .libPaths()))
suppressPackageStartupMessages({library(CanPAS); library(survival)})

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
bak <- file.path(root, "pipeline/backup/expr_pre_log2")
accs <- sub("\\.rds$", "", list.files(bak, pattern = "\\.rds$"))

# 论文工作示例用到的基因（GAPDH/MKI67/TP53/BIRC5/AURKA + 文中点名的两个对照 PPIA/TBP）
genes_paper <- c("GAPDH", "MKI67", "TP53", "BIRC5", "AURKA", "PPIA", "TBP")
paper_cohorts <- c("GSE31210", "GSE37745", "GSE14814")

surv_of <- function(acc) {
  f <- file.path(root, "data/processed/surv", paste0(acc, "_surv.rds"))
  if (!file.exists(f)) return(NULL)
  sv <- readRDS(f)
  data.frame(ID = as.character(rownames(sv)), sv, check.names = FALSE)
}

fit_one <- function(time, status, x) {
  ok <- is.finite(time) & !is.na(status) & is.finite(x)
  f <- tryCatch(coxph(Surv(time, status) ~ x, subset = ok), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- summary(f)$coefficients
  s2 <- summary(coxph(Surv(time, status) ~ scale(x), subset = ok))$coefficients
  list(n = sum(ok), events = sum(status[ok] == 1),
       hr = exp(s[1, 1]), lo = exp(s[1, 1] - 1.96 * s[1, 3]), hi = exp(s[1, 1] + 1.96 * s[1, 3]),
       p = s[1, 5], sd = sd(x[ok]),
       hr_sd = exp(s2[1, 1]), lo_sd = exp(s2[1, 1] - 1.96 * s2[1, 3]), hi_sd = exp(s2[1, 1] + 1.96 * s2[1, 3]),
       p_sd = s2[1, 5])
}

# ---- 取数：探针集来自 API（与 get_expr_data 完全一致），
#      转化前的值直接从备份文件里按同一批探针取 max 合并 ----------------
fetch_post <- function(acc, gene) {
  for (k in 1:6) {
    x <- try(get_expr_data(acc, gene), silent = TRUE)
    if (!inherits(x, "try-error")) return(x)
    Sys.sleep(4 * k)                       # API 会 429，退避重试
  }
  NULL
}
probe_max_pre <- function(acc, gene, probes) {
  f <- file.path(bak, paste0(acc, ".rds"))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "double"
  rownames(m) <- as.character(x[[1]])
  p <- intersect(probes, rownames(m)); if (!length(p)) return(NULL)
  v <- apply(m[p, , drop = FALSE], 2, max, na.rm = TRUE)
  data.frame(ID = names(v), value = as.numeric(v), stringsAsFactors = FALSE)
}

rows <- list()
for (a in accs) {
  for (g in c("GAPDH", if (a %in% paper_cohorts) setdiff(genes_paper, "GAPDH") else character(0))) {
    # API 会对密集请求返回 429，重试若干次并留出间隔
    e <- fetch_post(a, g)
    if (is.null(e)) { cat("x"); next }
    # ref_ids 是 data.frame（列：gene_id, 原始探针 id, symbol, ensembl），第 2 列是探针
    rid <- e$ref_ids
    probes <- if (is.data.frame(rid)) as.character(rid[[2]]) else
              as.character(unlist(rid[[length(rid)]]))
    pre_series <- probe_max_pre(a, g, probes)
    if (is.null(pre_series)) { cat("x"); next }
    post_series <- e$expr_data
    names(post_series)[1] <- "ID"; names(post_series)[2] <- g
    pre_series$ID <- as.character(pre_series$ID)
    post_series$ID <- as.character(post_series$ID)
    dd <- merge(post_series, surv_of(a), by = "ID")
    dd_pre <- merge(pre_series, surv_of(a), by = "ID")
    if (!nrow(dd) || !nrow(dd_pre)) { cat("x"); next }
    tok <- tryCatch(endpoint_resolve(a, "OS"), error = function(e) NA_character_)
    if (is.na(tok) || !all(c(paste0(tok, "_time"), paste0(tok, "_status"), g) %in% names(dd))) next
    tm <- dd[[paste0(tok, "_time")]]; st <- dd[[paste0(tok, "_status")]]
    # 按 ID 对齐后再比较
    k <- match(dd$ID, dd_pre$ID)
    y <- dd[[g]]; x <- dd_pre[[g]][k]
    pre  <- fit_one(tm, st, x)
    post <- fit_one(tm, st, y)
    if (is.null(pre) || is.null(post)) next
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = a, gene = g, endpoint = tok, n = pre$n, events = pre$events,
      sd_pre = pre$sd, sd_post = post$sd,
      hr_pre_unit = pre$hr, hr_post_unit = post$hr,
      hr_pre_sd = pre$hr_sd, lo_pre_sd = pre$lo_sd, hi_pre_sd = pre$hi_sd, p_pre = pre$p,
      hr_post_sd = post$hr_sd, lo_post_sd = post$lo_sd, hi_post_sd = post$hi_sd, p_post = post$p,
      stringsAsFactors = FALSE)
  }
  cat(".")
}
res <- do.call(rbind, rows)
dir.create(cpas_out_root(), showWarnings = FALSE)
write.csv(res, cpas_out("log2_impact.csv"), row.names = FALSE)
cat("\n")

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
md <- c("# log2 转化对结果的影响（同一批样本，反推对比）", "",
        paste0("- 时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0("- 覆盖: ", length(accs), " 个转化数据集；论文工作示例队列额外算 ",
               paste(setdiff(genes_paper, "GAPDH"), collapse = "/")),
        "- 每单位 HR 在转化前是「每 1 线性强度单位」，转化后是「每 1 log2 单位」；每 SD HR 两者都可比",
        "- p 值随尺度改变（log2 是非线性变换）", "",
        "## GAPDH（论文中的阴性对照基因）", "",
        "| dataset | n | events | 每单位 HR (前) | 每单位 HR (后) | 每 SD HR (前) | 每 SD HR (后) | p (前) | p (后) |",
        "|---|---|---|---|---|---|---|---|---|")
g1 <- res[res$gene == "GAPDH", ]
md <- c(md, sprintf("| %s | %d | %d | %s | %s | %s (%s-%s) | %s (%s-%s) | %.3g | %.3g |",
                    g1$dataset, g1$n, g1$events,
                    fmt(g1$hr_pre_unit, 4), fmt(g1$hr_post_unit, 4),
                    fmt(g1$hr_pre_sd), fmt(g1$lo_pre_sd), fmt(g1$hi_pre_sd),
                    fmt(g1$hr_post_sd), fmt(g1$lo_post_sd), fmt(g1$hi_post_sd),
                    g1$p_pre, g1$p_post))
for (a in intersect(paper_cohorts, unique(res$dataset))) {
  sub <- res[res$dataset == a & res$gene != "GAPDH", ]
  if (!nrow(sub)) next
  md <- c(md, "", paste0("## 论文工作示例队列 ", a, "（OS，", sub$n[1], " 例 / ", sub$events[1], " 事件）"), "",
          "| gene | 每 SD HR (前) | p (前) | 每 SD HR (后) | p (后) |",
          "|---|---|---|---|---|",
          sprintf("| %s | %s (%s-%s) | %.3g | %s (%s-%s) | %.3g |", sub$gene,
                  fmt(sub$hr_pre_sd), fmt(sub$lo_pre_sd), fmt(sub$hi_pre_sd), sub$p_pre,
                  fmt(sub$hr_post_sd), fmt(sub$lo_post_sd), fmt(sub$hi_post_sd), sub$p_post))
}
writeLines(md, cpas_out("REPORT_log2_impact.md"))
cat("已写出 pipeline/out/log2_impact.csv 与 REPORT_log2_impact.md\n")
cat("汇总：GAPDH 每 SD HR 中位变化倍数 =",
    round(median(g1$hr_post_sd / g1$hr_pre_sd, na.rm = TRUE), 3), "\n")
}
