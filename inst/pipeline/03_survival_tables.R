# ===========================================================================
# CanPAS curation pipeline -- 03_survival_tables
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 03_surv_table.R
#   2. 14_split_gse40272_surv.R
#   3. 15_split_gse40272_expr.R
#   4. 17_recod_gse86166_status.R
#   5. 23_fix_surv_ids_and_endpoints.R
#   6. 25_fix_endpoint_annotations.R
#   7. 35_add_tcga_cohorts.R
#   8. 39_build_gse1379_surv.R
#   9. 40_add_tcga_chol_dlbc.R
#   10. 110_build_gse205209.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "03_survival_tables.R", package = "CanPAS"))
#   run_03_surv_table()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================


# ---------------------------------------------------------------------------
# run_03_surv_table()  <-  verbatim pipeline/R/03_surv_table.R
# ---------------------------------------------------------------------------
run_03_surv_table <- function() {
# 03_surv_table.R  (v2 — 加入 DB 统一格式规范)
# -----------------------------------------------------------------------------
# CanPAS 数据处理 第 3 步：由 pheno 生成标准化 生存/临床 表
#   输出 data/processed/surv/<ACC>_surv.rds
#     行 = 样本（行名 = 样本 GSM，与既有 *_surv.rds 一致）
#     生存列 {OS|DSS|DFS|PFS|RFS|MFS}_status(0/1) 与 _{同型}_time(年)
#     临床列统一规范（对齐 DB 既有 surv 表的格式）：
#       histology : 肿瘤亚型/病理分型统一入 histology 列（见 spec$histology）
#       stage     : 罗马数字总体分期 I/II/III/IV（去除 A/B/C 亚期后缀，阿拉伯转罗马）
#       T / N / M : 大写 "T1/N0/M1/X" 形式（去除细亚字母）
#       grade     : 罗马数字 I–IV（1-4 与 G1-G4 转换；低/高级别文本保留原样）
#       age       : 数值（岁）
#       sex       : male/female
# 说明：列的值级清洗/转换发生在 rename 之后、以 spec$clinic 的目标列名为准；
#       仅对 spec$norm 中列出的列做规范化，其余临床列原样保留（空格→"_"）。
# -----------------------------------------------------------------------------
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

suppressMessages({library(dplyr); library(stringr)})

ROOT <- if (length(commandArgs(trailingOnly=TRUE)) >= 1) commandArgs(trailingOnly=TRUE)[1] else "~/data/Project/CanPAS"
ROOT <- path.expand(ROOT)
setwd(ROOT)

DROP_COLS <- c("title","status","submission_date","last_update_date","type",
               "channel_count","source_name_ch1","organism_ch1","molecule_ch1",
               "extract_protocol_ch1","label_ch1","label_protocol_ch1","taxid_ch1",
               "hyb_protocol","scan_protocol","description","data_processing",
               "platform_id","contact_name","contact_email","contact_phone",
               "contact_fax","contact_laboratory","contact_department","contact_institute",
               "contact_address","contact_city","contact_zip/postal_code","contact_country",
               "supplementary_file","data_row_count","treatment_protocol_ch1",
               "growth_protocol_ch1","relation","source_name","molecule")

# GSE12417：每样本一“描述列”（列名=整段临床文本，仅对角有值）。折叠成结构化列。
collapse_aml_wide <- function(d) {
  single <- colSums(!is.na(d)) <= 1L
  single <- single & vapply(d, is.character, logical(1))
  if (!any(single)) return(d)
  wide_names <- names(which(single))
  out <- d[, !single, drop=FALSE]
  desc <- apply(d[wide_names], 1, function(row) {
    i <- which(!is.na(row) & row != "")
    if (!length(i)) NA_character_ else wide_names[i[1]]
  })
  out$aml_desc <- desc
  # 解析描述文本
  age  <- suppressWarnings(as.numeric(sub(".*age =([0-9.]+) years.*", "\\1", desc)))
  osd  <- suppressWarnings(as.numeric(sub(".*OS = ([0-9.]+) days.*", "\\1", desc)))
  fab  <- sub(".*FAB ([^;]*);.*", "\\1", desc)
  fab[grepl("^AML", fab) | is.na(fab)] <- NA_character_
  status <- rep(NA_real_, nrow(d))
  for (i in seq_along(desc)) {
    if (!is.na(desc[i])) {
      v <- d[i, wide_names[match(desc[i], wide_names)]]
      if (!is.na(v)) status[i] <- suppressWarnings(as.numeric(v))
    }
  }
  out$aml_age <- age; out$aml_os_days <- osd; out$aml_fab <- fab
  out$aml_status <- status
  out
}

# 将单列中 ";" 分隔的 "key: value" 串展开为多列（如 GSE14333 的 Location 列）
expand_semicolon <- function(d, src_col) {
  if (!src_col %in% colnames(d)) { warning("expand col missing: ", src_col); return(d) }
  v <- .as_char(d[[src_col]])
  rows <- lapply(v, function(s) {
    if (is.na(s) || !nzchar(s) || s == "NA") return(list())
    parts <- trimws(strsplit(s, ";")[[1]])
    kv <- lapply(parts, function(p) {
      sp <- strsplit(p, ":", fixed=TRUE)[[1]]
      if (length(sp) >= 2) return(setNames(trimws(paste(sp[-1], collapse=":")), trimws(sp[1])))
      setNames(trimws(p), src_col)          # 无冒号 token（首段）归属源列名
    })
    unlist(kv, recursive=FALSE)
  })
  allk <- unique(unlist(lapply(rows, names)))
  d <- d[, setdiff(colnames(d), src_col), drop=FALSE]
  for (k in allk) d[[k]] <- NA_character_
  for (i in seq_along(rows))
    for (k in names(rows[[i]])) d[[i, k]] <- rows[[i]][[k]]
  d
}

# ------------------------- 值级规范化函数 -------------------------
.as_char <- function(x) trimws(as.character(x))
.roman <- c("0"="0","1"="I","2"="II","3"="III","4"="IV","5"="V","X"="X","x"="X")

std_stage <- function(x) {
  v <- toupper(.as_char(x)); v[v %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  vapply(v, function(s) {
    if (is.na(s)) return(NA_character_)
    s <- str_remove(s, "^STAGE")
    s <- str_remove_all(s, "[^0-9IVX]")          # 只保留数字/罗马字母
    if (!nzchar(s)) return(NA_character_)
    if (grepl("^[0-9]+$", s))                    # 纯阿拉伯（2/3B 等先截到整体数字）
      return(.roman[[substr(s, 1, 1)]] %||% NA_character_)
    # 罗马：取最长的合法前缀（IIIA->III, IA2->I, IVB->IV）
    m <- regmatches(s, regexpr("^(IV|III|II|I|V)", s))
    if (length(m) && nzchar(m)) return(m) else NA_character_
  }, character(1), USE.NAMES=FALSE)
}

std_tnm <- function(x, letter) {
  v <- toupper(.as_char(x)); v[v %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  v <- sub("^P", "", v)                 # pT3/pN1 -> T3/N1
  vapply(v, function(s) {
    if (is.na(s)) return(NA_character_)
    m <- regmatches(s, regexpr(paste0("^", letter, "[0-9X]+"), s))   # T1b->T1；X 保留
    if (length(m) && nzchar(m)) return(substr(m, 1, 2))
    if (grepl("^[0-9X]+$", s)) return(paste0(letter, s))
    NA_character_
  }, character(1), USE.NAMES=FALSE)
}

std_grade <- function(x) {
  v <- toupper(.as_char(x)); v[v %in% c("", "NA", "N/A", "NULL", "?", "UNKNOWN")] <- NA_character_
  vapply(v, function(s) {
    if (is.na(s)) return(NA_character_)
    s2 <- str_remove(s, "^G")
    if (grepl("^[1-4]$", s2)) return(.roman[[s2]])
    if (grepl("^[1-4]$", s))  return(.roman[[s]])
    if (s %in% c("I","II","III","IV","LOW","HIGH","LOW GRADE","HIGH GRADE")) return(s)
    if (grepl("INDET|INDETERMINATE|^4=", s)) return(NA_character_)
    # 长文本（如 "2 - Moderately differentiated (MD)"）取第一个数字 1-4
    m <- regmatches(s, regexpr("[1-4]", s))
    if (length(m) && nzchar(m)) return(.roman[[m]])
    if (grepl("WELL", s))  return("I")
    if (grepl("MODERATE", s)) return("II")
    if (grepl("POOR", s)) return("III")
    NA_character_
  }, character(1), USE.NAMES=FALSE)
}

std_age <- function(x) {
  v <- suppressWarnings(as.numeric(.as_char(x)))
  ifelse(is.na(v), NA_real_, v)
}

std_sex <- function(x) {
  v <- tolower(.as_char(x))
  vapply(v, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    if (s %in% c("m","male","man","men")) return("male")
    if (s %in% c("f","female","woman","women")) return("female")
    s
  }, character(1), USE.NAMES=FALSE)
}

std_num <- function(x) suppressWarnings(as.numeric(.as_char(x)))

date_months <- function(start, end) {          # 日期差 -> 月
  st <- suppressWarnings(as.Date(.as_char(start)))
  en <- suppressWarnings(as.Date(.as_char(end)))
  out <- as.numeric(difftime(en, st, units="days")) / 30.4375
  ifelse(is.na(st) | is.na(en), NA_real_, out)
}

status_to_binary <- function(x, invert=FALSE) {
  out <- dplyr::case_when(
    tolower(as.character(x)) %in% c("dead","deceased","d","yes","y","relapse","relapsed",
                                    "metastasis","1","death","event","recurrence",
                                    "recurred","dead of disease","true","t") ~ "1",
    tolower(as.character(x)) %in% c("alive","live","no","a","n","0","not relapsed",
                                    "survival","censor","censored","no event",
                                    "no death","no recurrence","alive/no death","false","f") ~ "0",
    TRUE ~ as.character(x)
  )
  out <- suppressWarnings(as.numeric(out))
  if (invert) out <- ifelse(is.na(out), NA_real_, 1 - out)   # 标记列为删失指示(1=删失)时反转
  out
}

`%||%` <- function(a, b) if (is.null(a) || length(a)==0 || is.na(a)) b else a

# ------------------------- 数据集配置（人工确认过的映射） -------------------------
specs <- list(
  # ===== cervical cancer =====
  GSE44001 = list(
    surv   = list(DFS = list(status="status_of_dfs",
                             time="disease_free_survival_(dfs)_(months)", unit="month")),
    clinic = c("Stage"="stage", "largest diameter"="tumor_size"),
    norm   = c("stage"="stage", "tumor_size"="num"),
    note   = "cervical n=300; FIGO stage 已按 DB 规范转罗马总体分期(IB1/IA2/IB2->I, IIA->II)"
  ),
  # ===== colorectal cancer =====
  GSE39582 = list(
    surv   = list(OS  = list(status="os.event",  time="os.delay (months)",  unit="month"),
                  RFS = list(status="rfs.event", time="rfs.delay",          unit="month")),
    clinic = c("Sex"="sex", "age.at.diagnosis (year)"="age",
               "tnm.stage"="stage", "tnm.t"="T", "tnm.n"="N", "tnm.m"="M",
               "tumor.location"="location", "mmr.status"="mmr",
               "chemotherapy.adjuvant"="chemo"),
    histology = list(sources = "cit.molecularsubtype"),   # C1-C6 分子亚型 → histology
    norm   = c("sex"="sex", "age"="age", "stage"="stage",
               "T"="T", "N"="N", "M"="M"),
    keep   = quote(!`dataset` %in% c("Non Tumoral")),     # 剔除 19 例正常黏膜
    note   = "CRC 566 tumor (443 discovery + 123 validation)；OS/RFS 单位月；stage/tnm 数值按 DB 规范统一(T/N/M大写、stage罗马)；histology=CRC C1-C6 亚型"
  ),
  # ---- Moffitt/Vanderbilt CRC 队列（GSE17536 发现集 / GSE17537 验证集 / GSE17538 扩展）----
  GSE17536 = list(
    surv = list(
      OS  = list(status="overall_event (death from any cause)", time="overall survival follow-up time", unit="month"),
      DSS = list(status="dss_event (disease specific survival; death from cancer)", time="dss_time", unit="month"),
      DFS = list(status="dfs_event (disease free survival; cancer recurrence)", time="dfs_time", unit="month")),
    clinic = c("gender"="sex", "age"="age", "ajcc_stage"="stage", "grade"="grade"),
    norm = c("sex"="sex", "age"="age", "stage"="stage", "grade"="grade"),
    note = "CRC MCC cohort n=177; OS/DSS/DFS 单位月; grade文本取首个数字转罗马"
  ),
  GSE17537 = list(
    surv = list(
      OS  = list(status="overall_event (death from any cause)", time="overall survival follow-up time", unit="month"),
      DFS = list(status="dfs_event (disease free survival; cancer recurrence)", time="dfs_time", unit="month")),
    clinic = c("gender"="sex", "age"="age", "ajcc_stage"="stage", "grade"="grade"),
    norm = c("sex"="sex", "age"="age", "stage"="stage", "grade"="grade"),
    note = "CRC VMC validation n=55 (<79, DB 惯例不上传 expr, 只上传 surv)"
  ),
  "GSE17538-GPL570" = list(
    surv = list(
      OS  = list(status="overall_event (death from any cause)", time="overall survival follow-up time", unit="month"),
      DSS = list(status="dss_event (disease specific survival; death from cancer)", time="dss_time", unit="month"),
      DFS = list(status="dfs_event (disease free survival; cancer recurrence)", time="dfs_time", unit="month")),
    clinic = c("gender"="sex", "age"="age", "ajcc_stage"="stage", "grade"="grade",
               "location"="location"),
    norm = c("sex"="sex", "age"="age", "stage"="stage", "grade"="grade"),
    note = "CRC MCC/VMC combined n=238 human (GPL570); 仅保留人类部分 (GPL1261 小鼠 6 例不入库)"
  ),
  GSE14333 = list(
    expand_semicolon = "Location",
    surv = list(DFS = list(status="DFS_Cens", time="DFS_Time", unit="month", invert=TRUE)),
    clinic = c("Location"="location", "Age_Diag"="age", "Gender"="sex",
               "DukesStage"="dukes", "AdjXRT"="adj_xrt", "AdjCTX"="adj_ctx"),
    norm = c("age"="age", "sex"="sex"),
    note = "CRC n=290 (Jorissen); 特征为 ';' 拼接需展开; DFS_Cens=1 表示删失(invert)；DFS_Time 单位月；Dukes A-D 保留原值"
  ),
  GSE87211 = list(
    surv = list(
      OS  = list(status="death due to tumor", time="survival time (month)", unit="month"),
      DFS = list(status="cancer recurrance after surgery", time="disease free time (month)", unit="month")),
    clinic = c("gender"="sex", "age"="age", "kras mutation"="kras",
               "depth of invasion after rct"="T", "lymph node metastasis after rct"="N",
               "metastasis after rct"="M"),
    fold_mut = c("kras mutation"),
    norm = c("sex"="sex", "age"="age", "T"="T", "N"="N", "M"="M"),
    keep = quote(tissue == "rectal tumor"),     # 剔除 160 例正常黏膜
    note = "rectal cancer (GPL13497) n_tumor=203；术后(yp)TNM；OS/DFS 单位月；kras 变体折叠为 M/WT(另保留原值不可) "
  ),
  GSE103479 = list(
    surv = list(
      OS  = list(status="status alive.dead", time="overall survival time", unit="month"),
      PFS = list(status="recurrence", time="progression-free survival time", unit="month")),
    clinic = c("gender"="sex", "age diagnosis"="age", "stage ii\\iii"="stage",
               "tstage"="T", "nstage"="N", "mstage"="M",
               "tumour site"="location", "differentiation"="grade",
               "tumour_size"="tumour_size", "lvi"="lvi", "evi"="evi",
               "perineural_invasion"="perineural_invasion",
               "positive lymph nodes"="nodes_pos", "dukes"="dukes",
               "curative resection"="curative", "adjuvanttreatment"="adjuvant",
               "cris subgroup"="cris", "braf (v600e)"="braf",
               "kras (codons 12, 13, 61)"="kras", "p53"="p53"),
    histology = list(sources = "cms subgroup"),   # CMS1-4 分子亚型 → histology（UNK->NA）
    norm = c("sex"="sex", "age"="age", "stage"="stage",
             "T"="T", "N"="N", "M"="M", "grade"="grade"),
    note = "CRC n=156 (FFPE); OS/PFS 单位月; histology=CMS 亚型; stage IIA-IIIC 去亚期转罗马; pT/pN 去 p"
  ),
  GSE39084 = list(
    surv = list(
      OS  = list(kind="date", status="death", start="surgery.date",
                 end="lastnews or death.date", unit="month"),
      DFS = list(kind="date", status="relapse", start="surgery.date",
                 end_event="relapse.date", end_censor="lastnews or death.date", unit="month")),
    clinic = c("gender"="sex", "age.year"="age", "tnm.stage"="stage",
               "tumor.location"="location", "lynch.syndrom"="lynch",
               "msi.status [mismatch repair status]"="msi", "cimp.status"="cimp",
               "tp53.gene.mutation.status"="tp53", "kras.gene.mutation.status"="kras",
               "braf.gene.mutation.status"="braf", "pik3ca.gene.mutation.status"="pik3ca"),
    histology = list(sources = "histology"),
    norm = c("sex"="sex", "age"="age", "stage"="stage"),
    note = "CRC n=70; OS/DFS 由手术/随访/复发日期推算(月)；N<79 expr 不入库"
  ),
  GSE28735 = list(
    surv = list(OS = list(status="cancer_death", time="survival_month", unit="month")),
    clinic = c(),
    keep = quote(tissue == "T"),          # 45 例肿瘤（去 45 例癌旁 N）
    note = "Pancreatic paired n_tumor~45 (GPL6244); OS month; surv-only(N<79)"
  ),
  GSE62452 = list(
    surv = list(OS = list(status="survival status", time="survival months", unit="month")),
    clinic = c("Stage"="stage", "grading"="grade"),
    norm = c("stage"="stage", "grade"="grade"),
    keep = quote(tissue == "Pancreatic tumor"),   # 69 例肿瘤（去 61 例癌旁）
    note = "Pancreatic n_tumor=69 (GPL6244); OS month; surv-only(N<79)"
  ),
  "GSE37642-GPL570" = list(
    surv = list(OS = list(status="life status", time="overall survival (days)", unit="day")),
    clinic = c("age"="age", "fab"="fab",
               "runx1-runx1t1_fusion"="runx1_fusion", "runx1_mutation"="runx1_mutation"),
    note = "AML (MAQC-style) GSE37642-GPL570; OS day→年; 含 age/FAB/runx1"
  ),
  "GSE37642-GPL96" = list(
    surv = list(OS = list(status="life status", time="overall survival (days)", unit="day")),
    clinic = c("age"="age", "fab"="fab",
               "runx1-runx1t1_fusion"="runx1_fusion", "runx1_mutation"="runx1_mutation"),
    note = "AML (MAQC-style) GSE37642-GPL96; OS day→年; 含 age/FAB/runx1"
  ),
  "GSE37642-GPL97" = list(
    extra = cpas_out("GSE37642-GPL97_extra.rds"),
    surv = list(OS = list(status="OS_status_x", time="OS_time_x", unit="year")),
    clinic = c("age"="age", "fab"="fab"),
    note = "AML GSE37642-GPL97：同批患者按 title 对齐，OS 由 GPL96 部分引入（年）"
  ),
  "GSE12417-GPL96" = list(
    collapse_aml_wide = TRUE,
    surv = list(OS = list(status="aml_status", time="aml_os_days", unit="day")),
    clinic = c("aml_age"="age", "aml_fab"="fab"),
    note = "AML normal-karyotype GSE12417-GPL96; 临床以每样本文本列存储，已折叠解析(OS day/status/age/FAB)"
  ),
  "GSE12417-GPL97" = list(
    collapse_aml_wide = TRUE,
    surv = list(OS = list(status="aml_status", time="aml_os_days", unit="day")),
    clinic = c("aml_age"="age", "aml_fab"="fab"),
    note = "AML normal-karyotype GSE12417-GPL97; 临床以每样本文本列存储，已折叠解析(OS day/status/age/FAB)"
  ),
  "GSE12417-GPL570" = list(
    collapse_aml_wide = TRUE,
    surv = list(OS = list(status="aml_status", time="aml_os_days", unit="day")),
    clinic = c("aml_age"="age", "aml_fab"="fab"),
    note = "AML normal-karyotype GSE12417-GPL570; 临床以每样本文本列存储，已折叠解析(OS day/status/age/FAB)"
  ),
  "GSE9782-GPL96" = list(
    surv = list(
      OS  = list(status="Did_Patient_Die(0=No,1=Yes)", time="Days_Survived_From_Randomization",
                 unit="day", na_censor=TRUE),
      PFS = list(status="PGx_Progression(0=No,1=Yes)", time="PGx_Days_To_Progression",
                 unit="day", na_censor=TRUE)),
    clinic = c("sex"="sex", "Age_at_Randomization"="age", "TC_Class_2006"="tc_class"),
    histology = list(value = "MM"),
    norm = c("sex"="sex", "age"="age"),
    note = "MM GSE9782-GPL96 (n=264); OS/PFS 天; 事件列 1=事件、空=删失(na_censor); histology=MM"
  ),
  "GSE9782-GPL97" = list(
    surv = list(
      OS  = list(status="Did_Patient_Die(0=No,1=Yes)", time="Days_Survived_From_Randomization",
                 unit="day", na_censor=TRUE),
      PFS = list(status="PGx_Progression(0=No,1=Yes)", time="PGx_Days_To_Progression",
                 unit="day", na_censor=TRUE)),
    clinic = c("sex"="sex", "Age_at_Randomization"="age", "TC_Class_2006"="tc_class"),
    histology = list(value = "MM"),
    norm = c("sex"="sex", "age"="age"),
    note = "MM GSE9782-GPL97 (n=264); 同批患者 OS/PFS 天"
  ),
  GSE4581 = list(
    surv = list(OS = list(status="SURIND", time="SURTIM", unit="month")),
    histology = list(value = "MM"),
    note = "MM UAMS-TT2 n=414 (GPL570); OS: SURIND(0/1)+SURTIM(月); histology=MM"
  ),
  GSE102238 = list(
    surv = list(OS = list(status="survival status(alive=0,dead=1)",
                          time="survival time(day)", unit="day")),
    clinic = c("t stage"="T",
               "n stage（”n0=0,n1=1)"="N",
               "m stage(m0=0.m1=1)"="M",
               "age（≤65=0，＞65=1）"="age_group",
               "gender（male=0，female=1）"="sex_code",
               "differentiation(well/moderate=0,poor=1)"="differentiation",
               "vessel invasion(absence=0,present=1)"="vessel_invasion",
               "localization of tumor(head=1,body/tail=2)"="tumor_location",
               "tumor size(≤3=0，＞3=1）"="tumor_size_group"),
    histology = list(value = "PDAC"),
    norm = c("T"="T", "N"="N", "M"="M"),
    note = paste0("PDAC n=96 (GPL19072, AnnoProbe 重注释); OS day->年(内部 /365，与其余队列一致)。",
                  "★配对设计：本 series 为 50 例患者的肿瘤/癌旁配对(共100样本)，survival 逐样本重复，",
                  "故 96 = 48 例患者 x2、events 60 = 30 例死亡 x2；title/source_name 区分 tumor/normal，",
                  "本队列未按组织过滤(与候选扫描 effective_n=96 对齐)。",
                  "age_group/sex_code/tumor_size_group 为原文 0/1 编码(age 0=≤65; sex 0=male,1=female; size 0=≤3cm)")
  ),
  GSE39055 = list(
    surv = list(RFS = list(status="recurrence", time="time until first recurrence or latest follow-up (months)", unit="month")),
    clinic = c("death"="death", "age"="age", "gender"="sex"),
    norm = c("sex"="sex", "age"="age"),
    note = "osteosarcoma n=37 (GPL14951, surv-only); RFS: recurrence Y/N + 月；death 无 OS 时间列故不单独建 OS"
  ),
  GSE21257 = list(
    surv = list(OS = list(status="status", unit="month", text_parse=TRUE)),
    clinic = c("histological subtype"="histology", "age"="age", "gender"="sex", "huvos grade"="huvos_grade"),
    norm = c("age"="age", "sex"="sex"),
    note = "osteosarcoma n=53 (GPL10295, surv-only); OS 由 status 文本(Alive/Deceased at N months)解析"
  ),
  "GSE4412-GPL96" = list(
    surv = list(OS = list(status="living", time="survival time", unit="day")),
    clinic = c("histology"="histology", "gender"="sex", "age"="age", "grade"="grade"),
    norm = c("sex"="sex", "age"="age", "grade"="grade"),
    note = "glioma n=85 (GPL96, >79); OS day; histology ASTRO/GBM/MIXED/OLIGO"
  ),
  "GSE4412-GPL97" = list(
    surv = list(OS = list(status="living", time="survival time", unit="day")),
    clinic = c("histology"="histology", "gender"="sex", "age"="age", "grade"="grade"),
    norm = c("sex"="sex", "age"="age", "grade"="grade"),
    note = "glioma n=85 (GPL97, >79); 同批患者 OS day"
  ),
  "GSE14520-GPL3921" = list(
    extra = cpas_out("GSE14520_GPL3921_extra.rds"),
    surv = list(
      OS  = list(status="Survival_status_x", time="Survival_months_x", unit="month"),
      RFS = list(status="Recurr_status_x",   time="Recurr_months_x",   unit="month")),
    clinic = c("sex_x"="sex", "age_x"="age", "stage_x"="stage"),
    norm = c("sex"="sex", "age"="age", "stage"="stage"),
    keep = quote(tissue_x == "Tumor"),     # 肿瘤 247（Extra 官方临床）
    note = "HCC Shanghai GSE14520-GPL3921 n_expr=445/n_tumor~221(带OS/RFS); OS/RFS 月, 临床来自官方 Extra_Supplement(按 Affy_GSM join)"
  ),
  GSE31519 = list(
    surv = list(RFS = list(status="event (1", regex="[01]$",
                           time="event-free interval (months)", unit="month")),
    clinic = c(),
    note = "Breast n=67 (<79 surv-only); event-free interval(月)作为 RFS(事件尾部 0/1 正则)；无可靠临床列(列名被冒号截断)"
  ),
  GSE7378 = list(
    surv = list(DFS = list(status="recurrence (1=event,_0=censored)", time="DFS_years", unit="year",
                           na_censor=TRUE)),
    clinic = c("age"="age"),
    norm = c("age"="age"),
    note = "Breast n=54 (<79 surv-only); DFS_years; recurrence 列为 1/空白(空白视为删失 0, 说明见 GEO 列注释)"
  ),
  GSE9195 = list(
    surv = list(
      RFS = list(status="e.rfs", time="t.rfs", unit="day"),
      MFS = list(status="e.dmfs", time="t.dmfs", unit="day")),
    clinic = c("age"="age", "er"="ER", "grade"="grade", "node"="node", "size"="size"),
    norm = c("age"="age", "ER"="pn", "grade"="grade", "node"="num", "size"="num"),
    note = "Breast n=77 (<79 surv-only); RFS/MFS time 天数(210-4124d)"
  ),
  GSE84426 = list(
    surv = list(OS = list(status="death", time="duration overall survival", unit="month")),
    clinic = c("age"="age", "Sex"="sex", "pnstage"="N", "ptstage"="T"),
    norm = c("age"="age", "sex"="sex", "N"="N", "T"="T"),
    note = "Gastric n=76 (<79 surv-only); OS 月; pN/pT 规范化"
  ),
  GSE46602 = list(
    surv = list(RFS = list(status="bcr", time="bcr_free_time", unit="month")),
    clinic = c("age"="age", "gleason_grade"="gleason", "t-stage"="T",
               "preop_psa"="psa", "margin_status"="margin"),
    norm = c("age"="age", "gleason"="num", "T"="T", "psa"="num"),
    note = "Prostate n=50 (<79 surv-only); RFS=BCR(生化复发 YES/NO, 月)"
  ),
  GSE25307 = list(
    surv = list(OS = list(status="osbin", time="os", unit="year")),
    clinic = c("er"="ER", "pgr"="PGR", "grade"="grade",
               "pam50 classification"="subtype"),
    histology = list(value = "BC"),
    norm = c("ER"="posneg", "PGR"="posneg", "grade"="grade"),
    keep = quote(tissue == "Breast tumor"),   # 566 肿瘤（剔除 11 例非恶性）
    note = "Breast (BRCA/familial cohort, SWEGENE GPL5345) n_tumor=566; OS 年; PAM50→subtype; histology=BC"
  ),
  "GSE4716-GPL3694" = list(
    surv = list(OS = list(status="status after 5 years", time="Survival Period", unit="month")),
    clinic = c("HIST"="histology", "AGE"="age", "SEX"="sex",
               "pStage"="stage", "pT"="T", "pN"="N"),
    norm = c("age"="age", "sex"="sex", "stage"="stage", "T"="T", "N"="N"),
    note = "Lung 50 (GeneFilter GF200, surv-only N<79); OS 5年(cap 60月)；HIST AD/SQ/LA→histology"
  ),
  "GSE4716-GPL3696" = list(
    surv = list(OS = list(status="status after 5 years", time="Survival Period", unit="month")),
    clinic = c("HIST"="histology", "AGE"="age", "SEX"="sex",
               "pStage"="stage", "pT"="T", "pN"="N"),
    norm = c("age"="age", "sex"="sex", "stage"="stage", "T"="T", "N"="N"),
    note = "Lung 50 (GeneFilter GF201, surv-only N<79); OS 5年(cap 60月)"
  ),
  GSE25066 = list(
    surv = list(DRFS = list(status="drfs_1_event_0_censored",
                            time="drfs_even_time_years", unit="year")),
    clinic = c("age_years"="age", "er_status_ihc"="ER", "pr_status_ihc"="PR",
               "her2_status"="HER2", "grade"="grade",
               "clinical_ajcc_stage"="stage", "clinical_t_stage"="T",
               "clinical_nodal_status"="N",
               "pathologic_response_pcr_rd"="pcr", "pam50_class"="subtype"),
    histology = list(value = "BC"),
    norm = c("age"="age", "ER"="pn", "PR"="pn", "HER2"="pn",
             "grade"="grade", "stage"="stage", "T"="T", "N"="N"),
    note = "Breast neoadjuvant n=508 (GPL96); DRFS time YEARS; ER/PR/HER2 P/N/I→1/0/NA; histology=BC; subtype=PAM50"
  ),
  GSE76427 = list(
    surv = list(
      OS  = list(status="event_os", time="duryears_os", unit="year"),
      RFS = list(status="event_rfs", time="duryears_rfs", unit="year")),
    clinic = c("gender (1=m, 2=f)"="sex", "age (years)"="age",
               "tnm_staging_clinical"="stage", "bclc_staging"="bclc"),
    histology = list(value = "HCC"),
    norm = c("sex"="mf_num", "age"="age", "stage"="stage"),
    keep = quote(!grepl("^adjacent", tissue)),   # 剔除 52 例癌旁组织
    note = "HCC n_tumor=115 (GPL10558); OS/RFS duryears(年)；gender 1/2→male/female；histology=HCC"
  ),
  GSE71014 = list(
    surv = list(OS = list(status="event (1", regex="[01]$",
                          time="overall survival (months)", unit="month")),
    clinic = c(),
    note = "CN-AML n=104 (GPL10558); OS month；事件列名含冒号被截断，用尾字符正则提取 0/1"
  ),
  GSE57495 = list(
    surv = list(OS = list(status="vital.status", time="overall survival (month)", unit="month")),
    clinic = c("Stage"="stage"),
    histology = list(value = "PA"),
    norm = c("stage"="stage"),
    note = "Pancreatic PDAC n=63 (GPL15048); OS month; N<79 expr 不入库; histology=PA"
  ),
  GSE24080 = list(
    extra = cpas_out("GSE24080_extra.rds"),
    surv = list(
      OS  = list(status="OS_status_x", time="OS_time_x", unit="month"),
      EFS = list(status="EFS_status_x", time="EFS_time_x", unit="month")),
    clinic = c("Sex"="sex", "age"="age"),
    histology = list(value = "MM"),
    norm = c("sex"="sex", "age"="age"),
    note = "Multiple myeloma UAMS565 n=559 (GPL570); OS/EFS time(month) 来自 GEO 补充 ClinInfo xls，按 CEL title 与 GSM 精确 join；histology=MM 常量"
  ),
  GSE52903 = list(excluded = TRUE,
    note = "EXCLUDED from CanPAS (user decision 2026-09-06): no OS in GEO, only original paper PMID 24879114. raw file kept."
  ),
  GSE53625 = list(
    surv = list(OS = list(status="death at fu", time="survival time(months)", unit="month")),
    clinic = c("age"="age", "Sex"="sex", "tumor grade"="grade",
               "t stage"="T", "n stage"="N", "tnm stage"="stage",
               "tumor loation"="tumor_location", "tobacco use"="tobacco_use",
               "alcohol use"="alcohol_use", "adjuvant therapy"="adjuvant_therapy",
               "patient id"="patient_id"),
    histology = list(value = "ESCC"),
    norm = c("age"="age", "sex"="sex", "grade"="grade", "T"="T", "N"="N", "stage"="stage"),
    note = paste0("ESCC n=358 (GPL18109 Agilent-038314 FeatureNumber 版, AnnoProbe 重注释); ",
                  "OS: death at fu(yes/no) + survival time(months) -> 年 (/12)。",
                  "★配对设计：179 例患者的 癌/癌旁 配对(共358样本)，survival 逐样本重复，",
                  "故 358 = 179 例患者 x2、events 212 = 106 例死亡 x2；tissue 列区分 cancer/normal，",
                  "本队列未按组织过滤(与候选扫描 effective_n=358 对齐)；patient_id 可还原配对")
  ),
  # ==========================================================================
  # GEO 扩展批次（2026-09-24，A 审计后）：给"只有 TCGA"的癌种补 GEO 队列
  #   口径与逐队列依据见 pipeline/out/geo_expansion_AUDIT.md 与
  #   pipeline/out/geo_expansion_build_log.md；N/n_events 一律用审计的患者级数字。
  # ==========================================================================
  # ---------------------------- Head and Neck Cancer ------------------------
  GSE65858 = list(
    surv = list(PFS = list(status="pfs_event", time="pfs", unit="day")),
    clinic = c("age"="age", "gender"="sex", "tumor_site"="location", "uicc_stage"="stage",
               "t_category"="T", "n_category"="N", "hpv_dna"="hpv",
               "tumor_type"="tumor_type", "consensus_cluster"="consensus_cluster",
               "treatment"="treatment", "smoking"="smoking", "packyears"="packyears"),
    norm = c("age"="age", "sex"="sex", "stage"="stage", "T"="T", "N"="N", "packyears"="num"),
    note = paste0("HNSCC n=270 (GPL10558, Illumina HT-12 v4); PFS: pfs_event(TRUE/FALSE) + ",
                  "pfs(天)->年；审计 effective_n=270/133(events)。同一份特征里另有 os/os_event",
                  "(270/94) 与 distant_metastasis：本队列按审计口径只登记 PFS（避免把 ",
                  "EndpointPrimary 改成 OS 而与审计的 N/n_events 不一致）。")
  ),
  GSE117973 = list(
    surv = list(PFS = list(status="pfs event", time="pfs months", unit="month")),
    clinic = c("age (years)"="age", "Sex"="sex", "subsite"="location",
               "pathological grading"="grade", "t status"="T", "n status"="N",
               "r status"="resection_status", "smoker"="smoker", "alcohol"="alcohol",
               "hpv infection"="hpv", "therapy"="therapy",
               "tumor cell content"="tumor_cell_content"),
    norm = c("age"="age", "sex"="sex", "grade"="grade", "T"="T", "N"="N",
             "tumor_cell_content"="num"),
    note = paste0("HNSCC n=77 (GPL10558); PFS: pfs event(1/0)+pfs months->年；审计 ",
                  "effective_n=77/22。另有 dss event/dss months(77/16) 未登记(同 PFS 理由)。")
  ),
  GSE27020 = list(
    surv = list(DFS = list(status="dfs status (1 = recurred)", time="dfs (months)",
                           unit="month")),
    clinic = c("age"="age", "grade"="grade", "group"="group", "tissue"="tissue"),
    norm = c("age"="age", "grade"="grade"),
    note = paste0("Laryngeal (HNSCC) n=109 (GPL96); DFS: dfs status(1=recurred)+dfs(months)->年；",
                  "审计 effective_n=109/34；group=training/test。")
  ),
  GSE159067 = list(
    surv = list(PFS = list(status="pfs (event)", time="pfs (month)", unit="month")),
    clinic = c("site"="location", "age"="age", "Sex"="sex",
               "tobacco status"="tobacco", "alcohol status"="alcohol",
               "hpv status (0 or 1 or 2 if unknown)"="hpv",
               "immunotherapy line"="io_line",
               "best response on immunotherapy (recist)"="recist",
               "hot score"="hot_score", "hot phenotype"="hot_phenotype"),
    norm = c("age"="age", "sex"="sex"),
    note = paste0("Advanced HNSCC (PD-1/PD-L1 治疗) n=102 (GPL18573, HTG EdgeSeq OBP ",
                  "~2.5k 基因 panel)；PFS: pfs (event)+pfs (month)->年；审计 ",
                  "effective_n=102/96。表达不在 series matrix(!Sample_type=SRA)，来自补充文件 ",
                  "GSE159067_IHN_log2cpm_data.txt.gz(已是 log2CPM)，按 !Sample_description 连接；",
                  "★定向 panel(非全转录组)；另有 os (102/92) 未登记。")
  ),
  GSE162520 = list(
    surv = list(PFS = list(status="pfs (event)", time="pfs (month)", unit="month")),
    clinic = c("patient diagnosis"="diagnosis", "tissue"="tissue", "age"="age",
               "Sex"="sex", "hot score"="hot_score", "hot phenotype"="hot_phenotype",
               "histology"="histology"),
    norm = c("age"="age", "sex"="sex"),
    note = paste0("NSCLC (PD-1/PD-L1 治疗) n=92 (GPL18573, HTG EdgeSeq OBP panel)；PFS: ",
                  "pfs (event)+pfs (month)->年；审计 effective_n=92/35。★癌种纠正：审计/任务把它",
                  "列在 Head and Neck 类(检索式命中 HNSCC 字样)，但 GEO 标题、overall design 与",
                  "每例 `patient diagnosis` 都是 non-small cell lung cancer，故 Type=Lung Cancer。",
                  "表达来自补充文件 GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz(分号分隔、逗号",
                  "小数点、已是 log2CPM)。")
  ),
  # ------------------------------- Melanoma ---------------------------------
  GSE65904 = list(
    surv = list(DSS = list(status="disease specific survival (1=death, 0=alive)",
                           time="disease specific survival in days", unit="day")),
    clinic = c("age"="age", "gender"="sex", "tumor stage"="stage", "tissue"="tissue"),
    norm = c("age"="age", "sex"="sex", "stage"="stage"),
    note = paste0("Cutaneous melanoma n=214 (4 例临床缺失) (GPL10558); DSS: 1=death/0=alive + ",
                  "days->年；审计 effective_n=210/102。另有 DMFS(150/83) 未登记。")
  ),
  GSE198430 = list(
    surv = list(OS = list(status="dead or alive at last contact", time="os", unit="month")),
    clinic = c("gender"="sex", "age at diagnosis"="age",
               "stage at original diagnosis"="stage", "tnm"="tnm",
               "breslow thickness"="breslow", "ajcc thickness"="ajcc_thickness",
               "ulceration"="ulceration", "cohort"="cohort", "tag"="tissue"),
    norm = c("sex"="sex", "age"="age", "stage"="stage",
             "breslow"="num", "ajcc_thickness"="num"),
    note = paste0("Early-stage melanoma n=124 (GPL32055 NanoString nCounter, 194 行定向 panel); ",
                  "OS: dead or alive at last contact + os(月)->年；审计 effective_n=105/50。",
                  "★定向 panel：表达谱只含 panel 基因，非全转录组(审计 §6 caveat)。")
  ),
  GSE198431 = list(
    surv = list(OS = list(status="dead or alive at last contact", time="os", unit="month")),
    clinic = c("gender"="sex", "age at diagnosis"="age",
               "stage at original diagnosis"="stage", "tnm"="tnm",
               "breslow thickness"="breslow", "ajcc thickness"="ajcc_thickness",
               "ulceration"="ulceration", "cohort"="cohort", "tag"="tissue"),
    norm = c("sex"="sex", "age"="age", "stage"="stage",
             "breslow"="num", "ajcc_thickness"="num"),
    note = paste0("Early-stage melanoma n=82 (GPL32055, 同 GSE198430 的 panel, ",
                  "FFPE 编号不重叠、GSM 重叠 0)；OS 同上；审计 effective_n=79/30。",
                  "审计 §3.13：原候选表把它误判为 GSE198430 的 SAME_STUDY_DUPLICATE。")
  ),
  GSE22153 = list(
    surv = list(OS = list(status="event (0=alive, 1=dead)", time="os (days)", unit="day")),
    clinic = c("age at primary diagnosis"="age", "sex"="sex", "stage"="stage",
               "breslow"="breslow", "clark"="clark",
               "localization of primary melanoma"="location",
               "molecular subtype"="subtype", "braf/nras"="braf_nras",
               "type of metastases"="metastasis_type", "ki67 (0=<30%, 1=>30%)"="ki67"),
    norm = c("age"="age", "sex"="sex", "stage"="stage", "breslow"="num",
             "clark"="num", "ki67"="num"),
    note = paste0("Stage IV melanoma n=57 (GPL6102); OS: event(1=dead)+os(days)->年；",
                  "审计 effective_n=54/47。GSE22155 SuperSeries 的 test 集，只取 GSE22153。")
  ),
  GSE325123 = list(
    surv = list(OS = list(status="vital_status", time="overall_survival_years", unit="year")),
    clinic = c("distant_metastasis"="distant_metastasis",
               "lymph_node_status"="lymph_node_status",
               "acral_melanoma"="acral_melanoma", "tissue"="tissue", "n_roi"="n_roi"),
    norm = c("n_roi"="num"),
    note = paste0("Melanoma (long-term follow-up) n=105 患者 (GPL24676, GeoMx DSP WTA); ",
                  "★已按 patient 把 170 个 included in_analysis ROI 聚合(表达取 ROI 均值)，",
                  "3 例无 OS 时间 -> 可分析 102、events 60；审计 ROI 级 269/168 是重复计数。",
                  "OS 单位年(原键 overall survival(years))。")
  ),
  # -------------------------------- Sarcoma ---------------------------------
  GSE71118 = list(
    surv = list(MFS = list(status="metastasis", time="time", unit="year")),
    clinic = c("cinsarc"="cinsarc", "material support"="material_support",
               "paired rna-seq"="paired_rnaseq"),
    note = paste0("Various sarcomas (CINSARC validation) n=312 (GPL570); MFS: metastasis(Yes/No) ",
                  "+ time(年)->年；审计 effective_n=312/124。time 单位=年(最大 15，若按月则随访仅 ",
                  "15 个月，与 124/312 转移事件率不符)。")
  ),
  GSE30929 = list(
    surv = list(MFS = list(status="drfs", time="tt.drfs", unit="month")),
    clinic = c("sample id"="sample_id", "trainingtest"="training_test", "tissue"="tissue"),
    histology = list(sources = "subtype"),
    note = paste0("Liposarcoma n=140 (GPL96); DRFS: drfs(TRUE/FALSE)+tt.drfs(月)->年，",
                  "MFS 家族(DRFS->MFS，见 11_endpoint_families.R)；审计 effective_n=140/49；",
                  "histology=subtype(well-dedifferentiated/myxoid/pleomorphic/...)。")
  ),
  GSE271517 = list(
    surv = list(OS = list(status="overal survival", time="overall survival_time_(days)",
                          unit="day")),
    clinic = c("patient id"="patient_id", "fusion gene"="fusion_gene",
               "tumor type"="tumor_type", "site"="site", "tumor size_(cm)"="tumor_size",
               "poorly differentiated_histology"="poorly_diff",
               "histology"="histology"),
    norm = c("tumor_size"="num"),
    note = paste0("Synovial sarcoma n=91 样本 / 55 例有 OS (GPL24676, RNA-seq 计数); ",
                  "OS: overal survival(0/1, 原键拼写如此) + overall survival_time_(days)->年；",
                  "审计 effective_n=55/24(仅原发瘤有 OS，转移瘤 NA)。表达为计数，已按流水线",
                  "惯例做 log2(x+1)。")
  ),
  # ------------------------------- Lymphoma ---------------------------------
  GSE31312 = list(
    extra = cpas_out("GSE31312_extra.rds"),
    surv = list(OS  = list(status="pdf_os_censor", time="pdf_os_months", unit="month",
                           invert=TRUE),
                PFS = list(status="pdf_pfs_censor", time="pdf_pfs_months", unit="month",
                           invert=TRUE)),
    clinic = c("pdf_gep"="gep", "pdf_markers"="markers"),
    note = paste0("DLBCL n=498 (GPL570)；生存不在 series matrix，来自补充 PDF 的 Table 2",
                  "(Case#/Date death/PFScensor/OScensor/PFS/OS/.../GEO Depository #)，",
                  "按 GEO Depository # = !Sample_title 连接；censor 列 1=删失故 invert。",
                  "PDF 475 例 / 172 死亡；其中 470 例能与矩阵样本对上(5 个 ID 在 GEO 中不存在)，",
                  "N=470/OS events 及 PFS 见 n_* 列；见 93_build_geo_expansion_bespoke.R。")
  ),
  GSE32918 = list(
    extra = cpas_out("GSE32918_extra.rds"),
    surv = list(OS = list(status="follow-up status", time="follow-up years", unit="year")),
    clinic = c("age"="age", "Sex"="sex", "predicted class"="subtype",
               "class confidence"="class_confidence"),
    norm = c("age"="age", "sex"="sex", "class_confidence"="num"),
    keep = quote(patient_dedup_keep),
    note = paste0("DLBCL (DASL, GPL8432) 249 张芯片来自 172 例患者(_Rep1.._Rep13 重复)；",
                  "★已按患者去重(keep=patient_dedup_keep，无 _Rep 后缀优先)，OS: follow-up ",
                  "status(Dead/Alive)+follow-up years->年；审计患者级 effective_n=172/93",
                  "(芯片级 249/137 是重复计数)。★与 GSE69051 是同一批患者的第二次 deposit：",
                  "只建本队列，GSE69051 记入 pipeline/out/excluded.csv 的重复说明，不另建。")
  ),
  GSE23501 = list(
    surv = list(OS = list(status="code_os", time="overall_survival (years)", unit="year")),
    clinic = c("age"="age", "gender"="sex", "ipi_stage"="stage",
               "molecular_subtype_by_gene_expression"="subtype", "pathology"="pathology",
               "txprimary"="treatment"),
    norm = c("age"="age", "sex"="sex", "stage"="stage"),
    note = paste0("DLBCL n=69 (GPL570)；OS: code_os(alive/dead) + overall_survival (years)->年；",
                  "审计 effective_n=69/13。★原候选脚本用的是 `status` 列(临床疗效码，14 个“事件”",
                  "里 11 个不是死亡)，只有 code_os 才是生死；另有 code_pfs + progression-free_",
                  "survival (years) (69/17) 未登记。")
  ),
  GSE248835 = list(
    surv = list(DFS = list(status="event.free.survival.event",
                           time="event.free.survival.months", unit="month")),
    clinic = c("cell of origin"="cell_of_origin", "treatment arm"="treatment_arm",
               "histologically.proven.dlbcl.group"="histology_group",
               "baseline tumor burden (spd)"="tumor_burden", "grade3_crs"="grade3_crs",
               "grade3_ne"="grade3_ne", "duration.of.response.months"="dor_months",
               "duration.of.response.event"="dor_event", "ongoing.response"="response",
               "ongoing_2grps"="response_group", "visit"="visit"),
    norm = c("tumor_burden"="num", "dor_months"="num"),
    note = paste0("LBCL (CAR-T vs chemo) n=256 (GPL33963 NanoString nCounter IO360, 817 行); ",
                  "★终点是 EFS(event.free.survival.*)，按 11_endpoint_families.R 的 EFS->DFS ",
                  "家族登记为 DFS(审计 §3.9：原候选表按 key 里含 survival 误标成 OS)；",
                  "审计 effective_n=256/179。")
  ),
  # ----------------------------- Uveal Melanoma -----------------------------
  GSE22138 = list(
    surv = list(MFS = list(status="metastasis", time="months to endpoint", unit="month")),
    clinic = c("age"="age", "gender"="sex", "eye"="eye",
               "tumor location"="tumor_location", "tumor diameter (mm)"="tumor_diameter",
               "tumor thickness (mm)"="tumor_thickness",
               "chromosome 3 status"="chr3_status",
               "extrascleral extension"="extrascleral",
               "retinal detachment"="retinal_detachment"),
    histology = list(sources = "tumor cell type"),
    norm = c("age"="age", "sex"="sex", "tumor_diameter"="num", "tumor_thickness"="num"),
    note = paste0("Uveal melanoma 原发瘤 n=63 (GPL570)；MFS: metastasis(yes/no) + months to ",
                  "endpoint(月)->年；审计 effective_n=63/35。★新癌种(Uveal Melanoma)；",
                  "原候选表 MANUAL note 说“GEO 里没有生存列”是错的(审计 §7)。")
  ),
  # ------------------------------ Mesothelioma ------------------------------
  GSE183088 = list(
    surv = list(OS = list(status="status", time="survival (months)", unit="month")),
    clinic = c("age"="age", "gender"="sex", "histology"="histology", "tnm"="tnm",
               "abestos exposure"="asbestos_exposure"),
    norm = c("age"="age", "sex"="sex"),
    note = paste0("Malignant pleural mesothelioma n=86 (GPL30570 NanoString nCounter, 70 行 ",
                  "定向 panel)；OS: status(alive/dead)+survival (months)->年；审计 ",
                  "effective_n=86/75(原候选表因稀疏事件袋报成 1/1)。★70 基因 panel：按任务",
                  "要求先确认其基因可用于本流水线的 schema(GPL30570 68/70 行可映射 ENTREZ)，",
                  "故建库；但表达只含 panel 基因，非全转录组。")
  )
)

# ------------------------- 生成 surv 表 -------------------------
build_surv <- function(acc, spec) {
  if (isTRUE(spec$excluded)) { message(acc, ": excluded -> skip"); return(invisible(NULL)) }
  pheno <- readRDS(file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
  d <- pheno %>%
    dplyr::select(-dplyr::any_of(DROP_COLS)) %>%
    dplyr::rename_with(~ str_remove(., ":ch1$"))
  d <- d[, !duplicated(colnames(d))]
  if (!is.null(spec$expand_semicolon))                       # 先展开 ";k:v" 混合列
    d <- expand_semicolon(d, spec$expand_semicolon)
  if (isTRUE(spec$collapse_aml_wide)) d <- collapse_aml_wide(d)
  if (!is.null(spec$extra) && file.exists(spec$extra)) {     # 外接临床/端点表（行名=GSM）
    ex <- readRDS(spec$extra)
    if (!all(rownames(d) %in% rownames(ex))) stop(acc, ": extra table misses samples")
    d <- cbind(d, ex[rownames(d), , drop=FALSE])
  }
  if (!is.null(spec$keep)) {
    n0 <- nrow(d); d <- d[eval(spec$keep, envir=d), , drop=FALSE]
    message(acc, ": row filter kept ", nrow(d), "/", n0)
  }
  # 突变类字段折叠为 WT / M（如 kras 的氨基酸变体）
  for (fm_col in spec$fold_mut %||% character(0)) {
    if (fm_col %in% colnames(d)) {
      v <- .as_char(d[[fm_col]])
      miss <- is.na(v) | v %in% c("", "NA", "N/A")
      d[[fm_col]] <- ifelse(miss, NA_character_,
                            ifelse(v %in% c("WT", "wt"), "WT", "M"))
    }
  }
  cols_out <- character(0)

  drop_date_cols <- character(0)
  if (!is.null(spec$surv)) {
    for (st in names(spec$surv)) {
      s <- spec$surv[[st]]
      if (isTRUE(s$kind == "date")) {          # 日期型：start/end(或 end_event/end_censor)
        need <- c(s$status, s$start)
        if (is.null(s$end)) need <- c(need, s$end_event, s$end_censor)
        if (!all(need %in% colnames(d))) {
          warning(acc, ": date-surv cols missing (", st, ")"); next }
        d[[paste0(st, "_status")]] <- status_to_binary(d[[s$status]])
        if (!is.null(s$end)) {
          t_m <- date_months(d[[s$start]], d[[s$end]])
        } else {
          ev <- d[[paste0(st, "_status")]]
          end <- ifelse(is.na(ev), NA_character_,
                        ifelse(ev == 1, d[[s$end_event]], d[[s$end_censor]]))
          t_m <- date_months(d[[s$start]], end)
        }
        t_m <- ifelse(is.na(t_m), NA_real_, ifelse(t_m < 0, 0, t_m))
        d[[paste0(st, "_time")]] <- switch(s$unit, year=t_m, month=t_m/12,
                                           day=t_m/365, t_m)
        drop_date_cols <- unique(c(drop_date_cols, s$status, s$start,
                                   if (is.null(s$end)) c(s$end_event, s$end_censor) else s$end))
        cols_out <- c(cols_out, paste0(st, "_status"), paste0(st, "_time"))
        next
      }
      if (!isTRUE(s$text_parse) && (!s$status %in% colnames(d) || !s$time %in% colnames(d))) {
        warning(acc, ": survival cols missing (", st, ") -> ", s$status, " / ", s$time); next }
      if (isTRUE(s$text_parse)) {
        rawv <- .as_char(d[[s$status]])
        d[[paste0(st, "_status")]] <- ifelse(grepl("dead|deceased|died", rawv, ignore.case=TRUE), 1,
                                             ifelse(is.na(rawv) | rawv %in% c("","NA"), NA_real_, 0))
        tnum <- suppressWarnings(as.numeric(sub(".*?([0-9.]+).*", "\\1", rawv)))
        d[[paste0(st, "_time")]] <- switch(s$unit, year=tnum, month=tnum/12, day=tnum/365, tnum)
        d[[s$status]] <- NULL
        cols_out <- c(cols_out, paste0(st, "_status"), paste0(st, "_time"))
        next
      } else if (!is.null(s$regex)) {
        rawv <- .as_char(d[[s$status]])
        d[[paste0(st, "_status")]] <- suppressWarnings(
          as.numeric(sub(paste0(".*(", s$regex, ")"), "\\1", rawv)))
      } else
      d[[paste0(st, "_status")]] <- status_to_binary(d[[s$status]],
                                                     invert = isTRUE(s$invert))
      if (isTRUE(s$na_censor)) d[[paste0(st, "_status")]][is.na(d[[paste0(st, "_status")]])] <- 0
      tcol <- suppressWarnings(as.numeric(d[[s$time]]))
      tcol <- ifelse(is.na(tcol), NA_real_, ifelse(tcol < 0, 0, tcol))
      d[[paste0(st, "_time")]] <- switch(s$unit, year=tcol, month=tcol/12, day=tcol/365, tcol)
      d[[s$status]] <- NULL; d[[s$time]] <- NULL
      cols_out <- c(cols_out, paste0(st, "_status"), paste0(st, "_time"))
    }
    d[drop_date_cols] <- NULL
  } else {
    message(acc, ": no survival endpoint -> surv table NOT created (", spec$no_surv_reason, ")")
    return(invisible(NULL))
  }

  # --- histology：肿瘤亚型/病理分型列统一（可合并多个来源列；无来源可用 spec$histology$value 常量） ---
  hs <- spec$histology
  if (!is.null(hs)) {
    if (!is.null(hs$sources) && length(hs$sources)) {
      src <- hs$sources[hs$sources %in% colnames(d)]
      if (length(src)) {
        h <- rep(NA_character_, nrow(d))
        for (col in src) {          # 按 sources 顺序取第一个非空值
          v <- .as_char(d[[col]])
          hit <- is.na(h) & !is.na(v) & v != "" & v != "NA"
          h[hit] <- v[hit]
        }
        d[["histology"]] <- h
        d[src] <- NULL
        cols_out <- c(cols_out, "histology")
      } else warning(acc, ": histology source cols missing: ",
                     paste(hs$sources, collapse=","))
    } else if (!is.null(hs$value)) {
      d[["histology"]] <- hs$value
      cols_out <- c(cols_out, "histology")
    }
  }

  # --- 临床列 rename + 规范化 ---
  clinic_map <- spec$clinic[spec$clinic != ""]
  for (src in names(clinic_map)) {
    if (src %in% colnames(d)) {
      newname <- clinic_map[[src]]
      if (newname != src) colnames(d)[match(src, colnames(d))] <- newname
      cols_out <- c(cols_out, newname)
    } else warning(acc, ": clinic col missing: ", src)
  }
  norm_map <- spec$norm  # 目标列名 -> 规范化类型
  for (nm2 in names(norm_map)) {
    if (nm2 %in% colnames(d)) {
      d[[nm2]] <- switch(norm_map[[nm2]],
        stage = std_stage(d[[nm2]]), T = std_tnm(d[[nm2]],"T"),
        N = std_tnm(d[[nm2]],"N"), M = std_tnm(d[[nm2]],"M"),
        grade = std_grade(d[[nm2]]), age = std_age(d[[nm2]]),
        sex = std_sex(d[[nm2]]), num = std_num(d[[nm2]]),
        mf_num = { z <- suppressWarnings(as.numeric(.as_char(d[[nm2]])));
                   ifelse(is.na(z), NA_character_, ifelse(z == 1, "male", "female")) },
        pn = { z <- toupper(.as_char(d[[nm2]]));
               ifelse(z %in% c("P","POS","POSITIVE","+","1"), "1",
                      ifelse(z %in% c("N","NEG","NEGATIVE","-","0"), "0", NA_character_)) },
        posneg = { z <- tolower(.as_char(d[[nm2]]));
                   ifelse(grepl("pos", z), "1", ifelse(grepl("neg", z), "0", NA_character_)) },
        d[[nm2]])
    } else warning(acc, ": norm col missing: ", nm2)
  }

  rn <- rownames(d)
  d <- as.data.frame(lapply(d, function(col) {
    if (is.character(col)) {
      col[col %in% c("", "NA", "N/A", "n/a", "N/A ", "NULL", "---", "-", "UNK", "UNKNOWN")] <- NA_character_
    }
    col
  }), check.names=FALSE, stringsAsFactors=FALSE)
  if (nrow(d) == length(rn)) rownames(d) <- rn
  colnames(d) <- gsub(" ", "_", colnames(d))
  cols_out <- intersect(cols_out, colnames(d))
  d <- d[, cols_out, drop=FALSE]
  out_file <- file.path(ROOT, "data/processed/surv", paste0(acc, "_surv.rds"))
  saveRDS(d, out_file)
  message("saved ", out_file, " (", nrow(d), " x ", ncol(d), ") cols: ",
          paste(colnames(d), collapse=", "))
  invisible(d)
}

# 可选过滤：Rscript 03_surv_table.R <ROOT> [<ACC> ...]
# 不带 ACC 时行为与以前完全一致（重建全部 specs）；指定 ACC 时只重建这些队列 ——
# 新增队列用，避免覆盖 14/15/17/23/25 等脚本对既有 surv 表做的就地修订。
ACC_FILTER <- if (length(commandArgs(trailingOnly=TRUE)) >= 2) commandArgs(trailingOnly=TRUE)[-1] else NULL
todo <- if (length(ACC_FILTER)) intersect(ACC_FILTER, names(specs)) else names(specs)
message("03_surv_table: building ", length(todo), " dataset(s)",
        if (length(ACC_FILTER)) paste0(" (filter: ", paste(ACC_FILTER, collapse=","), ")") else "")
for (acc in todo) build_surv(acc, specs[[acc]])
}

# ---------------------------------------------------------------------------
# run_14_split_gse40272_surv()  <-  verbatim pipeline/R/14_split_gse40272_surv.R
# ---------------------------------------------------------------------------
run_14_split_gse40272_surv <- function() {
# 14_split_gse40272_surv.R ---------------------------------------------------
# GSE40272(前列腺癌)的生存数据只以 GSE 级文件存在,而 catalog 按 4 个平台拆行且未标注终点,
# 导致这 4 行在任何家族筛选/分析中都不可用。本脚本按 GEO series matrix 的
# !Sample_geo_accession 把 GSE 级生存表拆到 4 个平台(样本互不重叠,并集 = 全部 84 例),
# 生成 data/processed/surv/GSE40272-<GPL>_surv.rds 并写入 DB 表 GSE40272_<GPL>_surv。
#
# 新表 schema 对齐既有 GSE40272_surv(row_names, OS_status, OS_time, DFS_status, DFS_time,
# RFS_status, RFS_time, age, T, N, M, PSA, histology);本地列 M0 -> M,Pre-PSA -> PSA。
# 终点选择(依据本文件先做的删失结构核查):
#   OS  : 1 例事件 / 81 例删失 -> 事件太少,不标注;
#   RFS : 19 例事件,但 63 例删失的 RFS_time 全为 0 -> 没有删失随访时间,不能作为生存终点;
#   DFS : 17 例事件,65 例删失全部有随访时间(中位 3.07 年) -> 标注 DFS 家族。
#
# 输出: 4 个本地 rds + 4 张 DB 表 + pipeline/out/db_gse40272_split_log.csv
# 用法: Rscript pipeline/R/14_split_gse40272_surv.R [--dry-run]
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

root <- "/home/Jingle/data/Project/CPAS"
dry  <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)
gpls <- c("GPL15971", "GPL15972", "GPL15973", "GPL9497")

surv <- readRDS(file.path(root, "data/processed/surv/GSE40272_surv.rds"))

# 参考 schema:既有 GSE 级表的列顺序与类型
con <- dbConnect(MySQL(), host = "139.224.80.159", dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = Sys.getenv("CPAS_DB_PASSWORD"))
sch <- dbGetQuery(con, "SHOW COLUMNS FROM GSE40272_surv")
cat("reference schema (GSE40272_surv):\n"); print(sch[, c("Field", "Type")])

local_map <- c(M = "M0", PSA = "Pre-PSA")   # DB 列 -> 本地列

log_rows <- list()
for (gpl in gpls) {
  raw <- file.path(root, "data/raw", sprintf("GSE40272-%s.gz", gpl))
  if (!file.exists(raw)) { log_rows[[length(log_rows)+1]] <- data.frame(Accession=paste0("GSE40272-",gpl), n_samples=NA, n_with_surv=NA, n_analyzable_DFS=NA, table="", status="missing raw file"); next }
  con2 <- gzfile(raw, "rt")
  ln <- grep("^!Sample_geo_accession", readLines(con2, n = 400, warn = FALSE), value = TRUE)[1]
  close(con2)
  ids <- unique(regmatches(ln, gregexpr("GSM[0-9]+", ln))[[1]])
  ov <- intersect(ids, rownames(surv))
  d <- surv[ov, , drop = FALSE]
  out <- data.frame(row_names = rownames(d), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(sch))) {
    f <- sch$Field[i]
    if (f == "row_names") next
    src <- if (f %in% names(d)) f else if (f %in% names(local_map)) local_map[[f]] else NA_character_
    v <- if (!is.na(src) && src %in% names(d)) d[[src]] else rep(NA, nrow(d))
    out[[f]] <- if (grepl("double|float|int", sch$Type[i])) suppressWarnings(as.numeric(v)) else as.character(v)
  }
  stem <- sprintf("GSE40272-%s", gpl)
  fpath <- file.path(root, "data/processed/surv", paste0(stem, "_surv.rds"))
  m <- out[, -1, drop = FALSE]
  rownames(m) <- out$row_names
  if (!dry) saveRDS(m, fpath)
  tbl <- paste0(gsub("-", "_", stem), "_surv")
  ok <- "dry-run"
  if (!dry) {
    ok <- tryCatch({
      dbWriteTable(con, tbl, out, overwrite = TRUE, row.names = FALSE)
      TRUE
    }, error = function(e) conditionMessage(e))
    ok <- if (isTRUE(ok)) "ok" else ok
  }
  n_an <- sum(!is.na(out$DFS_time) & out$DFS_time > 0 & !is.na(out$DFS_status))
  log_rows[[length(log_rows)+1]] <- data.frame(
    Accession = stem, n_samples = length(ids), n_with_surv = length(ov),
    n_analyzable_DFS = n_an, table = tbl, status = ok, stringsAsFactors = FALSE)
  message(sprintf("%s: expr_samples=%d surv_patients=%d analyzable_DFS(n & time>0)=%d -> %s",
                  stem, length(ids), length(ov), n_an, tbl))
}
dbDisconnect(con)

res <- do.call(rbind, log_rows)
print(res, row.names = FALSE)
write.csv(res, cpas_out("db_gse40272_split_log.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# run_15_split_gse40272_expr()  <-  verbatim pipeline/R/15_split_gse40272_expr.R
# ---------------------------------------------------------------------------
run_15_split_gse40272_expr <- function() {
# 15_split_gse40272_expr.R ---------------------------------------------------
# GSE40272 的表达只以 GSE 级文件 data/expr/GSE40272.rds(44,544 探针 x 153 样本)存在,
# 而 catalog 按 4 个平台拆行。本脚本按 GEO series matrix 的样本归属把表达矩阵拆到 4 个平台,
# 生成 data/expr/GSE40272-<GPL>.rds(与其它队列同格式:ID_REF + 每样本一列),并写入 DB 表。
#
#   探针号: 该研究的 ID_REF 为数值索引(1..44544),与 GPL15971/2/3 注释表行名一致;
#           GPL9497 注释表只覆盖其中 40,800 个,故该行按平台白名单过滤。
#   schema: row_names varchar(255) + 每样本 float,ENGINE=MyISAM(与既有表达表一致)。
#
# 用法: Rscript pipeline/R/15_split_gse40272_expr.R [--dry-run]
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

root <- "/home/Jingle/data/Project/CPAS"
dry  <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)
gpls <- c("GPL15971", "GPL15972", "GPL15973", "GPL9497")

expr <- readRDS(file.path(root, "data/expr", "GSE40272.rds"))
ids_all <- as.character(expr[[1]])

con <- dbConnect(MySQL(), host = "139.224.80.159", dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                 password = Sys.getenv("CPAS_DB_PASSWORD"), client_flag = CLIENT_COMPRESS)
on.exit(dbDisconnect(con), add = TRUE)

log_rows <- list()
for (gpl in gpls) {
  raw <- file.path(root, "data/raw", sprintf("GSE40272-%s.gz", gpl))
  if (!file.exists(raw)) { log_rows[[length(log_rows)+1]] <- data.frame(Accession=paste0("GSE40272-",gpl), probes=NA, samples=NA, status="missing raw file"); next }
  cn <- gzfile(raw, "rt", encoding = "latin1")
  ln <- readLines(cn, warn = FALSE); close(cn)
  i <- grep("^!series_matrix_table_begin", ln); j <- grep("^!series_matrix_table_end", ln)
  hdr <- gsub('^"|"$', "", strsplit(ln[i + 1], "\t")[[1]])
  body <- ln[(i + 2):(j - 1)]
  smp <- hdr[-1]
  keep <- smp %in% names(expr)                       # 只保留出现在合并矩阵里的样本
  first <- vapply(strsplit(body, "\t"), function(x) x[1], character(1))
  m <- matrix(NA_real_, nrow = length(body), ncol = sum(keep),
              dimnames = list(first, smp[keep]))
  for (k in seq_along(body)) {
    v <- strsplit(body[k], "\t")[[1]][-1]
    m[k, ] <- suppressWarnings(as.numeric(v[keep]))
  }
  gm <- readRDS(file.path(root, "data/processed/gpl", paste0(gpl, ".rds")))
  m <- m[rownames(m) %in% rownames(gm), , drop = FALSE]
  m <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]
  out <- data.frame(ID_REF = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  stem <- sprintf("GSE40272-%s", gpl)
  fout <- file.path(root, "data/expr", paste0(stem, ".rds"))
  if (!dry) saveRDS(out, fout)
  tbl <- gsub("-", "_", stem, fixed = TRUE)
  ok <- "dry-run"
  if (!dry) {
    dd <- m; storage.mode(dd) <- "double"
    dd <- as.data.frame(dd, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(dd) <- rownames(m)
    ft <- rep("FLOAT", ncol(dd)); names(ft) <- colnames(dd)
    ft <- c(row_names = "VARCHAR(255)", ft)
    ok <- tryCatch({
      dbWriteTable(con, tbl, dd, overwrite = TRUE, row.names = TRUE, field.types = ft)
      dbExecute(con, sprintf("ALTER TABLE `%s` ENGINE=MyISAM", tbl))
      TRUE
    }, error = function(e) conditionMessage(e))
    ok <- if (isTRUE(ok)) "ok" else ok
  }
  log_rows[[length(log_rows)+1]] <- data.frame(Accession = stem, probes = nrow(m),
                                               samples = ncol(m), status = ok, stringsAsFactors = FALSE)
  message(sprintf("%s: %d probes x %d samples -> %s (%s)", stem, nrow(m), ncol(m), tbl, ok))
}

res <- do.call(rbind, log_rows)
print(res, row.names = FALSE)
write.csv(res, cpas_out("db_gse40272_expr_log.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# run_17_recod_gse86166_status()  <-  verbatim pipeline/R/17_recod_gse86166_status.R
# ---------------------------------------------------------------------------
run_17_recod_gse86166_status <- function() {
# 17_recod_gse86166_status.R -------------------------------------------------
# GSE86166(乳腺癌,366 例)的 OS/RFS 状态在源数据中编码为 1/2,不是 0/1:
#   vital status : 1 = alive, 2 = dead   -> OS_status 必须反转
#   recurrence   : 1 = recurrence, 2 = none -> RFS_status 保持方向
#
# 证据(可核对):
#   * Prabhakaran et al. 2017 (PMID 28629479 / PMC5477261) 明确写道
#     "71.9% were alive at time of data collection" —— 366 x 71.9% = 263 例,
#     正是数据中 vital status = 1 的例数(263/366 = 71.86%);其余 103 例为死亡。
#   * "83 patients (23.4%) had disease recurrence";数据中 recurrence = 1 共 71 例。
#   * 时间关系:recurrence = 1 的 71 例中 69 例 rfs < os;recurrence = 2 的 284 例
#     全部 rfs == os(无复发者 RFS 在末次随访删失)——故 recurrence 1 = 事件。
#   * 时间单位已是年(0.21-18.0 年 = 2.5-216 月,与文献中位随访 66.3 月一致)。
#
# 输出: 就地更新 data/processed/surv/GSE86166_surv.rds + 同步 DB 表 GSE86166_surv
# 备份: pipeline/backup/surv_pre_status_recoding/GSE86166_surv.rds
# 用法: CPAS_DB_PASSWORD=... Rscript pipeline/R/17_recod_gse86166_status.R
# ----------------------------------------------------------------------------
suppressPackageStartupMessages(library(RMySQL))
root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
f <- file.path(root, "data/processed/surv/GSE86166_surv.rds")
s <- readRDS(f)
stopifnot(all(c("OS_status", "RFS_status") %in% names(s)))
cat("before: OS", paste(names(table(s$OS_status)), table(s$OS_status), sep=":", collapse=" "),
    "| RFS", paste(names(table(s$RFS_status)), table(s$RFS_status), sep=":", collapse=" "), "\n")
s$OS_status  <- ifelse(s$OS_status  == 1, 0, ifelse(s$OS_status  == 2, 1, NA))
s$RFS_status <- ifelse(s$RFS_status == 1, 1, ifelse(s$RFS_status == 2, 0, NA))
cat("after : OS", paste(names(table(s$OS_status)), table(s$OS_status), sep=":", collapse=" "),
    "| RFS", paste(names(table(s$RFS_status)), table(s$RFS_status), sep=":", collapse=" "), "\n")
saveRDS(s, f)
cat("local rds updated:", f, "\n")

pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
if (!nzchar(pw)) { cat("CPAS_DB_PASSWORD not set: DB not updated.\n"); quit(status = 0) }
con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                 dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"), password = pw)
sch <- dbGetQuery(con, "SHOW COLUMNS FROM GSE86166_surv")
id <- rownames(s); if (is.null(id)) id <- as.character(s[[1]])
out <- data.frame(row_names = id, stringsAsFactors = FALSE)
for (i in seq_len(nrow(sch))) {
  cc <- sch$Field[i]; if (cc == "row_names") next
  v <- if (cc %in% names(s)) s[[cc]] else rep(NA, nrow(s))
  out[[cc]] <- if (grepl("double|float|int", sch$Type[i])) suppressWarnings(as.numeric(v)) else as.character(v)
}
dbBegin(con); dbExecute(con, "DELETE FROM GSE86166_surv")
dbWriteTable(con, "GSE86166_surv", out, append = TRUE, row.names = FALSE); dbCommit(con)
chk <- dbGetQuery(con, "SELECT OS_status v, COUNT(*) n FROM GSE86166_surv GROUP BY OS_status ORDER BY v")
cat("DB after update:\n"); print(chk)
dbDisconnect(con)
}

# ---------------------------------------------------------------------------
# run_23_fix_surv_ids_and_endpoints()  <-  verbatim pipeline/R/23_fix_surv_ids_and_endpoints.R
# ---------------------------------------------------------------------------
run_23_fix_surv_ids_and_endpoints <- function() {
# 23_fix_surv_ids_and_endpoints.R -------------------------------------------
# 修复三类已确认的生存数据/标注缺陷。
#
# A. 生存表样本 ID 丢失（row_names 是行号 1,2,3…，与表达表的 GSM 无交集）
#      GSE4573        pheno 的 geo_accession 129 个；OS 68 事件（数据本身完好）
#      GSE3494_GPL96  pheno 的 geo_accession 179 个；DSS 36 事件
#      GSE3494_GPL97  同上（另一平台，样本集与 GPL96 不重叠）
#    修法：按行对齐校验通过后，把 rownames 换成 pheno$geo_accession。
#
# B. GSE48075 的 OS_status 整列 NA（源里是 "os censor: uncensored/censored"，
#    解析器未识别，只保留了 dss censor 列）
#    修法：OS_status = (os censor == uncensored)；OS_time 不变（已由
#    "survival (mo)" / 12 换算，逐样本核对 73/73 一致）。
#    不额外标注 DSS：源只给一列生存时间，DSS 需要"肿瘤特异性死亡时间"，
#    与 OS 共用同一列会把非肿瘤死亡算成事件。
#
# C. GSE70768 主要终点 OS 仅 1 个事件（源里 OS 为 Y=1/N=56/UNKNOWN=128）
#    修法：主要终点改为 RFS（BCR 事件，29 事件 / 41 完整记录），并取消 OS 标注
#    （沿用 10_parse_cgga/11_endpoint_families 的既有原则：终点要有事件才标注）。
#    生存表里的 OS 列保留，不删源数据。
#
# 用法
#   Rscript pipeline/R/23_fix_surv_ids_and_endpoints.R              # dry-run
#   Rscript pipeline/R/23_fix_surv_ids_and_endpoints.R --apply      # 改本地文件
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/23_fix_surv_ids_and_endpoints.R --apply --db
#   （--db 重新上传 4 张被修的 *_surv 表）
# 备份：pipeline/backup/surv_pre_idfix/
# ----------------------------------------------------------------------------
.libPaths(c("/home/Jingle/R/library", .libPaths()))
suppressPackageStartupMessages(library(RMySQL))

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
apply_l <- "--apply" %in% args
use_db  <- "--db" %in% args
bak_dir <- file.path(root, "pipeline/backup/surv_pre_idfix")
dir.create(bak_dir, recursive = TRUE, showWarnings = FALSE)
note <- function(...) cat(paste0(...), "\n")

read_surv <- function(a) readRDS(file.path(root, "data/processed/surv", paste0(a, "_surv.rds")))
read_pheno <- function(a) readRDS(file.path(root, "data/pheno", paste0(a, ".rds")))
save_surv <- function(a, x) {
  f <- file.path(root, "data/processed/surv", paste0(a, "_surv.rds"))
  b <- file.path(bak_dir, paste0(a, "_surv.rds"))
  if (!file.exists(b)) file.copy(f, b, overwrite = FALSE)
  saveRDS(x, f)
}

# ---------------- A) 恢复 GSM 行名 ----------------
note("## A) 生存表样本 ID 恢复")
id_fix <- list(GSE4573 = "AGE", `GSE3494-GPL96` = "age", `GSE3494-GPL97` = "age")
reportA <- data.frame()
for (a in names(id_fix)) {
  ph <- read_pheno(a); sv <- read_surv(a)
  if (is.null(ph$geo_accession)) { note("  ", a, " 跳过：pheno 无 geo_accession"); next }
  if (nrow(ph) != nrow(sv)) { note("  ", a, " 跳过：行数不一致 ", nrow(ph), " vs ", nrow(sv)); next }
  # 行对齐校验：用同名临床列逐行比对
  col_ph <- id_fix[[a]]; col_sv <- if (col_ph == "AGE") "age" else "age"
  same <- sum(as.numeric(ph[[col_ph]]) == as.numeric(sv[[col_sv]]), na.rm = TRUE)
  aligned <- same == sum(!is.na(as.numeric(sv[[col_sv]])))
  note(sprintf("  %-16s 行对齐校验 %d/%d 通过=%s | 现 rownames 例: %s",
               a, same, nrow(sv), aligned, paste(head(rownames(sv), 2), collapse = ",")))
  if (!aligned) { note("  ", a, " 跳过：行对齐失败，拒绝猜测映射"); next }
  old <- rownames(sv)
  rownames(sv) <- as.character(ph$geo_accession)
  # 修复后核对映射是否有效：与表达表 ID 应有交集（拿本地表达文件核对）
  f_expr <- file.path(root, "data/expr", paste0(a, ".rds"))
  inter <- NA_integer_
  if (file.exists(f_expr)) {
    ex <- readRDS(f_expr); inter <- length(intersect(rownames(sv), names(ex)))
  }
  note(sprintf("     -> rownames 改为 GSM（如 %s）| 与本地表达表交集 = %s",
               paste(head(rownames(sv), 2), collapse = ","), as.character(inter)))
  reportA <- rbind(reportA, data.frame(dataset = a, rows = nrow(sv), matched_expression = inter))
  if (apply_l) save_surv(a, sv)
}

# ---------------- B) GSE48075 由源推导 OS_status ----------------
note("\n## B) GSE48075 OS_status 由 'os censor' 推导")
a <- "GSE48075"
ph <- read_pheno(a); sv <- read_surv(a)
cen <- ph[["os censor:ch1"]]; tim <- ph[["survival (mo):ch1"]]
st <- ifelse(cen == "uncensored", 1L, ifelse(cen == "censored", 0L, NA_integer_))
tt <- as.numeric(tim) / 12
names(st) <- rownames(ph)
k <- match(rownames(sv), rownames(ph))
note(sprintf("  源: censored=%d uncensored=%d NA=%d | 推导后完整记录=%d 事件=%d",
             sum(cen == "censored", na.rm = TRUE), sum(cen == "uncensored", na.rm = TRUE), sum(is.na(cen)),
             sum(!is.na(st[k]) & !is.na(tt[k])), sum(st[k] == 1, na.rm = TRUE)))
same_time <- sum(abs(tt[k] - sv$OS_time) < 1e-6, na.rm = TRUE)
note(sprintf("  OS_time 与源推导一致: %d/%d（未改动该列）", same_time, sum(!is.na(sv$OS_time))))
sv$OS_status <- as.integer(st[k])
if (apply_l) save_surv(a, sv)

# ---------------- C) GSE70768 主要终点改 RFS ----------------
note("\n## C) GSE70768 主要终点 OS -> RFS（OS 仅 1 事件）")
sv7 <- read_surv("GSE70768")
for (tok in c("OS", "RFS")) {
  st <- sv7[[paste0(tok, "_status")]]; tt <- sv7[[paste0(tok, "_time")]]
  note(sprintf("  %-3s 完整记录=%3d 事件=%3d", tok, sum(!is.na(st) & !is.na(tt)), sum(st == 1 & !is.na(tt), na.rm = TRUE)))
}
di <- get(load(file.path(root, "CanPAS/data/dataset_info.rda")))
i <- which(di$Accession == "GSE70768")
note(sprintf("  catalog 现: SurvivalTypes=%s EndpointFamilies=%s EndpointPrimary=%s EP_OS=%s EP_DFS=%s",
             di$SurvivalTypes[i], di$EndpointFamilies[i], di$EndpointPrimary[i], di$EP_OS[i], di$EP_DFS[i]))
if (apply_l) {
  di$SurvivalTypes[i] <- "RFS"; di$EndpointFamilies[i] <- "DFS"
  di$EndpointPrimary[i] <- "RFS"; di$EP_OS[i] <- NA_character_
  dataset_info <- di
  save(dataset_info, file = file.path(root, "CanPAS/data/dataset_info.rda"), version = 2)
  write.csv(di, file.path(root, cpas_data("dataset_info.csv")), row.names = FALSE)
  note("  已改：SurvivalTypes=RFS, EndpointFamilies=DFS, EndpointPrimary=RFS, EP_OS=NA")
}

if (!apply_l) { note("\n[dry-run] 未修改任何文件。加 --apply 执行。"); quit(save = "no") }

# ---------------- D) 重新上传被修的生存表 ----------------
if (use_db) {
  pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
  if (!nzchar(pw)) stop("CPAS_DB_PASSWORD 未设置。")
  con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                   dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                   password = pw, client_flag = CLIENT_COMPRESS)
  on.exit(dbDisconnect(con), add = TRUE)
  note("\n## D) 重新上传生存表")
  for (a in c("GSE4573", "GSE3494-GPL96", "GSE3494-GPL97", "GSE48075")) {
    sv <- read_surv(a)
    out <- data.frame(row_names = rownames(sv), sv, check.names = FALSE, stringsAsFactors = FALSE)
    t <- paste0(gsub("-", "_", a), "_surv")
    dbWriteTable(con, name = t, value = out, overwrite = TRUE, row.names = FALSE)
    n <- dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM `%s`", t))$n
    note(sprintf("  %-16s 已上传 %d 行（原 %d 行）", t, n, nrow(sv)))
  }
}
note("\n完成。后续请重跑 pipeline/R/19_fix_catalog_N.R --write 以更新 catalog 样本量。")
}

# ---------------------------------------------------------------------------
# run_25_fix_endpoint_annotations()  <-  verbatim pipeline/R/25_fix_endpoint_annotations.R
# ---------------------------------------------------------------------------
run_25_fix_endpoint_annotations <- function() {
# 25_fix_endpoint_annotations.R ----------------------------------------------
# 处理三件与终点标注/缺失数据有关的事。
#
# A) GSE40272 四个平台行的标注（按"终点要有事件"原则逐行判断）
#      GPL15973  DFS n=40 事件=10  -> 保留
#      GPL15971  DFS n=17 事件= 2  -> 取消标注
#      GPL15972  DFS n=10 事件= 2  -> 取消标注
#      GPL9497   DFS n=15 事件= 3  -> 取消标注
#    四个平台样本互不重叠（38+13+78+24 = 153 = 该研究总样本数），取消的是标注不是数据。
#
# B) GSE5327（58 例乳腺癌，Minn 2005 肺转移队列）：本地 pheno 里一直有生存字段
#    但从未生成生存表 -> 由 `metastasis free survival (yr)` + `metastasis` 生成 MFS
#    （11 事件），并上传表达与生存表。
#
# C) GSE7849（78 例乳腺癌）：本地 raw series matrix 里有 DFS（月）+ Recurrence，
#    但解析器只认 "key: value" 格式，这批 "Key = Value" 对角排布的字段被丢弃。
#    本脚本直接从 data/raw/GSE7849.gz 逐样本解析 -> 生成 DFS 生存表并上传。
#
# 用法
#   Rscript pipeline/R/25_fix_endpoint_annotations.R                 # dry-run
#   Rscript pipeline/R/25_fix_endpoint_annotations.R --apply         # 写本地文件 + catalog
#   CPAS_DB_PASSWORD=... Rscript pipeline/R/25_fix_endpoint_annotations.R --apply --db
# 备份：pipeline/backup/cohort_build/
# ----------------------------------------------------------------------------
.libPaths(c("/home/Jingle/R/library", .libPaths()))
suppressPackageStartupMessages(library(RMySQL))

root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
apply_l <- "--apply" %in% args
use_db  <- "--db" %in% args
bak <- file.path(root, "pipeline/backup/cohort_build")
dir.create(bak, recursive = TRUE, showWarnings = FALSE)
note <- function(...) cat(paste0(...), "\n")

di <- get(load(file.path(root, "CanPAS/data/dataset_info.rda")))
MIN_EVENTS <- 5L      # 与 cpas_meta()/cpas_meta_panel() 的 min_events 默认值一致

# ============ A) GSE40272 逐平台判断 =========================================
note("## A) GSE40272 平台行标注（事件数 >= ", MIN_EVENTS, " 才标注）")
con <- NULL
if (use_db || TRUE) {
  pw <- Sys.getenv("CPAS_DB_PASSWORD", unset = "")
  if (!nzchar(pw)) stop("CPAS_DB_PASSWORD 未设置（本脚本需要读镜像的事件数）。")
  con <- dbConnect(MySQL(), host = Sys.getenv("CPAS_DB_HOST", unset = "139.224.80.159"),
                   dbname = "cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CanPAS"),
                   password = pw, client_flag = CLIENT_COMPRESS)
}
rows40272 <- grep("^GSE40272", di$Accession, value = TRUE)
changed <- character(0)
for (a in rows40272) {
  adb <- gsub("-", "_", a); st <- paste0(adb, "_surv")
  if (!st %in% dbGetQuery(con, "SHOW TABLES")[[1]]) { note("  ", a, " 无生存表"); next }
  sv <- dbGetQuery(con, sprintf("SELECT * FROM `%s`", st))
  ex <- setdiff(dbListFields(con, adb), "row_names")
  ids <- as.character(sv[[1]])
  fam <- CanPAS::endpoint_family(as.character(di$EndpointPrimary[di$Accession == a])[1])
  tok <- CanPAS::endpoint_resolve(a, fam)
  v <- sv[[paste0(tok, "_status")]]; tt <- sv[[paste0(tok, "_time")]]
  keep <- !is.na(v) & !is.na(tt) & ids %in% ex
  ev <- sum(v[keep] == 1, na.rm = TRUE)
  note(sprintf("  %-18s %s: n=%d 事件=%d -> %s", a, tok, sum(keep), ev,
               if (ev >= MIN_EVENTS) "保留标注" else "取消标注"))
  if (ev < MIN_EVENTS) {
    i <- which(di$Accession == a)
    di$SurvivalTypes[i] <- NA_character_
    di$EndpointFamilies[i] <- NA_character_
    di$EndpointPrimary[i] <- NA_character_
    for (cc in c("EP_OS", "EP_DSS", "EP_DFS", "EP_PFS", "EP_MFS")) di[[cc]][i] <- NA_character_
    changed <- c(changed, a)
  }
}

# ============ B) GSE5327 -> MFS ==============================================
note("\n## B) GSE5327 生成 MFS 生存表（本地 pheno 已有字段，此前未建表）")
build_surv <- function(acc, time_col, status_col, token, extra = list()) {
  p <- readRDS(file.path(root, "data/pheno", paste0(acc, ".rds")))
  if (!all(c(time_col, status_col) %in% names(p)))
    stop("GSE ", acc, ": pheno 缺少 ", time_col, " / ", status_col)
  tt <- suppressWarnings(as.numeric(p[[time_col]]))
  st <- suppressWarnings(as.integer(p[[status_col]]))
  if (any(!st %in% c(0L, 1L, NA_integer_))) stop(acc, ": status 不是 0/1")
  out <- data.frame(row.names = rownames(p))
  out[[paste0(token, "_status")]] <- st
  out[[paste0(token, "_time")]] <- tt
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}
sv5327 <- build_surv("GSE5327", "metastasis free survival (yr):ch1", "metastasis:ch1", "MFS",
                     extra = list(cohort = as.character(readRDS(file.path(root, "data/pheno", "GSE5327.rds"))[["cohort:ch1"]]),
                                  lms_status = as.integer(readRDS(file.path(root, "data/pheno", "GSE5327.rds"))[["lms status:ch1"]])))
note(sprintf("  构建 MFS: n=%d 事件=%d（追踪 %s 年）",
             sum(!is.na(sv5327$MFS_status) & !is.na(sv5327$MFS_time)),
             sum(sv5327$MFS_status == 1, na.rm = TRUE),
             paste(round(range(sv5327$MFS_time, na.rm = TRUE), 2), collapse = "-")))
if (apply_l) {
  f <- file.path(root, "data/processed/surv", "GSE5327_surv.rds")
  if (file.exists(f)) file.copy(f, file.path(bak, "GSE5327_surv.rds"), overwrite = FALSE)
  saveRDS(sv5327, f); note("  已写出 data/processed/surv/GSE5327_surv.rds")
}

# ============ C) GSE7849 -> DFS（解析 raw series matrix 的 "Key = Value"）=====
note("\n## C) GSE7849 生成 DFS 生存表（raw 里是 Key = Value 对角格式）")
raw <- readLines(gzfile(file.path(root, "data/raw/GSE7849.gz")), warn = FALSE)
gsm <- strsplit(grep("^!Sample_geo_accession", raw, value = TRUE)[1], "\t")[[1]][-1]
gsm <- gsub('^"|"$', "", gsm)
ch  <- grep("^!Sample_characteristics", raw, value = TRUE)
kv  <- lapply(seq_along(gsm), function(j) {
  vals <- vapply(ch, function(l) { v <- strsplit(l, "\t")[[1]][-1]; gsub('^"|"$', "", v[j]) }, character(1))
  vals <- vals[!is.na(vals) & nzchar(vals) & grepl("=", vals)]
  k <- trimws(sub("=.*$", "", vals)); vv <- trimws(sub("^[^=]*=", "", vals))
  stats::setNames(vv, k)
})
pick <- function(key) vapply(kv, function(m) if (key %in% names(m)) m[[key]] else NA_character_, character(1))
dfs_m <- suppressWarnings(as.numeric(pick("DFS")))
rec   <- suppressWarnings(as.integer(pick("Recurrence")))
age   <- suppressWarnings(as.numeric(pick("Age_at_Dx")))
grade <- pick("Nuclear_Grade")
er    <- pick("ER"); pr <- pick("PR"); nodes <- suppressWarnings(as.integer(pick("Regional_Nodes_Positive")))
# raw 里 `Recurrence` 只在复发样本上出现（14 例），删失样本无该字段：
# 有随访时间而无 Recurrence 的样本按删失（0）处理。
rec_full <- ifelse(!is.na(rec) & rec == 1L, 1L,
                   ifelse(!is.na(dfs_m), 0L, NA_integer_))
note(sprintf("  解析到 DFS 时间 %d 例（月）、Recurrence 标记 %d 例（仅复发样本有值）；补齐删失后 事件=%d / 可分析=%d",
             sum(!is.na(dfs_m)), sum(!is.na(rec)), sum(rec_full == 1, na.rm = TRUE),
             sum(!is.na(rec_full) & !is.na(dfs_m))))
sv7849 <- data.frame(DFS_status = rec_full, DFS_time = dfs_m / 12,
                     age = age, grade = grade, ER = er, PR = pr, nodes_positive = nodes,
                     row.names = gsm, check.names = FALSE)
note(sprintf("  构建 DFS: n=%d 事件=%d（追踪 %s 年）",
             sum(!is.na(sv7849$DFS_status) & !is.na(sv7849$DFS_time)),
             sum(sv7849$DFS_status == 1, na.rm = TRUE),
             paste(round(range(sv7849$DFS_time, na.rm = TRUE), 2), collapse = "-")))
if (apply_l) {
  f <- file.path(root, "data/processed/surv", "GSE7849_surv.rds")
  if (file.exists(f)) file.copy(f, file.path(bak, "GSE7849_surv.rds"), overwrite = FALSE)
  saveRDS(sv7849, f); note("  已写出 data/processed/surv/GSE7849_surv.rds")
}

# ============ D) catalog 标注（新增两个可分析队列 + GSE40272 改动）==========
note("\n## D) catalog 标注")
# 注意：函数内直接改 di 只会改到副本，必须返回后重新赋值
set_anno <- function(x, acc, tokens, families, primary, ep) {
  i <- which(x$Accession == acc)
  if (!length(i)) stop("accession not in catalog: ", acc)
  x$SurvivalTypes[i] <- paste(tokens, collapse = ",")
  x$EndpointFamilies[i] <- paste(families, collapse = ",")
  x$EndpointPrimary[i] <- primary
  for (cc in c("EP_OS", "EP_DSS", "EP_DFS", "EP_PFS", "EP_MFS")) x[[cc]][i] <- NA_character_
  for (nm in names(ep)) x[[nm]][i] <- ep[[nm]]
  x
}
di <- set_anno(di, "GSE5327", "MFS", "MFS", "MFS", list(EP_MFS = "MFS"))
di <- set_anno(di, "GSE7849", "DFS", "DFS", "DFS", list(EP_DFS = "DFS"))
note("  GSE5327 -> MFS；GSE7849 -> DFS（均为新可分析队列）")
note("  GSE40272 取消标注的行: ", if (length(changed)) paste(changed, collapse = ", ") else "（无）")

if (!apply_l) { note("\n[dry-run] 未写文件。加 --apply 执行。"); dbDisconnect(con); quit(save = "no") }
dataset_info <- di
save(dataset_info, file = file.path(root, "CanPAS/data/dataset_info.rda"), version = 2)
write.csv(di, file.path(root, cpas_data("dataset_info.csv")), row.names = FALSE)
note("  已写回 CanPAS/data/dataset_info.rda 与 data/dataset_info.csv")

# ============ E) 上传两张新生存表与表达表 ===================================
if (use_db) {
  note("\n## E) 上传到镜像（06_upload_db.R，强制覆盖）")
  for (a in c("GSE5327", "GSE7849")) {
    r <- system2("Rscript", c("pipeline/R/06_upload_db.R", a, root, "overwrite"),
                 stdout = TRUE, stderr = TRUE, env = paste0("CPAS_DB_PASSWORD=", 
                 Sys.getenv("CPAS_DB_PASSWORD")))
    note(paste0("  ", a, ": ", paste(tail(r, 3), collapse = " | ")))
  }
}
dbDisconnect(con)
note("\n完成。请接着运行 pipeline/R/19_fix_catalog_N.R --write 刷新样本量。")
}

# ---------------------------------------------------------------------------
# run_35_add_tcga_cohorts()  <-  verbatim pipeline/R/35_add_tcga_cohorts.R
# ---------------------------------------------------------------------------
run_35_add_tcga_cohorts <- function() {
# 35_add_tcga_cohorts.R ------------------------------------------------------
# 目的
#   把 16 个"已核验、尚未编目"的 TCGA 项目登记进 CanPAS：
#     1) 生成 data/processed/surv/TCGA-<PROJ>_surv.rds（每项目一份），列布局/命名/
#        单位与既有 15 份 TCGA-*_surv.rds 完全一致（不新造格式）；
#     2) 在 data/dataset_info.csv 追加 16 行（X 续号；不改任何既有行）；
#     3) 用既有 15 个项目做"生成器回归"：同一生成函数必须能逐格复现既有表，
#        从而证明新增表的列语义/单位/临床列映射与既有队列同源。
#
# 为什么
#   TCGA 队列按需从本地 rda 取数（data/tcga/tcga_surv.rda + tcga_clinical.rda），
#   镜像里刻意不存 TCGA 表，所以新增 TCGA 队列**不需要**任何数据库上传，
#   只需要：每项目一份 survival 表文件 + catalog 行（+ 重建包内 dataset_info.rda）。
#
# 数据源（只读）
#   data/tcga/tcga_surv.rda      object tcga_surv     : sample, OS, OS.time, DSS,
#                                DSS.time, DFI, DFI.time, PFI, PFI.time（status 0/1，time 天）
#   data/tcga/tcga_clinical.rda  object tcga_clinical : sample, patient, type,
#                                age_at_initial_pathologic_diagnosis, gender,
#                                ajcc_pathologic_tumor_stage, histological_type,
#                                histological_grade, ...
#
# 既有 15 份表的格式（2026-09-24 逐格核对所得，全部 15×13 列零差异）
#   rownames = TCGA 样本 barcode（如 TCGA-OR-A5J1-01），行序 = 样本名升序
#   列（13）：OS_status, OS_time, DSS_status, DSS_time, DFI_status, DFI_time,
#            PFI_status, PFI_time, age, sex, stage, histology, grade
#   OS/DSS/DFI/PFI_time = tcga_surv 的 <TOK>.time / 365.25（天 -> 年）
#   age       = age_at_initial_pathologic_diagnosis（数值，原样）
#   sex       = tolower(gender)（male/female）
#   stage     = ajcc_pathologic_tumor_stage：去掉前缀 "Stage "，只保留罗马主分期
#               I/II/III/IV（含 0）；"Stage X"/"[Discrepancy]" 等无法归约 -> NA
#   histology = histological_type（原样；整格以 "[" 开头的占位值如
#               "[Discrepancy]"/"[Not Available]" -> NA）
#   grade     = histological_grade：High Grade->G3、Low Grade->G2、G1-G4 原样，
#               其余（GX/GB/[Unknown]/[Discrepancy]）-> NA
#
# 端点 token 与 SurvivalTypes 约定（由既有 15 行反推 + 回归验证）
#   保留规则：某 token 在**该项目内**至少有 1 个 "status 非缺失且 time 非缺失" 的
#   样本且至少 1 个事件 -> 该 token 进入 SurvivalTypes；否则不登记。
#   验证：该规则恰好复现既有 15 行的 SurvivalTypes（LAML 只有 OS：其 DSS/DFI/PFI
#   在源表中全为 NA；GBM 的 DFI 只有 1 个事件但仍有 1 个可用样本，故保留）。
#   ⚠ 注意：任务书里写的"<5 事件则丢弃"与既有行不一致 —— GBM 的 DFI 仅 1 事件却
#     在既有 catalog 中；且 TGCT 的 OS 只有 4 个事件（任务书本身也要求登记
#     "TGCT 137/4" 这个 OS 队列），若按 <5 事件丢弃会把 OS 也丢掉。
#     因此本脚本默认 MIN_EVENTS = 1（= 既有约定），并在末尾打印 <5 事件的敏感性
#     对比，供你决定是否改口径。
#   ⚠ 既有 LAML 表**保留**了全 NA 的 DSS/DFI/PFI 列（只是 catalog 不登记这些 token），
#     所以"跳过 <5 事件的端点列"在本仓库的真实做法是"不进 SurvivalTypes"，
#     而不是物理删列；为保持 13 列布局一致，新表同样写满 8 个端点列。
#
# 与既有表的唯一差别（有意为之，见任务书）
#   新表只保留"至少有一个可用端点(status+time 都在)的样本"（drop_uninformative=TRUE）；
#   既有 15 份表里还有少量全 NA 行（如 LAML 173 行中 12 行无任何端点）。
#   回归检查因此以 drop_uninformative=FALSE 运行（=既有文件的真实行集）。
#
# 产物
#   data/processed/surv/TCGA-<PROJ>_surv.rds   ×16
#   data/dataset_info.csv                       （追加 16 行）
#   pipeline/out/tcga16_add_report.md           （回归 + 行数核对报告）
# 后续必须执行
#   Rscript pipeline/R/11_endpoint_families.R <ROOT>   # 重算家族列并重建 rda
# 用法
#   Rscript pipeline/R/35_add_tcga_cohorts.R [ROOT]
# 说明：不重建 tarball、不连数据库、不 git 提交、不动 manuscript。
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

ROOT <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(ROOT)) ROOT[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
stopifnot(dir.exists(file.path(ROOT, "data")))
setwd(ROOT)

MIN_EVENTS <- 1L              # 既有约定：>=1 事件即登记该 token（见头部说明）
TOKENS     <- c("OS", "DSS", "DFI", "PFI")   # 源表 token（顺序即列顺序）
CLIN_COLS  <- c("age", "sex", "stage", "histology", "grade")

# ---------------------------------------------------------------- 16 个新项目
# 有效 OS 样本 / OS 事件 来自 tcga_clinical$type ⋈ tcga_surv（sample 连接）
NEW <- data.frame(
  project = c("KIRC","THCA","HNSC","SKCM","KIRP","SARC","ESCA","UCEC","PCPG",
              "TGCT","THYM","KICH","MESO","UVM","ACC","UCS"),
  type = c("Kidney Cancer","Thyroid Cancer","Head and Neck Cancer","Melanoma",
           "Kidney Cancer","Sarcoma","Esophageal Cancer","Endometrial Cancer",
           "Pheochromocytoma","Testicular Cancer","Thymoma","Kidney Cancer",
           "Mesothelioma","Uveal Melanoma","Adrenocortical Cancer",
           "Uterine Carcinosarcoma"),
  stringsAsFactors = FALSE)
# Kidney Cancer 同时含 KIRC/KIRP/KICH，与既有 Type 归类习惯一致（如 Glioma Cancer
# 同时含 GBM/LGG，Colorectal Cancer 含 COAD/READ）。

EXISTING15 <- c("BLCA","BRCA","CESC","COAD","GBM","LAML","LGG","LIHC","LUAD",
                "LUSC","OV","PAAD","PRAD","READ","STAD")

# ---------------------------------------------------------------- 临床列映射
norm_stage <- function(v) {
  v <- trimws(as.character(v))
  v <- sub("^Stage[[:space:]]+", "", v)          # "Stage IIIA" -> "IIIA"
  out <- ifelse(grepl("^(IV|III|II|I|0)", v), sub("^(IV|III|II|I|0).*$", "\\1", v), NA_character_)
  out[v %in% c("", "NA")] <- NA_character_
  out
}
norm_histology <- function(v) {
  v <- trimws(as.character(v))
  v[grepl("^\\[", v)] <- NA_character_
  v[v %in% c("", "NA")] <- NA_character_
  v
}
norm_grade <- function(v) {
  v <- trimws(as.character(v))
  v[!is.na(v) & v == "High Grade"] <- "G3"
  v[!is.na(v) & v == "Low Grade"]  <- "G2"
  v[!is.na(v) & !grepl("^G[1-4]$", v)] <- NA_character_
  v
}

# ---------------------------------------------------------------- 生成器
# drop_uninformative = TRUE : 只保留至少一个可用端点（新表，任务书要求）
#                    = FALSE: 保留项目内全部 tcga_surv 行（既有 15 份表的真实行集）
build_tcga_surv <- function(project, clinical, surv, drop_uninformative = TRUE) {
  cl <- clinical[clinical$type == project, , drop = FALSE]
  sv <- surv[surv$sample %in% cl$sample, , drop = FALSE]
  if (!nrow(sv)) stop("no survival rows for project ", project)
  cl <- cl[match(sv$sample, cl$sample), , drop = FALSE]

  tab <- data.frame(row.names = sv$sample, check.names = FALSE)
  for (tok in TOKENS) {
    tab[[paste0(tok, "_status")]] <- suppressWarnings(as.numeric(sv[[tok]]))
    tab[[paste0(tok, "_time")]]   <- suppressWarnings(as.numeric(sv[[paste0(tok, ".time")]])) / 365.25
  }
  tab$age       <- suppressWarnings(as.numeric(cl$age_at_initial_pathologic_diagnosis))
  tab$sex       <- tolower(as.character(cl$gender))
  tab$stage     <- norm_stage(cl$ajcc_pathologic_tumor_stage)
  tab$histology <- norm_histology(cl$histological_type)
  tab$grade     <- norm_grade(cl$histological_grade)
  tab <- tab[order(rownames(tab)), , drop = FALSE]

  ep <- lapply(TOKENS, function(tok) {
    s <- tab[[paste0(tok, "_status")]]; t <- tab[[paste0(tok, "_time")]]
    ok <- !is.na(s) & !is.na(t)
    c(n = sum(ok), events = sum(s[ok] == 1))
  })
  names(ep) <- TOKENS

  if (drop_uninformative) {
    keep <- rep(FALSE, nrow(tab))
    for (tok in TOKENS) keep <- keep | (!is.na(tab[[paste0(tok, "_status")]]) & !is.na(tab[[paste0(tok, "_time")]]))
    tab <- tab[keep, , drop = FALSE]
  }
  list(table = tab, ep = ep)
}

kept_tokens <- function(ep, min_events = MIN_EVENTS)
  TOKENS[vapply(TOKENS, function(tk) unname(ep[[tk]]["events"]) >= min_events, logical(1))]

# ---------------------------------------------------------------- 端点家族映射
# 与 11_endpoint_families.R 的 FAMILY_OF 保持一致（该脚本会重算并覆盖这些列）
FAMILY_OF <- c(OS = "OS", DSS = "DSS", CSS = "DSS", BCSS = "DSS",
               DFS = "DFS", RFS = "DFS", EFS = "DFS", DFI = "DFS",
               PFS = "PFS", PFI = "PFS", MFS = "MFS", DRFS = "MFS")
BROWSE_OF <- c(OS = "OS", DSS = "DSS", DFS = "DFS", PFS = "PFS", MFS = "PFS")
POOL_ORDER <- c("OS", "DSS", "DFS", "PFS", "MFS")
DERIVED <- c("DFI", "PFI")

family_cols <- function(tokens) {
  fams <- unname(FAMILY_OF[tokens])
  ep <- setNames(rep(NA_character_, 5), POOL_ORDER)
  for (f in POOL_ORDER) { hit <- tokens[fams == f]; if (length(hit)) ep[[f]] <- hit[1] }
  avail <- POOL_ORDER[!is.na(ep)]
  list(EndpointFamilies = if (length(avail)) paste(unique(unname(BROWSE_OF[avail])), collapse = ",") else NA_character_,
       EP = as.list(ep),
       EndpointPrimary = if (length(avail)) ep[[avail[1]]] else NA_character_,
       EndpointDerived = { d <- tokens[tokens %in% DERIVED]; if (length(d)) paste(d, collapse = ",") else NA_character_ })
}

# ---------------------------------------------------------------- 载入数据
e <- new.env(); load(file.path(ROOT, "data/tcga/tcga_surv.rda"), envir = e)
load(file.path(ROOT, "data/tcga/tcga_clinical.rda"), envir = e)
tcga_surv <- get("tcga_surv", envir = e); tcga_clinical <- get("tcga_clinical", envir = e)

di <- read.csv(file.path(ROOT, cpas_data("dataset_info.csv")), check.names = FALSE,
               stringsAsFactors = FALSE)
if (nrow(di) != 136L) stop("expected the pre-change 136-row catalog, got ", nrow(di))
if (any(paste0("TCGA-", NEW$project) %in% di$Accession))
  stop("some of the 16 projects are already catalogued; aborting to stay idempotent-safe")

report <- c("# TCGA 16 队列登记报告", sprintf("日期: %s", format(Sys.time())),
            "", "## 1. 生成器回归（既有 15 项目，逐格比对）", "")

# ---------------------------------------------------------------- 回归检查
for (p in EXISTING15) {
  f <- file.path(ROOT, "data/processed/surv", paste0("TCGA-", p, "_surv.rds"))
  ex <- readRDS(f)
  gen <- build_tcga_surv(p, tcga_clinical, tcga_surv, drop_uninformative = FALSE)$table
  ex <- ex[order(rownames(ex)), , drop = FALSE]
  ok_cols <- identical(colnames(ex), colnames(gen))
  ok_rows <- identical(rownames(ex), rownames(gen))
  nd <- 0L
  if (ok_cols && ok_rows) {
    for (cn in colnames(ex)) {
      a <- ex[[cn]]; b <- gen[[cn]]
      nd <- nd + if (is.numeric(a))
        sum(!((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-9)))
      else sum(!((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)))
    }
  } else nd <- NA_integer_
  row <- di[di$Accession == paste0("TCGA-", p), ]
  st_gen <- paste(kept_tokens(build_tcga_surv(p, tcga_clinical, tcga_surv, FALSE)$ep), collapse = ",")
  st_ex <- as.character(row$SurvivalTypes[1])
  report <- c(report, sprintf("- TCGA-%s: rows=%d/%d cols=%s rownames=%s cells_diff=%s SurvivalTypes=%s (existing=%s) %s",
                              p, nrow(gen), nrow(ex), ok_cols, ok_rows, nd, st_gen, st_ex,
                              if (identical(st_gen, st_ex)) "OK" else "MISMATCH"))
}
# 15 个既有表都通过（rows/cols/cells 全 0 差异）时回归算通过
n_bad <- sum(grepl("MISMATCH", report)) +
  sum(!grepl("cells_diff=0 ", report[grepl("^- TCGA-", report)]))
cat(sprintf("[regression] existing 15 TCGA tables: %s\n",
            if (n_bad == 0) "ALL MATCH (cell-by-cell, 0 differences)" else paste(n_bad, "problem(s) - see report")))

# ---------------------------------------------------------------- 生成 16 份新表
report <- c(report, "", "## 2. 新增 16 项目", "",
            "| project | type | rows(n_surv) | OS eff/ev | DSS eff/ev | DFI eff/ev | PFI eff/ev | SurvivalTypes |",
            "|---|---|---|---|---|---|---|---|")
new_rows <- list()
for (i in seq_len(nrow(NEW))) {
  p <- NEW$project[i]
  g <- build_tcga_surv(p, tcga_clinical, tcga_surv, drop_uninformative = TRUE)
  tab <- g$table; ep <- g$ep
  out <- file.path(ROOT, "data/processed/surv", paste0("TCGA-", p, "_surv.rds"))
  saveRDS(tab, out)
  tk <- kept_tokens(ep)
  fc <- family_cols(tk)
  os <- ep[["OS"]]
  new_rows[[i]] <- data.frame(
    X = NA_integer_,
    Type = NEW$type[i],
    Accession = paste0("TCGA-", p),
    SurvivalTypes = paste(tk, collapse = ","),
    EndpointFamilies = fc$EndpointFamilies,
    EP_OS = fc$EP$OS, EP_DSS = fc$EP$DSS, EP_DFS = fc$EP$DFS,
    EP_PFS = fc$EP$PFS, EP_MFS = fc$EP$MFS,
    EndpointPrimary = fc$EndpointPrimary, EndpointDerived = fc$EndpointDerived,
    GPL = "TCGA_HiSeqV2", N = unname(os["n"]), method = "TCGA-RNAseq",
    n_expr = NA_integer_, n_surv = nrow(tab), n_events = unname(os["events"]),
    n_OS = unname(ep[["OS"]]["n"]), n_DSS = unname(ep[["DSS"]]["n"]),
    n_DFS = unname(ep[["DFI"]]["n"]), n_PFS = unname(ep[["PFI"]]["n"]),
    n_MFS = 0L, expr_in_mirror = FALSE, CohortGroup = NA_character_,
    Note = NA_character_, stringsAsFactors = FALSE)
  report <- c(report, sprintf("| TCGA-%s | %s | %d | %d/%d | %d/%d | %d/%d | %d/%d | %s |",
                              p, NEW$type[i], nrow(tab),
                              ep[["OS"]]["n"], ep[["OS"]]["events"],
                              ep[["DSS"]]["n"], ep[["DSS"]]["events"],
                              ep[["DFI"]]["n"], ep[["DFI"]]["events"],
                              ep[["PFI"]]["n"], ep[["PFI"]]["events"],
                              paste(tk, collapse = ",")))
  cat(sprintf("[new] TCGA-%-5s rows=%-4d OS=%d/%d SurvivalTypes=%s -> %s\n",
              p, nrow(tab), os["n"], os["events"], paste(tk, collapse = ","), basename(out)))
}
new_df <- do.call(rbind, new_rows)

# ---------------------------------------------------------------- 写 catalog
new_df$X <- (max(di$X) + 1L):(max(di$X) + nrow(new_df))   # 续号（catalog 的 X 有历史空洞，不能重排既有 X）
new_df <- new_df[, colnames(di)]
di2 <- rbind(di, new_df)
if (nrow(di2) != 152L) stop("expected 152 catalog rows, got ", nrow(di2))
if (!identical(colnames(di2), colnames(di))) stop("column set changed")
if (any(duplicated(di2$Accession))) stop("duplicated accession after append")
write.csv(di2, file.path(ROOT, cpas_data("dataset_info.csv")), row.names = FALSE)

# ---------------------------------------------------------------- 行数核对
report <- c(report, "", "## 3. 行数与任务书 OS 有效样本对照", "",
            "| project | table rows | OS effective (task) | diff | 说明 |", "|---|---|---|---|---|")
TASK_OS <- c(KIRC = 603, THCA = 571, HNSC = 563, SKCM = 455, KIRP = 320, SARC = 264,
             ESCA = 195, UCEC = 193, PCPG = 185, TGCT = 137, THYM = 120, KICH = 89,
             MESO = 86, UVM = 79, ACC = 77, UCS = 57)
for (i in seq_len(nrow(NEW))) {
  p <- NEW$project[i]; n <- new_df$n_surv[i]; d <- n - TASK_OS[[p]]
  note <- if (d == 0) "与 OS 有效样本一致" else
    sprintf("多 %d 行（仅 DSS/DFI/PFI 可用，无 OS 对），按任务书保留", d)
  report <- c(report, sprintf("| TCGA-%s | %d | %d | %+d | %s |", p, n, TASK_OS[[p]], d, note))
}

# ---------------------------------------------------------------- 敏感性
report <- c(report, "", "## 4. 端点 token 门槛敏感性（MIN_EVENTS = 1 vs 5）", "",
            "既有 15 行的 SurvivalTypes 只有在 `>=1 事件` 规则下才能完全复现：",
            "GBM 的 DFI 只有 1 个事件却在既有 catalog 中；若用 `>=5 事件`，15 行中会有 1 行不一致。",
            "对新增项目，`>=5 事件` 还会额外丢掉 PCPG 的 DFI(4 事件)、TGCT 的 DSS(3 事件)、",
            "KICH 的 DFI(4 事件) —— 以及 TGCT 的 OS（只有 4 个事件，与任务书要求的 TGCT 137/4 冲突）。", "",
            "| project | tokens @>=1 | tokens @>=5 |", "|---|---|---|")
for (i in seq_len(nrow(NEW))) {
  p <- NEW$project[i]
  ep <- build_tcga_surv(p, tcga_clinical, tcga_surv, TRUE)$ep
  report <- c(report, sprintf("| TCGA-%s | %s | %s |", p,
                              paste(kept_tokens(ep, 1), collapse = ","),
                              paste(kept_tokens(ep, 5), collapse = ",")))
}
writeLines(report, cpas_out("tcga16_add_report.md"))
cat("\nwrote data/dataset_info.csv rows:", nrow(di2), "| new tables: 16 |",
    "report: pipeline/out/tcga16_add_report.md\n")
cat("NEXT: Rscript pipeline/R/11_endpoint_families.R", ROOT, "\n")
}

# ---------------------------------------------------------------------------
# run_39_build_gse1379_surv()  <-  verbatim pipeline/R/39_build_gse1379_surv.R
# ---------------------------------------------------------------------------
run_39_build_gse1379_surv <- function() {
# 39_build_gse1379_surv.R ----------------------------------------------------
# 目的
#   重建 GSE1379（breast, GPL1223, whole-tissue sections）的交付生存表
#   data/processed/surv/GSE1379_surv.rds —— 该队列曾经编目、后被作者要求移除，
#   现在作者要求恢复（60 例 / 28 DFS 事件）。
#
# 为什么单独写
#   GSE1379 的随访信息**不在** `!Sample_characteristics_ch1`（该系列一行都没有），
#   而在 series matrix 的 `!Sample_description` 自由文本里，形如
#     Clinical information: Sample.ID=25;Tumor.type=D;Size=3.5;Grade=2;Nodes=0/9;
#     ER=Neg;PR=Pos;HER2=Neg;Age=62;DFS=75;Status=recur
#
# 解析逻辑来源
#   任务书指到的 pipeline/R/33_parse_missed_endpoints.R **在本 checkout 中已不存在**
#   （只剩 run/CanPAS-tool-paper/analysis/db_outcome_audit_20260923/
#    REPORT_missed_endpoints_recovered.md）。同一批审计里保留下来、且仍在仓库中的
#   等价实现是 pipeline/R/36_screen_geo_candidates.R —— 本脚本复用它的三条关键约定：
#     1) strip_quotes()：先 trimws() 去 CRLF 的 \r，再剥单元格两端引号
#        （不要自造 whitespace= 正则，否则 `$` 锚点会静默失效）；
#     2) 否定词先于肯定词：`non-recur` 必须先判为 0，否则 grepl("recur") 会把它
#        误判成事件（这是当初的真实 bug，把 32 个删失都算成事件）；
#     3) 只按顶层 `;` 切段、段内取第一个 `=`（括号内不切）。
#   同一解析家族当时产出的表（GSE10885-GPL1390_surv.rds / GSE20624-GPL1390_surv.rds）
#   的格式即本表格式：行名 = GSM，时间列 = **年**，`<TOKEN>_status` / `<TOKEN>_time`
#   在前、临床列在后。
#
# 单位约定（沿用该队列首次建表时的口径）
#   源字段 DFS 的单位是**月**；交付表统一用**年**：DFS_time = DFS / 12。
#
# 产物
#   data/processed/surv/GSE1379_surv.rds
#   列：DFS_status, DFS_time, tumor_type, size, grade, nodes, ER, PR, HER2, age
#   （DFS 时间外的临床列 = 该自由文本格式提供的字段，原样转录；
#     Size/Grade/Age 为数值，Nodes/ER/PR/HER2 保留源文本，"ND" -> NA）
#
# 用法
#   Rscript pipeline/R/39_build_gse1379_surv.R [ROOT]
# 说明：不连数据库、不写 catalog、不动 tarball/git/manuscript。
# ----------------------------------------------------------------------------

ROOT <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(ROOT)) ROOT[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
stopifnot(dir.exists(file.path(ROOT, "data")))
RAW <- file.path(ROOT, "data/raw/GSE1379.gz")
stopifnot(file.exists(RAW))

# ---------------------------------------------------------------- 基础工具（同 36）
strip_quotes <- function(x) {
  x <- trimws(x)                                  # 去掉 CRLF 的 \r（关键）
  x <- sub('^"', '', x); x <- sub('"$', '', x)
  x <- gsub('""', '"', x, fixed = TRUE)
  trimws(x)
}
split_top <- function(s, seps = ";") {             # 顶层（括号外）切分
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

# ---------------------------------------------------------------- 读 header
h     <- read_header(RAW)
tag   <- sub("\t.*$", "", h)
body  <- sub("^[^\t]*\t?", "", h)
split_vals <- function(i) vapply(strsplit(body[i], "\t", fixed = TRUE)[[1]],
                                 strip_quotes, character(1), USE.NAMES = FALSE)
idx <- function(p) grep(p, tag, ignore.case = TRUE)

gi <- idx("^!Sample_geo_accession$"); stopifnot(length(gi) >= 1L)
gsm <- split_vals(gi[1])

di_idx <- idx("^!Sample_description$"); stopifnot(length(di_idx) >= 1L)
clin_row <- NA_integer_; cells <- character(0)
for (i in di_idx) {
  v <- split_vals(i)
  if (length(v) == length(gsm) && sum(grepl("^Clinical information", v)) == length(gsm)) {
    clin_row <- i; cells <- v; break
  }
}
stopifnot(!is.na(clin_row))
cat(sprintf("[read] GSE1379: %d samples, clinical info in !Sample_description row %d\n",
            length(gsm), clin_row))

# ---------------------------------------------------------------- 解析 key=value
parse_cell <- function(cell) {
  cell <- strip_quotes(cell)
  cell <- sub("^Clinical information:\\s*", "", cell)
  k <- character(0); v <- character(0)
  for (s in split_top(cell)) {
    pos <- regexpr("=", s, fixed = TRUE)[1]
    if (pos < 2L || pos >= nchar(s)) next
    key <- trimws(substr(s, 1L, pos - 1L))
    val <- trimws(substr(s, pos + 1L, nchar(s)))
    if (nzchar(key)) { k <- c(k, key); v <- c(v, val) }
  }
  if (!length(k)) NULL else stats::setNames(v, k)
}
fields <- lapply(cells, parse_cell)
stopifnot(all(vapply(fields, function(z) !is.null(z), logical(1))))

getf <- function(key) vapply(fields, function(z) {
  if (is.null(z) || !key %in% names(z)) NA_character_ else z[[key]]
}, character(1))

# ---------------------------------------------------------------- 事件编码（否定词优先）
# 关键修复：`non-recur` 必须先判 0，否则 grepl("recur") 会误判为事件。
NEG_RE <- "^(non[ _-]?recur|non[ _-]?relaps|no[ _-]?recur|not[ _-]?recur|alive|ned|disease[ _-]?free|no evidence)"
POS_RE <- "^(recur|relaps|progress|dead|death|deceased)"
parse_status <- function(v) {
  v <- tolower(trimws(v))
  out <- rep(NA_real_, length(v))
  out[grepl(NEG_RE, v)] <- 0
  out[is.na(out) & grepl(POS_RE, v)] <- 1
  out
}
num <- function(v) suppressWarnings(as.numeric(trimws(v)))
chr_na <- function(v) { v <- trimws(v); v[tolower(v) %in% c("nd", "na", "n/a", "", "unknown")] <- NA_character_; v }

st  <- parse_status(getf("Status"))
cat("[parse] Status values:", paste(sort(unique(tolower(trimws(getf("Status"))))), collapse = " | "), "\n")
stopifnot(!anyNA(st))                       # 任何无法解码的 Status 都必须暴露出来
dmo <- num(getf("DFS")); stopifnot(!anyNA(dmo))

tab <- data.frame(
  DFS_status = as.numeric(st),
  DFS_time   = dmo / 12,                     # 月 -> 年（首次建表的口径）
  tumor_type = chr_na(getf("Tumor.type")),
  size       = num(gsub("^[^0-9.]*", "", trimws(getf("Size")))),   # ">2.0" -> 2.0
  grade      = num(getf("Grade")),
  nodes      = chr_na(getf("Nodes")),
  ER         = chr_na(getf("ER")),
  PR         = chr_na(getf("PR")),
  HER2       = chr_na(getf("HER2")),
  age        = num(getf("Age")),
  row.names  = gsm, check.names = FALSE, stringsAsFactors = FALSE
)

# ---------------------------------------------------------------- 校验
stopifnot(nrow(tab) == 60L, !anyDuplicated(rownames(tab)))
stopifnot(sum(tab$DFS_status == 1) == 28L)                  # 旧构建实测 60 / 28
stopifnot(all(is.finite(tab$DFS_time)), all(tab$DFS_time > 0))
cat(sprintf("[check] rows=%d | DFS events=%d | DFS censored=%d | DFS_time %.3f-%.2f y\n",
            nrow(tab), sum(tab$DFS_status == 1), sum(tab$DFS_status == 0),
            min(tab$DFS_time), max(tab$DFS_time)))

# 表达 ∩ 生存 必须 = 60（表达表 data/expr/GSE1379.rds 首列 ID_REF，其余列为 GSM）
ex <- readRDS(file.path(ROOT, "data/expr/GSE1379.rds"))
ex_samples <- setdiff(colnames(ex), "ID_REF")
inter <- intersect(ex_samples, rownames(tab))
cat(sprintf("[check] expr samples=%d | surv rows=%d | intersection=%d\n",
            length(ex_samples), nrow(tab), length(inter)))
if (length(ex_samples) != 60L || !identical(sort(ex_samples), sort(rownames(tab))))
  stop("expression and survival sample sets differ (expr=", length(ex_samples),
       " surv=", nrow(tab), " intersection=", length(inter), ")")
cat("[check] expression sample set == survival sample set (60), intersection = 60\n")

out <- file.path(ROOT, "data/processed/surv/GSE1379_surv.rds")
saveRDS(tab, out)
cat("[write]", out, "\n")
cat("NEXT: Rscript pipeline/R/04_qc_report.R", ROOT, "GSE1379\n")
}

# ---------------------------------------------------------------------------
# run_40_add_tcga_chol_dlbc()  <-  verbatim pipeline/R/40_add_tcga_chol_dlbc.R
# ---------------------------------------------------------------------------
run_40_add_tcga_chol_dlbc <- function() {
# 40_add_tcga_chol_dlbc.R ----------------------------------------------------
# 目的
#   把两个此前因"低于尺寸门槛"被搁置的 TCGA 项目补进交付层：
#     TCGA-CHOL（Liver Cancer，45 例 / 45 个可用 OS 对 / 23 事件）
#     TCGA-DLBC（Lymphoma，47 / 47 / 9 事件，OS 时间分布含 1 行 time = 0）
#   生成 data/processed/surv/TCGA-<PROJ>_surv.rds，格式与既有 31 份 TCGA 表**逐格一致**。
#
# 为什么这样做
#   TCGA 队列按需从本地 rda 取数（data/tcga/tcga_surv.rda + tcga_clinical.rda），镜像里
#   刻意不存 TCGA 表，所以新增 TCGA 队列**不需要**任何数据库上传，只需要：
#     每项目一份 survival 表文件 + catalog 行（+ 重建包内 dataset_info.rda）。
#   生成器与 pipeline/R/35_add_tcga_cohorts.R 完全同一套（13 列；时间 = 天 / 365.25；
#   sex/stage/histology/grade 归一化口径一致），本脚本原样复制该函数以便逐格回归。
#
# 保真度证明（本脚本会打印）
#   对既有 31 份表逐一重建并逐格比对：
#     * 35 号脚本当年写入的 15 份旧表 -> drop_uninformative = FALSE（保留全 NA 行）
#     * 35 号脚本当年写入的 16 份新表 -> drop_uninformative = TRUE（丢全 NA 行）
#   两份都至少复现其一，TCGA-UCEC 必须 0 差异（任务书指定的回归锚点）。
#
# 产物
#   data/processed/surv/TCGA-CHOL_surv.rds
#   data/processed/surv/TCGA-DLBC_surv.rds
#   pipeline/out/tcga_chol_dlbc_report.md   （回归 + 新表统计 + catalog 行取值）
# 用法
#   Rscript pipeline/R/40_add_tcga_chol_dlbc.R [ROOT]
# 说明：不连数据库、不写 catalog（catalog 行由 41_register_final3_rows.R 追加）、
#       不动 tarball/git/manuscript。
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

ROOT <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(ROOT)) ROOT[1] else Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(ROOT); stopifnot(dir.exists("data/tcga"))

MIN_EVENTS <- 1L                              # 与 35 号脚本同一约定（既有 catalog 的实际口径）
TOKENS     <- c("OS", "DSS", "DFI", "PFI")
CLIN_COLS  <- c("age", "sex", "stage", "histology", "grade")
NEW_PROJECTS <- c("CHOL", "DLBC")
OLD15 <- c("BLCA","BRCA","CESC","COAD","GBM","LAML","LGG","LIHC","LUAD",
           "LUSC","OV","PAAD","PRAD","READ","STAD")
NEW16 <- c("KIRC","THCA","HNSC","SKCM","KIRP","SARC","ESCA","UCEC","PCPG",
           "TGCT","THYM","KICH","MESO","UVM","ACC","UCS")

# ---------------------------------------------------------------- 临床列映射（同 35）
norm_stage <- function(v) {
  v <- trimws(as.character(v))
  v <- sub("^Stage[[:space:]]+", "", v)
  out <- ifelse(grepl("^(IV|III|II|I|0)", v), sub("^(IV|III|II|I|0).*$", "\\1", v), NA_character_)
  out[v %in% c("", "NA")] <- NA_character_
  out
}
norm_histology <- function(v) {
  v <- trimws(as.character(v)); v[grepl("^\\[", v)] <- NA_character_
  v[v %in% c("", "NA")] <- NA_character_; v
}
norm_grade <- function(v) {
  v <- trimws(as.character(v))
  v[!is.na(v) & v == "High Grade"] <- "G3"
  v[!is.na(v) & v == "Low Grade"]  <- "G2"
  v[!is.na(v) & !grepl("^G[1-4]$", v)] <- NA_character_
  v
}

# ---------------------------------------------------------------- 生成器（原样复制自 35）
build_tcga_surv <- function(project, clinical, surv, drop_uninformative = TRUE) {
  cl <- clinical[clinical$type == project, , drop = FALSE]
  sv <- surv[surv$sample %in% cl$sample, , drop = FALSE]
  if (!nrow(sv)) stop("no survival rows for project ", project)
  cl <- cl[match(sv$sample, cl$sample), , drop = FALSE]

  tab <- data.frame(row.names = sv$sample, check.names = FALSE)
  for (tok in TOKENS) {
    tab[[paste0(tok, "_status")]] <- suppressWarnings(as.numeric(sv[[tok]]))
    tab[[paste0(tok, "_time")]]   <- suppressWarnings(as.numeric(sv[[paste0(tok, ".time")]])) / 365.25
  }
  tab$age       <- suppressWarnings(as.numeric(cl$age_at_initial_pathologic_diagnosis))
  tab$sex       <- tolower(as.character(cl$gender))
  tab$stage     <- norm_stage(cl$ajcc_pathologic_tumor_stage)
  tab$histology <- norm_histology(cl$histological_type)
  tab$grade     <- norm_grade(cl$histological_grade)
  tab <- tab[order(rownames(tab)), , drop = FALSE]

  ep <- lapply(TOKENS, function(tok) {
    s <- tab[[paste0(tok, "_status")]]; t <- tab[[paste0(tok, "_time")]]
    ok <- !is.na(s) & !is.na(t)
    c(n = sum(ok), events = sum(s[ok] == 1))
  })
  names(ep) <- TOKENS

  if (drop_uninformative) {
    keep <- rep(FALSE, nrow(tab))
    for (tok in TOKENS) keep <- keep | (!is.na(tab[[paste0(tok, "_status")]]) & !is.na(tab[[paste0(tok, "_time")]]))
    tab <- tab[keep, , drop = FALSE]
  }
  list(table = tab, ep = ep)
}
kept_tokens <- function(ep, min_events = MIN_EVENTS)
  TOKENS[vapply(TOKENS, function(tk) unname(ep[[tk]]["events"]) >= min_events, logical(1))]

# 逐格比对：返回不一致单元格数（列集/行集不同则 NA）
cells_diff <- function(a, b) {
  a <- a[order(rownames(a)), , drop = FALSE]
  if (!identical(colnames(a), colnames(b)) || !identical(rownames(a), rownames(b))) return(NA_integer_)
  nd <- 0L
  for (cn in colnames(a)) {
    x <- a[[cn]]; y <- b[[cn]]
    nd <- nd + if (is.numeric(x))
      sum(!((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) < 1e-9)))
    else sum(!((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)))
  }
  nd
}

# ---------------------------------------------------------------- 载入
e <- new.env(); load("data/tcga/tcga_surv.rda", envir = e); load("data/tcga/tcga_clinical.rda", envir = e)
tcga_surv <- get("tcga_surv", envir = e); tcga_clinical <- get("tcga_clinical", envir = e)

di <- read.csv(cpas_data("dataset_info.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cat(sprintf("[state] catalog rows = %d\n", nrow(di)))
if (any(paste0("TCGA-", NEW_PROJECTS) %in% di$Accession))
  stop("TCGA-CHOL / TCGA-DLBC already catalogued; aborting to stay idempotent-safe")

rep_lines <- c("# TCGA-CHOL / TCGA-DLBC 交付表生成报告",
               sprintf("日期: %s", format(Sys.time())), "",
               "## 1. 生成器回归（既有 31 份 TCGA 表，逐格比对）", "",
               "| project | batch | rows gen/existing | cells_diff | SurvivalTypes |",
               "|---|---|---|---|---|")

# ---------------------------------------------------------------- 回归
check_one <- function(p, du) {
  f <- file.path("data/processed/surv", paste0("TCGA-", p, "_surv.rds"))
  if (!file.exists(f)) return(NULL)
  ex <- readRDS(f)
  gen <- build_tcga_surv(p, tcga_clinical, tcga_surv, drop_uninformative = du)$table
  nd <- cells_diff(ex, gen)
  list(rows = c(nrow(gen), nrow(ex)), diff = nd,
       st = paste(kept_tokens(build_tcga_surv(p, tcga_clinical, tcga_surv, FALSE)$ep), collapse = ","),
       st_ex = as.character(di$SurvivalTypes[di$Accession == paste0("TCGA-", p)])[1])
}
reg <- list(); n_bad <- 0L
for (p in OLD15) {
  r <- check_one(p, FALSE)
  if (is.null(r)) { reg[[p]] <- "missing"; n_bad <- n_bad + 1L; next }
  ok <- identical(r$diff, 0L) && identical(r$st, r$st_ex)
  reg[[p]] <- r
  rep_lines <- c(rep_lines, sprintf("| TCGA-%s | 15-old | %d/%d | %s | %s (existing %s) %s |",
                                    p, r$rows[1], r$rows[2], r$diff, r$st, r$st_ex,
                                    if (ok) "OK" else "MISMATCH"))
  if (!ok) n_bad <- n_bad + 1L
}
for (p in NEW16) {
  r <- check_one(p, TRUE)
  if (is.null(r)) { reg[[p]] <- "missing"; n_bad <- n_bad + 1L; next }
  ok <- identical(r$diff, 0L) && identical(r$st, r$st_ex)
  reg[[p]] <- r
  rep_lines <- c(rep_lines, sprintf("| TCGA-%s | 16-new | %d/%d | %s | %s (existing %s) %s |",
                                    p, r$rows[1], r$rows[2], r$diff, r$st, r$st_ex,
                                    if (ok) "OK" else "MISMATCH"))
  if (!ok) n_bad <- n_bad + 1L
}
cat(sprintf("[regression] existing 31 TCGA tables: %s\n",
            if (n_bad == 0L) "ALL MATCH (cell-by-cell, 0 differences)" else paste(n_bad, "problem(s)")))
u <- reg[["UCEC"]]
cat(sprintf("[regression] anchor TCGA-UCEC: rows %d/%d, cells_diff=%s, SurvivalTypes=%s (existing %s)\n",
            u$rows[1], u$rows[2], u$diff, u$st, u$st_ex))
if (!identical(u$diff, 0L)) stop("regression anchor TCGA-UCEC did not reproduce cell-by-cell")

# ---------------------------------------------------------------- 生成新表
rep_lines <- c(rep_lines, "", "## 2. 新增两个项目", "",
               "| project | rows(n_surv) | OS | DSS | DFI | PFI | OS time (y) | OS time==0 |",
               "|---|---|---|---|---|---|---|---|")
new_rows <- list()
for (p in NEW_PROJECTS) {
  g <- build_tcga_surv(p, tcga_clinical, tcga_surv, drop_uninformative = TRUE)
  tab <- g$table; ep <- g$ep
  out <- file.path("data/processed/surv", paste0("TCGA-", p, "_surv.rds"))
  saveRDS(tab, out)
  tk <- kept_tokens(ep)
  os <- ep[["OS"]]
  os_t <- tab$OS_time[!is.na(tab$OS_status) & !is.na(tab$OS_time)]
  z <- sum(os_t == 0)
  cat(sprintf("[new] TCGA-%-4s rows=%-3d OS=%d/%d DSS=%d/%d DFI=%d/%d PFI=%d/%d tokens=%s zero-time=%d -> %s\n",
              p, nrow(tab), ep$OS["n"], ep$OS["events"], ep$DSS["n"], ep$DSS["events"],
              ep$DFI["n"], ep$DFI["events"], ep$PFI["n"], ep$PFI["events"],
              paste(tk, collapse = ","), z, basename(out)))
  if (p == "CHOL") stopifnot(nrow(tab) == 45L, unname(os["n"]) == 45L, unname(os["events"]) == 23L,
                             unname(ep$DSS["n"]) == 43L, unname(ep$DFI["n"]) == 32L,
                             unname(ep$PFI["n"]) == 45L, z == 0L,
                             abs(min(os_t) - 0.0274) < 0.01, abs(max(os_t) - 5.41) < 0.01)
  if (p == "DLBC") stopifnot(nrow(tab) == 47L, unname(os["n"]) == 47L, unname(os["events"]) == 9L,
                             unname(ep$DSS["n"]) == 47L, unname(ep$DFI["n"]) == 27L,
                             unname(ep$PFI["n"]) == 47L, z == 1L,
                             max(os_t) > 17.5)
  new_rows[[p]] <- list(tab = tab, ep = ep, tk = tk, zero = z)
  rep_lines <- c(rep_lines, sprintf("| TCGA-%s | %d | %d/%d | %d/%d | %d/%d | %d/%d | %.3f-%.2f | %d |",
                                    p, nrow(tab), ep$OS["n"], ep$OS["events"], ep$DSS["n"], ep$DSS["events"],
                                    ep$DFI["n"], ep$DFI["events"], ep$PFI["n"], ep$PFI["events"],
                                    min(os_t), max(os_t), z))
}
rep_lines <- c(rep_lines, "", "## 3. 回归汇总", "",
               sprintf("- 既有 31 份 TCGA 表逐格比对：%s",
                       if (n_bad == 0L) "全部 0 差异" else sprintf("%d 份不一致", n_bad)),
               "- 回归锚点 TCGA-UCEC：0 差异（任务书指定）",
               "- 新表口径：drop_uninformative = TRUE（与 35 号脚本写入的 16 份新表一致）",
               "", "## 4. 新表的 catalog 取值（由 41_register_final3_rows.R 写入）", "",
               "| Accession | Type | SurvivalTypes | EndpointFamilies | EP_OS | EP_DSS | EP_DFS | EP_PFS | EndpointPrimary | EndpointDerived | N | n_surv | n_events | n_OS | n_DSS | n_DFS | n_PFS | n_MFS |",
               "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (p in NEW_PROJECTS) {
  nr <- new_rows[[p]]; ep <- nr$ep; tk <- nr$tk
  fam <- c(OS = "OS", DSS = "DSS", DFI = "DFS", PFI = "PFS")[tk]
  rep_lines <- c(rep_lines, sprintf("| TCGA-%s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d | %d | %d | %d | %d | %d | %d | 0 |",
    p, if (p == "CHOL") "Liver Cancer" else "Lymphoma", paste(tk, collapse = ","),
    paste(unique(unname(fam)), collapse = ","),
    if ("OS" %in% tk) "OS" else NA, if ("DSS" %in% tk) "DSS" else NA,
    if ("DFI" %in% tk) "DFI" else NA, if ("PFI" %in% tk) "PFI" else NA,
    if (length(tk)) tk[1] else NA,
    paste(tk[tk %in% c("DFI", "PFI")], collapse = ","),
    ep$OS["n"], nrow(nr$tab), ep$OS["events"], ep$OS["n"], ep$DSS["n"], ep$DFI["n"], ep$PFI["n"]))
}
writeLines(rep_lines, cpas_out("tcga_chol_dlbc_report.md"))
cat("\nwrote pipeline/out/tcga_chol_dlbc_report.md\n")
cat("NEXT: Rscript pipeline/R/41_register_final3_rows.R", ROOT, "\n")
}

# ---------------------------------------------------------------------------
# run_110_build_gse205209()  <-  verbatim pipeline/R/110_build_gse205209.R
# ---------------------------------------------------------------------------
run_110_build_gse205209 <- function() {
# 110_build_gse205209.R -------------------------------------------------------
# 构建 GSE205209（子宫内膜癌 / 子宫浆液性癌 USC, NanoString PanCancer IO 360）
#   1) data/expr/GSE205209.rds                  ID_REF + 每 GSM 一列
#   2) data/processed/gpl/GPL27956.rds          探针 -> gene_id(Entrez)
#   3) data/processed/surv/GSE205209_surv.rds   29 患者行（配对 primary/met 去重）
#
# 授权: 作者显式决定纳入，尽管 N=29 低于 relaxed gate (>=30)。
# 用法: Rscript pipeline/R/110_build_gse205209.R [--apply]
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

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("CPAS_DATA_ROOT", unset = "/home/Jingle/data/Project/CPAS")
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
apply_l <- "--apply" %in% args
ACC <- "GSE205209"; GPL <- "GPL27956"
RAW <- file.path(root, "data/raw", paste0(ACC, ".gz"))
WORK <- cpas_out("gse205209_raw")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)

# ---- 1) 解析 series matrix --------------------------------------------------
gz <- gzfile(RAW, "rt")
lines <- readLines(gz, warn = FALSE); close(gz)
hdr <- lines[startsWith(lines, "!")]
tab_b <- which(lines == "!series_matrix_table_begin"); tab_e <- which(lines == "!series_matrix_table_end")

field <- function(key) {
  l <- hdr[startsWith(hdr, key)]
  if (!length(l)) return(NULL)
  v <- strsplit(l[1], "\t")[[1]][-1]
  gsub('^"|"$', "", v)
}
gsm   <- field("!Sample_geo_accession")
title <- field("!Sample_title")
chars <- lapply(hdr[startsWith(hdr, "!Sample_characteristics_ch1")],
                function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
n_s <- length(gsm)
# 值 = 转义冒号还原后最后一个非空 ':' 分段的 trim；这样 "tumor type: M-Metastatic,
# D-Distant, P-Primary: M" -> "M"，而 "time to recurrence (ttr in months) \: 5.2"
# -> "5.2"。
strip_key <- function(x) {
  y <- gsub("\\\\:", ":", x)
  vapply(strsplit(y, ":", fixed = TRUE), function(p) {
    p <- trimws(p); p <- p[nzchar(p)]
    if (length(p)) p[length(p)] else NA_character_
  }, character(1))
}
chr_row <- function(pat) {
  for (v in chars) if (length(v) && grepl(pat, v[1], ignore.case = TRUE)) return(v)
  rep(NA_character_, n_s)
}
# 精确行定位: chars[[i]] 的 key 部分（第一个未转义 ':' 之前）必须恰好等于 key。
# 不能用 grepl(key, ...) —— "tumor purity..." 行里含 "plk3 est. cn" 之后的单元，
# 早期的 grepl(pat) 会把 "tp53 adj. est cn (for tumor purity)" 与
# "tumor purity estimate from wes analysis" 混起来。
chr_row_exact <- function(key) {
  for (v in chars) {
    if (!length(v) || is.na(v[1])) next
    k <- trimws(sub("^([^:]*):.*$", "\\1", gsub("\\\\:", "\u0001", v[1])))
    if (identical(k, key)) return(v)
  }
  rep(NA_character_, n_s)
}
num_of <- function(v) {
  s <- sub("%", "", strip_key(v))
  suppressWarnings(as.numeric(sub("^\\s*([-+]?[0-9]*\\.?[0-9]+).*$", "\\1", s)))
}
# 特征行按 index 取（series matrix 的 !Sample_characteristics_ch1 顺序固定，
# 见上方 chr_row_exact 注释；index 由 chars 的原始顺序确定，与 GEO 一致）。
keymap <- c(tissue = 1, clin_outcome = 2, age = 3, ancestry = 4, ethnicity = 5,
            stage = 6, date_chemo = 7, date_recurrence = 8, date_death = 9,
            ttr_months = 10, os_months = 11, ttr_cens = 12, os_cens = 13,
            site_metastasis = 14, site_recurrence = 15, sample_name = 16,
            tumor_type = 17, wes = 18, purity = 19, immune_class = 20,
            plk3_est_cn = 21, plk3_cna = 22, tp53_cna = 23, log2_plk3 = 24,
            plk3_adj_cn = 25, tp53_adj_cn = 26, tcell_dna = 27)
stopifnot(max(keymap) <= length(chars))
G <- function(nm) chars[[keymap[[nm]]]]
strip1 <- function(nm) strip_key(G(nm))
num1   <- function(nm) num_of(G(nm))
sample_name <- strip1("sample_name")
tumor_type  <- strip1("tumor_type")
tissue      <- strip1("tissue")
os_m   <- num1("os_months")
os_c   <- num1("os_cens")
ttr_m  <- num1("ttr_months")
ttr_c  <- num1("ttr_cens")
age    <- num1("age")
clin_out <- strip1("clin_outcome")
stage    <- strip1("stage")
ethn     <- strip1("ethnicity")
anc      <- strip1("ancestry")
site_met <- strip1("site_metastasis")
date_death <- strip1("date_death")
purity   <- strip1("purity")
tp53_cna <- strip1("tp53_cna")
plk3_cna <- strip1("plk3_cna")
clinical_out <- strip1("immune_class")
# 独立交叉校验：按 key 逐行核对上面的 index 没有错位
stopifnot(grepl("^tissue", G("tissue")[1], ignore.case = TRUE))
stopifnot(grepl("^time to recurrence \\(ttr", G("ttr_months")[1], ignore.case = TRUE))
stopifnot(grepl("^overall survival \\(os", G("os_months")[1], ignore.case = TRUE))
stopifnot(grepl("^overall survival \\(months\\) - censored", G("os_cens")[1], ignore.case = TRUE))
stopifnot(grepl("^tumor purity estimate", G("purity")[1], ignore.case = TRUE))
stopifnot(grepl("^favorable immune response", G("immune_class")[1], ignore.case = TRUE))
stopifnot(grepl("^plk3 cna call", G("plk3_cna")[1], ignore.case = TRUE))
stopifnot(grepl("^tp53 cna call", G("tp53_cna")[1], ignore.case = TRUE))
stopifnot(grepl("^sample_name", G("sample_name")[1], ignore.case = TRUE))
stopifnot(grepl("^tumor type", G("tumor_type")[1], ignore.case = TRUE))
desc <- field("!Sample_description")

# 受试者 id: USC<n> from sample_name (USC7MD / USC20MD are 2nd metastases of USC7/USC20)
subject <- sub("^(USC[0-9]+).*$", "\\1", sample_name)
ty1 <- toupper(substr(tumor_type, 1, 1))          # M / D / P
stopifnot(length(subject) == n_s, !any(is.na(subject)))

# 表达矩阵
d <- fread(RAW, skip = "ID_REF", header = TRUE, na.strings = c("null", "NA", ""),
           data.table = FALSE)
d <- d[!is.na(d$ID_REF), ]
ids <- as.character(d$ID_REF)
m <- as.matrix(d[, setdiff(colnames(d), "ID_REF"), drop = FALSE])
storage.mode(m) <- "double"
empty <- colnames(m)[colSums(is.na(m)) == nrow(m)]
keep_s <- setdiff(colnames(m), empty)
cat(sprintf("series matrix: %d probes x %d samples; %d all-NA columns dropped: %s\n",
            nrow(m), ncol(m), length(empty), paste(empty, collapse = ",")))
# row.names = NULL 很关键：与既有 data/expr/*.rds 一致（.row_names_info < 0 自动行名），
# 否则 06_upload_db.R 的 tibble::column_to_rownames() 会因为显式行名报错。
expr <- data.frame(ID_REF = ids, m[, keep_s, drop = FALSE], check.names = FALSE,
                   row.names = NULL)

# ---- 2) 探针 -> Entrez -----------------------------------------------------
# 平台 SOFT: ID(=SYMBOL) / ORF(=SYMBOL) / SPOT_ID(空) -> 无 Entrez 列；
# 本地 data/gpl/<GPL>.rds 不存在 -> AnnoProbe::idmap(GPL27956) 试一下；
# 仍无则用 org.Hs.eg.db 的 SYMBOL->ENTREZID（与 92/95/96/102 的 symbol2entrez 同口径）。
symbol2entrez <- function(sym) {
  sym <- unique(as.character(sym))
  out <- setNames(rep(NA_character_, length(sym)), sym)
  db <- tryCatch(getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"), error = function(e) NULL)
  if (is.null(db)) return(out)
  mm <- tryCatch(suppressWarnings(AnnotationDbi::select(db, keys = sym,
        columns = "ENTREZID", keytype = "SYMBOL")), error = function(e) NULL)
  if (is.null(mm)) return(out)
  mm <- mm[!is.na(mm$ENTREZID), , drop = FALSE]
  if (nrow(mm)) {
    agg <- tapply(mm$ENTREZID, mm$SYMBOL, function(z) paste(sort(unique(z)), collapse = " /// "))
    out[names(agg)] <- as.character(agg)
  }
  out
}
dest_gpl <- file.path(root, "data/processed/gpl", paste0(GPL, ".rds"))
dir.create(dirname(dest_gpl), showWarnings = FALSE, recursive = TRUE)
sym <- toupper(ids)
gid <- unname(symbol2entrez(sym))
route <- "org.Hs.eg.db SYMBOL->ENTREZID"
if (sum(!is.na(gid)) == 0 && requireNamespace("AnnoProbe", quietly = TRUE)) {
  res <- tryCatch(suppressWarnings(AnnoProbe::idmap(GPL, type = "pipe")), error = function(e) NULL)
  if (!is.null(res) && nrow(res)) {
    symcol <- intersect(c("symbol", "Symbol", "GENE_SYMBOL"), colnames(res))[1]
    s2e <- symbol2entrez(res[[symcol]])
    gid <- unname(s2e[sym])
    route <- "AnnoProbe::idmap(pipe) SYMBOL->ENTREZID"
  }
}
gm <- data.frame(gene_id = gid, row.names = ids, stringsAsFactors = FALSE)
gm <- gm[!duplicated(rownames(gm)), , drop = FALSE]
cov <- 100 * sum(!is.na(gm$gene_id)) / nrow(gm)
cat(sprintf("GPL map (%s): %d probes, %d with ENTREZ = %.2f%%\n",
            route, nrow(gm), sum(!is.na(gm$gene_id)), cov))

# ---- 3) 患者级生存表 -------------------------------------------------------
# 去重规则: 每位受试者优先取 **在表达矩阵里可用** 的 primary (P) 阵列；该受试者
# 没有可用 primary 时，退到第一个可用的 metastatic (M/D) 阵列。
usable_s <- gsm %in% keep_s
pick <- vapply(sort(unique(subject)), function(s) {
  i <- which(subject == s & usable_s)
  if (!length(i)) return(NA_integer_)
  p <- i[ty1[i] == "P"]
  if (length(p)) p[1] else i[1]
}, integer(1))
names(pick) <- sort(unique(subject))
# USC9 的两个阵列（GSM6208643 Met / GSM6208644 Primary）都是全 null —— 正是论文
# Materials and Methods 点名的"patient ID #9 Met and Primary"不可用样本。该受试者
# 保留在生存表里（OS 有值），但没有表达阵列 -> sample_id/array 为 NA，
# expr_available = FALSE。其余 28 位受试者必须有可用阵列。
no_expr <- names(pick)[is.na(pick)]
stopifnot(identical(no_expr, "USC9"))
# 每个受试者第一个阵列（不管有无表达）—— 用于取受试者级临床/终点字段
picksub <- vapply(names(pick), function(s) which(subject == s)[1], integer(1))
# 该受试者被丢弃的配对阵列（含表达不可用的那些）
disc <- setdiff(seq_len(n_s), pick)
pat <- data.frame(
  subject      = names(pick),
  sample_id    = ifelse(is.na(pick), NA_character_, gsm[pick]),
  array        = ifelse(is.na(pick), NA_character_, sample_name[pick]),
  tumor_type   = ifelse(is.na(pick), NA_character_, ty1[pick]),
  discarded    = vapply(names(pick), function(s) {
    psel <- if (is.na(pick[s])) character(0) else sample_name[pick[s]]
    dd <- setdiff(sample_name[subject == s], psel)
    if (length(dd)) paste(dd, collapse = ";") else NA_character_
  }, character(1)),
  # 临床/终点是受试者级的，不依赖阵列是否可用 —— USC9 没有可用表达阵列，但它
  # 的 OS/TTR 值仍然存在（picksub 指向该受试者的第一个阵列）。
  OS_time_mo   = os_m[picksub],
  OS_cens      = os_c[picksub],
  TTR_time_mo  = ttr_m[picksub],
  TTR_cens     = ttr_c[picksub],
  age          = age[picksub],
  stage        = stage[picksub],
  platinum     = clin_out[picksub],
  ethnicity    = ethn[picksub],
  ancestry     = anc[picksub],
  site_metastasis = site_met[picksub],
  date_of_death = date_death[picksub],
  tumor_purity_pct = num1("purity")[picksub],
  tp53_cna_call = tp53_cna[picksub],
  plk3_cna_call = plk3_cna[picksub],
  immune_response_class = clinical_out[picksub],
  stringsAsFactors = FALSE)

stopifnot(nrow(pat) == 29L)
# 单位: GEO 是月 -> 年
pat$OS_time    <- pat$OS_time_mo / 12
pat$OS_status  <- ifelse(is.na(pat$OS_cens), NA_integer_, ifelse(pat$OS_cens == 1, 0L, 1L))
# TTR: 3 例 (USC7/USC12/USC31) 的时间是 0 个月且事件标记为 0 —— 0 不是合法的生存
# 时间，置 NA（不进 N）。其余 26 例 (os_cens=0 -> 事件, =1 -> 删失) 可用。
pat$TTR_time   <- ifelse(is.na(pat$TTR_time_mo) | pat$TTR_time_mo <= 0, NA_real_,
                         pat$TTR_time_mo / 12)
pat$TTR_status <- ifelse(is.na(pat$TTR_time), NA_integer_,
                         ifelse(pat$TTR_cens == 1, 0L, 1L))
sv <- pat[, c("OS_status", "OS_time", "TTR_status", "TTR_time", "age", "stage",
              "platinum", "ethnicity", "ancestry", "site_metastasis",
              "date_of_death", "tumor_purity_pct", "tp53_cna_call", "plk3_cna_call",
              "immune_response_class", "sample_id", "array", "tumor_type", "discarded")]
sv$age <- suppressWarnings(as.numeric(sv$age))
# 0/1 而不是 logical：全部 200 个既有 _surv.rds 都没有 logical 列，RMySQL 对
# logical 的导入/回读会变成 0/TRUE 之类的不一致写法。
sv$expr_available <- as.integer(!is.na(sv$sample_id))
# 行名 = 保留阵列的 GSM（受试者 id 放进 subject 列）。16_verify_catalog_mirror.R 的
# 第 (a) 类检查要求 survival 表首列与 expression 表列名至少有一个交集；受试者级行名
# (USC1...) 与表达表的 GSM 列名交集为 0，会被判 ERROR "shares no sample ID"。
# 行数仍是 29（每位受试者一行）：USC9 没有可用阵列，其 sample_id 为 NA，行名退化为
# 受试者 id "USC9"，该行因此不参与连接 —— 与"该受试者无表达数据"的事实一致。
rn <- ifelse(is.na(pat$sample_id), pat$subject, pat$sample_id)
stopifnot(length(rn) == 29L, !any(duplicated(rn)))
sv$subject <- pat$subject
sv <- sv[, c("subject", setdiff(colnames(sv), "subject"))]
rownames(sv) <- rn

cat(sprintf("\npatients=%d  OS pairs=%d  OS events=%d  censored=%d\n",
            nrow(sv), sum(!is.na(sv$OS_time) & !is.na(sv$OS_status)),
            sum(sv$OS_status == 1, na.rm = TRUE), sum(sv$OS_status == 0, na.rm = TRUE)))
cat(sprintf("OS months range = %.1f .. %.1f  -> years %.3f .. %.3f\n",
            min(pat$OS_time_mo), max(pat$OS_time_mo), min(sv$OS_time), max(sv$OS_time)))
cat(sprintf("TTR pairs=%d  events=%d  censored=%d  (dropped zero-time: %s)\n",
            sum(!is.na(sv$TTR_time)), sum(sv$TTR_status == 1, na.rm = TRUE),
            sum(sv$TTR_status == 0, na.rm = TRUE),
            paste(pat$subject[is.na(sv$TTR_time)], collapse = ",")))
cat("\nexpr usable samples:", ncol(expr) - 1L, " from subjects:",
    length(unique(subject[subject %in% sub("^(USC[0-9]+)$", "\\1", names(pick))])), "\n")
expr_pat <- unique(subject[colnames(expr)[-1] %in% gsm])
cat("expression subjects:", length(expr_pat), "of 29; subjects with NO usable array:",
    paste(sort(setdiff(names(pick), expr_pat)), collapse = ","), "\n")
cat("surv rows without expression:", paste(rownames(sv)[!rownames(sv) %in% expr_pat], collapse = ","), "\n")
cat("\n-- unit cross-check (raw series matrix vs table, years) --\n")
for (s in c("USC1", "USC6", "USC34")) {
  i <- which(gsm == pat$sample_id[pat$subject == s])
  cat(sprintf("  %s: raw os_months=%.4f  -> %.6f yr ; table OS_time=%.6f ; raw cens=%s -> status %d\n",
              s, os_m[i], os_m[i] / 12, sv$OS_time[rownames(sv) == s],
              as.character(os_c[i]), sv$OS_status[rownames(sv) == s]))
}
print(pat[, c("subject","sample_id","array","tumor_type","discarded","OS_time_mo","OS_cens","TTR_time_mo","TTR_cens")])

if (apply_l) {
  saveRDS(expr, file.path(root, "data/expr", paste0(ACC, ".rds")))
  saveRDS(gm, dest_gpl)
  saveRDS(sv, file.path(root, "data/processed/surv", paste0(ACC, "_surv.rds")))
  write.csv(pat, file.path(WORK, "patient_table.csv"), row.names = FALSE)
  cat("\n[applied] wrote data/expr/GSE205209.rds,", dest_gpl,
      ", data/processed/surv/GSE205209_surv.rds\n")
} else cat("\n[dry-run] nothing written\n")
}
