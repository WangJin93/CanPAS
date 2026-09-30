# ===========================================================================
# CanPAS curation pipeline -- 01_ingest_parse
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 01_parse_gse.R
#   2. 10_parse_cgga.R
#   3. 90_build_gse108474_suppl.R
#   4. 91_complete_gse14520_surv.R
#   5. 92_build_geo_expansion_expr_pheno.R
#   6. 93_build_geo_expansion_bespoke.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "01_ingest_parse.R", package = "CanPAS"))
#   run_01_parse_gse()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================


# ---------------------------------------------------------------------------
# run_01_parse_gse()  <-  verbatim pipeline/R/01_parse_gse.R
# ---------------------------------------------------------------------------
run_01_parse_gse <- function() {
# 01_parse_gse.R
# -----------------------------------------------------------------------------
# CanPAS 数据处理 第 1 步：下载 GEO series matrix 并解析为
#   data/raw/<ACC>.gz       原始 series matrix（保留）
#   data/expr/<ACC>.rds     表达谱 data.frame（首列 ID_REF = 探针，其余列为样本）
#   data/pheno/<ACC>.rds    临床信息 data.frame（含 ":ch1" 后缀的特征列，
#                           去除了 !Sample_ 通用元数据中的 characteristics 合并列）
# 并登记一行到 data/dataset_info.csv（Type, Accession, GPL, N, method）
#
# 解析逻辑与 parseGEO.R::parseGSEMatrix 完全一致（去掉 GEOquery getGPL 部分，
# 平台注释由 pipeline/R/02_gpl_map.R 单独处理，与既有 data/processed/gpl 一致）。
#
# 用法: Rscript 01_parse_gse.R <GSE-id> <CancerType> [<ProjectRoot>]
# 例  : Rscript 01_parse_gse.R GSE44001 "cervical cancer" ~/data/Project/CanPAS
# -----------------------------------------------------------------------------

suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(curl); library(Biobase)
})

args <- commandArgs(trailingOnly = TRUE)
ID        <- args[1]
Type      <- args[2]
ROOT      <- if (length(args) >= 3) args[3] else "~/data/Project/CanPAS"
ROOT      <- path.expand(ROOT)
destdir   <- ROOT
setwd(ROOT)

# ---- 1. parseGSEMatrix: 直接复用 parseGEO.R 中的解析函数（getGPL 分支关闭）----
parseGSEMatrix <- function(local_file, AnnotGPL=FALSE, getGPL=TRUE,
                           parseCharacteristics=TRUE) {
  dat <- readr::read_lines(local_file)
  series_header_row_count <- sum(grepl("^!Series_", dat))
  sample_header_start <- grep("^!Sample_", dat)[1]
  samples_header_row_count <- sum(grepl("^!Sample_", dat))
  series_table_begin_line <- grep("^!series_matrix_table_begin", dat)
  if(length(series_table_begin_line) != 1)
    stop("parsing failed--expected only one '!series_matrix_table_begin'")

  header <- read.table(local_file, sep="\t", header=FALSE, nrows=series_header_row_count)
  tmpdat <- read.table(local_file, sep="\t", header=FALSE,
                       nrows=samples_header_row_count, skip=sample_header_start-1)
  headertmp <- t(header); headerdata <- rbind(data.frame(), headertmp[-1,])
  colnames(headerdata) <- sub('!Series_','',as.character(header[,1]))

  tmptmp <- t(tmpdat)
  sampledat <- rbind(data.frame(), tmptmp[-1,])
  colnames(sampledat) <- make.unique(sub('!Sample_','',as.character(tmpdat[,1])))
  sampledat[['geo_accession']] <- as.character(sampledat[['geo_accession']])
  rownames(sampledat) <- sampledat[['geo_accession']]

  if(length(grep('characteristics_ch',colnames(sampledat)))>0 && parseCharacteristics) {
    pd <- sampledat %>%
      dplyr::select(dplyr::contains('characteristics_ch')) %>%
      dplyr::mutate(accession = rownames(.)) %>%
      mutate_all(as.character) %>%
      tidyr::gather(characteristics, kvpair, -accession) %>%
      dplyr::mutate(kvpair = dplyr::case_when(
        grepl(':', kvpair) ~ kvpair,
        grepl(' = ', kvpair, fixed=TRUE) ~ sub(' = ', ':', kvpair, fixed=TRUE),   # "k = v"
        grepl('^\\[', trimws(kvpair)) ~ {  # "[k=v (desc)]" 形式
          x <- trimws(kvpair); x <- sub('^\\[', '', x); x <- sub('\\]$', '', x)
          k <- sub('=.*$', '', x)
          v <- sub('^[^=]*=\\s*', '', x)
          v <- sub('\\s*\\(.*$', '', v)
          paste0(trimws(k), ':', trimws(v)) },
        TRUE ~ kvpair)) %>%
      dplyr::filter(grepl(':',kvpair) & !is.na(kvpair))
    if(nrow(pd)) {
      pd <- pd %>%
        dplyr::mutate(characteristics = ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
        tidyr::separate(kvpair, into=c('k','v'), sep=":", fill='right', extra='merge') %>%
        dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
        dplyr::select(-characteristics) %>%
        dplyr::filter(!is.na(v)) %>%
        dplyr::group_by(accession,k) %>%
        dplyr::mutate(v = paste0(trimws(v), collapse=";")) %>%
        unique() %>%
        tidyr::spread(k,v)
    } else {
      pd <- pd %>% dplyr::select(accession)
    }
    sampledat <- sampledat %>% dplyr::left_join(pd, by=c('geo_accession'='accession'))
  }

  Sys.setenv("VROOM_CONNECTION_SIZE" = 2147483647)
  datamat <- read_tsv(local_file, quote='"',
                      na=c('NA','null','NULL','Null'),
                      skip = series_table_begin_line,
                      comment = '!series_matrix_table_end',
                      skip_empty_rows = FALSE,
                      show_col_types = FALSE)
  tmprownames <- datamat[[1]]
  datamat <- as.matrix(datamat[!is.na(tmprownames),-1])
  rownames(datamat) <- tmprownames[!is.na(tmprownames)]
  datamat <- as.matrix(datamat)
  rownames(sampledat) <- colnames(datamat)
  if(is.null(nrow(datamat))) {
    tmpnames <- names(datamat)
    rownames(sampledat) <- tmpnames
    datamat <- matrix(nrow=0, ncol=nrow(sampledat)); colnames(datamat) <- tmpnames
  }
  eset <- new('ExpressionSet', phenoData=as(sampledat,'AnnotatedDataFrame'),
              annotation=as.character(sampledat[1, grep('platform_id',colnames(sampledat),
                                                        ignore.case=TRUE)]),
              exprs=as.matrix(datamat))
  return(list(GPL=as.character(sampledat[1, grep('platform_id',colnames(sampledat),
                                                 ignore.case=TRUE)]),
              eset=eset))
}

# ---- 2. 从 FTP 目录取 series matrix 文件名 ----
matrix_files <- function(GSE_id, base="https://ftp.ncbi.nlm.nih.gov/geo/series/") {
  series <- str_remove(GSE_id, "-GPL[0-9]+$")          # 平台分文件：父 series 目录
  prefix <- substr(series, 1, nchar(series)-3)
  url <- paste0(base, prefix, "nnn/", series, "/matrix/")
  html <- tryCatch(readLines(url, warn=FALSE), error=function(e) stop("FTP listing failed: ", url))
  fnames <- unique(str_extract(html, "GSE[0-9]+(-GPL[0-9]+)?_series_matrix\\.txt\\.gz"))
  fnames <- fnames[!is.na(fnames)]
  if (GSE_id != series) fnames <- fnames[grepl(paste0("^", GSE_id, "_series"), fnames)]
  if (!length(fnames)) stop("no series matrix file found at ", url)
  message(sprintf("Found %d file(s)", length(fnames))); fnames
}

# ---- 3. 登记 dataset_info.csv ----
append_dataset_info <- function(Type, acc, gpl, n, method) {
  csv <- file.path(ROOT, "data/dataset_info.csv")
  dd <- read.csv(csv, row.names=1, stringsAsFactors=FALSE)
  row <- data.frame(Type=Type, Accession=acc, GPL=gpl, N=n, method=method,
                    stringsAsFactors=FALSE)
  dd <- rbind(dd, row)
  dd <- dd[!duplicated(dd$Accession), ]
  write.csv(dd, csv)
  invisible(acc)
}

# ---- 4. 单个数据集处理 ----
process_data <- function(ID, Type, force=FALSE) {
  fnames <- matrix_files(ID)
  for (fname in fnames) {
    series <- str_remove(ID, "-GPL[0-9]+$")
    prefix <- substr(series, 1, nchar(series)-3)
    matrix_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                         prefix, "nnn/", series, "/matrix/", fname)
    base <- str_remove(fname, "_series_matrix.txt.gz")       # GSE52903[-GPLxxx]
    local_file <- file.path(ROOT, "data/raw", paste0(base, ".gz"))
    if (!file.exists(local_file) || force) {
      message("downloading ", matrix_url)
      curl::curl_download(matrix_url, local_file)
    } else message("raw file exists: ", local_file)

    Sys.setenv("VROOM_CONNECTION_SIZE" = 2147483647)
    res <- parseGSEMatrix(local_file, AnnotGPL=FALSE, getGPL=FALSE,
                          parseCharacteristics=TRUE)
    eset <- res[["eset"]]
    pd <- Biobase::pData(eset)
    type <- if ("type" %in% colnames(pd)) as.character(pd[["type"]][1]) else NA
    append_dataset_info(Type, base,
                        Biobase::annotation(eset),
                        length(pd[["geo_accession"]]), type)
    if (!is.na(type) && stringr::str_detect(type, "SRA")) {
      message("SRA-type dataset skipped for expr/pheno saving: ", base); next
    }
    expr_data <- as.data.frame(Biobase::exprs(eset)) %>%
      tibble::rownames_to_column("ID_REF")
    saveRDS(expr_data, file=file.path(ROOT, "data/expr", paste0(base, ".rds")))
    pheno_data <- pd %>% dplyr::select(-contains("characteristics"))
    saveRDS(pheno_data, file=file.path(ROOT, "data/pheno", paste0(base, ".rds")))
    message("saved expr/pheno for ", base,
            "  (samples=", ncol(expr_data)-1, ", probes=", nrow(expr_data),
            ", GPL=", Biobase::annotation(eset), ", type=", type, ")")
  }
  invisible(ID)
}

process_data(ID, Type)
}

# ---------------------------------------------------------------------------
# run_10_parse_cgga()  <-  verbatim pipeline/R/10_parse_cgga.R
# ---------------------------------------------------------------------------
run_10_parse_cgga <- function() {
# 10_parse_cgga.R ------------------------------------------------------------
# 解析 CGGA.zip 中的三个数据集(mRNA-array_301 / mRNAseq_325 / mRNAseq_693),
# 按 CanPAS 既有 GEO 清洗策略产出可直接上传的三件套:
#   data/expr/<ACC>.rds                  # ID_REF(基因) + 样本列
#   data/processed/gpl/<GPL>.rds         # rownames=基因, gene_id=Entrez(未匹配='---')
#   data/processed/surv/<ACC>_surv.rds   # rownames=样本ID, 终点+临床(标准化)
# 规范化规则与 GEO 一致: 缺失值->NA; sex male/female; age 数值;
#                       grade -> G2/G3/G4(WHO II/III/IV); 时间->年; 状态 0/1
# RNA-seq: log2(RSEM+1); array: 保留原有标准化值
# ----------------------------------------------------------------------------
suppressMessages({library(dplyr)})
root <- "/home/Jingle/data/Project/CPAS"
raw <- file.path(root, "data/raw/cgga")
setwd(root)

SPECS <- list(
  list(acc = "CGGA_693", gpl = "CGGA_693_PLAT", method = "RNA-seq",
       dir = "mRNAseq_693", expr_pat = "RSEM-genes", clin_pat = "clinical",
       transform = "log2"),
  list(acc = "CGGA_325", gpl = "CGGA_325_PLAT", method = "RNA-seq",
       dir = "mRNAseq_325", expr_pat = "RSEM-genes", clin_pat = "clinical",
       transform = "log2"),
  list(acc = "CGGA_301", gpl = "CGGA_301_PLAT", method = "array",
       dir = "mRNA-array_301", expr_pat = "gene_level", clin_pat = "clinical",
       transform = "none")
)

MISS <- c("", "NA", "N/A", "N.A.", "#NA", "#N/A", "NULL", "NaN", "-", "--",
          "—", "–", ".", "?", "UNKNOWN", "UNK", "NOT AVAILABLE", "NOT REPORTED",
          "NOT SPECIFIED", "NOT APPLICABLE", "MISSING", "NONE")
clean_chr <- function(x) {
  x <- as.character(x)
  t <- trimws(x)
  x[tolower(t) %in% tolower(MISS) | grepl("^[#\\s]*n/?a[#\\s]*$", t, ignore.case = TRUE)] <- NA
  x <- gsub("[[:space:]]+", " ", trimws(x))
  x[trimws(x) == ""] <- NA_character_
  x
}
read_txt <- function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  txt <- gsub("\r\n", "\n", txt); txt <- gsub("\r", "\n", txt)
  read.delim(text = txt, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "")
}

ID_map <- local({
  e <- new.env(); utils::data("ID_map", package = "CanPAS", envir = e)
  get("ID_map", e)
})
sym2entrez <- function(sym) {
  hit <- ID_map$gene_id[match(sym, ID_map$Symbol)]
  hit[is.na(hit)] <- "---"
  as.character(hit)
}

report <- c("# CGGA 解析与清洗报告", sprintf("日期: %s", Sys.time()), "")
for (sp in SPECS) {
  d <- file.path(raw, sp$dir)
  ef <- list.files(d, pattern = sp$expr_pat, full.names = TRUE)[1]
  cf <- list.files(d, pattern = sp$clin_pat, full.names = TRUE)[1]
  message("=== ", sp$acc, "\n  expr: ", basename(ef), "\n  clin: ", basename(cf))

  # ---------- clinical ----------
  cl_raw <- read_txt(cf)
  nms <- names(cl_raw)
  find_col <- function(pats) {
    for (p in pats) { i <- grep(p, nms, ignore.case = TRUE)[1]
      if (!is.na(i)) return(nms[i]) }
    NA_character_
  }
  c_id    <- find_col(c("^CGGA_ID$"))
  c_prs   <- find_col(c("^PRS_type$"))
  c_hist  <- find_col(c("^Histology$"))
  c_grade <- find_col(c("^Grade$"))
  c_sex   <- find_col(c("^Gender$"))
  c_age   <- find_col(c("^Age$"))
  c_os    <- find_col(c("^OS$"))
  c_cen   <- find_col(c("^Censor"))
  c_rad   <- find_col(c("^Radio_status"))
  c_che   <- find_col(c("^Chemo_status"))
  c_idh   <- find_col(c("^IDH_mutation_status$"))
  c_cod   <- find_col(c("1p19q.*odeletion"))
  c_mgmt  <- find_col(c("^MGMTp_methylation"))
  c_sub   <- find_col(c("TCGA_subtypes"))
  cl <- data.frame(ID = cl_raw[[c_id]], stringsAsFactors = FALSE)
  cl$prs_type   <- tolower(clean_chr(cl_raw[[c_prs]]))
  h <- clean_chr(cl_raw[[c_hist]])
  h <- ifelse(grepl("^r[A-Z]", h), sub("^r", "", h), h)   # rGBM -> GBM
  cl$histology  <- ifelse(toupper(h) == "SGBM", "GBM", h)  # sGBM -> GBM
  g <- toupper(clean_chr(cl_raw[[c_grade]]))
  cl$grade <- ifelse(g == "WHO II", "G2",
              ifelse(g == "WHO III", "G3",
              ifelse(g == "WHO IV", "G4", NA)))
  sx <- tolower(clean_chr(cl_raw[[c_sex]]))
  cl$sex <- ifelse(sx == "male", "male", ifelse(sx == "female", "female", NA))
  cl$age <- suppressWarnings(as.numeric(clean_chr(cl_raw[[c_age]])))
  cl$OS_time <- suppressWarnings(as.numeric(clean_chr(cl_raw[[c_os]]))) / 365.25
  cl$OS_status <- suppressWarnings(as.numeric(clean_chr(cl_raw[[c_cen]])))
  cl$radio_status <- suppressWarnings(as.numeric(clean_chr(cl_raw[[c_rad]])))
  cl$chemo_status <- suppressWarnings(as.numeric(clean_chr(cl_raw[[c_che]])))
  cl$idh_mutation <- tolower(clean_chr(cl_raw[[c_idh]]))
  cl$codeletion_1p19q <- tolower(clean_chr(cl_raw[[c_cod]]))
  mg <- tolower(clean_chr(cl_raw[[c_mgmt]]))
  cl$mgmt_methylation <- ifelse(mg %in% c("un-methylated","un-methylatated","unmethylated"),
                                "unmethylated",
                         ifelse(mg == "methylated", "methylated", NA))
  if (!is.na(c_sub)) cl$tcga_subtype <- tolower(clean_chr(cl_raw[[c_sub]]))
  cl$ID <- clean_chr(cl$ID)
  keep <- c("OS_status","OS_time","age","sex","grade","histology","prs_type",
            "idh_mutation","codeletion_1p19q","mgmt_methylation",
            "chemo_status","radio_status",
            if ("tcga_subtype" %in% colnames(cl)) "tcga_subtype")
  sv <- cl[, keep, drop = FALSE]
  rownames(sv) <- cl$ID
  n_clin <- nrow(cl)

  # ---------- expression ----------
  ex <- read_txt(ef)
  names(ex)[1] <- "Gene_Name"
  ex$Gene_Name <- clean_chr(ex$Gene_Name)
  ex <- ex[!is.na(ex$Gene_Name), ]
  samples <- setdiff(colnames(ex), "Gene_Name")
  ex[samples] <- lapply(ex[samples], function(x) suppressWarnings(as.numeric(x)))
  # 同基因多行 -> 取最大值(与 get_expr_data 默认 max 一致)
  ex <- ex %>%
    group_by(Gene_Name) %>%
    summarise(across(all_of(samples), ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)),
              .groups = "drop")
  genes <- ex$Gene_Name
  m <- as.matrix(ex[, samples, drop = FALSE]); rownames(m) <- genes
  if (sp$transform == "log2") m <- log2(m + 1)

  # 样本取交集(表达 ∩ 临床)
  common <- intersect(colnames(m), rownames(sv))
  m <- m[, common, drop = FALSE]
  sv <- sv[common, , drop = FALSE]
  cat(sprintf("%s: genes=%d, expr samples=%d, clinical=%d, intersection=%d\n",
              sp$acc, nrow(m), length(samples), nrow(cl), length(common)))

  expr_df <- data.frame(ID_REF = rownames(m), m, check.names = FALSE,
                        stringsAsFactors = FALSE)
  rownames(expr_df) <- NULL        # 与 GEO expr rds 一致: 仅用自动行名
  saveRDS(expr_df, file.path(root, "data/expr", paste0(sp$acc, ".rds")))
  gpl_df <- data.frame(gene_id = sym2entrez(rownames(m)),
                       row.names = rownames(m), stringsAsFactors = FALSE)
  saveRDS(gpl_df, file.path(root, "data/processed/gpl", paste0(sp$gpl, ".rds")))
  saveRDS(sv, file.path(root, "data/processed/surv", paste0(sp$acc, "_surv.rds")))

  report <- c(report,
    sprintf("## %s (%s)", sp$acc, sp$method),
    sprintf("- 表达: %d 基因 × %d 样本 (变化: %s)", nrow(m), ncol(m), sp$transform),
    sprintf("- 临床: %d 行; 交集样本 %d", n_clin, length(common)),
    sprintf("- 事件(OS): %d ; 中位随访(年): %.2f",
            sum(sv$OS_status == 1, na.rm = TRUE),
            stats::median(sv$OS_time, na.rm = TRUE)),
    sprintf("- grade: %s", paste(names(table(sv$grade)), table(sv$grade),
                                 sep = "=", collapse = ", ")),
    sprintf("- histology: %s", paste(names(table(sv$histology)),
                                     table(sv$histology), sep = "=", collapse = ", ")),
    sprintf("- sex: %s", paste(names(table(sv$sex)), table(sv$sex),
                               sep = "=", collapse = ", ")),
    "")
}
writeLines(report, file.path(root, "pipeline/out/REPORT_CGGA_parse.md"))
cat("done. report -> pipeline/out/REPORT_CGGA_parse.md\n")
}

# ---------------------------------------------------------------------------
# run_90_build_gse108474_suppl()  <-  verbatim pipeline/R/90_build_gse108474_suppl.R
# ---------------------------------------------------------------------------
run_90_build_gse108474_suppl <- function() {
# 90_build_gse108474_suppl.R -------------------------------------------------
# GSE108474 (REMBRANDT, glioma, GPL570) — build pheno / expr / surv from the
# series matrix + the series' own supplementary clinical files.
#
# 背景
#   GSE108474 的 series matrix 里 **没有** 生存列（!Sample_characteristics_ch1 只有
#   provider / disease / tissue / tumor grade / extract name / assay name），
#   事件与时间都在 supplementary：
#     GSE108474_REMBRANDT_clinical.data.txt.gz          (SUBJECT_ID 级 EVENT_OS / OVERALL_SURVIVAL_MONTHS)
#     GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz
#
# join key
#   !Sample_title  ==  clinical$SUBJECT_ID        (541/550 精确)
#   9 例 title 形如 "<SUBJECT>_duplicate[ B| C_T| C_NT]" 的重复切片样本，
#   去掉 "_duplicate*" 后缀后 550/550 全部命中 -> 采用该规则。
#
# 输出（全部沿用既有格式）
#   data/pheno/GSE108474.rds                     (与 01_parse_gse.R 完全同构)
#   data/expr/GSE108474.rds                      ID_REF + 每 GSM 一列；探针**不折叠**
#                                                (与 GSE12417-GPL570 / GSE39582 等一致)
#   pipeline/out/GSE108474_extra.rds             行名 = GSM 的补充临床表
#                                                (对应 GSE14520_GPL3921_extra.rds 的用法)
#   data/processed/surv/GSE108474_surv.rds       OS_status(0/1) + OS_time(年) + 临床列
#                                                (遵循 03_surv_table.R 约定)
#
# 用法: Rscript pipeline/R/90_build_gse108474_suppl.R [ROOT]      (默认 CPAS root)
# 注意: 本脚本 **不** 写 data/dataset_info.csv / .rda。
# ----------------------------------------------------------------------------
suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(Biobase)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT)
setwd(ROOT)

ACC   <- "GSE108474"
RAW   <- file.path(ROOT, "data/raw", paste0(ACC, ".gz"))
SUPPL <- file.path(ROOT, "data/suppl")
CLIN  <- file.path(SUPPL, "GSE108474_REMBRANDT_clinical.data.txt.gz")
BMAP  <- file.path(SUPPL, "GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz")
stopifnot(file.exists(RAW), file.exists(CLIN), file.exists(BMAP))

# ---- 1. parseGSEMatrix（与 01_parse_gse.R 中的实现逐字一致，去掉 GPL/csv 分支）----
parseGSEMatrix <- function(local_file, AnnotGPL=FALSE, getGPL=TRUE,
                           parseCharacteristics=TRUE) {
  dat <- readr::read_lines(local_file)
  series_header_row_count <- sum(grepl("^!Series_", dat))
  sample_header_start <- grep("^!Sample_", dat)[1]
  samples_header_row_count <- sum(grepl("^!Sample_", dat))
  series_table_begin_line <- grep("^!series_matrix_table_begin", dat)
  if(length(series_table_begin_line) != 1)
    stop("parsing failed--expected only one '!series_matrix_table_begin'")

  header <- read.table(local_file, sep="\t", header=FALSE, nrows=series_header_row_count)
  tmpdat <- read.table(local_file, sep="\t", header=FALSE,
                       nrows=samples_header_row_count, skip=sample_header_start-1)
  headertmp <- t(header); headerdata <- rbind(data.frame(), headertmp[-1,])
  colnames(headerdata) <- sub('!Series_','',as.character(header[,1]))

  tmptmp <- t(tmpdat)
  sampledat <- rbind(data.frame(), tmptmp[-1,])
  colnames(sampledat) <- make.unique(sub('!Sample_','',as.character(tmpdat[,1])))
  sampledat[['geo_accession']] <- as.character(sampledat[['geo_accession']])
  rownames(sampledat) <- sampledat[['geo_accession']]

  if(length(grep('characteristics_ch',colnames(sampledat)))>0 && parseCharacteristics) {
    pd <- sampledat %>%
      dplyr::select(dplyr::contains('characteristics_ch')) %>%
      dplyr::mutate(accession = rownames(.)) %>%
      mutate_all(as.character) %>%
      tidyr::gather(characteristics, kvpair, -accession) %>%
      dplyr::mutate(kvpair = dplyr::case_when(
        grepl(':', kvpair) ~ kvpair,
        grepl(' = ', kvpair, fixed=TRUE) ~ sub(' = ', ':', kvpair, fixed=TRUE),
        grepl('^\\[', trimws(kvpair)) ~ {
          x <- trimws(kvpair); x <- sub('^\\[', '', x); x <- sub('\\]$', '', x)
          k <- sub('=.*$', '', x)
          v <- sub('^[^=]*=\\s*', '', x)
          v <- sub('\\s*\\(.*$', '', v)
          paste0(trimws(k), ':', trimws(v)) },
        TRUE ~ kvpair)) %>%
      dplyr::filter(grepl(':',kvpair) & !is.na(kvpair))
    if(nrow(pd)) {
      pd <- pd %>%
        dplyr::mutate(characteristics = ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
        tidyr::separate(kvpair, into=c('k','v'), sep=":", fill='right', extra='merge') %>%
        dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
        dplyr::select(-characteristics) %>%
        dplyr::filter(!is.na(v)) %>%
        dplyr::group_by(accession,k) %>%
        dplyr::mutate(v = paste0(trimws(v), collapse=";")) %>%
        unique() %>%
        tidyr::spread(k,v)
    } else {
      pd <- pd %>% dplyr::select(accession)
    }
    sampledat <- sampledat %>% dplyr::left_join(pd, by=c('geo_accession'='accession'))
  }

  Sys.setenv("VROOM_CONNECTION_SIZE" = 2147483647)
  datamat <- read_tsv(local_file, quote='"',
                      na=c('NA','null','NULL','Null'),
                      skip = series_table_begin_line,
                      comment = '!series_matrix_table_end',
                      skip_empty_rows = FALSE,
                      show_col_types = FALSE)
  tmprownames <- datamat[[1]]
  datamat <- as.matrix(datamat[!is.na(tmprownames),-1])
  rownames(datamat) <- tmprownames[!is.na(tmprownames)]
  datamat <- as.matrix(datamat)
  rownames(sampledat) <- colnames(datamat)
  if(is.null(nrow(datamat))) {
    tmpnames <- names(datamat)
    rownames(sampledat) <- tmpnames
    datamat <- matrix(nrow=0, ncol=nrow(sampledat)); colnames(datamat) <- tmpnames
  }
  eset <- new('ExpressionSet', phenoData=as(sampledat,'AnnotatedDataFrame'),
              annotation=as.character(sampledat[1, grep('platform_id',colnames(sampledat),
                                                        ignore.case=TRUE)]),
              exprs=as.matrix(datamat))
  return(list(GPL=as.character(sampledat[1, grep('platform_id',colnames(sampledat),
                                                 ignore.case=TRUE)]),
              eset=eset))
}

message("parsing ", RAW, " ...")
res  <- parseGSEMatrix(RAW, AnnotGPL=FALSE, getGPL=FALSE, parseCharacteristics=TRUE)
eset <- res$eset
pd   <- Biobase::pData(eset)
message("GPL=", Biobase::annotation(eset), "  samples=", nrow(pd),
        "  probes=", nrow(Biobase::exprs(eset)))

# ---- 2. expr / pheno（与 01_parse_gse.R 输出同构；探针不折叠）----
expr_data <- as.data.frame(Biobase::exprs(eset)) %>% tibble::rownames_to_column("ID_REF")
# 该 series matrix 里 8 个 "*_duplicate*" GSM 的整列为空（data_row_count = 0，
# 即 GEO 上没有任何表达值）。既有 146 个 expr 文件都没有「整列全空」的先例，
# 保留会破坏镜像表的隐含不变量，故剔除这 8 列（均为无数据的重复切片样本）。
empty_cols <- names(expr_data)[-1][colSums(is.na(expr_data[, -1, drop = FALSE])) == nrow(expr_data)]
if (length(empty_cols)) {
  message("dropping ", length(empty_cols), " all-NA sample columns: ",
          paste(empty_cols, collapse = ", "))
  expr_data <- expr_data[, setdiff(colnames(expr_data), empty_cols), drop = FALSE]
}
saveRDS(expr_data, file = file.path(ROOT, "data/expr", paste0(ACC, ".rds")))
pheno_data <- pd %>% dplyr::select(-contains("characteristics"))
saveRDS(pheno_data, file = file.path(ROOT, "data/pheno", paste0(ACC, ".rds")))
message("saved data/expr/", ACC, ".rds (", nrow(expr_data), " x ", ncol(expr_data)-1,
        ") and data/pheno/", ACC, ".rds (", nrow(pheno_data), ")")

# ---- 3. 补充临床：SUBJECT_ID -> EVENT_OS / OVERALL_SURVIVAL_MONTHS ----
clin <- read.delim(CLIN, check.names=FALSE, stringsAsFactors=FALSE,
                   na.strings=c("", "NA"))
bmap <- read.delim(BMAP, check.names=FALSE, stringsAsFactors=FALSE,
                   na.strings=c("", "NA"))
stopifnot(!anyDuplicated(clin$SUBJECT_ID))

gsm    <- rownames(pd)
titles <- trimws(as.character(pd$title))
subj   <- sub("_duplicate.*$", "", titles)          # 重复切片样本还原到患者
n_exact <- sum(titles %in% clin$SUBJECT_ID)
n_after <- sum(subj   %in% clin$SUBJECT_ID)
message(sprintf("join key: title==SUBJECT_ID exact %d/%d ; after stripping '_duplicate*' %d/%d",
                n_exact, length(titles), n_after, length(subj)))
if (n_after != length(subj)) stop("unmatched samples remain: ",
                                  paste(titles[!subj %in% clin$SUBJECT_ID], collapse=","))

m  <- match(subj, clin$SUBJECT_ID)
ex <- data.frame(
  subject_id          = subj,
  EVENT_OS            = clin$EVENT_OS[m],
  OVERALL_SURVIVAL_MONTHS = clin$OVERALL_SURVIVAL_MONTHS[m],
  DISEASE_TYPE        = clin$DISEASE_TYPE[m],
  WHO_GRADE           = clin$WHO_GRADE[m],
  GENDER              = clin$GENDER[m],
  AGE_RANGE           = clin$AGE_RANGE[m],
  RACE                = clin$RACE[m],
  INSTITUTION_NAME    = clin$INSTITUTION_NAME[m],
  KARNOFSKY           = clin$KARNOFSKY[m],
  row.names = gsm, stringsAsFactors = FALSE)
attr(ex, "join_key") <- "!Sample_title (stripped of '_duplicate*') == clinical.data SUBJECT_ID"
attr(ex, "match_rate") <- sprintf("%d/%d exact, %d/%d after strip", n_exact,
                                  length(titles), n_after, length(subj))
saveRDS(ex, file.path(ROOT, "pipeline/out/GSE108474_extra.rds"))
message("saved pipeline/out/GSE108474_extra.rds (", nrow(ex), " x ", ncol(ex), ")")

# ---- 4. surv 表（03_surv_table.R 约定：tumor-only 行 / OS 单位月 -> 年 / 临床规范化）----
clinic <- pd[, c("title", "disease:ch1", "tumor grade:ch1", "extract name:ch1"),
             drop = FALSE]
d <- cbind(as.data.frame(ex, stringsAsFactors = FALSE),
           histology_geo = as.character(clinic[["disease:ch1"]]),
           grade_geo     = as.character(clinic[["tumor grade:ch1"]]))
stopifnot(identical(rownames(d), gsm))

# 肿瘤样本：clinical DISEASE_TYPE == "NON_TUMOR" 为癌旁/正常组织 -> 剔除
n0 <- nrow(d)
keep <- !is.na(d$DISEASE_TYPE) & toupper(d$DISEASE_TYPE) != "NON_TUMOR"
d <- d[keep, , drop = FALSE]
message("row filter (DISEASE_TYPE != NON_TUMOR) kept ", nrow(d), "/", n0)

# OS 事件/时间：EVENT_OS 已是 0/1；时间单位 = 月 -> 年
d$OS_status <- suppressWarnings(as.numeric(d$EVENT_OS))
d$OS_time   <- suppressWarnings(as.numeric(d$OVERALL_SURVIVAL_MONTHS)) / 12
d$OS_status[!d$OS_status %in% c(0, 1)] <- NA_real_
d$OS_time[!is.na(d$OS_time) & d$OS_time < 0] <- 0

# 临床规范化（复刻 03_surv_table.R 的 std_* 规则）
.roman <- c("0"="0","1"="I","2"="II","3"="III","4"="IV","5"="V","X"="X","x"="X")
std_grade <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v[v %in% c("", "NA", "N/A", "NULL", "?", "UNKNOWN")] <- NA_character_
  vapply(v, function(s) {
    if (is.na(s)) return(NA_character_)
    s2 <- sub("^G", "", s)
    if (grepl("^[1-4]$", s2)) return(.roman[[s2]])
    if (grepl("^[1-4]$", s))  return(.roman[[s]])
    if (s %in% c("I","II","III","IV","LOW","HIGH","LOW GRADE","HIGH GRADE")) return(s)
    m <- regmatches(s, regexpr("[1-4]", s))
    if (length(m) && nzchar(m)) return(.roman[[m]])
    if (grepl("WELL", s)) return("I")
    if (grepl("MODERATE", s)) return("II")
    if (grepl("POOR", s)) return("III")
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}
std_sex <- function(x) {
  v <- tolower(trimws(as.character(x)))
  vapply(v, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    if (s %in% c("m","male","man","men")) return("male")
    if (s %in% c("f","female","woman","women")) return("female")
    s
  }, character(1), USE.NAMES = FALSE)
}
# AGE_RANGE 为 5 岁分箱文本（如 "50-54"）-> 数值取组中值；原始分箱另存 age_range
age_mid <- function(x) {
  v <- trimws(as.character(x))
  lo <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", v)))
  hi <- suppressWarnings(as.numeric(sub("^[0-9]+-([0-9]+)$", "\\1", v)))
  out <- ifelse(is.na(lo) | is.na(hi), NA_real_, (lo + hi) / 2)
  out
}

d$age       <- age_mid(d$AGE_RANGE)
d$age_range <- d$AGE_RANGE
d$sex       <- std_sex(d$GENDER)
# grade: WHO_GRADE 优先（已是罗马数字），GEO "grade 2" 兜底
g <- ifelse(!is.na(d$WHO_GRADE) & nzchar(as.character(d$WHO_GRADE)),
            as.character(d$WHO_GRADE), as.character(d$grade_geo))
d$grade <- std_grade(g)
# histology: DISEASE_TYPE，UNKNOWN/UNCLASSIFIED -> NA
h <- toupper(trimws(as.character(d$DISEASE_TYPE)))
h[h %in% c("", "NA", "UNKNOWN", "UNCLASSIFIED")] <- NA_character_
d$histology <- h

out <- d[, c("OS_status", "OS_time", "age", "age_range", "sex", "grade", "histology")]
out <- as.data.frame(lapply(out, function(col) {
  if (is.character(col)) col[col %in% c("", "NA", "N/A", "NULL", "UNK", "UNKNOWN")] <- NA_character_
  col
}), check.names = FALSE, stringsAsFactors = FALSE)
rownames(out) <- rownames(d)
rownames(out) <- sub("^X", "", rownames(out)); rownames(out) <- gsm[keep]

out_file <- file.path(ROOT, "data/processed/surv", paste0(ACC, "_surv.rds"))
saveRDS(out, out_file)

n_ok <- sum(!is.na(out$OS_status) & !is.na(out$OS_time))
expr_cols <- setdiff(colnames(expr_data), "ID_REF")
n_analysable <- sum(!is.na(out$OS_status) & !is.na(out$OS_time) & rownames(out) %in% expr_cols)
subj_ok <- unique(as.character(d$subject_id[!is.na(out$OS_status) & !is.na(out$OS_time)]))
dup_subj <- names(which(table(as.character(d$subject_id[!is.na(out$OS_status) & !is.na(out$OS_time)])) > 1))
cat("\n================ GSE108474 build summary ================\n")
cat("expr  : ", nrow(expr_data), " probes x ", ncol(expr_data)-1, " samples (",
    length(empty_cols), " all-NA columns dropped)\n", sep="")
cat("pheno : ", nrow(pheno_data), "\n", sep="")
cat("surv  : ", nrow(out), " rows (tumor only)\n", sep="")
cat("effective n (OS status & time non-missing): ", n_ok, "\n", sep="")
cat("N = expr ∩ usable OS: ", n_analysable, "\n", sep="")
cat("unique patients among usable: ", length(subj_ok),
    " ; patients contributing >1 usable sample: ", length(dup_subj),
    " (", sum(table(as.character(d$subject_id[!is.na(out$OS_status) & !is.na(out$OS_time)]))[dup_subj]),
    " samples )\n", sep="")
cat("OS events: ", sum(out$OS_status == 1, na.rm=TRUE),
    "  censored: ", sum(out$OS_status == 0, na.rm=TRUE),
    "  NA: ", sum(is.na(out$OS_status)), "\n", sep="")
cat("OS_time (years) range: ",
    paste(round(range(out$OS_time, na.rm=TRUE), 3), collapse=" - "), "\n", sep="")
cat("columns: ", paste(colnames(out), collapse=", "), "\n", sep="")
cat("saved: ", out_file, "\n", sep="")
cat("========================================================\n")
}

# ---------------------------------------------------------------------------
# run_91_complete_gse14520_surv()  <-  verbatim pipeline/R/91_complete_gse14520_surv.R
# ---------------------------------------------------------------------------
run_91_complete_gse14520_surv <- function() {
# 91_complete_gse14520_surv.R -------------------------------------------------
# GSE14520 (HCC, GPL3921) — 用 series 自带的补充临床文件补全 surv 表。
#
# 现状
#   data/processed/surv/GSE14520-GPL3921_surv.rds = 225 行（GPL3921 的 225 例 Tumor）
#   列 = OS_status / OS_time / RFS_status / RFS_time / sex / age / stage
#   其中 221 例有完整 OS+RFS 对，4 例（GSM363228, GSM363334, GSM712532, GSM712534）
#   在补充文件里本来就没有生存信息。
#
# 补充文件
#   data/suppl/GSE14520_Extra_Supplement.txt.gz （GEO GSE14520 的官方 Extra_Supplement）
#   488 行 = GPL3921 445 + GPL571 43；列含 LCS ID / Affy_GSM / Agilent_GSM /
#   Tissue Type / Survival status+months / Recurr status+months / Gender / Age /
#   HBV / ALT / Main Tumor Size / Multinodular / Cirrhosis / TNM / BCLC / CLIP / AFP
#
# join key : supp$Affy_GSM  <->  surv 行名 (GSM)       匹配率 225/225
#
# 输出（**不覆盖**原文件）
#   data/processed/surv/GSE14520-GPL3921_surv_completed.rds
#     = 原 7 列原样保留 + 补充文件的临床列
#     （若有 GPL3921 Tumor 且带可用生存对、但原表缺失的样本，则追加为新行）
#
# 约束：GPL571 的队列已被作者从工具中移除，本脚本**不**读取/构建 GPL571 队列，
#       也不写任何 GSE14520-GPL571 产物（仅用于核对计数）。
#
# 用法: Rscript pipeline/R/91_complete_gse14520_surv.R [ROOT]    (不写 dataset_info)
# ----------------------------------------------------------------------------
suppressMessages({library(dplyr)})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT)
setwd(ROOT)

ACC_3921 <- "GSE14520-GPL3921"
SUPPL    <- file.path(ROOT, "data/suppl/GSE14520_Extra_Supplement.txt.gz")
OLD      <- file.path(ROOT, "data/processed/surv", paste0(ACC_3921, "_surv.rds"))
OUT      <- file.path(ROOT, "data/processed/surv",
                      paste0(ACC_3921, "_surv_completed.rds"))
stopifnot(file.exists(SUPPL), file.exists(OLD))

supp <- read.delim(SUPPL, check.names=FALSE, stringsAsFactors=FALSE,
                   na.strings=c("", "NA"))
supp$Affy_GSM <- trimws(as.character(supp$Affy_GSM))
supp$`Tissue Type` <- trimws(as.character(supp$`Tissue Type`))
if (anyDuplicated(supp$Affy_GSM)) stop("supplement Affy_GSM not unique")

old <- readRDS(OLD)
p3921 <- rownames(readRDS(file.path(ROOT, "data/pheno", paste0(ACC_3921, ".rds"))))

# ---- 匹配率 ----
m <- match(rownames(old), supp$Affy_GSM)
cat(sprintf("join key: supp$Affy_GSM == rownames(old surv)   matched %d/%d (%.1f%%)\n",
            sum(!is.na(m)), nrow(old), 100*sum(!is.na(m))/nrow(old)))
stopifnot(!anyNA(m))

# ---- (b) 补充文件中 GPL3921 的 Tumor 行（带可用生存对）中，原表缺失的样本 ----
supp3921_tum <- supp[supp$Affy_GSM %in% p3921 & supp$`Tissue Type` == "Tumor", , drop=FALSE]
usable <- !is.na(supp3921_tum$`Survival status`) & !is.na(supp3921_tum$`Survival months`)
missing_new <- setdiff(supp3921_tum$Affy_GSM[usable], rownames(old))
cat(sprintf("GPL3921 Tumor rows in suppl: %d ; with usable OS pair: %d ; already in old table: %d ; NEW: %d\n",
            nrow(supp3921_tum), sum(usable), sum(supp3921_tum$Affy_GSM[usable] %in% rownames(old)),
            length(missing_new)))

# ---- (c) 组装：旧列原样 + 新临床列 ----
new_cols <- data.frame(
  lcs_id             = supp$`LCS ID`[m],
  sample_id          = supp$ID[m],
  tissue             = supp$`Tissue Type`[m],
  risk_signature     = supp$`Predicted risk Metastasis Signature`[m],
  cgh_survival_group = supp$`CGH_survival_group`[m],
  hbv_status         = supp$`HBV viral status`[m],
  alt                = supp$`ALT(>/<=50U/L)`[m],
  tumor_size         = supp$`Main Tumor Size (>/<=5 cm)`[m],
  multinodular       = supp$Multinodular[m],
  cirrhosis          = supp$Cirrhosis[m],
  tnm_stage          = supp$`TNM staging`[m],
  bclc               = supp$`BCLC staging`[m],
  clip               = supp$`CLIP staging`[m],
  afp                = supp$`AFP (>/<=300ng/ml)`[m],
  row.names = rownames(old), stringsAsFactors = FALSE, check.names = FALSE)

out <- cbind(old, new_cols)

# 补充文件用 "." 表示缺失（HBV/TNM/BCLC/CLIP/AFP/tumor size 共 74 处）-> NA
new_names <- colnames(new_cols)
out[new_names] <- lapply(out[new_names], function(x) {
  x <- as.character(x)
  x[trimws(x) %in% c("", ".", "NA", "N/A", "NULL", "UNK", "UNKNOWN")] <- NA_character_
  x
})

if (length(missing_new)) {
  mm <- match(missing_new, supp$Affy_GSM)
  extra <- data.frame(
    OS_status  = suppressWarnings(as.numeric(supp$`Survival status`[mm])),
    OS_time    = suppressWarnings(as.numeric(supp$`Survival months`[mm])) / 12,
    RFS_status = suppressWarnings(as.numeric(supp$`Recurr status`[mm])),
    RFS_time   = suppressWarnings(as.numeric(supp$`Recurr months`[mm])) / 12,
    sex        = supp$Gender[mm], age = supp$Age[mm],
    stage      = supp$`TNM staging`[mm],
    lcs_id = supp$`LCS ID`[mm], sample_id = supp$ID[mm],
    tissue = supp$`Tissue Type`[mm],
    risk_signature = supp$`Predicted risk Metastasis Signature`[mm],
    cgh_survival_group = supp$`CGH_survival_group`[mm],
    hbv_status = supp$`HBV viral status`[mm], alt = supp$`ALT(>/<=50U/L)`[mm],
    tumor_size = supp$`Main Tumor Size (>/<=5 cm)`[mm],
    multinodular = supp$Multinodular[mm], cirrhosis = supp$Cirrhosis[mm],
    tnm_stage = supp$`TNM staging`[mm], bclc = supp$`BCLC staging`[mm],
    clip = supp$`CLIP staging`[mm], afp = supp$`AFP (>/<=300ng/ml)`[mm],
    row.names = missing_new, stringsAsFactors = FALSE, check.names = FALSE)
  out <- rbind(out, extra[, colnames(out)])
}

saveRDS(out, OUT)

# ---- 报告 ----
ev <- function(x) sum(x == 1, na.rm = TRUE)
cat("\n================ GSE14520-GPL3921 completion summary ================\n")
cat("rows           : ", nrow(old), " -> ", nrow(out), "\n", sep="")
cat("OS  events     : ", ev(old$OS_status), " -> ", ev(out$OS_status),
    "   (censored ", sum(old$OS_status==0,na.rm=TRUE), " -> ", sum(out$OS_status==0,na.rm=TRUE),
    ", NA ", sum(is.na(old$OS_status)), " -> ", sum(is.na(out$OS_status)), ")\n", sep="")
cat("RFS events     : ", ev(old$RFS_status), " -> ", ev(out$RFS_status),
    "   (censored ", sum(old$RFS_status==0,na.rm=TRUE), " -> ", sum(out$RFS_status==0,na.rm=TRUE),
    ", NA ", sum(is.na(old$RFS_status)), " -> ", sum(is.na(out$RFS_status)), ")\n", sep="")
cat("usable OS pair : ", sum(!is.na(old$OS_status) & !is.na(old$OS_time)), " -> ",
    sum(!is.na(out$OS_status) & !is.na(out$OS_time)), "\n", sep="")
cat("usable RFS pair: ", sum(!is.na(old$RFS_status) & !is.na(old$RFS_time)), " -> ",
    sum(!is.na(out$RFS_status) & !is.na(out$RFS_time)), "\n", sep="")
cat("columns        : ", paste(colnames(old), collapse=", "), "\n", sep="")
cat("              -> ", paste(colnames(out), collapse=", "), "\n", sep="")
cat("\nnew-column non-NA counts:\n")
print(sapply(new_cols, function(x) sum(!is.na(x))))
cat("\nsaved: ", OUT, "\n", sep="")
cat("====================================================================\n")
}

# ---------------------------------------------------------------------------
# run_92_build_geo_expansion_expr_pheno()  <-  verbatim pipeline/R/92_build_geo_expansion_expr_pheno.R
# ---------------------------------------------------------------------------
run_92_build_geo_expansion_expr_pheno <- function() {
# 92_build_geo_expansion_expr_pheno.R ----------------------------------------
# GEO 扩展批次（A 审计后）第 1 步：把每个新队列的
#   data/raw/<ACC>.gz        series matrix（本地缺失时下载）
#   data/pheno/<ACC>.rds     与 01_parse_gse.R 同构的临床表
#   data/expr/<ACC>.rds      ID_REF + 每 GSM 一列（探针**不折叠**，与既有文件一致）
#   data/processed/gpl/<GPL>.rds   探针 -> ENTREZ gene_id（缺失时构建）
# 全部沿用既有脚本的格式约定（01_parse_gse.R / 90_build_gse108474_suppl.R /
# 02_gpl_map.R）；本脚本 **不** 写 data/dataset_info.csv / .rda，也不建 surv 表。
#
# 表达式来源有两类：
#   expr_mode = "matrix"  —— 表达值就在 series matrix 里（01/90 的老路径）
#   expr_mode = "suppl"   —— GEO 的 series matrix 是空表（!Sample_type = SRA），
#                            表达值在 supplementary 文件里，需按各自的连接键装配：
#       GSE159067  GSE159067_IHN_log2cpm_data.txt.gz   列名 = !Sample_description
#       GSE162520  GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz 列名 = !Sample_description
#                  （分号分隔 + 逗号小数点，需 as.numeric(sub(",", ".", x))）
#       GSE271517  GSE271517_Sample_Counts.csv.gz      列名 = T<编号>，与 !Sample_title
#                  的编号对齐；行名 ENSG -> ENTREZ（REF_ID/www/biomart_table_Hs.Rds）
#
# 用法: Rscript pipeline/R/92_build_geo_expansion_expr_pheno.R [ROOT] [ACC ...]
# ----------------------------------------------------------------------------
suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(Biobase)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
ONLY <- if (length(args) >= 2) args[-1] else NULL

SUPPL <- file.path(ROOT, "data/suppl")
dir.create(SUPPL, showWarnings = FALSE, recursive = TRUE)

# ---- 队列清单（acc, GPL, 表达式来源）--------------------------------------
SPEC <- list(
  list(acc = "GSE65858",  gpl = "GPL10558", mode = "matrix"),
  list(acc = "GSE117973", gpl = "GPL10558", mode = "matrix"),
  list(acc = "GSE27020",  gpl = "GPL96",    mode = "matrix"),
  list(acc = "GSE65904",  gpl = "GPL10558", mode = "matrix"),
  list(acc = "GSE22153",  gpl = "GPL6102",  mode = "matrix"),
  list(acc = "GSE71118",  gpl = "GPL570",   mode = "matrix"),
  list(acc = "GSE30929",  gpl = "GPL96",    mode = "matrix"),
  list(acc = "GSE22138",  gpl = "GPL570",   mode = "matrix"),
  list(acc = "GSE23501",  gpl = "GPL570",   mode = "matrix"),
  list(acc = "GSE32918",  gpl = "GPL8432",  mode = "matrix"),
  list(acc = "GSE31312",  gpl = "GPL570",   mode = "matrix"),
  list(acc = "GSE198430", gpl = "GPL32055", mode = "matrix"),
  list(acc = "GSE198431", gpl = "GPL32055", mode = "matrix"),
  list(acc = "GSE248835", gpl = "GPL33963", mode = "suppl"),
  # GeoMx：supplementary 是 ROI 级表达，必须按 patient 聚合（审计 §3.6：ROI != 患者）
  list(acc = "GSE325123", gpl = "GPL24676", mode = "roi_patient"),
  list(acc = "GSE159067", gpl = "GPL18573", mode = "suppl"),
  list(acc = "GSE162520", gpl = "GPL18573", mode = "suppl"),
  list(acc = "GSE271517", gpl = "GPL24676", mode = "suppl"),
  list(acc = "GSE183088", gpl = "GPL30570", mode = "matrix")
)

# ---- 与 01_parse_gse.R 逐字一致的 series matrix 解析器 ---------------------
parseGSEMatrix <- function(local_file, AnnotGPL = FALSE, getGPL = TRUE,
                           parseCharacteristics = TRUE) {
  dat <- readr::read_lines(local_file)
  series_header_row_count <- sum(grepl("^!Series_", dat))
  sample_header_start <- grep("^!Sample_", dat)[1]
  samples_header_row_count <- sum(grepl("^!Sample_", dat))
  series_table_begin_line <- grep("^!series_matrix_table_begin", dat)
  if (length(series_table_begin_line) != 1)
    stop("parsing failed--expected only one '!series_matrix_table_begin'")

  header <- read.table(local_file, sep = "\t", header = FALSE, nrows = series_header_row_count)
  tmpdat <- read.table(local_file, sep = "\t", header = FALSE,
                       nrows = samples_header_row_count, skip = sample_header_start - 1)
  headertmp <- t(header); headerdata <- rbind(data.frame(), headertmp[-1, ])
  colnames(headerdata) <- sub('!Series_', '', as.character(header[, 1]))

  tmptmp <- t(tmpdat)
  sampledat <- rbind(data.frame(), tmptmp[-1, ])
  colnames(sampledat) <- make.unique(sub('!Sample_', '', as.character(tmpdat[, 1])))
  sampledat[['geo_accession']] <- as.character(sampledat[['geo_accession']])
  rownames(sampledat) <- sampledat[['geo_accession']]

  if (length(grep('characteristics_ch', colnames(sampledat))) > 0 && parseCharacteristics) {
    pd <- sampledat %>%
      dplyr::select(dplyr::contains('characteristics_ch')) %>%
      dplyr::mutate(accession = rownames(.)) %>%
      mutate_all(as.character) %>%
      tidyr::gather(characteristics, kvpair, -accession) %>%
      dplyr::mutate(kvpair = dplyr::case_when(
        grepl(':', kvpair) ~ kvpair,
        grepl(' = ', kvpair, fixed = TRUE) ~ sub(' = ', ':', kvpair, fixed = TRUE),
        grepl('^\\[', trimws(kvpair)) ~ {
          x <- trimws(kvpair); x <- sub('^\\[', '', x); x <- sub('\\]$', '', x)
          k <- sub('=.*$', '', x)
          v <- sub('^[^=]*=\\s*', '', x)
          v <- sub('\\s*\\(.*$', '', v)
          paste0(trimws(k), ':', trimws(v)) },
        TRUE ~ kvpair)) %>%
      dplyr::filter(grepl(':', kvpair) & !is.na(kvpair))
    if (nrow(pd)) {
      pd <- pd %>%
        dplyr::mutate(characteristics = ifelse(grepl('_ch2', characteristics), 'ch2', 'ch1')) %>%
        tidyr::separate(kvpair, into = c('k', 'v'), sep = ":", fill = 'right', extra = 'merge') %>%
        dplyr::mutate(k = paste(k, characteristics, sep = ":")) %>%
        dplyr::select(-characteristics) %>%
        dplyr::filter(!is.na(v)) %>%
        dplyr::group_by(accession, k) %>%
        dplyr::mutate(v = paste0(trimws(v), collapse = ";")) %>%
        unique() %>%
        tidyr::spread(k, v)
    } else {
      pd <- pd %>% dplyr::select(accession)
    }
    sampledat <- sampledat %>% dplyr::left_join(pd, by = c('geo_accession' = 'accession'))
  }

  # 注意：这里把 !Sample_description 等行原样保留在 sampledat 里；SRA 型 series
  # （矩阵体为空）的样本标识（如 GSE159067 的 "125-AN14015199-_1"）就在
  # !Sample_description 里，是后续与 supplementary 表达表连接的唯一键。
  Sys.setenv("VROOM_CONNECTION_SIZE" = 2147483647)
  datamat <- read_tsv(local_file, quote = '"',
                      na = c('NA', 'null', 'NULL', 'Null'),
                      skip = series_table_begin_line,
                      comment = '!series_matrix_table_end',
                      skip_empty_rows = FALSE, show_col_types = FALSE)
  tmprownames <- datamat[[1]]
  datamat <- as.matrix(datamat[!is.na(tmprownames), -1])
  rownames(datamat) <- tmprownames[!is.na(tmprownames)]
  datamat <- as.matrix(datamat)
  if (ncol(datamat)) rownames(sampledat) <- colnames(datamat)
  eset <- new('ExpressionSet', phenoData = as(sampledat, 'AnnotatedDataFrame'),
              annotation = as.character(sampledat[1, grep('platform_id', colnames(sampledat),
                                                          ignore.case = TRUE)]),
              exprs = as.matrix(datamat))
  list(GPL = as.character(sampledat[1, grep('platform_id', colnames(sampledat),
                                            ignore.case = TRUE)]), eset = eset)
}

# ---- 下载（本地缺失时）-----------------------------------------------------
fetch_matrix <- function(acc) {
  dest <- file.path(ROOT, "data/raw", paste0(acc, ".gz"))
  if (file.exists(dest) && file.size(dest) > 1000) { message("raw exists: ", dest); return(dest) }
  pre <- sub("[0-9]{3}$", "nnn", acc)
  dir <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", pre, acc)
  html <- readLines(dir, warn = FALSE)
  fn <- unique(str_extract(html, "GSE[0-9]+_series_matrix\\.txt\\.gz"))
  fn <- fn[!is.na(fn)]
  if (!length(fn)) stop("no series matrix listed at ", dir)
  url <- paste0(dir, fn[1])
  message("downloading ", url)
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

fetch_suppl <- function(acc, file) {
  dest <- file.path(SUPPL, file)
  if (file.exists(dest) && file.size(dest) > 100) { return(dest) }
  pre <- sub("[0-9]{3}$", "nnn", acc)
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s", pre, acc, file)
  message("downloading ", url)
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

# ---- supplementary 表达表装配 ---------------------------------------------
expr_from_suppl <- function(acc, pd) {
  if (acc == "GSE159067") {
    f <- fetch_suppl(acc, "GSE159067_IHN_log2cpm_data.txt.gz")
    x <- read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
    ids <- x[[1]]; m <- as.data.frame(x[, -1, drop = FALSE], check.names = FALSE)
    colnames(m) <- colnames(x)[-1]
    key <- trimws(as.character(pd[["description"]]))
    idx <- match(key, colnames(m))
    if (anyNA(idx)) stop(acc, ": suppl columns do not match !Sample_description (",
                         sum(is.na(idx)), " unmatched)")
    m <- m[, idx, drop = FALSE]; colnames(m) <- rownames(pd)
    return(data.frame(ID_REF = ids, m, check.names = FALSE))
  }
  if (acc == "GSE162520") {
    f <- fetch_suppl(acc, "GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz")
    # 分号分隔、逗号小数点、首列无列名
    ln <- readLines(gzfile(f), warn = FALSE)
    hdr <- strsplit(ln[1], ";", fixed = TRUE)[[1]]
    cols <- gsub('"', '', hdr)[-1]
    body <- strsplit(ln[-1], ";", fixed = TRUE)
    ids <- vapply(body, function(z) gsub('"', '', z[1]), character(1))
    mat <- matrix(NA_real_, nrow = length(ids), ncol = length(cols))
    for (i in seq_along(body)) {
      v <- gsub('"', '', body[[i]][-1])
      mat[i, ] <- suppressWarnings(as.numeric(sub(",", ".", v, fixed = TRUE)))
    }
    m <- as.data.frame(mat, check.names = FALSE); colnames(m) <- cols
    key <- trimws(as.character(pd[["description"]]))
    idx <- match(key, cols)
    if (anyNA(idx)) stop(acc, ": suppl columns do not match !Sample_description (",
                         sum(is.na(idx)), " unmatched)")
    m <- m[, idx, drop = FALSE]; colnames(m) <- rownames(pd)
    return(data.frame(ID_REF = ids, m, check.names = FALSE))
  }
  if (acc == "GSE271517") {
    f <- fetch_suppl(acc, "GSE271517_Sample_Counts.csv.gz")
    x <- read.csv(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
    ids <- x[[1]]; m <- as.data.frame(x[, -1, drop = FALSE], check.names = FALSE)
    # 列名 = T<编号>；!Sample_title 就是 T<编号>
    key <- str_extract(trimws(as.character(pd[["title"]])), "^T[0-9]+")
    idx <- match(key, colnames(m))
    if (anyNA(idx)) stop(acc, ": suppl columns do not match !Sample_title (",
                         sum(is.na(idx)), " unmatched)")
    m <- m[, idx, drop = FALSE]; colnames(m) <- rownames(pd)
    return(data.frame(ID_REF = ids, m, check.names = FALSE))
  }
  if (acc == "GSE248835") {
    # series matrix 的行标识是 1..817 的探针序号，supplementary 的同名表才是
    # 基因名（nCounter IO360 panel）。用 supplementary 的列名（= !Sample_title）
    # 连接，行名去掉 ".IO360" 后缀作为基因符号，与 GPL33963 注释表一致。
    f <- fetch_suppl(acc, "GSE248835_Normalized_data.csv.gz")
    x <- read.csv(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
    # 保留原始行名（含 ".IO360" 后缀）作为 ID_REF —— 去掉后缀会产生 6 个重名，
    # 而 06_upload_db.R 的 column_to_rownames("ID_REF") 不接受重复行名。
    # 后缀只在 GPL 注释阶段（build_gpl_map）剥掉用于 SYMBOL 查询。
    ids <- as.character(x[[1]])
    m <- as.data.frame(x[, -1, drop = FALSE], check.names = FALSE)
    key <- trimws(as.character(pd[["title"]]))
    idx <- match(key, colnames(m))
    if (anyNA(idx)) stop(acc, ": suppl columns do not match !Sample_title (",
                         sum(is.na(idx)), " unmatched)")
    m <- m[, idx, drop = FALSE]; colnames(m) <- rownames(pd)
    return(data.frame(ID_REF = ids, m, check.names = FALSE))
  }
  stop("no supplementary expression loader for ", acc)
}

# ---- 探针 -> ENTREZ 注释 ---------------------------------------------------
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
ensg2entrez <- function(ens) {
  bm <- readRDS(file.path(ROOT, "REF_ID/www/biomart_table_Hs.Rds"))
  bm_map <- stats::setNames(as.character(bm$Entrez), bm$Ensembl_gene)
  ens <- unique(sub("\\..*$", "", as.character(ens)))
  out <- setNames(rep(NA_character_, length(ens)), ens)
  hit <- intersect(ens, names(bm_map))
  out[hit] <- unname(bm_map[hit])
  out
}
# RefSeq mRNA 访问号（nCounter panel 的常见行标识，如 GPL32055 的 NM_000045.3）
refseq2entrez <- function(ids) {
  ids0 <- sub("\\..*$", "", as.character(ids))
  out <- setNames(rep(NA_character_, length(unique(ids0))), unique(ids0))
  db <- tryCatch(getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"), error = function(e) NULL)
  if (is.null(db)) return(unname(out[ids0]))
  m <- tryCatch(suppressWarnings(AnnotationDbi::select(db, keys = unique(ids0),
                                                      columns = "ENTREZID", keytype = "REFSEQ")),
                error = function(e) NULL)
  if (is.null(m)) return(unname(out[ids0]))
  m <- m[!is.na(m$ENTREZID), , drop = FALSE]
  if (nrow(m)) {
    agg <- tapply(m$ENTREZID, m$REFSEQ, function(z) paste(sort(unique(z)), collapse = " /// "))
    out[names(agg)] <- as.character(agg)
  }
  unname(out[ids0])
}
# 行标识 -> ENTREZ 的统一分派（GPL 注释表缺失时的兜底路线）。
# 同一个矩阵里可能混着多种标识（如 GPL32055 的 nCounter panel = NM_* RefSeq +
# ERCC_* spike-in），因此逐行标识判型，各走各的映射，未识别的留给下一步（NA）。
ids2entrez <- function(ids) {
  ids <- as.character(ids)
  if (!length(ids)) return(character(0))
  kind_of <- function(z) ifelse(grepl("^ENSG", z), "ENSG",
                         ifelse(grepl("^[0-9]+$", z), "ENTREZ",
                         ifelse(grepl("^N[MR]_", z), "REFSEQ", "SYMBOL")))
  u <- unique(ids); ku <- kind_of(u)
  lut <- setNames(rep(NA_character_, length(u)), u)
  if (any(ku == "ENSG"))   { m <- ensg2entrez(u[ku == "ENSG"]); lut[names(m)] <- m }
  if (any(ku == "ENTREZ"))   lut[u[ku == "ENTREZ"]] <- u[ku == "ENTREZ"]
  if (any(ku == "REFSEQ"))   lut[u[ku == "REFSEQ"]] <- refseq2entrez(u[ku == "REFSEQ"])
  if (any(ku == "SYMBOL")) { m <- symbol2entrez(u[ku == "SYMBOL"]); lut[u[ku == "SYMBOL"]] <- m }
  unname(lut[ids])
}

# 平台注释：优先 AnnoProbe::idmap(<GPL>, type="pipe")（审计指定的回退路线），
# 否则用项目根的 <GPL>.soft.gz（02_gpl_map.R 的路径）；两者都拿不到时，若矩阵行
# 本身就是基因标识（RNA-seq / nCounter），直接合成（ID_REF -> ENTREZ）。
build_gpl_map <- function(gpl, probe_ids) {
  dest <- file.path(ROOT, "data/processed/gpl", paste0(gpl, ".rds"))
  ids <- unique(as.character(probe_ids))
  existing <- if (file.exists(dest)) readRDS(dest) else NULL
  # 同一个 GPL 可能被不同队列用不同的行标识空间使用（如 GPL24676 同时被
  # GSE271517 的 ENSG 计数矩阵和 GSE325123 的基因名矩阵使用），因此这里做
  # 「按需扩充」：现有表能覆盖 >50% 的当前行标识就直接用，否则新建后并集合并。
  if (!is.null(existing) && length(ids) && mean(ids %in% rownames(existing)) > 0.5) {
    miss <- setdiff(ids, rownames(existing))
    if (!length(miss)) return(existing)
    message("  ", gpl, ": extending existing map with ", length(miss), " new probe ids")
    add <- data.frame(gene_id = ids2entrez(miss), row.names = miss, stringsAsFactors = FALSE)
    out <- rbind(existing, add); out <- out[!duplicated(rownames(out)), , drop = FALSE]
    saveRDS(out, dest)
    return(out)
  }
  merge_map <- function(new) {
    if (is.null(existing)) return(new)
    add <- new[setdiff(rownames(new), rownames(existing)), , drop = FALSE]
    out <- rbind(existing, add)
    out[!duplicated(rownames(out)), , drop = FALSE]
  }
  soft <- file.path(ROOT, paste0(gpl, ".soft.gz"))
  if (file.exists(soft)) {
    message("  ", gpl, ": building from local SOFT via the 02_gpl_map.R route")
    system2("Rscript", c("pipeline/R/02_gpl_map.R", gpl, ROOT))
    if (file.exists(dest)) return(readRDS(dest))
  }
  if (requireNamespace("AnnoProbe", quietly = TRUE)) {
    res <- tryCatch(suppressWarnings(AnnoProbe::idmap(gpl, type = "pipe")),
                    error = function(e) NULL)
    if (!is.null(res) && nrow(res)) {
      message("  ", gpl, ": AnnoProbe pipe map ", nrow(res), " rows, cols=",
              paste(colnames(res), collapse = ","))
      symcol <- intersect(c("symbol", "Symbol", "GENE_SYMBOL", "gene_symbol"), colnames(res))[1]
      prbcol <- colnames(res)[1]
      s2e <- symbol2entrez(res[[symcol]])
      gid <- unname(s2e[as.character(res[[symcol]])])
      out <- data.frame(gene_id = gid, row.names = as.character(res[[prbcol]]),
                        stringsAsFactors = FALSE)
      out <- out[!duplicated(rownames(out)), , drop = FALSE]
      out <- merge_map(out)
      saveRDS(out, dest)
      return(out)
    }
  }
  # 基因级矩阵：行名即基因标识（nCounter panel 的 ".IO360" 后缀只剥掉做 SYMBOL 查询，
  # 行名本身保留原样，避免重名）
  gid <- ids2entrez(sub("\\.IO360$", "", ids))
  out <- data.frame(gene_id = gid, row.names = ids, stringsAsFactors = FALSE)
  out <- out[!duplicated(rownames(out)), , drop = FALSE]
  message("  ", gpl, ": synthesized gene-level map from matrix row ids (",
          sum(!is.na(out$gene_id)), "/", nrow(out), " mapped)")
  out <- merge_map(out)
  saveRDS(out, dest)
  out
}

# ---- GeoMx ROI -> 患者聚合（GSE325123）-------------------------------------
# supplementary 表达表只含 included in_analysis == 1 的 170 个 ROI（列名 = aoi）；
# 审计 §3.6 指出 series matrix 的 272 列是 ROI 不是患者（114 个 patient 值，其中
# 3 个是 control）。此处按 patient 聚合：表达取该患者 ROI 的均值，临床/终点取该
# 患者 ROI 上的唯一非空值；control 患者因 included in_analysis == 0 已被剔除。
roi_patient_tables <- function(acc, pd) {
  stopifnot(acc == "GSE325123")
  # parseGSEMatrix 的特征列名是 "<key>:ch1"，这里统一取列
  getc <- function(nm) {
    v <- pd[[nm]]
    if (is.null(v)) v <- pd[[paste0(nm, ":ch1")]]
    if (is.null(v)) stop(acc, ": pheno column not found: ", nm)
    as.character(v)
  }
  f <- fetch_suppl(acc, "GSE325123_ITM_log2CPM.csv.gz")
  x <- read.csv(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(x[[1]])
  m <- as.data.frame(x[, -1, drop = FALSE], check.names = FALSE)
  pt <- getc("patient"); aoi <- getc("aoi")
  sel <- which(aoi %in% colnames(m))
  if (length(sel) != ncol(m))
    stop(acc, ": supplementary columns (", ncol(m), ") do not match aoi (",
         length(sel), " matched)")
  keep_pat <- unique(pt[sel])
  mm <- matrix(NA_real_, nrow = length(ids), ncol = length(keep_pat),
               dimnames = list(ids, keep_pat))
  for (p in keep_pat) {
    j <- sel[pt[sel] == p]
    mm[, p] <- rowMeans(as.matrix(m[, aoi[j], drop = FALSE]), na.rm = TRUE)
  }
  mmdf <- as.data.frame(mm, check.names = FALSE); rownames(mmdf) <- NULL
  expr <- data.frame(ID_REF = ids, mmdf, check.names = FALSE)
  pick <- function(col) {                 # 患者内取唯一非空值（ROI 间应一致）
    base <- getc(col)
    vapply(keep_pat, function(p) {
      v <- base[sel[pt[sel] == p]]
      v <- v[!is.na(v) & nzchar(v) & v != "NA"]
      if (!length(v)) NA_character_ else v[1]
    }, character(1))
  }
  nroi <- vapply(keep_pat, function(p) sum(pt[sel] == p), integer(1))
  pheno <- data.frame(
    patient                = keep_pat,
    overall_survival_years = pick("overall survival(years)"),
    vital_status           = pick("vital status"),
    distant_metastasis     = pick("distant metastasis_(1=yes,_0=no)"),
    lymph_node_status      = pick("lymph node_status_(1=yes,_0=no)"),
    acral_melanoma         = pick("acral melanoma_(1=yes,_0=no)"),
    tissue                 = pick("tissue"),
    n_roi                  = nroi,
    row.names = keep_pat, stringsAsFactors = FALSE)
  list(expr = expr, pheno = pheno, n_roi_in = length(sel))
}

# ---- 表达尺度归一（与 18_log2_transform.R 同规则；只作用于本脚本刚建的文件）----
# 18_log2_transform.R 的判定规则：min<0 或 (median>0 且 max/median<=50) 或
# (median<=0 且 max<=25) => 已是 log；否则 log2(x+1)。
# 这里**不能**直接跑 `18 --apply`：它的规则对既有队列 GSE91061（已 log2，max/median=52.8）
# 会误判成 needs_log2，重跑会把既有数据二次转化。因此本脚本只对本批次新建的
# expr 文件套用同一规则，原件备份到 pipeline/backup/expr_pre_log2/。
fix_scale <- function(acc) {
  f <- file.path(ROOT, "data/expr", paste0(acc, ".rds"))
  if (!file.exists(f)) return(NA_character_)
  x <- readRDS(f); ids <- x[[1]]
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "double"
  v <- m[is.finite(m)]
  verdict <- if (min(v) < 0) "already_log" else
             if (median(v) > 0 && max(v) / median(v) <= 50) "already_log" else
             if (median(v) <= 0 && max(v) <= 25) "already_log" else "needs_log2"
  if (verdict == "needs_log2") {
    bak <- file.path(ROOT, "pipeline/backup/expr_pre_log2", paste0(acc, ".rds"))
    dir.create(dirname(bak), showWarnings = FALSE, recursive = TRUE)
    if (!file.exists(bak)) file.copy(f, bak)
    m2 <- log2(m + 1)
    out <- data.frame(ID_REF = ids, m2, check.names = FALSE); names(out) <- names(x)
    saveRDS(out, f)
    message("  log2(x+1) applied: min=", round(min(m2, na.rm = TRUE), 3),
            " median=", round(median(m2, na.rm = TRUE), 3),
            " max=", round(max(m2, na.rm = TRUE), 3), " (原件备份 expr_pre_log2/)")
  }
  verdict
}

# ---- 主循环 ---------------------------------------------------------------
if (length(ONLY)) SPEC <- Filter(function(s) s$acc %in% ONLY, SPEC)
summ <- list()
for (s in SPEC) {
  acc <- s$acc; gpl <- s$gpl
  message("=== ", acc, " (", gpl, ", expr_mode=", s$mode, ")")
  raw <- fetch_matrix(acc)
  res <- parseGSEMatrix(raw, AnnotGPL = FALSE, getGPL = FALSE, parseCharacteristics = TRUE)
  eset <- res$eset; pd <- Biobase::pData(eset)
  n_expr_row <- nrow(Biobase::exprs(eset))
  message("  parsed: ", nrow(pd), " samples, matrix rows = ", n_expr_row,
          ", GPL(annotation) = ", Biobase::annotation(eset))

  expr_data <- if (s$mode == "matrix" && n_expr_row > 0) {
    as.data.frame(Biobase::exprs(eset)) %>% tibble::rownames_to_column("ID_REF")
  } else if (s$mode == "suppl") {
    expr_from_suppl(acc, pd)
  } else if (s$mode == "roi_patient") {
    rp <- roi_patient_tables(acc, pd)
    message("  ROI -> patient: ", rp$n_roi_in, " in-analysis ROIs -> ",
            nrow(rp$pheno), " patients")
    pd <- rp$pheno
    rp$expr
  } else if (s$mode == "matrix" && n_expr_row == 0) {
    stop(acc, ": empty series matrix but mode=matrix (expression must come from a supplement)")
  } else NULL

  if (!is.null(expr_data)) {
    # 与 90 相同的「整列全空」清理（GEO 里 data_row_count = 0 的占位样本）
    empty_cols <- names(expr_data)[-1][
      colSums(is.na(expr_data[, -1, drop = FALSE])) == nrow(expr_data)]
    if (length(empty_cols)) {
      message("  dropping ", length(empty_cols), " all-NA sample columns: ",
              paste(empty_cols, collapse = ","))
      expr_data <- expr_data[, setdiff(names(expr_data), empty_cols), drop = FALSE]
    }
    expr_data[-1] <- lapply(expr_data[-1], function(x) suppressWarnings(as.numeric(x)))
    saveRDS(expr_data, file.path(ROOT, "data/expr", paste0(acc, ".rds")))
  }
  pheno_data <- pd %>% dplyr::select(-dplyr::contains("characteristics"))
  saveRDS(pheno_data, file.path(ROOT, "data/pheno", paste0(acc, ".rds")))

  gm <- build_gpl_map(gpl, if (!is.null(expr_data)) expr_data$ID_REF else character(0))
  n_in <- if (!is.null(expr_data)) sum(as.character(expr_data$ID_REF) %in% rownames(gm)) else NA_integer_
  n_gene <- sum(!is.na(gm$gene_id))
  message("  GPL map: ", nrow(gm), " probes, ", n_gene, " with ENTREZ; matrix probes in map = ",
          n_in, if (!is.na(n_in)) sprintf(" (%.1f%% of %d)", 100 * n_in / max(1, nrow(expr_data)), nrow(expr_data)) else "")

  summ[[length(summ) + 1L]] <- data.frame(
    acc = acc, gpl = gpl, mode = s$mode, n_pheno = nrow(pd),
    expr_rows = if (!is.null(expr_data)) nrow(expr_data) else NA_integer_,
    expr_cols = if (!is.null(expr_data)) ncol(expr_data) - 1L else NA_integer_,
    matrix_rows = n_expr_row, gpl_rows = nrow(gm), gpl_mapped = n_gene,
    probes_in_map = n_in,
    map_rate = if (!is.null(expr_data) && nrow(expr_data) > 0)
      round(100 * n_in / nrow(expr_data), 2) else NA_real_,
    gene_map_rate = round(100 * n_gene / nrow(gm), 2),
    scale_verdict = fix_scale(acc), stringsAsFactors = FALSE)
}
out <- do.call(rbind, summ)
write.csv(out, file.path(ROOT, "pipeline/out/geo_expansion_build_expr_pheno.csv"), row.names = FALSE)
cat("\n=== 92 build summary ===\n"); print(out, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# run_93_build_geo_expansion_bespoke()  <-  verbatim pipeline/R/93_build_geo_expansion_bespoke.R
# ---------------------------------------------------------------------------
run_93_build_geo_expansion_bespoke <- function() {
# 93_build_geo_expansion_bespoke.R -------------------------------------------
# GEO 扩展批次（A 审计后）第 2 步：两个需要"额外一步"的队列，产物沿用既有
# pipeline/out/<ACC>_extra.rds 约定（行名 = GSM，对应 03_surv_table.R 的 spec$extra，
# 与 GSE37642-GPL97 / GSE14520-GPL3921 / GSE24080 用法一致）：
#
#   GSE31312 (Lymphoma, DLBCL, GPL570, n=498)
#     series matrix 里**没有**生存列；临床表是补充 PDF
#     GSE31312_Microarray_and_clinical_data_DLBCL_475_cases_PMID_22437443.pdf.gz
#     的 Table 2（列: Case # | Date death | PFScensor | OScensor | PFS | OS | ... |
#     GEO Depository #）。用 pdftotext -layout 抽出后按 `GEO Depository #`
#     （= !Sample_title）与 498 个样本连接。OScensor/PFScensor 为**删失**指示
#     （1=删失，0=事件）-> 03 spec 用 invert=TRUE。
#
#   GSE32918 (Lymphoma, DLBCL, GPL8432, n=249)
#     GEO design = 172 例患者，但矩阵有 249 张芯片（_Rep1.._Rep13 重复）。审计
#     §3.5 指出按芯片计数会把 172 例患者重复计成 249。这里为每张芯片登记
#     patient_id / array_rep / patient_dedup_keep（每位患者保留一张：无 _Rep
#     后缀者优先，否则编号最小者），03 spec 用 keep = quote(patient_dedup_keep)。
#
# 用法: Rscript pipeline/R/93_build_geo_expansion_bespoke.R [ROOT]
# ----------------------------------------------------------------------------
suppressMessages({library(dplyr)})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
SUPPL <- file.path(ROOT, "data/suppl")
OUT   <- file.path(ROOT, "pipeline/out")
summ <- list()

# ---------------------------------------------------------------- GSE31312 PDF
acc <- "GSE31312"
pd  <- readRDS(file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
gsm <- rownames(pd); titles <- trimws(as.character(pd$title))
pdfgz <- file.path(SUPPL, "GSE31312_Microarray_and_clinical_data_DLBCL_475_cases_PMID_22437443.pdf.gz")
stopifnot(file.exists(pdfgz))
tmp_pdf <- tempfile(fileext = ".pdf"); tmp_txt <- tempfile(fileext = ".txt")
system2("zcat", c(shQuote(pdfgz)), stdout = tmp_pdf)
system2("pdftotext", c("-layout", tmp_pdf, tmp_txt))
lines <- readLines(tmp_txt, warn = FALSE)
unlink(c(tmp_pdf, tmp_txt))

RE <- paste0("^\\s*(\\d+)\\s+(?:(\\d+/\\d+/\\d+)\\s+)?([01])\\s+([01])\\s+",
             "([0-9.]+)\\s+([0-9.]+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+",
             "(\\S+)\\s+(\\S+)\\s+(DLBCL GEP data \\S+)\\s*$")
m <- regmatches(lines, regexec(RE, lines, perl = TRUE))
m <- m[vapply(m, length, integer(1)) == 15L]
message(acc, ": parsed ", length(m), " clinical rows from the supplementary PDF")
if (length(m) != 475L) warning(acc, ": expected 475 rows, got ", length(m))
pdf_tab <- do.call(rbind, lapply(m, function(z) data.frame(
  pdf_case = as.integer(z[2]), pdf_date_death = z[3],
  pdf_pfs_censor = as.numeric(z[4]), pdf_os_censor = as.numeric(z[5]),
  pdf_pfs_months = as.numeric(z[6]), pdf_os_months = as.numeric(z[7]),
  pdf_cd10 = as.numeric(z[8]), pdf_bcl6 = as.numeric(z[9]), pdf_foxp1 = as.numeric(z[10]),
  pdf_gcet = as.numeric(z[11]), pdf_mum1 = as.numeric(z[12]),
  pdf_markers = z[13], pdf_gep = z[14], title = z[15], stringsAsFactors = FALSE)))
stopifnot(!anyDuplicated(pdf_tab$title), !anyDuplicated(pdf_tab$pdf_case))
ex <- pdf_tab[match(titles, pdf_tab$title), , drop = FALSE]
ex$title <- NULL; rownames(ex) <- gsm
n_join <- sum(!is.na(ex$pdf_case))
message(acc, ": joined ", n_join, "/", length(gsm), " samples to the PDF table",
        " (OS non-missing ", sum(!is.na(ex$pdf_os_months)),
        ", deaths ", sum(ex$pdf_os_censor == 0, na.rm = TRUE), ")")
saveRDS(ex, file.path(OUT, paste0(acc, "_extra.rds")))
summ[[length(summ) + 1L]] <- data.frame(
  acc = acc, step = "pdf_clinical_join", n_rows = nrow(ex), n_joined = n_join,
  detail = sprintf("PDF rows=%d; join by GEO Depository # == !Sample_title; OS deaths=%d of %d joined",
                   nrow(pdf_tab), sum(ex$pdf_os_censor == 0, na.rm = TRUE), n_join),
  stringsAsFactors = FALSE)

# ------------------------------------------------------- GSE32918 患者级去重
acc <- "GSE32918"
pd  <- readRDS(file.path(ROOT, "data/pheno", paste0(acc, ".rds")))
titles <- trimws(as.character(pd$title))
base   <- sub("_Rep[0-9]+$", "", titles)
repn   <- ifelse(grepl("_Rep[0-9]+$", titles),
                 suppressWarnings(as.integer(sub(".*_Rep([0-9]+)$", "\\1", titles))), 0L)
ord <- order(base, repn)
keep <- logical(length(base)); seen <- character(0)
for (i in ord) if (!base[i] %in% seen) { seen <- c(seen, base[i]); keep[i] <- TRUE }
ex <- data.frame(patient_id = base, array_rep = repn, patient_dedup_keep = keep,
                 row.names = rownames(pd), stringsAsFactors = FALSE)
# pheno 的特征列带 ":ch1" 后缀（01/92 的 parseGSEMatrix 约定）
st <- as.character(pd[["follow-up status:ch1"]])
yr <- suppressWarnings(as.numeric(pd[["follow-up years:ch1"]]))
k <- keep & !is.na(yr) & st %in% c("Dead", "Alive")
message(acc, ": ", nrow(ex), " arrays -> ", sum(keep), " patients; usable ",
        sum(k), ", deaths ", sum(st[k] == "Dead"))
saveRDS(ex, file.path(OUT, paste0(acc, "_extra.rds")))
summ[[length(summ) + 1L]] <- data.frame(
  acc = acc, step = "patient_dedup", n_rows = nrow(ex), n_joined = sum(keep),
  detail = sprintf("%d arrays -> %d patients (%d with >1 array); usable %d / deaths %d",
                   nrow(ex), sum(keep), sum(table(base) > 1), sum(k), sum(st[k] == "Dead")),
  stringsAsFactors = FALSE)

res <- do.call(rbind, summ)
write.csv(res, file.path(OUT, "geo_expansion_bespoke_log.csv"), row.names = FALSE)
cat("\n=== 93 bespoke build summary ===\n"); print(res, row.names = FALSE)
}
