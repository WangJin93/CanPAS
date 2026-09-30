# ===========================================================================
# CanPAS curation pipeline -- 02_platform_maps
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. 02_gpl_map.R
#   2. 12_gpl_map_symbol.R
#   3. 38_geo_expansion_annotation_check.R
#   4. 102_extend_gpl4133_agilent_name.R
#   5. 103_build_embl_gpl_maps.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "02_platform_maps.R", package = "CanPAS"))
#   run_02_gpl_map()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================


# ---------------------------------------------------------------------------
# run_02_gpl_map()  <-  verbatim pipeline/R/02_gpl_map.R
# ---------------------------------------------------------------------------
run_02_gpl_map <- function() {
# 02_gpl_map.R
# -----------------------------------------------------------------------------
# CanPAS 数据处理 第 2 步：生成平台注释“探针 -> ENTREZ gene_id”表
#   输出 data/processed/gpl/<GPL>.rds
#     行名 = 探针 ID；唯一列 gene_id = ENTREZ（多基因探针以 " /// " 连接，如既有文件）
# 输入：项目根目录下的 <GPL>.soft.gz（与既有 GPL*.soft.gz 相同：纯文本 SOFT 或 gzip 均可）。
#
# 优先使用 SOFT platform_table 中的 Entrez 列
#   （ENTREZ_GENE_ID / Entrez_Gene_ID / GENE_ID / Gene ID）
# 若平台表没有 Entrez 列（例如 Affymetrix Gene-ST：GPL6244），则解析 mrna_assignment
#   中的 "gene:ENSGxxxx" 并用 REF_ID/www/biomart_table_Hs.Rds 转换为 ENTREZ。
#
# 用法: Rscript 02_gpl_map.R <GPL-id> [<ProjectRoot>]
# 例  : Rscript 02_gpl_map.R GPL6244 ~/data/Project/CanPAS
# -----------------------------------------------------------------------------

suppressMessages({library(readr); library(dplyr); library(stringr)})

args  <- commandArgs(trailingOnly=TRUE)
GPLid <- args[1]
ROOT  <- if (length(args) >= 2) args[2] else "~/data/Project/CanPAS"
ROOT  <- path.expand(ROOT)
setwd(ROOT)

soft <- file.path(ROOT, paste0(GPLid, ".soft.gz"))
if (!file.exists(soft)) stop("SOFT file not found: ", soft, " (put it in the project root first)")

is_gz <- function(f) {
  con <- file(f, "rb"); m <- readBin(con, "raw", n=2); close(con)
  identical(as.integer(m), as.integer(c(0x1f, 0x8b)))
}

read_soft_lines <- function(f) {
  if (is_gz(f)) readLines(gzfile(f), warn=FALSE, encoding="latin1")
  else readLines(f, warn=FALSE, encoding="latin1")
}

# 只读取 platform_table 区域（family soft 里平台表之后还有大量 SAMPLE 记录，不整读）
read_platform_table <- function(f) {
  con <- if (is_gz(f)) gzfile(f, "rt", encoding="latin1") else file(f, "rt", encoding="latin1")
  started <- FALSE
  parts <- list()
  repeat {
    chunk <- readLines(con, n=10000, warn=FALSE)
    if (!length(chunk)) break
    if (!started) {
      i <- grep("^!platform_table_begin", chunk)
      if (!length(i)) next
      started <- TRUE
      chunk <- chunk[(i[1]+1):length(chunk)]
    }
    j <- grep("^!platform_table_end", chunk)
    if (length(j)) {
      chunk <- chunk[1:(j[1]-1)]
      parts[[length(parts)+1]] <- sub("\r$", "", chunk)
      close(con); return(unlist(parts, use.names=FALSE))
    }
    parts[[length(parts)+1]] <- sub("\r$", "", chunk)
  }
  close(con)
  unlist(parts, use.names=FALSE)
}

tab_lines <- read_platform_table(soft)
if (!length(tab_lines)) stop("no platform_table in ", soft)
tab_lines <- tab_lines[nzchar(trimws(tab_lines))]

hdr <- strsplit(tab_lines[1], "\t")[[1]]
hdr <- sub("\r$", "", hdr)
body <- read.delim(text=paste(tab_lines[-1], collapse="\n"), sep="\t",
                   header=FALSE, quote="", stringsAsFactors=FALSE,
                   col.names=hdr, check.names=FALSE)
message(GPLid, ": table ", nrow(body), " x ", ncol(body),
        "  cols: ", paste(hdr, collapse=", "))
body <- body[!is.na(body[[1]]) & body[[1]] != "" & body[[1]] != "---", , drop=FALSE]

probe_col <- hdr[1]
entrez_candidates <- c("ENTREZ_GENE_ID", "Entrez_Gene_ID", "GENE_ID",
                       "Gene ID", "GENE", "ENTREZID", "gene_id",
                       "Entrez.Gene", "Entrez Gene ID", "Entrez.ID")
use_col <- intersect(hdr, entrez_candidates)[1]

gene_vec <- if (!is.na(use_col)) {
  message("using column: ", use_col)
  v <- as.character(body[[use_col]])
  v <- ifelse(v %in% c("---", "", "NA"), NA_character_, v)
  v
} else {
  bm <- readRDS(file.path(ROOT, "REF_ID/www/biomart_table_Hs.Rds"))
  bm_map <- stats::setNames(as.character(bm$Entrez), bm$Ensembl_gene)
  ensg_cols <- intersect(hdr, c("OligoSet_ensemblGene", "Ensembl", "Ensembl_gene",
                                "ENSEMBL", "ENSEMBL_ID", "ensembl"))
  if (length(ensg_cols)) {
    message("using ENSG column: ", ensg_cols[1], " -> biomart Entrez")
    vapply(body[[ensg_cols[1]]], function(s) {
      ens <- unique(regmatches(s, gregexpr("ENSG[0-9]+", s))[[1]])
      ent <- unique(stats::na.omit(bm_map[ens]))
      if (!length(ent)) NA_character_ else paste(sort(ent), collapse=" /// ")
    }, character(1), USE.NAMES=FALSE)
  } else {
    message("no Entrez/ENSG column -> parse mrna_assignment gene:ENSG..., map via biomart_table_Hs")
    ensg_of <- function(s) {
      if (is.na(s) || !nzchar(s)) return(character(0))
      m <- regmatches(s, gregexpr("gene:ENSG[0-9]+", s))[[1]]
      unique(sub("gene:", "", m))
    }
    vapply(body[[which(hdr=="mrna_assignment")[1]]], function(s) {
      ens <- ensg_of(s)
      ent <- unique(stats::na.omit(bm_map[ens]))
      if (!length(ent)) NA_character_ else paste(sort(ent), collapse=" /// ")
    }, character(1), USE.NAMES=FALSE)
  }
}

out <- data.frame(gene_id = gene_vec, stringsAsFactors=FALSE, row.names = body[[1]])
out <- out[!is.na(rownames(out)) & rownames(out) != "" & rownames(out) != "---", , drop=FALSE]
out <- out[!duplicated(rownames(out)), , drop=FALSE]
# 与既有文件一致：把空字符串统一为 ""（不上传），NA 保留（gpl_app 上传时 na.omit 会去掉）
saveRDS(out, file=file.path(ROOT, "data/processed/gpl", paste0(GPLid, ".rds")))

n_na   <- sum(is.na(out$gene_id) | out$gene_id=="")
n_multi <- sum(grepl(" /// ", out$gene_id), na.rm=TRUE)
message("saved data/processed/gpl/", GPLid, ".rds : ", nrow(out), " probes, ",
        n_na, " unmapped, ", n_multi, " multi-gene probes")
}

# ---------------------------------------------------------------------------
# run_12_gpl_map_symbol()  <-  verbatim pipeline/R/12_gpl_map_symbol.R
# ---------------------------------------------------------------------------
run_12_gpl_map_symbol <- function() {
# 12_gpl_map_symbol.R --------------------------------------------------------
# 为「平台表没有 Entrez 列、只有基因符号列」的 GEO 平台生成探针 -> ENTREZ 表。
# 与 02_gpl_map.R 输出格式一致(行名 = 探针 ID,唯一列 gene_id = ENTREZ,多基因以 " /// " 连接),
# 区别只在于映射来源:此处用 CanPAS 包内置 ID_map(Symbol -> Entrez)而非平台表的 Entrez 列。
#
# 输入: 项目根目录下的 <GPL>.soft.gz / <GPL>.annot.gz(需含 !platform_table_begin 平台表)
# 输出: data/processed/gpl/<GPL>.rds
# 用法: Rscript pipeline/R/12_gpl_map_symbol.R <GPL-id> [<ProjectRoot>]
# 例  : Rscript pipeline/R/12_gpl_map_symbol.R GPL10295 ~/data/Project/CanPAS
# ----------------------------------------------------------------------------
suppressMessages({library(readr)})

args  <- commandArgs(trailingOnly = TRUE)
GPLid <- args[1]
ROOT  <- path.expand(if (length(args) >= 2) args[2] else "~/data/Project/CanPAS")

is_gz <- function(f) {
  con <- file(f, "rb"); m <- readBin(con, "raw", n = 2); close(con)
  identical(as.integer(m), as.integer(c(0x1f, 0x8b)))
}
soft <- file.path(ROOT, paste0(GPLid, ".soft.gz"))
if (!file.exists(soft)) soft <- file.path(ROOT, paste0(GPLid, ".annot.gz"))
if (!file.exists(soft)) stop("platform file not found for ", GPLid)

read_platform_table <- function(f) {
  con <- if (is_gz(f)) gzfile(f, "rt", encoding = "latin1") else file(f, "rt", encoding = "latin1")
  on.exit(close(con), add = TRUE)
  started <- FALSE; parts <- list()
  repeat {
    chunk <- readLines(con, n = 20000, warn = FALSE)
    if (!length(chunk)) break
    if (!started) {
      i <- grep("^!platform_table_begin", chunk)
      if (!length(i)) next
      started <- TRUE
      chunk <- chunk[(i[1] + 1):length(chunk)]
    }
    j <- grep("^!platform_table_end", chunk)
    if (length(j)) { parts[[length(parts) + 1L]] <- sub("\r$", "", chunk[1:(j[1] - 1)]); break }
    parts[[length(parts) + 1L]] <- sub("\r$", "", chunk)
  }
  unlist(parts, use.names = FALSE)
}

tl <- read_platform_table(soft)
if (!length(tl)) stop("no platform_table in ", soft)
tl <- tl[nzchar(trimws(tl))]
hdr <- strsplit(tl[1], "\t")[[1]]
tab <- read.delim(text = paste(tl[-1], collapse = "\n"), sep = "\t", header = FALSE,
                  quote = "", col.names = hdr, check.names = FALSE,
                  stringsAsFactors = FALSE)
message(GPLid, ": platform table ", nrow(tab), " x ", ncol(tab), " cols: ", paste(hdr, collapse = ", "))

sym_col <- intersect(c("Symbol", "GENE_SYMBOL", "Gene symbol", "ILMN_Gene", "Gene", "GENE"),
                     hdr)[1]
if (is.na(sym_col))
  stop("no gene-symbol column in platform table of ", GPLid, "; use 02_gpl_map.R instead")

ID_map <- get(load(file.path(ROOT, "CanPAS/data/ID_map.rda")))
smap <- ID_map$gene_id
names(smap) <- ID_map$Symbol
smap <- smap[!is.na(smap) & nzchar(names(smap))]
smap <- smap[!duplicated(names(smap))]

v <- as.character(tab[[sym_col]])
gene_vec <- vapply(v, function(s) {
  if (is.na(s) || !nzchar(s)) return(NA_character_)
  syms <- unique(trimws(unlist(strsplit(s, "///|//|;|,"))))
  syms <- syms[nzchar(syms)]
  ent <- unique(unname(smap[syms]))
  ent <- ent[!is.na(ent)]
  if (!length(ent)) NA_character_ else paste(ent, collapse = " /// ")
}, character(1), USE.NAMES = FALSE)

out <- data.frame(gene_id = gene_vec, stringsAsFactors = FALSE, row.names = tab[[1]])
out <- out[!is.na(rownames(out)) & rownames(out) != "" & rownames(out) != "---", , drop = FALSE]
out <- out[!duplicated(rownames(out)), , drop = FALSE]
saveRDS(out, file = file.path(ROOT, "data/processed/gpl", paste0(GPLid, ".rds")))
message("saved data/processed/gpl/", GPLid, ".rds : ", nrow(out), " probes, ",
        sum(is.na(out$gene_id)), " unmapped, ",
        sum(grepl(" /// ", out$gene_id), na.rm = TRUE), " multi-gene probes ",
        "(mapped from column '", sym_col, "')")
}

# ---------------------------------------------------------------------------
# run_38_geo_expansion_annotation_check()  <-  verbatim pipeline/R/38_geo_expansion_annotation_check.R
# ---------------------------------------------------------------------------
run_38_geo_expansion_annotation_check <- function() {
# 38_geo_expansion_annotation_check.R ----------------------------------------
# 目的
#   对 pipeline/out/geo_expansion_candidates.csv 里通过门槛的候选所用平台，判断
#   "探针 -> 基因"注释是否可得，并把结论写回
#   pipeline/out/geo_expansion_annotation_check.csv（由 37 脚本读成 annotation_status 列）。
#
# 三级判定
#   1) local_rds                      : 项目里已有 data/gpl/<GPL>.rds
#   2) annoprobe_pipe                 : AnnoProbe::idmap(<GPL>, type="pipe") 能返回
#   3) gene_level_no_probe_map_needed : RNA-seq / NanoString 这类平台的矩阵行本身就是基因，
#                                       不需要探针映射（GPL 表也可以没有注释）
#   4) geo_annot_symbol:<cols>        : AnnoProbe 没有该平台，但 GEO 自带 GPLxxx.annot.gz，
#                                       其表头含基因列（GF 用 HTTP Range 只取前 200KB 再部分解压，
#                                       不下载完整注释表）
#   5) geo_annot_no_symbol:<cols>     : GEO annot 表头没有任何基因列（如纯序列平台）
#   6) idmap_failed:<msg>             : 以上都不成立
#
# 用法: Rscript pipeline/R/38_geo_expansion_annotation_check.R [ROOT]
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
OUT  <- cpas_out_root()
CAND <- file.path(OUT, "geo_expansion_candidates.csv")
RES  <- file.path(OUT, "geo_expansion_annotation_check.csv")
if (!file.exists(CAND)) stop("run 37_geo_expansion_screen.R first: ", CAND)
cd <- read.csv(CAND, check.names = FALSE, stringsAsFactors = FALSE)
pass <- cd[cd$passes_gate %in% c(TRUE, "TRUE"), ]
gpls <- sort(unique(pass$platform))
cls  <- tapply(pass$platform_class, pass$platform, function(z) z[1])
cat("platforms used by passing candidates:", length(gpls), "\n")

# ---- GEO annot 表头（HTTP Range + 部分 gzip 解压，只读前 200KB） --------------
GEO_ANNOT_COLS <- c("gene symbol", "genesymbol", "gene_symbol", "symbol", "orf",
                    "entrez_gene_id", "entrez gene id", "gene id", "gene_id",
                    "gene assignment", "gene name", "gb_acc", "genbank accession",
                    "refseq", "ensembl", "unigene", "gene")
geo_table_header <- function(gpl, sizes = c(200000L, 1000000L, 4000000L)) {
  num <- as.integer(sub("^GPL", "", gpl))
  dir <- sprintf("GPL%dnnn", num %/% 1000)
  base <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s", dir, gpl)
  urls <- c(sprintf("%s/annot/%s.annot.gz", base, gpl),
            sprintf("%s/soft/%s_family.soft.gz", base, gpl))
  for (url in urls) for (nb in sizes) {
    tmp <- tempfile(fileext = ".gz")
    ok <- tryCatch({ system2("curl", c("-s", "-f", "-r", sprintf("0-%d", nb - 1L), "-o", tmp, url),
                              stdout = FALSE, stderr = FALSE)
                     file.exists(tmp) && file.size(tmp) > 1000 }, error = function(e) FALSE)
    if (!ok) { unlink(tmp); next }
    raw <- readBin(tmp, "raw", n = file.size(tmp)); unlink(tmp)
    con <- tryCatch(gzcon(rawConnection(raw, "rb")), error = function(e) NULL)
    if (is.null(con)) next
    txt <- tryCatch(readLines(con, n = 400000L, warn = FALSE), error = function(e) character(0))
    try(close(con), silent = TRUE)
    i <- grep("^!platform_table_begin", txt)
    if (length(i) && i[1] < length(txt)) {
      hdr <- txt[(i[1] + 1L):min(i[1] + 3L, length(txt))]
      hdr <- hdr[nzchar(hdr)]
      if (length(hdr)) {
        cols <- trimws(strsplit(hdr[1], "\t")[[1]])
        if (length(cols) < 2) cols <- trimws(strsplit(hdr[1], ",")[[1]])
        return(list(cols = cols, src = basename(url)))
      }
    }
    if (!grepl("\\.gz$", url)) break
  }
  NULL
}

rows <- list()
for (g in gpls) {
  cl <- if (!is.na(cls[g])) as.character(cls[g]) else "unknown"
  rds <- file.path(ROOT, "data/gpl", paste0(g, ".rds"))
  if (file.exists(rds)) {
    rows[[length(rows) + 1L]] <- data.frame(gpl = g, platform_class = cl, status = "local_rds",
      detail = sprintf("data/gpl/%s.rds present (%d bytes)", g, file.info(rds)$size),
      checked_at = format(Sys.time()), stringsAsFactors = FALSE); next
  }
  if (cl %in% c("rna_seq", "ncounter")) {
    rows[[length(rows) + 1L]] <- data.frame(gpl = g, platform_class = cl,
      status = "gene_level_no_probe_map_needed",
      detail = "matrix rows are gene-level identifiers; no probe->gene map required",
      checked_at = format(Sys.time()), stringsAsFactors = FALSE); next
  }
  if (requireNamespace("AnnoProbe", quietly = TRUE)) {
    st <- "annoprobe_pipe"; dt <- ""
    ok <- tryCatch({ x <- suppressWarnings(AnnoProbe::idmap(g, type = "pipe"))
                     if (is.null(x) || !nrow(x)) stop("idmap returned 0 rows")
                     dt <- sprintf("%d probe->gene rows; cols=%s", nrow(x), paste(colnames(x), collapse = ",")); TRUE },
                   error = function(e) { st <<- paste0("idmap_failed:", conditionMessage(e)); FALSE })
    if (ok) {
      rows[[length(rows) + 1L]] <- data.frame(gpl = g, platform_class = cl, status = st, detail = dt,
        checked_at = format(Sys.time()), stringsAsFactors = FALSE); next
    }
  } else st <- "idmap_failed:AnnoProbe_not_installed"
  gh <- geo_table_header(g)
  cols <- if (is.null(gh)) NULL else gh$cols
  if (!is.null(cols)) {
    low <- tolower(cols)
    hit <- intersect(low, GEO_ANNOT_COLS)
    rows[[length(rows) + 1L]] <- data.frame(gpl = g, platform_class = cl,
      status = if (length(hit)) "geo_annot_symbol" else "geo_annot_no_symbol",
      detail = sprintf("GEO %s columns: %s", gh$src, paste(cols, collapse = ",")),
      checked_at = format(Sys.time()), stringsAsFactors = FALSE); next
  }
  rows[[length(rows) + 1L]] <- data.frame(gpl = g, platform_class = cl, status = st,
    detail = "no AnnoProbe pipe map and no GEO platform-table header could be read",
    checked_at = format(Sys.time()), stringsAsFactors = FALSE)
  cat(sprintf("  %s -> %s\n", g, rows[[length(rows)]]$status))
}
out <- do.call(rbind, rows)
write.csv(out, RES, row.names = FALSE)
cat("wrote", RES, "\n"); print(table(out$status))
}

# ---------------------------------------------------------------------------
# run_102_extend_gpl4133_agilent_name()  <-  verbatim pipeline/R/102_extend_gpl4133_agilent_name.R
# ---------------------------------------------------------------------------
run_102_extend_gpl4133_agilent_name <- function() {
# 102_extend_gpl4133_agilent_name.R ------------------------------------------
# GPL4133（Agilent-014850 Whole Human Genome 4x44K G4112F）的既有映射表
# data/processed/gpl/GPL4133.rds 是以平台表的 **ID 列（1..45220 数字行号）** 为行名建的，
# 而 ArrayExpress E-MTAB-1980 的加工矩阵用的是 **NAME 列（A_23_Pxxxxxx）** 作探针名，
# 两者 0% 交集 -> 无法经 ID_map → GPL4133 → expr 路径检索。
#
# 本脚本**纯追加**：保留原有数字行名不动，把 NAME（A_23_P…）→ ENTREZ 作为新行加入，
# 绝不更新/删除既有行（与 99_extend_gpl_db.R 的“只追加”口径一致）。
# 维护后 data/processed/gpl/GPL4133.rds 同时支持两种探针命名。
#
# 用法: Rscript pipeline/R/102_extend_gpl4133_agilent_name.R [ROOT] [--write]
# ----------------------------------------------------------------------------
suppressMessages({ library(readr) })
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else "/home/Jingle/data/Project/CPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
WRITE <- "--write" %in% args

gpl <- "GPL4133"
dest <- file.path("data/processed/gpl", paste0(gpl, ".rds"))
soft <- file.path(ROOT, paste0(gpl, ".soft.gz"))
if (!file.exists(soft)) stop("platform file not found: ", soft)
existing <- readRDS(dest)
message(gpl, ": existing map ", nrow(existing), " rows; rowname sample: ",
        paste(head(rownames(existing), 3), collapse = ", "))

# 只读 platform_table 区域
con <- file(soft, "rt", encoding = "latin1")
on.exit(close(con), add = TRUE)
started <- FALSE; parts <- list()
repeat {
  chunk <- readLines(con, n = 20000, warn = FALSE)
  if (!length(chunk)) break
  if (!started) {
    i <- grep("^!platform_table_begin", chunk)
    if (!length(i)) next
    started <- TRUE
    chunk <- chunk[(i[1] + 1):length(chunk)]
  }
  j <- grep("^!platform_table_end", chunk)
  if (length(j)) { parts[[length(parts) + 1L]] <- sub("\r$", "", chunk[1:(j[1] - 1)]); break }
  parts[[length(parts) + 1L]] <- sub("\r$", "", chunk)
}
tl <- unlist(parts, use.names = FALSE); tl <- tl[nzchar(trimws(tl))]
hdr <- strsplit(tl[1], "\t")[[1]]
tab <- read.delim(text = paste(tl[-1], collapse = "\n"), sep = "\t", header = FALSE,
                  quote = "", col.names = hdr, check.names = FALSE, stringsAsFactors = FALSE)
message("  platform table ", nrow(tab), " x ", ncol(tab), " cols; using NAME + GENE_SYMBOL")
nm <- as.character(tab[["NAME"]]); sym <- as.character(tab[["GENE_SYMBOL"]])
keep <- !is.na(nm) & nzchar(nm) & nm != "---"
nm <- nm[keep]; sym <- sym[keep]
# 同名探针保留第一次出现（保持 1:1 行名）
dup <- duplicated(nm)
message("  NAME rows=", length(nm), " duplicated NAME=", sum(dup))
nm <- nm[!dup]; sym <- sym[!dup]

miss <- setdiff(nm, rownames(existing))
message("  NAME ids not yet in map: ", length(miss), " / ", length(nm))
if (!length(miss)) { message("nothing to add"); quit(save = "no") }
sym_miss <- sym[match(miss, nm)]

# SYMBOL -> ENTREZ（与 96 的 extend_gpl 同口径，用 org.Hs.eg.db）
db <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
gid <- setNames(rep(NA_character_, length(miss)), miss)
keys <- unique(sym_miss[!is.na(sym_miss) & nzchar(sym_miss)])
m <- suppressWarnings(AnnotationDbi::select(db, keys = keys, columns = "ENTREZID", keytype = "SYMBOL"))
m <- m[!is.na(m$ENTREZID), , drop = FALSE]
if (nrow(m)) {
  agg <- tapply(m$ENTREZID, m$SYMBOL, function(z) paste(sort(unique(z)), collapse = " /// "))
  for (i in seq_along(miss)) {
    s <- sym_miss[i]
    if (!is.na(s) && s %in% names(agg)) gid[i] <- as.character(agg[[s]])
  }
}
add_gid <- unname(gid)
rn <- c(rownames(existing), names(gid))
out <- data.frame(gene_id = c(as.character(existing$gene_id), add_gid), stringsAsFactors = FALSE)
attr(out, "row.names") <- rn
out <- out[!duplicated(rownames(out)), , drop = FALSE]
# 防御：调用点若传错 ids，这里必须报错而不是静默写坏表（同 96）
if (!all(miss %in% rownames(out))) stop("extend: not all requested NAME ids became row names")
if (!all(rownames(existing) %in% rownames(out))) stop("extend: existing rows were lost")
message(sprintf("  would append %d NAME rows (%d with ENTREZ) -> %d rows total",
                length(gid), sum(!is.na(add_gid)), nrow(out)))
if (!WRITE) { cat("\n[dry-run] 加 --write 才写回 ", dest, "\n"); quit(save = "no") }
saveRDS(out, dest)
message("saved ", dest, " : ", nrow(out), " rows, ",
        sum(is.na(out$gene_id)), " unmapped, ",
        sum(grepl(" /// ", out$gene_id), na.rm = TRUE), " multi-gene")
}

# ---------------------------------------------------------------------------
# run_103_build_embl_gpl_maps()  <-  verbatim pipeline/R/103_build_embl_gpl_maps.R
# ---------------------------------------------------------------------------
run_103_build_embl_gpl_maps <- function() {
# 103_build_embl_gpl_maps.R ---------------------------------------------------
# Phase-0 of the EMBL-EBI build: build probe -> ENTREZ maps for the two platforms
# the audit proved AnnoProbe::idmap() cannot serve (GPL17585 HTA-2.0,
# GPL16686 HuGene-2.0-st), straight from the GEO platform SOFT files.
#
# Why not 02_gpl_map.R:
#   GPL17585's platform_table has NO Entrez column and NO mrna_assignment column
#   (cols: ID, probeset_id, SPOT_ID, transcript_cluster_ids, chromosome,
#    RANGE_*, transcript_cluster_*, probe_count, gene_symbols, descriptions),
#   so 02_gpl_map.R's Entrez/ENSG/mrna_assignment routes all miss.  The gene
#   column that IS present is `gene_symbols` (Affymetrix "A /// B" convention),
#   so this script maps Symbol -> Entrez via the project's own
#   REF_ID/www/biomart_table_Hs.Rds (the same table 02_gpl_map.R uses).
#
# Outputs (both conventions, so nothing downstream has to guess):
#   data/gpl/<GPL>.rds              data.frame(ID_REF, ENTREZ_GENE_ID)  <- neighbouring maps
#   data/processed/gpl/<GPL>.rds    data.frame(gene_id, row.names=probe) <- 92_build_gpl_map() reads this
#
# Usage: Rscript pipeline/R/103_build_embl_gpl_maps.R [<GPL> ...] [ROOT]
# ----------------------------------------------------------------------------
suppressMessages({ library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- c("GPL17585", "GPL16686")
known <- c("GPL17585", "GPL16686")
gpls  <- intersect(args, known)
rest  <- setdiff(args, known)
ROOT  <- if (length(rest)) rest[length(rest)] else "/home/Jingle/data/Project/CPAS"
ROOT  <- path.expand(ROOT); setwd(ROOT)
GENELEVEL_MODE <- any(grepl("^--genelevel", commandArgs(trailingOnly = TRUE)))
if (!length(gpls) && !GENELEVEL_MODE) stop("no known GPL id given (known: ", paste(known, collapse=", "), ")")

is_gz <- function(f) {
  con <- file(f, "rb"); m <- readBin(con, "raw", n = 2); close(con)
  identical(as.integer(m), as.integer(c(0x1f, 0x8b)))
}

# Stream only the platform_table region.  Tolerates a *truncated* .gz (the family
# SOFT puts the platform table before the SAMPLE records, so a bounded-prefix
# download is enough) - gzfile() warns on the truncated tail, which we silence.
read_platform_table <- function(f) {
  con <- if (is_gz(f)) gzfile(f, "rt", encoding = "latin1") else file(f, "rt", encoding = "latin1")
  on.exit(close(con), add = TRUE)
  started <- FALSE; parts <- list()
  repeat {
    chunk <- tryCatch(withCallingHandlers(
      readLines(con, n = 20000, warn = FALSE),
      warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) character(0))
    if (!length(chunk)) break
    if (!started) {
      i <- grep("^!platform_table_begin", chunk)
      if (!length(i)) next
      started <- TRUE; chunk <- chunk[(i[1] + 1):length(chunk)]
    }
    j <- grep("^!platform_table_end", chunk)
    if (length(j)) {
      parts[[length(parts) + 1]] <- sub("\r$", "", chunk[1:(j[1] - 1)])
      return(unlist(parts, use.names = FALSE))
    }
    parts[[length(parts) + 1]] <- sub("\r$", "", chunk)
  }
  unlist(parts, use.names = FALSE)
}

# Symbol -> Entrez via the project biomart table (a symbol may map to >1 Entrez;
# join with " /// " exactly like the existing maps).
bm <- readRDS(file.path(ROOT, "REF_ID/www/biomart_table_Hs.Rds"))
bm <- bm[!is.na(bm$Symbol) & nzchar(bm$Symbol) & !is.na(bm$Entrez), c("Symbol", "Entrez")]
sym2ent <- tapply(as.character(bm$Entrez), as.character(bm$Symbol), function(x)
  paste(sort(unique(x)), collapse = " /// "))
message("biomart: ", length(sym2ent), " symbols -> Entrez")

ENTREZ_CAND <- c("ENTREZ_GENE_ID", "Entrez_Gene_ID", "GENE_ID", "Gene ID", "GENE",
                 "ENTREZID", "gene_id", "Entrez.Gene", "Entrez Gene ID", "Entrez.ID")
SYMBOL_CAND <- c("gene_symbols", "Gene symbol", "GENE_SYMBOL", "GENE_SYMBOLS",
                 "Symbol", "symbol", "ILMN_Gene", "GeneSymbol")

res <- list()
for (GPLid in gpls) {  # nolint
  soft <- file.path(ROOT, paste0(GPLid, ".soft.gz"))
  if (!file.exists(soft)) { message("!! ", GPLid, ": SOFT not found: ", soft); next }
  cat("\n====", GPLid, "====\n")
  tab <- read_platform_table(soft)
  if (!length(tab)) { message("!! ", GPLid, ": no platform_table"); next }
  tab <- tab[nzchar(trimws(tab))]
  hdr <- strsplit(tab[1], "\t")[[1]]
  cat("platform_table rows (data):", length(tab) - 1, " cols:", length(hdr), "\n")
  cat("header:", paste(hdr, collapse = " | "), "\n")
  body <- read.delim(text = paste(tab[-1], collapse = "\n"), sep = "\t", header = FALSE,
                     quote = "", stringsAsFactors = FALSE, col.names = hdr,
                     check.names = FALSE)
  probe <- as.character(body[[1]])
  keep  <- !is.na(probe) & nzchar(probe) & probe != "---"
  body <- body[keep, , drop = FALSE]

  use_ent <- intersect(hdr, ENTREZ_CAND)[1]
  use_sym <- intersect(hdr, SYMBOL_CAND)[1]
  cat("entrez column:", if (is.na(use_ent)) "<none>" else use_ent,
      " | symbol column:", if (is.na(use_sym)) "<none>" else use_sym, "\n")

  gene_vec <- if (!is.na(use_ent)) {
    v <- as.character(body[[use_ent]]); ifelse(v %in% c("---", "", "NA"), NA_character_, v)
  } else if (!is.na(use_sym)) {
    raw <- as.character(body[[use_sym]]); raw[raw %in% c("---", "", "NA")] <- NA_character_
    vapply(strsplit(raw, "\\s*///\\s*"), function(ss) {
      ss <- toupper(trimws(ss)); ss <- ss[nzchar(ss)]
      ent <- unique(unlist(sym2ent[ss], use.names = FALSE))
      ent <- ent[!is.na(ent) & nzchar(ent)]
      if (!length(ent)) NA_character_ else paste(sort(ent), collapse = " /// ")
    }, character(1), USE.NAMES = FALSE)
  } else {
    # GPL16686's platform_table carries NO gene column at all (only genomic ranges and
    # GB_ACC), so the map comes from the matching Bioconductor transcript-cluster db.
    ANNO_DB <- c(GPL16686 = "hugene20sttranscriptcluster.db")
    pkg <- ANNO_DB[[GPLid]]
    if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE))
      stop(GPLid, ": neither an Entrez nor a symbol column exists in the platform table",
           " and no Bioconductor annotation db is available")
    message("no gene column in the SOFT -> using Bioconductor db ", pkg)
    suppressMessages(library(pkg, character.only = TRUE))
    db <- get(pkg, envir = asNamespace(pkg))
    kk <- AnnotationDbi::keys(db)
    mp <- suppressMessages(AnnotationDbi::select(db, keys = kk, columns = "ENTREZID",
                                                keytype = "PROBEID"))
    mp <- mp[!is.na(mp$ENTREZID), ]
    agg <- tapply(as.character(mp$ENTREZID), as.character(mp$PROBEID), function(x)
      paste(sort(unique(x)), collapse = " /// "))
    unname(agg[as.character(body[[1]])])
  }

  id <- as.character(body[[1]])
  out12 <- data.frame(ID_REF = id, ENTREZ_GENE_ID = gene_vec,
                      stringsAsFactors = FALSE, check.names = FALSE)
  out12 <- out12[!duplicated(out12$ID_REF), , drop = FALSE]
  outgene <- data.frame(gene_id = out12$ENTREZ_GENE_ID, row.names = out12$ID_REF,
                        stringsAsFactors = FALSE)
  dir.create(file.path(ROOT, "data/gpl"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(ROOT, "data/processed/gpl"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out12,  file.path(ROOT, "data/gpl", paste0(GPLid, ".rds")))
  saveRDS(outgene, file.path(ROOT, "data/processed/gpl", paste0(GPLid, ".rds")))
  n <- nrow(out12); nm <- sum(!is.na(gene_vec) & nzchar(gene_vec))
  cat(sprintf("%s: %d probes, %d mapped to ENTREZ (%.2f%%), %d unmapped\n",
              GPLid, n, nm, 100 * nm / n, n - nm))
  cat("  multi-gene probes:", sum(grepl(" /// ", gene_vec), na.rm = TRUE), "\n")
  cat("  saved data/gpl/", GPLid, ".rds + data/processed/gpl/", GPLid, ".rds\n", sep = "")
  res[[GPLid]] <- c(probes = n, mapped = nm, pct = round(100 * nm / n, 2))
}
if (length(res)) { cat("\n---- summary ----\n"); print(do.call(rbind, res)) }

# ---------------------------------------------------------------------------
# Gene-level cohorts: the matrix rows ARE gene identifiers, so no platform SOFT
# exists. Build the same two-column map from the matrix's own row ids.
#   GPL16422           E-MTAB-3580  NanoString nCounter gene symbols
#   E-MTAB-4321_PLAT   E-MTAB-4321  RNA-seq FPKM, ENSG (versioned) row ids
# (`<ACC>_PLAT` follows the existing CGGA_<n>_PLAT convention; 16_verify excludes
#  mirror tables ending in _PLAT from the orphan check.)
# ----------------------------------------------------------------------------
if ("--genelevel" %in% commandArgs(trailingOnly = TRUE) || length(grep("^--genelevel", args))) {
  GENELEVEL <- list(
    GPL16422         = list(acc = "E-MTAB-3580", kind = "symbol"),
    `E-MTAB-4321_PLAT` = list(acc = "E-MTAB-4321", kind = "ensg")
  )
  bm2 <- readRDS(file.path(ROOT, "REF_ID/www/biomart_table_Hs.Rds"))
  sym2e <- tapply(as.character(bm2$Entrez[!is.na(bm2$Entrez)]),
                  as.character(bm2$Symbol[!is.na(bm2$Entrez)]),
                  function(x) paste(sort(unique(x)), collapse = " /// "))
  ens2e <- tapply(as.character(bm2$Entrez[!is.na(bm2$Entrez)]),
                  as.character(bm2$Ensembl_gene[!is.na(bm2$Entrez)]),
                  function(x) paste(sort(unique(x)), collapse = " /// "))
  for (tok in names(GENELEVEL)) {
    z <- GENELEVEL[[tok]]
    ef <- file.path(ROOT, "data/expr", paste0(z$acc, ".rds"))
    if (!file.exists(ef)) { message("!! ", tok, ": ", ef, " missing"); next }
    ids <- as.character(readRDS(ef)$ID_REF)
    gv <- if (z$kind == "symbol") unname(sym2e[ids]) else unname(ens2e[sub("\\..*$", "", ids)])
    gv[is.na(gv) | !nzchar(gv)] <- NA_character_
    o12 <- data.frame(ID_REF = ids, ENTREZ_GENE_ID = gv, stringsAsFactors = FALSE,
                      check.names = FALSE)
    og <- data.frame(gene_id = gv, row.names = ids, stringsAsFactors = FALSE)
    saveRDS(o12, file.path(ROOT, "data/gpl", paste0(tok, ".rds")))
    saveRDS(og, file.path(ROOT, "data/processed/gpl", paste0(tok, ".rds")))
    cat(sprintf("%s (%s): %d row ids, %d mapped to ENTREZ (%.2f%%)\n",
                tok, z$acc, length(ids), sum(!is.na(gv)), 100 * mean(!is.na(gv))))
  }
}
}
