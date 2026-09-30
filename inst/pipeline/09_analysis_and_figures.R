# ===========================================================================
# CanPAS curation pipeline -- 09_analysis_and_figures
# ---------------------------------------------------------------------------
# Consolidated script file.  Original pipeline/R scripts merged here, in this
# order:
#   1. batch_integrate.R
#   2. batch_integrate2.R
#   3. demo_meta_lung.R
#   4. demo_tcga_integration.R
#   5. demo_tcga_ondemand.R
#   6. demo_unified_reader.R
#   7. test_cpas_dataset.R
#   8. test_cpas_GSE44001.R
#   9. validate_cpas.R
#
# Each block below is the VERBATIM text of the named original script, wrapped
# in a zero-argument runner function so that source()-ing this file defines
# functions only and has no side effects.  Run one step with, e.g.:
#   source(system.file("pipeline", "09_analysis_and_figures.R", package = "CanPAS"))
#   run_batch_integrate()
# Steps that read/write the project tree take the data root from
# CPAS_DATA_ROOT (or from their own defaults), exactly as the originals did.
#
# Shipped as scripts only: no data files are included in this package.
# ===========================================================================


# ---------------------------------------------------------------------------
# run_batch_integrate()  <-  verbatim pipeline/R/batch_integrate.R
# ---------------------------------------------------------------------------
run_batch_integrate <- function() {
# batch_integrate.R — 步骤1(签名/多基因) + 步骤2(分癌种整合) + 步骤3(多因素/LOO/Stouffer/方向一致性) + 整合KM
suppressMessages({library(CanPAS); library(survival); library(ggplot2); library(RMySQL); library(dplyr)})
ROOT <- "~/data/Project/CanPAS"; ROOT <- path.expand(ROOT)
OUT <- file.path(ROOT, "pipeline/out")
dataset_info <- read.csv(file.path(ROOT, "data/dataset_info.csv"), stringsAsFactors=FALSE)
con <- dbConnect(MySQL(), host="139.224.80.159", dbname="cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CPAS"), password = Sys.getenv("CPAS_DB_PASSWORD"))
tabs <- dbListTables(con)
surv_db <- sub("_surv$", "", tabs[grepl("_surv$", tabs)])
expr_db <- setdiff(tabs[grepl("^GSE", tabs)], paste0(surv_db, "_surv"))
dbDisconnect(con)

# 队列缓存（含 fetch retry）
cache <- new.env()
fetch_merged <- function(acc, genes, max_try=4) {
  key <- paste0(acc, "|", paste(genes, collapse=","))
  if (exists(key, cache)) return(get(key, cache))
  for (i in seq_len(max_try)) {
    r <- tryCatch({
      ex <- get_expr_data(acc, genes, process_duplicates="max")
      es <- merge_surv_expr(acc, ex)
      list(merged = es$merged_data, genes = setdiff(colnames(ex$expr_data), "ID"))
    }, error=function(e) e)
    if (!inherits(r,"error")) { assign(key, r, cache); return(r) }
    Sys.sleep(3 + 2*i)
  }
  NULL
}
endpoint_of <- function(merged) {
  cols <- colnames(merged)
  hit <- intersect(c("OS","DSS","DFS","RFS","PFS","MFS","DRFS","EFS"), sub("_(status|time)$","",cols))
  if (length(hit)) hit[1] else NA_character_
}

cohorts_of <- function(Type) {
  accs <- unique(dataset_info$Accession[dataset_info$Type == Type])
  accs[accs %in% expr_db & accs %in% surv_db]
}

save_meta <- function(res, tag) {
  write.csv(res$per_cohort, file.path(OUT, paste0("int_", tag, "_per_cohort.csv")), row.names=FALSE)
  write.csv(res$pooled,   file.path(OUT, paste0("int_", tag, "_pooled.csv")),    row.names=FALSE)
  if (nrow(res$per_cohort) >= 2) {
    png(file.path(OUT, paste0("int_", tag, "_forest.png")), width=1500, height=400+80*nrow(res$per_cohort), res=130)
    print(plot_meta_forest(res)); dev.off()
  }
  cat("SAVED", tag, "| k=", nrow(res$per_cohort), " pooledHR=", round(res$pooled$HR,3),
      " p=", format(res$pooled$p, digits=2), " I2=", round(100*res$pooled$I2,0), " stouffer_p=", format(res$pooled$p_stouffer,digits=2), "\n")
}

all_res <- list()

# ============ 步骤1：肺癌 多基因签名公式 ============
cat("\n### STEP1 lung signature 0.5*TP53+0.3*GAPDH (OS)\n")
lung <- cohorts_of("Lung cancer")
m_genes <- c("TP53","GAPDH")
lung_merged <- lapply(lung, fetch_merged, genes=m_genes)
names(lung_merged) <- lung
lung_merged <- lung_merged[!vapply(lung_merged, is.null, logical(1))]
r1 <- cpas_meta(tables=names(lung_merged), marker="0.5*TP53+0.3*GAPDH", type="OS",
                method="RE", merged=lapply(lung_merged, `[[`, "merged"))
save_meta(r1, "lung_OS_sig_TP53GAPDH"); all_res[["lung_OS_sig"]] <- r1

# ============ 步骤2：分癌种整合（GAPDH，每癌种选样本量最多的共同终点）============
cat("\n### STEP2 per-cancer integration (GAPDH)\n")
meta_cancers <- list("Breast cancer"=c("RFS","DFS","OS","MFS"), "colorectal cancer"=c("OS","RFS"),
                     "gastric cancer"="OS", "ovarian"="OS", "Lung cancer"="OS", "bladder cancer"="OS",
                     "Prostate cancer"="RFS", "Multiple myeloma"="OS")
for (Type in names(meta_cancers)) {
  accs <- cohorts_of(Type)
  if (length(accs) < 3) { cat("skip", Type, length(accs), "\n"); next }
  cand <- lapply(accs, fetch_merged, genes="GAPDH")
  names(cand) <- accs
  cand <- cand[!vapply(cand, is.null, logical(1))]
  ep_count <- table(unlist(lapply(cand, function(r) endpoint_of(r$merged))))
  prefer <- names(meta_cancers[[Type]])
  if (!length(prefer)) prefer <- names(ep_count)
  use_ep <- prefer[prefer %in% names(ep_count)][1]
  if (is.na(use_ep) || ep_count[use_ep] < 3) {
    use_ep <- names(sort(ep_count, decreasing=TRUE))[1]
    if (is.na(use_ep)) { cat("no endpoint for", Type, "\n"); next }
  }
  cohorts_ep <- names(cand)[vapply(cand, function(r) paste0(use_ep,"_time") %in% colnames(r$merged), logical(1))]
  if (length(cohorts_ep) < 3) { cat("skip", Type, use_ep, length(cohorts_ep), "\n"); next }
  r <- cpas_meta(tables=cohorts_ep, marker="GAPDH", type=use_ep, method="RE",
                 merged=lapply(cand[cohorts_ep], `[[`, "merged"))
  tag <- paste0(gsub(" |/","_",Type), "_", use_ep)
  save_meta(r, tag); all_res[[tag]] <- r
}

# ============ 步骤3：功能升级（以肺癌 OS GAPDH 为例）============
cat("\n### STEP3 upgrades on lung OS GAPDH\n")
r0 <- cpas_meta(tables=names(lung_merged), marker="GAPDH", type="OS", method="RE",
                merged=lapply(lung_merged, `[[`, "merged"))
write.csv(r0$pooled, file.path(OUT,"int_lung_OS_upgrades_pooled.csv"), row.names=FALSE)
# 3a 多因素(可用年龄/性别/分期时) → 每队列 adjusted
r_adj <- cpas_meta(tables=names(lung_merged), marker="GAPDH", type="OS", method="RE",
                   confounders=c("age","sex","stage"),
                   merged=lapply(lung_merged, `[[`, "merged"))
save_meta(r_adj, "lung_OS_GAPDH_adjusted"); all_res[["lung_OS_adj"]] <- r_adj
# 3b leave-one-out
loo <- loo_meta(r0); write.csv(loo, file.path(OUT,"int_lung_OS_loo.csv"), row.names=FALSE)
cat("LOO rows:", nrow(loo), " pooled HR range:", round(range(loo$HR),3), "\n")
# 3c Stouffer 已在 pooled
# 3d 方向一致性（队列 HR vs pooled 方向 + 二项检验）
dir_tab <- r0$per_cohort; pooled_dir <- r0$pooled$logHR > 0
concord <- sum((dir_tab$logHR > 0) == pooled_dir, na.rm=TRUE)
cat("direction concordance:", concord, "/", nrow(dir_tab),
    " binomial p:", round(binom.test(concord, nrow(dir_tab))$p.value,4), "\n")

# ============ 整合 KM（lung OS GAPDH + breast RFS GAPDH）============
cat("\n### INTEGRATED KM\n")
for (spec in list(c("lung_OS_km","Lung cancer","OS","GAPDH"),
                  c("breast_RFS_km","Breast cancer","RFS","GAPDH"),
                  c("crc_OS_km","colorectal cancer","OS","GAPDH"))) {
  tag <- spec[1]; Type <- spec[2]; type <- spec[3]; gene <- spec[4]
  accs <- cohorts_of(Type)
  mds <- lapply(accs, fetch_merged, genes=gene); names(mds) <- accs
  mds <- mds[!vapply(mds,is.null,logical(1))]
  km <- tryCatch(cpas_km_pooled(lapply(mds,`[[`,"merged"), marker=gene, type=type), error=function(e)e)
  if (inherits(km,"error")) { cat("KM fail", tag, km$message,"\n"); next }
  cat("KM", tag, "| cohorts:", length(km$cohorts), " pooled logrank p:", format(km$logrank_p,digits=3), "\n")
  png(file.path(OUT, paste0("int_", tag, "_km.png")), width=1500, height=1300, res=140)
  print(plot_cpas_km(km)); dev.off()
}
cat("\nALL DONE\n")
}

# ---------------------------------------------------------------------------
# run_batch_integrate2()  <-  verbatim pipeline/R/batch_integrate2.R
# ---------------------------------------------------------------------------
run_batch_integrate2 <- function() {
# batch_integrate2.R — 补跑：正确癌种名的 gastric/ovarian/bladder + MM 重试 + KM 计数确认
suppressMessages({library(CanPAS); library(survival); library(ggplot2); library(RMySQL); library(dplyr)})
ROOT <- "~/data/Project/CanPAS"; ROOT <- path.expand(ROOT); OUT <- file.path(ROOT,"pipeline/out")
dataset_info <- read.csv(file.path(ROOT,"data/dataset_info.csv"), stringsAsFactors=FALSE)
con <- dbConnect(MySQL(), host="139.224.80.159", dbname="cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CPAS"), password = Sys.getenv("CPAS_DB_PASSWORD"))
tabs <- dbListTables(con); surv_db <- sub("_surv$","",tabs[grepl("_surv$",tabs)])
expr_db <- setdiff(tabs[grepl("^GSE",tabs)], paste0(surv_db,"_surv")); dbDisconnect(con)
cache <- new.env()
fetch_merged <- function(acc, genes="GAPDH", max_try=5){
  key<-paste0(acc,"|",genes)
  if(exists(key,cache)) return(get(key,cache))
  for(i in 1:max_try){ r<-tryCatch({ex<-get_expr_data(acc,genes,process_duplicates="max");es<-merge_surv_expr(acc,ex)
      list(merged=es$merged_data, genes=setdiff(colnames(ex$expr_data),"ID"))},error=function(e)e)
    if(!inherits(r,"error")){assign(key,r,cache);return(r)}; Sys.sleep(4+2*i)}; NULL }
endpoint_of <- function(m) intersect(c("OS","DSS","DFS","RFS","PFS","MFS","DRFS","EFS"), sub("_(status|time)$","",colnames(m)))[1]
cohorts_of <- function(Type) unique(dataset_info$Accession[dataset_info$Type==Type])
save_meta <- function(res,tag){ write.csv(res$per_cohort,file.path(OUT,paste0("int_",tag,"_per_cohort.csv")),row.names=FALSE)
  write.csv(res$pooled,file.path(OUT,paste0("int_",tag,"_pooled.csv")),row.names=FALSE)
  cat("SAVED",tag,"k=",nrow(res$per_cohort),"HR=",round(res$pooled$HR,3),"p=",format(res$pooled$p,digits=2),"I2=",round(100*res$pooled$I2,0),"\n") }
for (Type in c("Gastric cancer","Ovarian cancer","Bladder cancer","Multiple myeloma")) {
  accs <- cohorts_of(Type)
  cand <- lapply(accs, fetch_merged); names(cand) <- accs
  cand <- cand[!vapply(cand,is.null,logical(1))]
  if(!length(cand)){cat("no data",Type,"\n");next}
  eps <- table(unlist(lapply(cand,function(r) endpoint_of(r$merged))))
  use_ep <- names(sort(eps,decreasing=TRUE))[1]
  tab <- names(cand)[vapply(cand,function(r) paste0(use_ep,"_time") %in% colnames(r$merged),logical(1))]
  if(length(tab)>=3){ r<-cpas_meta(tables=tab,marker="GAPDH",type=use_ep,method="RE",merged=lapply(cand[tab],`[[`,"merged"))
    save_meta(r,paste0(gsub(" |/","_",Type),"_",use_ep)) } else cat("skip",Type,"ep=",use_ep,"k=",length(tab),"\n")
}
# KM 计数确认（三种）
for (spec in list(c("lung_OS_km","Lung cancer","OS"),c("breast_RFS_km","Breast cancer","RFS"),c("crc_OS_km","colorectal cancer","OS"))) {
  accs <- cohorts_of(spec[2]); mds <- lapply(accs, fetch_merged); names(mds)<-accs
  mds<-mds[!vapply(mds,is.null,logical(1))]
  km <- cpas_km_pooled(lapply(mds,`[[`,"merged"), marker="GAPDH", type=spec[3])
  cat("KM",spec[1],"cohorts:",paste(km$cohorts,collapse=",")," n=",nrow(km$df)," logrank_p=",format(km$logrank_p,digits=3),"\n")
}
cat("DONE2\n")
}

# ---------------------------------------------------------------------------
# run_demo_meta_lung()  <-  verbatim pipeline/R/demo_meta_lung.R
# ---------------------------------------------------------------------------
run_demo_meta_lung <- function() {
# demo_meta_lung.R — 肺癌 OS 多队列整合分析演示（两阶段 meta）
suppressMessages({library(CanPAS); library(survival); library(ggplot2); library(RMySQL); library(dplyr)})
ROOT <- "~/data/Project/CanPAS"; ROOT <- path.expand(ROOT)
dataset_info <- read.csv(file.path(ROOT, "data/dataset_info.csv"), stringsAsFactors=FALSE)
con <- dbConnect(MySQL(), host="139.224.80.159", dbname="cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CPAS"), password = Sys.getenv("CPAS_DB_PASSWORD"))
tabs <- dbListTables(con)
surv_db <- sub("_surv$", "", tabs[grepl("_surv$", tabs)])
expr_db <- setdiff(tabs[grepl("^GSE", tabs)], paste0(surv_db, "_surv"))
dbDisconnect(con)

lung <- dataset_info$Accession[dataset_info$Type == "Lung cancer"]
lung <- unique(lung[lung %in% expr_db & lung %in% surv_db])   # 只留 表达+生存 齐全
cat("lung OS cohorts to integrate:", length(lung), "\n"); print(lung)

marker <- "GAPDH"; type <- "OS"
res <- cpas_meta(tables = lung, marker = marker, type = type, method = "RE", min_events = 5, max_try = 4)
cat("\n==== pooled (RE) ====\n"); print(res$pooled, row.names=FALSE)
cat("\n==== per cohort ====\n"); print(res$per_cohort[, c("table","n","events","HR","lower","upper","p")], row.names=FALSE)
cat("\ncohort errors (skipped):\n"); print(res$errors)

write.csv(res$per_cohort, file.path(ROOT,"pipeline/out","meta_lung_per_cohort.csv"), row.names=FALSE)
write.csv(res$pooled,  file.path(ROOT,"pipeline/out","meta_lung_pooled.csv"),  row.names=FALSE)
png(file.path(ROOT,"pipeline/out","meta_lung_forest_GAPDH.png"), width=1600, height=1500, res=150)
print(plot_meta_forest(res)); dev.off()
cat("saved: meta_lung_per_cohort.csv / meta_lung_pooled.csv / meta_lung_forest_GAPDH.png\n")
}

# ---------------------------------------------------------------------------
# run_demo_tcga_integration()  <-  verbatim pipeline/R/demo_tcga_integration.R
# ---------------------------------------------------------------------------
run_demo_tcga_integration <- function() {
suppressMessages({library(CanPAS); library(survival); library(dplyr)})
# 1) 保留癌种清单 + 每项目样本/OS 概况（离线 tcga_surv.rda）
rows <- lapply(names(CanPAS::tcga_retained), function(proj){
  st <- tryCatch(tcga_surv_table(proj), error=function(e) NULL)
  if (is.null(st)) return(NULL)
  data.frame(project=proj, canonical_type=canonical_type(proj),
             dataset_short=short_name(proj),
             n=sum(!is.na(st$OS_status)), events=sum(st$OS_status==1,na.rm=TRUE),
             has_OS=("OS_status" %in% colnames(st)), has_DFI=("DFI_status" %in% colnames(st)))
})
cat <- do.call(rbind, rows)
write.csv(cat, "pipeline/out/tcga_catalog.csv", row.names=FALSE)
print(cat, row.names=FALSE)
# 2) 端到端：LUAD 与 BRCA 用 tcga_merged -> 直接进 CanPAS COX 与 meta
d1 <- tcga_merged("LUAD", c("TP53","GAPDH"), "OS")
d2 <- tcga_merged("BRCA", c("TP53","GAPDH"), "OS")
cat("\nLUAD merged:", nrow(d1), " BRCA merged:", nrow(d2), "\n")
# meta（两队列 OS GAPDH，merged-cache 用法，与 GEO 队列完全一致）
r <- cpas_meta(tables=c("TCGA-LUAD","TCGA-BRCA"), marker="GAPDH", type="OS", method="RE",
               merged=list("TCGA-LUAD"=d1, "TCGA-BRCA"=d2))
print(r$per_cohort[,c("table","n","events","HR","lower","upper","p")], row.names=FALSE)
print(r$pooled)
cat("\nLUAD+BRCA TCGA meta done\n")
}

# ---------------------------------------------------------------------------
# run_demo_tcga_ondemand()  <-  verbatim pipeline/R/demo_tcga_ondemand.R
# ---------------------------------------------------------------------------
run_demo_tcga_ondemand <- function() {
suppressMessages({library(CanPAS); library(survival); library(dplyr)})
# TCGA 按需(genes)从 Xena 取表达 + 本地 rda 生存 —— 全链路测试
cat("== tcga_surv_table(LUAD) cols/rows:\n")
s <- tcga_surv_table("LUAD")
cat(nrow(s), "x", ncol(s), "|", paste(head(colnames(s), 8), collapse=","), "\n")
cat("OS events:", sum(s$OS_status==1,na.rm=TRUE), " age NA:", sum(is.na(s$age)), " stage sample:", paste(head(unique(s$stage),3),collapse=","), "\n")
cat("\n== tcga_merged(LUAD, TP53+GAPDH, OS) —— 表达按需从 Xena:\n")
d <- tcga_merged("LUAD", c("TP53","GAPDH"), "OS")
cat("merged:", nrow(d), "x", ncol(d), "\n")
cat("OS events:", sum(d$OS_status==1,na.rm=TRUE), " GAPDH nonNA:", sum(!is.na(d$GAPDH)), "\n")
cox <- COX_analysis(d, type="OS", cont_Variates=c("TP53","GAPDH"), cate_Variates=NULL, method="uni")
print(cox$results_table)
cat("\n== cohort_merged 统一读取: TCGA-BRCA vs GEO 同结构\n")
b <- cohort_merged("TCGA-BRCA", c("TP53"), "OS")
cat("TCGA-BRCA merged:", nrow(b), "x", ncol(b), " cols:", paste(head(colnames(b),6),collapse=","), "\n")
cat("DONE\n")
}

# ---------------------------------------------------------------------------
# run_demo_unified_reader()  <-  verbatim pipeline/R/demo_unified_reader.R
# ---------------------------------------------------------------------------
run_demo_unified_reader <- function() {
suppressMessages({library(CanPAS); library(survival); library(dplyr)})
ROOT <- "/home/Jingle/data/Project/CPAS"
# 1) 物化 15 个 TCGA 项目的 GEO 同构 surv（tcga_clinical + tcga_surv）
projects <- names(CanPAS::tcga_retained)
cat("materializing surv for", length(projects), "projects:\n")
for (pr in projects) {
  f <- tryCatch(tcga_surv_file(pr, root = ROOT), error = function(e) NA)
  if (is.na(f)) { cat("  ", pr, "FAIL\n"); next }
  s <- readRDS(f)
  cat(sprintf("  TCGA-%-4s n=%d OS_ev=%d cols=%s\n", pr, nrow(s),
              sum(s$OS_status == 1, na.rm = TRUE), paste(colnames(s), collapse = ",")))
}
# 2) 表达物化（示例基因 TP53/GAPDH，LUAD/LUSC/BRCA）
for (pr in c("LUAD","LUSC","BRCA")) {
  f <- tcga_expr_file(pr, genes = c("TP53","GAPDH"), root = ROOT)
  ex <- readRDS(f); cat("  expr file:", basename(f), nrow(ex), "genes x", ncol(ex)-1, "samples\n")
}
# 3) 统一读取器：TCGA 与 GEO 结构一致并可直接用包函数
d <- cohort_merged("TCGA-LUAD", genes = c("TP53","GAPDH"), type = "OS")
cat("\ncohort_merged(TCGA-LUAD):", nrow(d), "x", ncol(d), "\n")
cat("cols:", paste(head(colnames(d), 12), collapse=", "), "\n")
cox <- COX_analysis(d, type="OS", cont_Variates=c("TP53","GAPDH"), cate_Variates=NULL, method="uni")
print(cox$results_table)
# 4) 混合整合：TCGA-LUAD + TCGA-LUSC（本地）与 GEO 队列 GSE13213（API）同列表 cpas_meta
d2 <- cohort_merged("TCGA-LUSC", c("TP53","GAPDH"), "OS")
geo <- tryCatch(cohort_merged("GSE13213", c("TP53","GAPDH"), "OS"), error=function(e) NULL)
merged_list <- list("TCGA-LUAD"=d, "TCGA-LUSC"=d2)
if (!is.null(geo)) merged_list[["GSE13213(GEO)"]] <- geo
r <- cpas_meta(tables=names(merged_list), marker="GAPDH", type="OS", method="RE",
               merged=lapply(merged_list, function(x) x))
print(r$per_cohort[, c("table","n","events","HR","lower","upper","p")], row.names=FALSE)
print(r$pooled)
cat("MIXED TCGA+GEO meta OK\n")
}

# ---------------------------------------------------------------------------
# run_test_cpas_dataset()  <-  verbatim pipeline/R/test_cpas_dataset.R
# ---------------------------------------------------------------------------
run_test_cpas_dataset <- function() {
# test_cpas_dataset.R — 本地 CanPAS 端到端测试（任意已上传数据集）
# 用法: Rscript pipeline/R/test_cpas_dataset.R <ACC> <genes,comma> <type> <marker> [<ProjectRoot>] [<predict.time>]
# 例  : Rscript pipeline/R/test_cpas_dataset.R GSE39582 TP53,GAPDH,ACTB OS TP53 ~/data/Project/CanPAS 3
suppressMessages({
  library(CanPAS); library(survival); library(survminer)
  library(ggplot2); library(dplyr)
})
args <- commandArgs(trailingOnly=TRUE)
ACC  <- args[1]
genes <- strsplit(args[2], ",")[[1]]
TYPE <- args[3]; MARKER <- args[4]
ROOT <- if (length(args) >= 5) args[5] else "~/data/Project/CanPAS"
PT   <- if (length(args) >= 6) as.numeric(args[6]) else 3
ROOT <- path.expand(ROOT)
OUT  <- file.path(ROOT, "pipeline/out"); dir.create(OUT, showWarnings=FALSE)
data(ID_map)
if (file.exists(file.path(ROOT, "data/dataset_info.csv"))) {
  dataset_info <- read.csv(file.path(ROOT, "data/dataset_info.csv"),
                           stringsAsFactors=FALSE)   # 始终用最新登记表
} else data(dataset_info)
cat("== dataset_info rows:", nrow(dataset_info), "|", ACC, "registered:",
    ACC %in% dataset_info$Accession, "\n")
stopifnot(ACC %in% dataset_info$Accession)

expr_data <- get_expr_data(ACC, genes, process_duplicates="max")
cat("\n[get_expr_data] genes matched:", expr_data$metadata$genes_matched_platform,
    " samples:", expr_data$metadata$sample_count, "\n")
es <- merge_surv_expr(ACC, expr_data)
cat("[merge_surv_expr] merged:", es$metadata$samples_merged,
    " surv cols:", paste(es$metadata$survival_columns, collapse=","), "\n")
df <- es$merged_data
cat("cols:", paste(colnames(df), collapse=", "), " dim:", dim(df), "\n")
cat("\n-- head --\n"); print(head(df, 4))
ev_col <- paste0(TYPE, "_status"); tm_col <- paste0(TYPE, "_time")
tval <- suppressWarnings(as.numeric(df[[tm_col]]))
cat("\n-- ", TYPE, " events:", sum(df[[ev_col]]==1, na.rm=TRUE), "/",
    sum(!is.na(df[[ev_col]])), " time range:",
    paste(round(range(tval, na.rm=TRUE),3), collapse=" - "), "(yr)\n")

cox <- COX_analysis(df, type=TYPE, cont_Variates=genes,
                    cate_Variates=NULL, method="uni", precision=3)
cat("\n-- univariate COX (", TYPE, ") --\n")
print(cox$results_table)
write.csv(cox$results_table, file.path(OUT, paste0(ACC, "_COX_uni_", TYPE, ".csv")),
          row.names=FALSE)

p <- plot_km(df, type=TYPE, marker=MARKER, pval=TRUE, legend="bottom", risk.table=FALSE)
png(file.path(OUT, paste0(ACC, "_KM_", TYPE, "_", MARKER, ".png")),
    width=1400, height=1200, res=150); print(p); dev.off()
r <- plot_roc(df, type=TYPE, marker=MARKER, predict.time=PT)
png(file.path(OUT, paste0(ACC, "_ROC_", TYPE, "_", MARKER, ".png")),
    width=1400, height=1200, res=150); print(r); dev.off()
cat("\nplots saved: ", ACC, "_KM_", TYPE, "_", MARKER, ".png / _ROC_... (", PT, "y)\n", sep="")
}

# ---------------------------------------------------------------------------
# run_test_cpas_GSE44001()  <-  verbatim pipeline/R/test_cpas_GSE44001.R
# ---------------------------------------------------------------------------
run_test_cpas_GSE44001 <- function() {
# test_cpas_GSE44001.R — 本地使用 CanPAS 包对 GSE44001（cervical, GPL14951, DFS）做端到端测试
# 流程: get_expr_data(基因->平台探针) -> merge_surv_expr(与 GSE44001_surv 合并)
#       -> COX_analysis / plot_km / plot_roc (type="DFS")
# 用法: Rscript pipeline/R/test_cpas_GSE44001.R [<ProjectRoot>]
suppressMessages({
  library(CanPAS); library(survival); library(survminer)
  library(ggplot2); library(dplyr)
})
ROOT <- if (length(commandArgs(trailingOnly=TRUE)) >= 1)
  commandArgs(trailingOnly=TRUE)[1] else "~/data/Project/CanPAS"
ROOT <- path.expand(ROOT)
OUT  <- file.path(ROOT, "pipeline/out"); dir.create(OUT, showWarnings=FALSE)

data(ID_map)
data(dataset_info)
cat("== CanPAS pkg dataset_info 行数:", nrow(dataset_info),
    "| GSE44001 已登记:", "GSE44001" %in% dataset_info$Accession, "\n")
stopifnot("GSE44001" %in% dataset_info$Accession)

table <- "GSE44001"
genes <- c("TP53", "GAPDH", "ACTB")

cat("\n[1] get_expr_data(", table, ", genes:", paste(genes, collapse="/"), ")\n")
expr_data <- get_expr_data(table, genes, process_duplicates = "max")
print(expr_data$metadata)
cat("\n- expr head:\n"); print(head(expr_data$expr_data, 3))

cat("\n[2] merge_surv_expr(", table, ")\n")
es <- merge_surv_expr(table, expr_data)
print(es$metadata)
cat("\n- 合并表列:", paste(colnames(es$merged_data), collapse=", "), "\n")
cat("- 合并表维度:", dim(es$merged_data), "\n")
df <- es$merged_data
print(head(df, 5))
cat("- DFS_status 计数:\n"); print(table(df$DFS_status, useNA="ifany"))
cat("- DFS_time(年) 范围:", range(df$DFS_time, na.rm=TRUE), "\n")

cat("\n[3] COX_analysis (uni, type=DFS)\n")
cox <- COX_analysis(df, type="DFS", cont_Variates=genes,
                    cate_Variates=NULL, method="uni", precision=3)
print(cox$results_table)
write.csv(cox$results_table, file.path(OUT, "GSE44001_COX_uni_DFS.csv"), row.names=FALSE)

cat("\n[4] plot_km TP53 (DFS, median cut)\n")
p_km <- plot_km(df, type="DFS", marker="TP53", pval=TRUE, legend="bottom",
                risk.table=FALSE)
png(file.path(OUT, "GSE44001_KM_DFS_TP53.png"), width=1400, height=1200, res=150)
print(p_km); dev.off()
cat("  -> saved GSE44001_KM_DFS_TP53.png\n")

cat("\n[5] plot_roc TP53 (DFS, predict.time=3y)\n")
p_roc <- plot_roc(df, type="DFS", marker="TP53", predict.time=3)
png(file.path(OUT, "GSE44001_ROC_DFS_TP53.png"), width=1400, height=1200, res=150)
print(p_roc); dev.off()
cat("  -> saved GSE44001_ROC_DFS_TP53.png\n")

cat("\n== CanPAS local test on GSE44001 finished OK ==\n")
}

# ---------------------------------------------------------------------------
# run_validate_cpas()  <-  verbatim pipeline/R/validate_cpas.R
# ---------------------------------------------------------------------------
run_validate_cpas <- function() {
# validate_cpas.R — 用 CanPAS 包批量验证每个入库数据集：取表达→合并生存→单因素COX
# 用法: Rscript pipeline/R/validate_cpas.R [<ProjectRoot>]
# 输出: pipeline/out/CPAS_validation_<date>.csv  （逐数据集 PASS/FAIL + 诊断）
suppressMessages({library(CanPAS); library(survival); library(dplyr); library(RMySQL)})

ROOT <- if (length(commandArgs(trailingOnly=TRUE)) >= 1) commandArgs(trailingOnly=TRUE)[1] else "~/data/Project/CanPAS"
ROOT <- path.expand(ROOT); setwd(ROOT)
data(ID_map)
dataset_info <- read.csv("data/dataset_info.csv", stringsAsFactors=FALSE)   # 最新登记（Accession 与 DB 表名一致）

con <- dbConnect(MySQL(), host="139.224.80.159", dbname="cpas", user = Sys.getenv("CPAS_DB_USER", unset = "CPAS"), password = Sys.getenv("CPAS_DB_PASSWORD"))
tabs <- dbListTables(con)
surv_db <- sub("_surv$", "", tabs[grepl("_surv$", tabs)])
expr_db <- setdiff(tabs[grepl("^GSE", tabs)], paste0(surv_db, "_surv"))
dbDisconnect(con)

PRIORITY <- c("OS","DSS","DFS","RFS","PFS","MFS","DRFS","EFS")
GENES <- c("GAPDH","ACTB","TP53","RPLP0","B2M")

pick_endpoint <- function(cols) {
  hit <- intersect(PRIORITY, sub("_(status|time)$", "", cols))
  if (length(hit)) hit[1] else {
    t <- sub("_time$","", cols[grepl("_time$", cols)])
    if (length(t)) t[1] else NA_character_
  }
}

run_one <- function(acc) {
  has_expr <- acc %in% expr_db
  has_surv <- acc %in% surv_db
  status <- "OK"; detail <- ""
  if (!has_expr && !has_surv) return(data.frame(Accession=acc, type="not-in-db", status="SKIP", endpoint=NA, samples=NA, genes=NA, events=NA, detail="no DB expr/surv table"))
  typ <- if (has_expr && has_surv) "full" else if (has_expr) "expr-only" else "surv-only"

  # 1) 取数（表达谱 + 生存）
  expr <- tryCatch(get_expr_data(acc, GENES, process_duplicates="max"), error=function(e) e)
  if (inherits(expr, "error")) {
    if (typ == "surv-only") {   # 无表达谱表属预期
      s <- tryCatch(get_data(table=acc, action="surv_data"), error=function(e) e)
      if (inherits(s, "error")) return(data.frame(Accession=acc,type=typ,status="FAIL",endpoint=NA,samples=NA,genes=NA,events=NA,detail=paste0("surv fetch err: ", s$message)))
      n <- nrow(s$response)
      cols <- colnames(s$response)[-1]
      return(data.frame(Accession=acc,type=typ,status="OK",endpoint=pick_endpoint(cols),samples=n,genes=NA,events=NA,detail="surv-only (no expr table)"))
    }
    return(data.frame(Accession=acc,type=typ,status="FAIL",endpoint=NA,samples=NA,genes=NA,events=NA,detail=paste0("get_expr_data err: ", substr(expr$message,1,120))))
  }
  gmatch <- expr$metadata$genes_matched_platform
  nsamp   <- expr$metadata$sample_count
  if (is.null(gmatch) || gmatch == 0 || is.null(nsamp) || nsamp == 0)
    return(data.frame(Accession=acc,type=typ,status="FAIL",endpoint=NA,samples=nsamp,genes=gmatch,events=NA,detail="no genes matched / no samples"))

  # 2) 合并生存
  es <- tryCatch(merge_surv_expr(acc, expr), error=function(e) e)
  if (inherits(es, "error"))
    return(data.frame(Accession=acc,type=typ,status="FAIL",endpoint=NA,samples=nsamp,genes=gmatch,events=NA,detail=paste0("merge err: ", substr(es$message,1,120))))
  merged <- es$merged_data
  ep <- pick_endpoint(colnames(merged))
  if (is.na(ep))
    return(data.frame(Accession=acc,type=typ,status="FAIL",endpoint=ep,samples=nsamp,genes=gmatch,events=NA,detail="no *_status/_time endpoints in merged data"))
  ev <- sum(as.numeric(merged[[paste0(ep,"_status")]]) == 1, na.rm=TRUE)
  geneCols <- intersect(setdiff(colnames(expr$expr_data), "ID"), colnames(merged))
  if (!length(geneCols)) geneCols <- intersect(c("GAPDH","ACTB","TP53"), colnames(merged))
  if (!length(geneCols))
    return(data.frame(Accession=acc,type=typ,status="FAIL",endpoint=ep,samples=nsamp,genes=gmatch,events=ev,detail="no gene columns in merged data"))

  # 3) 单因素 COX
  cox <- tryCatch(COX_analysis(df=merged, type=ep, cont_Variates=geneCols,
                               cate_Variates=NULL, method="uni", precision=3), error=function(e) e)
  if (inherits(cox, "error"))
    return(data.frame(Accession=acc,type=typ,status="FAIL",endpoint=ep,samples=nsamp,genes=gmatch,events=ev,detail=paste0("cox err: ", substr(cox$message,1,120))))
  rows_ok <- !is.null(cox$results_table) && nrow(cox$results_table) > 0
  status <- if (rows_ok) "OK" else "WARN"
  data.frame(Accession=acc,type=typ,status=status,endpoint=ep,samples=nsamp,
             genes=gmatch,events=ev,
             detail=paste0("cox rows=", if (rows_ok) nrow(cox$results_table) else 0,
                           if (ev==0) " [0 events!]" else ""))
}

targets <- unique(c(expr_db, surv_db))
targets <- targets[grepl("^GSE", targets)]
cat("validating", length(targets), "DB datasets (full/surv-only/expr-only)\n")
run_retry <- function(acc) {
  r <- tryCatch(run_one(acc), error=function(e)
    data.frame(Accession=acc, type="full", status="FAIL", endpoint=NA, samples=NA,
               genes=NA, events=NA, detail=paste0("run err: ", substr(e$message,1,100))))
  attempt <- 1L
  while (r$status == "FAIL" && attempt <= 4 &&
         grepl("connection|json|timeout|lexical|500|curl|error", r$detail, ignore.case=TRUE)) {
    Sys.sleep(3 + 2*attempt)
    message("retry ", acc, " (", attempt, ")")
    r <- tryCatch(run_one(acc), error=function(e)
      data.frame(Accession=acc, type="full", status="FAIL", endpoint=NA, samples=NA,
                 genes=NA, events=NA, detail=paste0("run err: ", substr(e$message,1,100))))
    attempt <- attempt + 1L
  }
  Sys.sleep(0.6)
  r
}
res <- do.call(rbind, lapply(targets, run_retry))
outf <- file.path(ROOT, "pipeline/out", paste0("CPAS_validation_", format(Sys.Date(), "%Y%m%d"), ".csv"))
write.csv(res, outf, row.names=FALSE)
cat("\n==== summary ====\n")
print(table(res$status, res$type, useNA="ifany"))
cat("OK 全链路:", sum(res$status=="OK" & res$type=="full"), "\n")
cat("FAIL:", sum(res$status=="FAIL"), "\n")
fails <- res[res$status=="FAIL", ]
if (nrow(fails)) { cat("--- FAIL details ---\n"); print(fails[, c("Accession","type","detail")], row.names=FALSE) }
message("written: ", outf, " (", nrow(res), " rows)")
}
