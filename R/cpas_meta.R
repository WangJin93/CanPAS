#' 两阶段跨队列整合（meta）分析
#' @description
#' 对同一癌种的多个队列执行两阶段整合分析：
#'   阶段1 每队列 coxph(Surv(time,status) ~ marker) 得到 logHR ± SE；
#'   阶段2 逆方差加权合并（固定 FE / 随机 RE: DerSimonian-Laird），
#'         输出合并 HR(95%CI)、Z/p、异质性 Q/df/I²/tau²。
#' @param datasets 字符向量：数据集 accession（与 dataset_info 的 Accession 一致）。
#' @param marker 单基因符号（如 "TP53"）或签名公式字符串（如 "0.5*TP53+0.3*GAPDH"）。
#' @param type 终点家族 (OS/DSS/DFS/PFS/MFS) 或原始终点 token (OS/DSS/DFS/RFS/PFS/MFS/DFI/PFI/DRFS/EFS…)。
#'   家族会按队列解析为具体 token (见 \code{\link{endpoint_resolve}}), 并在 per_dataset$endpoint 中记录。
#' @param method "RE"（随机效应，默认）或 "FE"（固定效应）。
#' @param confounders 可选：多因素校正协变量名向量（须存在于 merged 数据，按队列可用集合取交集）。
#' @param min_events 队列纳入最低事件数（默认 5）。
#' @param max_try 每个队列取数失败时的重试次数（应对 API 抖动）。
#' @param merged 可选：命名列表（队列名 -> 已 merge 好的数据框）。提供后不再联网取数。

#' @return class "cpas_meta" 列表：per_dataset、pooled、input、errors。
#' @examples
#' \dontrun{
#'    ## Real cohorts, real endpoint: each cohort contributes the DFS-family
#'    ## token it actually has, which the function resolves and reports.
#'    m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
#'    m$per_dataset[, c("dataset", "endpoint", "n", "events", "HR")]
#'    m$pooled
#'    plot_meta_forest(m)
#'    loo_meta(m)          # leave-one-out sensitivity
#'    print(m)
#' }
#' @export
cpas_meta <- function(datasets, marker, type = "OS",
                      method = c("RE", "FE"),
                      confounders = NULL,
                      min_events = 5,
                      max_try = 3,
                      merged = NULL) {
  method <- match.arg(method)
  if (missing(datasets) || is.null(datasets))
    stop("'datasets' must contain at least one accession.", call. = FALSE)
  datasets <- as.character(datasets)
  if (!length(datasets)) stop("'datasets' must contain at least one accession.", call. = FALSE)
  genes <- if (grepl("[+*:]", marker)) {
    unique(unlist(regmatches(marker, gregexpr("[A-Za-z][A-Za-z0-9._]*", marker))))
  } else marker

  if (anyDuplicated(datasets))
    stop("'datasets' contains duplicated accessions: ",
         paste(unique(datasets[duplicated(datasets)]), collapse = ", "),
         ". Each dataset can be used once.", call. = FALSE)
  if (!is.null(merged) && anyDuplicated(names(merged)))
    stop("'merged' contains duplicated names: ",
         paste(unique(names(merged)[duplicated(names(merged))]), collapse = ", "),
         ". Use unique dataset names.", call. = FALSE)
  per <- list(); errors <- list()
  for (t in datasets) {
    r <- NULL
    if (!is.null(merged) && !is.null(merged[[t]])) {
      df <- merged[[t]]
      ok <- TRUE
    } else {
      ok <- FALSE; err <- NULL
      if (startsWith(t, "TCGA-")) {
        r <- tryCatch(list(df = tcga_merged(sub("^TCGA-", "", t), genes, type)),
                      error = function(e) e)
        if (!inherits(r, "error")) ok <- TRUE else err <- r
      } else {
        for (i in seq_len(max_try)) {
          r <- tryCatch({
            ex <- get_expr_data(t, genes, process_duplicates = "max")
            es <- merge_surv_expr(t, ex)
            list(df = es$merged_data, geneCols = setdiff(colnames(ex$expr_data), "ID"))
          }, error = function(e) e)
          if (!inherits(r, "error")) { ok <- TRUE; break }
          err <- r; Sys.sleep(2 + 1.5 * i)
        }
      }
      if (!ok) { errors[[t]] <- paste0("fetch/merge: ", substr(err$message, 1, 100)); next }
      df <- r$df
    }
    # 家族 (OS/DSS/DFS/PFS/MFS) -> 该队列可用的具体终点; 也兼容原始 token
    tok <- endpoint_resolve(t, type)
    if (is.na(tok)) {
      # 回退: 队列不在 catalog 中(如用户自备数据), 若数据本身含该终点列则按原始 token 处理
      if (all(c(paste0(type, "_time"), paste0(type, "_status")) %in% colnames(df)))
        tok <- as.character(type)[1]
      else { errors[[t]] <- paste0("no ", type, " endpoint"); next }
    }
    tc <- paste0(tok, "_time"); sc <- paste0(tok, "_status")
    if (!all(c(tc, sc) %in% colnames(df))) {
      errors[[t]] <- paste0("no ", type, " endpoint (resolved ", tok,
                            " not in merged data)"); next
    }
    df[[tc]] <- suppressWarnings(as.numeric(df[[tc]]))
    df[[sc]] <- suppressWarnings(as.numeric(df[[sc]]))
    for (gn in genes) if (gn %in% colnames(df)) df[[gn]] <- suppressWarnings(as.numeric(df[[gn]]))

    if (grepl("[+*:]", marker)) {
      gv <- intersect(genes, colnames(df))
      missing_sig <- setdiff(genes, colnames(df))
      if (!length(gv)) {
        errors[[t]] <- paste0("signature genes absent: ",
                              paste(missing_sig, collapse = ",")); next
      }
      score <- tryCatch(eval(parse(text = marker), envir = as.list(df[gv])),
                        error = function(e) NULL)
      if (is.null(score) || all(is.na(score))) {
        errors[[t]] <- paste0("signature eval failed (missing on platform: ",
                              paste(missing_sig, collapse = ","), ")"); next
      }
      df$marker <- score
    } else {
      if (!marker %in% colnames(df)) { errors[[t]] <- paste0("gene missing: ", marker); next }
      df$marker <- suppressWarnings(as.numeric(df[[marker]]))
    }

    # 分析样本 = time/status/marker 与所用协变量的完整个案（保证标准化与建模同一样本）
    conf_avail <- intersect(confounders, colnames(df))
    # 剔除该队列中恒定(单水平)的协变量, 避免 coxph 报错导致整队列被弃
    if (length(conf_avail)) {
      okv <- vapply(conf_avail, function(cn)
        length(unique(df[[cn]][!is.na(df[[cn]])])) > 1L, logical(1))
      conf_avail <- conf_avail[okv]
    }
    for (cn in conf_avail)
      if (is.character(df[[cn]])) {
        num <- suppressWarnings(as.numeric(df[[cn]]))
        df[[cn]] <- if (all(is.na(num) == is.na(df[[cn]]))) num else factor(df[[cn]])
      } else if (is.logical(df[[cn]])) df[[cn]] <- factor(df[[cn]])
    req <- c(tc, sc, "marker", conf_avail)
    keep <- stats::complete.cases(df[req]) & is.finite(df[[tc]]) &
      is.finite(df$marker) & df[[tc]] >= 0
    df <- df[keep, , drop = FALSE]
    if (nrow(df) < 10) { errors[[t]] <- sprintf("insufficient after cleaning (n=%d)", nrow(df)); next }
    sdx <- stats::sd(df$marker)
    if (is.na(sdx) || sdx <= 0) { errors[[t]] <- "marker constant"; next }
    df$marker <- (df$marker - mean(df$marker)) / sdx   # per-SD 标准化(跨平台可比)
    events <- sum(df[[sc]] == 1, na.rm = TRUE)
    if (events < min_events || length(unique(df[[sc]])) < 2) {
      errors[[t]] <- sprintf("insufficient (n=%d, events=%d)", nrow(df), events); next
    }
    vars <- c("marker", conf_avail)
    diag <- .cpas_COX_diag(stats::as.formula("survival::Surv(time, status) ~ ."),
                           data.frame(time = df[[tc]], status = df[[sc]],
                                      df[, vars, drop = FALSE]))
    if (!diag$ok) {
      errors[[t]] <- paste0("cox: ", substr(diag$reason, 1, 120))
      next
    }
    fit <- diag$fit
    b <- stats::coef(fit)[["marker"]]; se <- sqrt(diag(stats::vcov(fit)))[["marker"]]
    per[[t]] <- data.frame(dataset = t, endpoint = tok, n = nrow(df), events = events,
                           HR = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
                           logHR = b, se = se,
                           p = 2 * stats::pnorm(-abs(b / se)), stringsAsFactors = FALSE)
  }

  if (!length(per)) {
    # keep the per-dataset reasons: the caller needs to know whether a gene was
    # missing on the platform, the marker was constant or the events were few
    why <- if (length(errors))
      paste(sprintf("%s: %s", names(errors), unlist(errors)), collapse = "; ") else "no reason recorded"
    stop("no dataset produced estimable results. Reasons: ", why, call. = FALSE)
  }
  if (length(per) == 1L)
    message("Only one dataset produced estimable results; the 'pooled' row is that dataset's estimate.")
  pc <- do.call(rbind, per)
  tok_used <- unique(pc$endpoint)
  if (length(tok_used) > 1L)
    warning("Cohorts contribute different endpoint tokens within family '", type, "': ",
            paste(sprintf("%s=%s", pc$dataset, pc$endpoint), collapse = ", "),
            ". Tokens of one family are pooled by design; the per-cohort token is kept in ",
            "per_dataset$endpoint and should be reported.", call. = FALSE)
  pooled <- meta_pool(pc$logHR, pc$se, p = pc$p, method = method)
  pooled$total_n <- sum(pc$n)
  pooled$total_events <- sum(pc$events)
  out <- list(input = list(datasets = datasets, marker = marker, type = type,
                           method = method, confounders = confounders,
                           time = Sys.time()),
              per_dataset = pc, pooled = pooled, errors = errors)
  class(out) <- "cpas_meta"
  out
}

meta_pool <- function(b, se, p = NULL, method = c("RE", "FE")) {
  method <- match.arg(method)
  k <- length(b)
  if (k == 0L) stop("meta_pool(): no dataset estimates were supplied.")
  w <- 1 / se^2
  bfe <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bfe)^2)
  df <- k - 1
  zi <- if (!is.null(p) && length(p) == k) sign(b) * stats::qnorm(1 - p / 2) else NULL
  pi_lower <- NA_real_; pi_upper <- NA_real_
  if (k == 1L) {
    # A single cohort is not a meta-analysis: report that cohort's own estimate
    # and leave the heterogeneity statistics undefined (NA, never NaN).
    return(data.frame(
      method = method, k = 1L,
      total_n = NA_integer_, total_events = NA_integer_,
      HR = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
      logHR = b, se = se, p = 2 * stats::pnorm(-abs(b / se)),
      Q = 0, df = 0L, p_heterogeneity = NA_real_, I2 = NA_real_, tau2 = NA_real_,
      pi_lower = NA_real_, pi_upper = NA_real_,
      z_stouffer = if (is.null(zi)) NA_real_ else sum(zi),
      p_stouffer = if (is.null(zi)) NA_real_ else 2 * stats::pnorm(-abs(sum(zi)))))
  }
  tau2 <- if (method == "RE") max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w))) else 0
  w2 <- 1 / (se^2 + tau2)
  bm <- sum(w2 * b) / sum(w2)
  seM <- sqrt(1 / sum(w2))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) else 0
  # prediction interval: where the true effect of a NEW cohort is expected to lie
  # (t distribution with k - 2 df; undefined for k < 3)
  if (k >= 3L) {
    tq <- stats::qt(0.975, df = k - 2)
    se_pi <- sqrt(seM^2 + tau2)
    pi_lower <- exp(bm - tq * se_pi); pi_upper <- exp(bm + tq * se_pi)
  }
  zs <- NA_real_; ps <- NA_real_
  if (!is.null(zi)) {
    zs <- sum(zi) / sqrt(k)
    ps <- 2 * stats::pnorm(-abs(zs))
  }
  data.frame(method = method, k = k, total_n = NA_integer_, total_events = NA_integer_,
             HR = exp(bm), lower = exp(bm - 1.96 * seM), upper = exp(bm + 1.96 * seM),
             logHR = bm, se = seM, p = 2 * stats::pnorm(-abs(bm / seM)),
             Q = Q, df = df, p_heterogeneity = stats::pchisq(Q, df, lower.tail = FALSE),
             I2 = I2, tau2 = tau2,
             pi_lower = pi_lower, pi_upper = pi_upper,
             z_stouffer = zs, p_stouffer = ps)
}

#' @title Print method for cpas_meta objects
#' @description Prints a concise summary of a \code{cpas_meta} analysis result.
#' @param x An object of class \code{cpas_meta}.
#' @param ... Unused, kept for S3 compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.cpas_meta <- function(x, ...) {
  cat("CanPAS integrative (meta) analysis\n")
  cat("marker:", x$input$marker, "| type:", x$input$type,
      "| method:", x$input$method,
      "| datasets:", nrow(x$per_dataset), "(failed:", length(x$errors), ")\n")
  print(x$pooled, row.names = FALSE)
  if (!is.null(x$pooled$pi_lower) && is.finite(x$pooled$pi_lower))
    cat(sprintf("95%% prediction interval for a new dataset: [%.3f, %.3f]\n",
                x$pooled$pi_lower, x$pooled$pi_upper))
  cat("\nPer-dataset:\n"); print(x$per_dataset[, c("dataset","n","events","HR","lower","upper","p")], row.names = FALSE)
  invisible(x)
}

#' @title 整合分析森林图
#' @description 绘制 \code{\link{cpas_meta}} 的逐队列 HR(95\% CI) 与合并 HR 森林图。
#' @param x 一个 \code{cpas_meta} 对象。
#' @param digits 图中标注保留的小数位(默认 4)。
#' @param show_stars 是否在面板左缘绘制逐队列显著性星号(\code{*} p<0.05、\code{**} p<0.01、\code{***} p<0.001)，默认 \code{TRUE}。
#' @param label_size 合并 HR / 95\% PI 标注的字号(默认 4.2)。
#' @param star_size 星号与 \code{Overall} 文字的字号，默认 \code{0.85 * label_size}。
#' @param y_headroom y 轴顶部留白，避免 Overall 菱形及其 CI 被面板裁切(默认 3.6)。
#' @param x_frac 面板左缘起的横向比例位置，用于锚定星号与合并标注(默认 0.02)。
#' @param label_lift 合并标注相对 Overall 行的纵向偏移(默认 1)。
#' @param label_where 合并标注位置：\code{"inside"} 置于面板左上角内，\code{"subtitle"} 置于面板上方(默认 \code{"inside"})。
#' @param ... 保留，暂未使用。
#' @return ggplot 对象。
#' @details
#' 自本版起，星号与合并标注锚定在由固定坐标轴留白反推得到的\strong{有限}数据值上
#' (不再用 \code{x = -Inf}：那会经 \code{log10} 变成 \code{NaN}，文本图层被静默丢弃并告警)。
#' 逐队列显著性星号位于面板内左缘，合并 HR / 95\% PI 标注默认置于面板左上角内。
#' @examples
#' \dontrun{
#'    m <- cpas_meta(c("GSE31210", "GSE37745"), marker = "TP53", type = "DFS")
#'    plot_meta_forest(m, digits = 4)
#' }
#' @export
plot_meta_forest <- function(x, digits = 4,
                             show_stars = TRUE,
                             label_size = 4.2,
                             star_size  = NULL,
                             y_headroom = 3.6,
                             x_frac = 0.02,
                             label_lift = 1.00,
                             label_where = c("inside", "subtitle"), ...) {
  label_where <- match.arg(label_where)
  pc <- x$per_dataset
  po <- x$pooled
  pc$se <- ifelse(pc$se <= 0, NA, pc$se)
  w <- ifelse(is.na(pc$se), 0, 1 / pc$se^2)
  pc$w <- w
  order_t <- pc$dataset
  pc$dataset <- factor(pc$dataset, levels = rev(order_t))
  k <- nrow(pc)
  y_overall <- k + 1                      # Overall 行所在的离散 y
  if (is.null(star_size)) star_size <- 0.85 * label_size
  PT <- 2.845                             # ggplot size -> pt
  ## 面板 x 边界（对数空间 5% 比例留白）：lo/hi 为面板左右缘对应的数据值
  lo_d <- min(pc$lower, na.rm = TRUE); hi_d <- max(pc$upper, na.rm = TRUE)
  Rlog <- log10(hi_d) - log10(lo_d)
  lo <- 10^(log10(lo_d) - 0.05 * Rlog); hi <- 10^(log10(hi_d) + 0.05 * Rlog)
  frac_x <- function(f) 10^(log10(lo) + f * (log10(hi) - log10(lo)))   # 面板内固定比例 -> 有限数据值
  x_anchor <- frac_x(x_frac)

  ## 合并标注文本（两行；内容与包内一致）
  pooled_lab <- if (is.finite(po$pi_lower))
    sprintf(paste0("Pooled HR %.", digits, "f [%.", digits, "f, %.", digits, "f] | I2=%.0f%%\n",
                   "95%% PI [%.", digits, "f, %.", digits, "f]"),
            po$HR, po$lower, po$upper, 100 * po$I2, po$pi_lower, po$pi_upper)
  else sprintf(paste0("Pooled HR %.", digits, "f [%.", digits, "f, %.", digits, "f], I2=%.0f%%"),
               po$HR, po$lower, po$upper, 100 * po$I2)

  ## ③ 星号文本（按行；x 由面板比例决定，与 HR 数值无关）
  st <- rep("", k)
  if (isTRUE(show_stars) && "p" %in% names(pc)) {
    pv <- suppressWarnings(as.numeric(pc$p))
    st <- ifelse(is.na(pv), "",
          ifelse(pv < 0.001, "***",
          ifelse(pv < 0.01,  "**",
          ifelse(pv < 0.05,  "*", ""))))
  }
  ## 行对齐：直接复用 pc$dataset 这个**同一个因子**作 y 美学，星号便与它自己的队列行
  ## 一同落在离散轴的同一位置上（v7 原用 y = seq_len(k) 的数值 1..k，而 levels = rev(order_t)
  ## 使数据第 i 行位于面板第 k-i+1 位，于是星号整列被上下镜像到错误队列上）。
  star_df <- data.frame(dataset = pc$dataset, st = st, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(pc, ggplot2::aes(x = HR, y = dataset)) +
    ggplot2::geom_point(size = 3 * sqrt(pc$w) / max(sqrt(pc$w)) + 1, color = "steelblue") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
                           orientation = "y", width = 0.25, color = "grey30") +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
    ## Overall 行（红菱形 + CI），靠顶部留白保证在框内
    ggplot2::geom_point(data = data.frame(HR = po$HR), ggplot2::aes(x = HR, y = y_overall),
                        shape = 18, size = 6, color = "firebrick", inherit.aes = FALSE) +
    ggplot2::geom_errorbar(data = data.frame(HR = po$HR, lower = po$lower, upper = po$upper),
                           ggplot2::aes(xmin = lower, xmax = upper, y = y_overall),
                           orientation = "y", width = 0.15, color = "firebrick", inherit.aes = FALSE) +
    ## x 轴只留右侧少量空白（左侧不再需要为文本让位）
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.6, y_headroom))) +
    ggplot2::labs(title = paste0("Meta forest: ", x$input$marker, " (", x$input$type, ", ",
                                 x$input$method, ")"),
                  x = "Hazard ratio (95% CI, log scale)", y = NULL) +
    ggplot2::theme_bw() +
    ## ③ 星号（x = 面板左缘起固定比例，y = 各队列行）
    ggplot2::geom_text(data = star_df, ggplot2::aes(x = x_anchor, y = dataset, label = st),
                       inherit.aes = FALSE, hjust = 0, size = star_size,
                       fontface = "bold", colour = "grey20") +
    ## Overall 文字（同一列）
    ggplot2::annotate("text", x = x_anchor, y = y_overall, label = "Overall",
                      hjust = 0, size = star_size, fontface = "bold", colour = "firebrick") +
    ## ② 合并标注：inside = 面板内（同一列、Overall 行之上）；subtitle = 面板上方
    { if (label_where == "inside")
        ggplot2::annotate("text", x = x_anchor, y = y_overall + label_lift, label = pooled_lab,
                          hjust = 0, vjust = 0, size = label_size, colour = "grey15") } +
    { if (label_where == "subtitle")
        ggplot2::labs(subtitle = pooled_lab) } +
    { if (label_where == "subtitle")
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = label_size,
                                                             hjust = 0, colour = "grey15",
                                                             margin = ggplot2::margin(b = 4))) }
  p
}

#' @title 敏感性分析：Leave-one-out
#' @description 对每个队列逐一剔除后重新计算合并估计，评估单一队列对整体结果的影响。
#' @param x 一个 \code{cpas_meta} 对象。
#' @return 数据框：\code{left_out} 为被剔除的队列，其余列为重新合并的汇总。
#' @examples
#' \dontrun{
#'    m <- cpas_meta(c("GSE31210", "GSE37745", "GSE42127"), marker = "TP53", type = "DFS")
#'    loo_meta(m)
#' }
#' @export
loo_meta <- function(x) {
  pc <- x$per_dataset
  if (nrow(pc) < 2L)
    stop("loo_meta(): leave-one-out sensitivity analysis needs at least 2 datasets.")
  out <- lapply(seq_len(nrow(pc)), function(i) {
    d <- pc[-i, , drop = FALSE]
    # with 2 cohorts the remaining estimate is that cohort's own result (k = 1)
    po <- meta_pool(d$logHR, d$se, p = d$p, method = x$input$method)
    data.frame(left_out = pc$dataset[i], po)
  })
  do.call(rbind, out)
}

#' @title 整合 KM（各队列内中位分组）
#' @description 每个队列 marker 中位切 High/Low，然后按 method 合并：
#'   "ipd"  直接合并患者做 pooled survfit（log-rank 同时给 pooled 与按队列分层两个 p）；
#'   "meta" 两阶段：对每个时点各队列 KM 的生存概率 S(t) 用 log(-log S) 转换后逆方差合并(RE/FE)，
#'           得到合并生存曲线与 1/3/5 年表；
#'   "both" 同时输出两者。
#' @param merged 命名列表：各队列 merge_surv_expr 的 merged data.frame；
#'   名称需为 catalog 中的 accession, 以便按家族解析具体终点。
#'   返回对象包含 \code{dataset_endpoints} (各数据集实际使用的终点)。
#' @param marker 基因列名。
#' @param type 终点（OS/RFS/...）。
#' @param method "ipd"/"meta"/"both"。
#' @param landmarks 输出时点(年)。
#' @param meta_method RE = 随机效应(DerSimonian-Laird,默认),FE = 固定效应。
#' @param cut 高/低分组规则,三种:\code{"median"}(默认)= 各队列内部取前 50%,
#'   即 High = marker 高于**该队列自身**的中位数;\code{"top_pct"}= 各队列内部按
#'   表达**从高到低排序取前 \code{top_pct}\%** 为 High(阈值 = 该队列的
#'   \eqn{100 - top_pct} 百分位,\code{top_pct = 25} 即前 25% 为高表达);
#'   \code{"custom"}= 用绝对阈值 \code{cut_value},大于阈值为 High、小于等于为 Low。
#'   绝对阈值在这里是有意义的,因为镜像提供的所有矩阵都在 log2 尺度上(论文 2.3 节),
#'   同一个数值在每个队列里含义相同;但按百分位分组才是"高表达 vs 低表达"的常规读法,
#'   应用界面只提供 median 与 top_pct 两种。这里刻意**不**提供"逐队列搜索最佳切点":
#'   那会把论文 3.4 节量化的切点搜索膨胀在每个队列各做一遍,合并后误差被叠加。
#' @param top_pct \code{cut = "top_pct"} 时的百分比:1-99 之间的单个有限数值,表达
#'   最高的 \code{top_pct}\% 患者记为 High。阈值用 \code{stats::quantile()} 的
#'   默认插值(type 7)计算;当多名患者的表达值恰好等于阈值时,分入 High 的比例可能
#'   略高于 \code{top_pct}\%,各队列实际阈值与人数都会返回(\code{cohort_thresholds}、
#'   \code{cutpoint})。若某队列因此一侧为空,该队列记入 \code{empty_cohorts} 并附
#'   原因,不参与合并。
#' @param cut_value \code{cut = "custom"} 时的阈值:单个有限数值,作用于 log2 表达
#'   (或签名得分)。若某队列在阈值一侧没有病人,该队列记入 \code{empty_cohorts}
#'   并附原因,不参与合并。
#' @return 普通 \code{list}(没有 class 属性;不是 S3 对象,没有对应的方法分派):
#'   \item{\code{df}:}{合并后的分析数据(time/status/marker/dataset/group)}
#'   \item{\code{datasets}, \code{n_high}, \code{n_low}, \code{method}, \code{cutpoint}:}{纳入的数据集与分组规模}
#'   \item{\code{dataset_endpoints}:}{各数据集实际使用的终点 token(命名向量)}
#'   \item{\code{cut}, \code{top_pct}, \code{cut_value}, \code{cutpoint}:}{实际使用的
#'     分组规则、其参数或阈值,以及可直接打印的规则说明;\code{cohort_thresholds}
#'     给出 \code{cut = "top_pct"} 时各队列的实际阈值}
#'   \item{\code{n_dropped}, \code{empty_cohorts}, \code{empty_reasons}:}{分组规则对纳入的
#'     队列是穷尽的(\code{n_dropped} 恒为 0);阈值/百分位规则在某队列里把一侧取空时,
#'     该队列记入 \code{empty_cohorts} 并附原因,不参与合并}
#'   \item{\code{skipped_cohorts}, \code{skipped_reasons}:}{因数据不足被排队的队列及其原因
#'     ——终点无法解析、缺少列、完整(time, status, marker)三元组少于 10 行,或 marker
#'     在该队列中只有一个取值。这些队列不参与合并,但**不会被静默丢弃**}
#'   \item{\code{fit}, \code{logrank_p}, \code{logrank_p_stratified}:}{method 为 ipd/both 时的 survfit 与 log-rank p(合并 / 按队列分层)}
#'   \item{\code{meta_landmarks}, \code{meta_curve}:}{method 为 meta/both 时的时点合并生存表与细网格曲线}
#' @details
#' 每队列内部按自身 marker 中位数分组(High = marker > median),因此比较的是"高于本队列中位"
#' 与"低于本队列中位"的风险,而不是跨队列的绝对阈值;需要绝对阈值时请用 \code{plot_km()} 逐队列分析。
#' IPD 合并同时给出未分层与按队列分层的 log-rank p,报告时建议给出后者。
#' 时点合并中,超出某队列最长随访的时点会沿用该队列最后一次观察到的 S(t)(\code{extend = TRUE}),
#' 每个时点实际贡献的队列数记录在 \code{meta_landmarks$k} 中。
#' @examples
#' \dontrun{
#'    ## Pooled Kaplan-Meier from an explicitly built list of merged cohorts
#'    ## (names must be catalog accessions so the endpoint is resolved per cohort).
#'    merged <- list(GSE31210 = cohort_merged("GSE31210", "GAPDH", type = "RFS"),
#'                    GSE37745 = cohort_merged("GSE37745", "GAPDH", type = "RFS"))
#'    km <- cpas_km_pooled(merged, marker = "GAPDH", type = "RFS", method = "both")
#'    km$datasets; km$dataset_endpoints; km$logrank_p
#' 
#'    plot_cpas_km(km)                       # pooled curve + landmark table
#'    plot_cpas_km_perdataset(km, ncol = 2)  # one panel per cohort
#' }
#' @export
cpas_km_pooled <- function(merged, marker, type = "OS",
                           method = c("both", "ipd", "meta"),
                           landmarks = c(1, 3, 5),
                           meta_method = "RE",
                           cut = c("median", "top_pct", "custom"),
                           top_pct = NULL,
                           cut_value = NULL) {
  method <- match.arg(method)
  cut <- match.arg(cut)
  empty <- list()                      # cohorts a threshold left one-sided
  skipped <- list()                    # cohorts dropped for insufficient/constant data
  med <- list()                        # per-cohort marker median (reported on failure)
  thr <- list()                        # per-cohort top-percent threshold
  if (identical(cut, "top_pct") &&
      (is.null(top_pct) || length(top_pct) != 1L || !is.finite(top_pct) ||
       top_pct <= 0 || top_pct >= 100))
    stop("cut = \"top_pct\" needs a single finite 'top_pct' between 1 and 99 ",
         "(the highest top_pct% of each cohort is High; 25 means the top 25%).",
         call. = FALSE)
  if (identical(cut, "custom") &&
      (is.null(cut_value) || length(cut_value) != 1L || !is.finite(cut_value)))
    stop("cut = \"custom\" needs a single finite 'cut_value' (a threshold on ",
         "the log2 expression scale, or on the signature score).", call. = FALSE)
  ep_used <- character(0)
  parts <- lapply(names(merged), function(t) {
    d <- merged[[t]]
    tok <- endpoint_resolve(t, type)
    if (is.na(tok)) {
      if (all(c(paste0(type, "_time"), paste0(type, "_status")) %in% colnames(d)))
        tok <- as.character(type)[1]
      else {
        skipped[[t]] <<- sprintf("no %s endpoint for this cohort and no %s_time/%s_status columns",
                                 type, type, type)
        return(NULL)
      }
    }
    tc <- paste0(tok, "_time"); sc <- paste0(tok, "_status")
    if (!all(c(tc, sc, marker) %in% colnames(d))) {
      skipped[[t]] <<- sprintf("column(s) absent: %s",
                               paste(setdiff(c(tc, sc, marker), colnames(d)), collapse = ", "))
      return(NULL)
    }
    ep_used[[t]] <<- tok
    st <- suppressWarnings(as.numeric(d[[sc]]))
    tm <- suppressWarnings(as.numeric(d[[tc]]))
    if (any(!is.na(st) & !st %in% c(0, 1)))
      stop("Cohort ", t, ": status column ", sc,
           " must be coded 0 (censored) / 1 (event).", call. = FALSE)
    if (any(!is.na(tm) & tm < 0))
      stop("Cohort ", t, ": negative survival time in ", tc, ".", call. = FALSE)
    v <- suppressWarnings(as.numeric(d[[marker]]))
    keep <- !is.na(v) & !is.na(suppressWarnings(as.numeric(d[[tc]]))) &
      !is.na(suppressWarnings(as.numeric(d[[sc]])))
    if (sum(keep) < 10) {
      skipped[[t]] <<- sprintf("only %d rows with a complete (time, status, marker) triple (minimum 10)",
                               sum(keep))
      return(NULL)
    }
    if (length(unique(v[keep])) < 2) {
      skipped[[t]] <<- sprintf("the marker takes a single value (%g) in this cohort, so no split is possible",
                               unique(v[keep])[1])
      return(NULL)
    }
    dd <- data.frame(time = suppressWarnings(as.numeric(d[[tc]][keep])),
                     status = suppressWarnings(as.numeric(d[[sc]][keep])),
                     marker = v[keep], dataset = t)
    ## Three rules:
    ##   "median"   the top 50% of THIS cohort (High = marker > cohort median);
    ##   "top_pct"  the top top_pct% of THIS cohort by expression, i.e. High =
    ##              marker > this cohort's (100 - top_pct)th percentile;
    ##   "custom"   an absolute threshold: High = marker > cut_value, Low otherwise.
    ## A custom absolute value is meaningful here because every matrix served by
    ## the mirror is on the log2 scale (Section 2.3 of the accompanying paper), so
    ## one threshold means the same thing in every cohort; it is also the rule a
    ## reader can reproduce without recomputing percentiles. A per-cohort search
    ## for the best cut point is deliberately not offered: that search inflates
    ## the p-value (quantified in Section 3.4) and pooling it over k cohorts
    ## compounds the inflation.
    med[[t]] <<- stats::median(dd$marker)
    if (identical(cut, "top_pct")) {
      thr <- stats::quantile(dd$marker, probs = 1 - top_pct / 100, names = FALSE)
      thr[[t]] <<- thr
      dd$group <- ifelse(dd$marker > thr, "High", "Low")
      if (length(unique(dd$group)) < 2L) {
        ## the top-percent rule can still leave one side empty when every value
        ## in the cohort is identical (or the ties span the whole cohort)
        empty[[t]] <<- sprintf("a top-%g%% split left no patient %s the cut (%g) (n = %d)",
                               top_pct,
                               if (all(dd$marker > thr)) "at or below" else "above",
                               thr, nrow(dd))
        return(NULL)
      }
    } else if (identical(cut, "custom")) {
      dd$group <- ifelse(dd$marker > cut_value, "High", "Low")
      ## one side missing entirely (min(table()) < 1 never fires: it means "no
      ## rows at all", not "no rows on that side")
      if (length(unique(dd$group)) < 2L) {
        ## the threshold leaves one side empty in THIS cohort: record it with the
        ## reason and leave the cohort out of the pool instead of failing later
        empty[[t]] <<- sprintf("no patient %s the threshold %g (n = %d)",
                               if (all(dd$marker > cut_value)) "at or below" else "above",
                               cut_value, nrow(dd))
        return(NULL)
      }
    } else {
      dd$group <- ifelse(dd$marker > stats::median(dd$marker), "High", "Low")
    }
    dd
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    if (length(empty)) {
      rule_desc <- switch(cut,
        top_pct = sprintf("A top-%g%% split", top_pct),
        custom  = sprintf("A custom cut of %g", cut_value),
        sprintf("Cut rule \"%s\"", cut))
      stop(rule_desc, " left every cohort one-sided (",
           paste(sprintf("%s: %s", names(empty), unlist(empty)), collapse = "; "),
           "). Pick a threshold inside the marker's observed range; the cohort ",
           "medians were ", paste(sprintf("%s = %.3g", names(med), unlist(med)),
                                  collapse = ", "), ".",
           if (length(skipped))
             paste0(" Cohorts skipped for insufficient or constant data: ",
                    paste(sprintf("%s (%s)", names(skipped), unlist(skipped)),
                          collapse = "; "), ".")
           else "", call. = FALSE)
    }
    stop("No cohort provided usable (time, status, marker) data for endpoint ",
         type, ".", call. = FALSE)
  }
  df <- do.call(rbind, parts)
  n_dropped <- sum(is.na(df$group))          # kept for callers: both rules are exhaustive
  df$group <- factor(df$group, levels = c("Low", "High"))
  df <- df[stats::complete.cases(df[c("time", "status", "marker", "group")]), , drop = FALSE]
  if (!sum(df$status == 1))
    stop("No event was observed in any cohort for endpoint ", type,
         "; a log-rank test / pooled KM cannot be estimated.", call. = FALSE)
  n_hi <- sum(df$group == "High"); n_lo <- sum(df$group == "Low")
  if (min(n_hi, n_lo) < 1L)
    stop("Median split produced an empty group.")
  cut_label <- switch(cut,
    median = "50% split within each cohort: High = marker > that cohort's median",
    top_pct = sprintf(paste0("top %g%% of each cohort by expression: ",
                             "High = marker > that cohort's %gth percentile"),
                      top_pct, 100 - top_pct),
    custom = sprintf(paste0("custom absolute threshold: High = marker > %g, ",
                            "Low = marker <= %g (log2 scale, the same value in every cohort)"),
                     cut_value, cut_value))
  out <- list(df = df, datasets = unique(df$dataset), method = method,
              n_high = n_hi, n_low = n_lo,
              dataset_endpoints = ep_used, cut = cut, top_pct = top_pct,
              cut_value = cut_value,
              n_dropped = n_dropped, cohort_medians = unlist(med),
              cohort_thresholds = unlist(thr),
              empty_cohorts = names(empty),
              empty_reasons = unlist(empty),
              skipped_cohorts = names(skipped),
              skipped_reasons = unlist(skipped), cutpoint = cut_label)
  if (method %in% c("ipd", "both")) {
    fit <- survival::survfit(survival::Surv(time, status) ~ group, data = df)
    lr  <- survival::survdiff(survival::Surv(time, status) ~ group, data = df)
    p   <- 1 - stats::pchisq(lr$chisq, length(lr$n) - 1)
    lrs <- tryCatch(survival::survdiff(survival::Surv(time, status) ~ group + strata(dataset), data = df),
                    error = function(e) NULL)
    p_strat <- if (!is.null(lrs)) 1 - stats::pchisq(lrs$chisq, length(lrs$n) - 1) else NA_real_
    out$fit <- fit; out$logrank_p <- p; out$logrank_p_stratified <- p_strat
  }
  if (method %in% c("meta", "both")) {
    ms <- km_surv_meta(df, times = landmarks, meta_method = meta_method)
    out$meta_landmarks <- ms$landmarks
    out$meta_curve <- ms$curve
  }
  invisible(out)
}

# 时点生存概率 meta：按 group 对每队列在 times 上的 S(t) 做 log(-log) 逆方差合并
km_surv_meta <- function(df, times, meta_method = "RE") {
  grp <- sort(unique(df$group))
  res_l <- lapply(grp, function(g) {
    dg <- df[df$group == g, ]
    cg <- unique(dg$dataset)
    fits <- lapply(cg, function(cf) {
      dd <- dg[dg$dataset == cf, ]
      s <- survival::survfit(survival::Surv(time, status) ~ 1, data = dd)
      sm <- summary(s, times = times, extend = TRUE)
      data.frame(dataset = cf, time = sm$time, S = sm$surv, SE = sm$std.err)
    })
    m <- do.call(rbind, fits); m <- m[is.finite(m$SE) & m$SE > 0 & m$S > 0 & m$S < 1, ]
    theta <- log(-log(m$S)); se_t <- m$SE / (m$S * abs(log(m$S)))
    tab <- lapply(unique(m$time), function(tm) {
      mm <- m[m$time == tm, ]
      if (nrow(mm) < 1) return(NULL)
      po <- meta_pool(theta[m$time == tm], se_t[m$time == tm], method = meta_method)
      data.frame(group = g, time = tm, k = nrow(mm),
                 S = exp(-exp(po$logHR)),
                 lower = exp(-exp(po$logHR + 1.96 * po$se)),
                 upper = exp(-exp(po$logHR - 1.96 * po$se)),
                 I2 = po$I2, p_het = po$p_heterogeneity)
    })
    do.call(rbind, tab)
  })
  land <- do.call(rbind, res_l)
  # 细网格曲线（每 0.5 年一步）
  grid_t <- seq(0.5, max(times) + 0.5, 0.5)
  res_c <- lapply(grp, function(g) {
    dg <- df[df$group == g, ]
    fits <- lapply(unique(dg$dataset), function(cf) {
      dd <- dg[dg$dataset == cf, ]
      s <- survival::survfit(survival::Surv(time, status) ~ 1, data = dd)
      sm <- summary(s, times = grid_t, extend = TRUE)
      data.frame(dataset = cf, time = sm$time, S = sm$surv, SE = sm$std.err)
    })
    m <- do.call(rbind, fits)
    m <- m[is.finite(m$SE) & m$SE > 0 & m$S > 0 & m$S < 1, ]
    tab <- lapply(unique(m$time), function(tm) {
      mm <- m[m$time == tm, ]; if (nrow(mm) < 1) return(NULL)
      po <- meta_pool(log(-log(mm$S)), mm$SE / (mm$S * abs(log(mm$S))), method = meta_method)
      data.frame(group = g, time = tm, S = exp(-exp(po$logHR)),
                 lower = exp(-exp(po$logHR + 1.96 * po$se)),
                 upper = exp(-exp(po$logHR - 1.96 * po$se)))
    })
    do.call(rbind, tab)
  })
  list(landmarks = land, curve = do.call(rbind, res_c))
}

#' @title 各队列 KM 小图网格（每个队列内中位分组的独立 KM）
#' @description 将 \code{\link{cpas_km_pooled}} 结果按数据集拆分为独立的 KM 小图并排布为网格。
#' @param km 一个 \code{cpas_km_pooled} 对象。
#' @param ncol 网格列数。
#' @param draw \code{TRUE}（默认）立即用 \code{gridExtra::grid.arrange} 画到当前设备；
#'   \code{FALSE} 时改为返回一个可组合的 \code{patchwork} 对象，供调用方与其它图
#'   拼版（把 \code{grid.arrange} 的返回值交给 patchwork 不会报错但会被静默丢弃，
#'   因此拼接场景请用 \code{draw = FALSE}）。
#' @return \code{draw = TRUE} 时不可见返回 \code{NULL}（图已画出）；
#'   \code{draw = FALSE} 时返回 \code{patchwork} 对象。
#' @export
plot_cpas_km_perdataset <- function(km, ncol = 4, draw = TRUE) {
  df <- km$df
  plist <- lapply(sort(unique(df$dataset)), function(cf) {
    d <- df[df$dataset == cf, ]
    f <- tryCatch(survival::survfit(survival::Surv(time, status) ~ group, data = d),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    sp <- survminer::ggsurvplot(f, data = d, pval = TRUE, legend = "none",
                                title = cf, xlab = "Time (years)", conf.int = FALSE,
                                palette = unname(.cpas_group_colours()))
    sp$plot
  })
  plist <- plist[!vapply(plist, is.null, logical(1))]
  if (isTRUE(draw)) return(do.call(gridExtra::grid.arrange, c(plist, list(ncol = ncol))))
  patchwork::wrap_plots(plist, ncol = ncol)
}

#' @title Plot pooled / meta Kaplan-Meier result
#' @description Visualizes a \code{\link{cpas_km_pooled}} result: when
#' method \code{"ipd"}/\code{"both"} was used a classical stratified
#' Kaplan-Meier plot is drawn; otherwise (method \code{"meta"}) the merged
#' meta survival curve with confidence band is plotted.
#' @param km An object returned by \code{\link{cpas_km_pooled}}.
#' @param which Which route to draw: \code{"auto"} (default) draws the
#'   individual-patient-data curve when that route was run and the time-point meta
#'   curve otherwise, exactly as before; \code{"ipd"} or \code{"meta"} selects a
#'   route explicitly. A \code{method = "both"} result holds both, so a caller that
#'   wants to show them together (as the Shiny application does) can draw each in
#'   turn instead of only ever getting the IPD curve.
#' @return A \code{ggsurvplot} object for the IPD route (with its risk table) or a
#' ggplot object for the meta route, carrying the requested landmark estimates and
#' their confidence intervals as points on the curve.
#' @export
plot_cpas_km <- function(km, which = c("auto", "ipd", "meta")) {
  which <- match.arg(which)
  has_ipd <- !is.null(km$fit)
  has_meta <- !is.null(km$meta_curve)
  if (which == "auto") which <- if (has_ipd) "ipd" else "meta"
  if (which == "ipd") {
    if (!has_ipd)
      stop("This result has no individual-patient-data route (method = \"meta\"). ",
           "Run with method = \"ipd\" or \"both\".", call. = FALSE)
    return(survminer::ggsurvplot(km$fit, data = km$df, pval = TRUE, risk.table = TRUE,
                                 legend.title = "marker",
                                 legend.labs = levels(km$df$group),
                                 palette = unname(.cpas_group_colours()),
                                 xlab = "Time (years)", conf.int = TRUE))
  }
  if (!has_meta)
    stop("This result has no time-point meta route (method = \"ipd\"). ",
         "Run with method = \"meta\" or \"both\".", call. = FALSE)
  ## the meta curve, with the requested landmark estimates drawn on it so the
  ## 1/3/5-year values that the landmark table reports are visible on the figure
  d <- km$meta_curve
  p <- ggplot(d, aes(x = time, y = S, color = group, fill = group)) +
    geom_step(linewidth = 0.7) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    scale_colour_manual(values = .cpas_group_colours()) +
    scale_fill_manual(values = .cpas_group_colours()) +
    scale_y_continuous(limits = c(0, 1)) + theme_bw() +
    labs(x = "Time (years)", y = "Survival probability",
         title = "Time-point meta curve (S(t) pooled across cohorts)")
  lm <- km$meta_landmarks
  if (!is.null(lm) && nrow(lm)) {
    lm <- lm[stats::complete.cases(lm[, c("time", "S", "lower", "upper")]), , drop = FALSE]
    if (nrow(lm)) {
      p <- p +
        geom_point(data = lm, aes(x = time, y = S, colour = group),
                   size = 2.1, inherit.aes = FALSE) +
        geom_errorbar(data = lm, aes(x = time, ymin = lower, ymax = upper, colour = group),
                      width = 0.12, linewidth = 0.45, inherit.aes = FALSE) +
        geom_text(data = lm, aes(x = time, y = lower, colour = group,
                                label = sprintf("%s\n%.3f", group, S)),
                  vjust = 1.5, size = 2.6, show.legend = FALSE, inherit.aes = FALSE)
    }
  }
  p
}