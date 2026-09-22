test_that("plot_km returns a ggsurvplot", {
  d <- mk_cohort(200, seed = 20)
  kp <- plot_km(d, type = "OS", marker = "GAPDH")
  expect_true(inherits(kp, "ggsurvplot"))
})

test_that("plot_km custom cutpoint creates correct group sizes", {
  d <- mk_cohort(120, seed = 22)
  kp <- plot_km(d, type = "OS", marker = "GAPDH", cutpoint = 13)
  ds <- kp$data.survplot
  at_risk <- ds[!duplicated(ds$group), c("group", "n.risk")]
  names(at_risk)[2] <- "n0"
  expect_equal(sort(at_risk$n0),
               sort(as.integer(table(ifelse(d$GAPDH > 13, "High", "Low")))))
})

test_that("plot_roc returns ggplot with AUC attribute", {
  d <- mk_cohort(250, seed = 21)
  g <- plot_roc(d, type = "OS", marker = "GAPDH", predict.time = 2)
  expect_s3_class(g, "ggplot")
  expect_true(is.finite(attr(g, "AUC")))
})

# --- forest_plot(): HR order must never break a label from its estimate -----
# A categorical covariate contributes one row per level, and a bare level name
# ("G3") is only meaningful next to its covariate name. The old implementation
# sorted by HR and then re-attached the covariate headers by moving whole
# position ranges, which could place a level under another covariate's header.

fp_table <- function(n = 300, seed = 4) {
  set.seed(seed)
  mk <- stats::rnorm(n, 10, 2); age <- stats::rnorm(n, 62, 9)
  grade <- factor(sample(c("G1", "G2", "G3"), n, TRUE))
  node <- factor(sample(c("N0", "N1", "N2"), n, TRUE))
  sex <- factor(sample(c("female", "male"), n, TRUE))
  hr <- exp(0.30 * (mk - 10) + 0.02 * (age - 62) - 0.8 * (grade == "G3") +
              0.5 * (node == "N2") + 0.1 * (sex == "male"))
  tt <- stats::rexp(n, hr / 10); st <- stats::rbinom(n, 1, 0.7)
  d <- data.frame(ID = sprintf("S%03d", seq_len(n)), OS_time = tt, OS_status = st,
                  marker = mk, age = age, grade = grade, N = node, sex = sex)
  COX_analysis(d, type = "OS", cont_Variates = c("marker", "age"),
               cate_Variates = c("grade", "N", "sex"), method = "uni")$results_table
}

test_that("sorted rows are monotone in HR and every label names its covariate", {
  tb <- fp_table()
  p <- forest_plot(tb)                    # default: HR_order = "decrease"
  dd <- p$data[order(p$data$row), , drop = FALSE]
  expect_true(all(diff(dd$HR) <= 0))                        # strictly sorted
  expect_false(any(dd$Type == "cate_header"))               # headers replaced by labels
  lvl <- dd[dd$Type == "cate_level", , drop = FALSE]
  expect_gt(nrow(lvl), 0L)
  expect_true(all(startsWith(lvl$label, paste0(lvl$Variates, ": "))))
  # no estimate is lost: every estimable row of the input is plotted exactly once
  est <- tb[!is.na(tb$HR), , drop = FALSE]
  expect_equal(nrow(dd), nrow(est))
  expect_setequal(dd$HR, est$HR)

  up <- forest_plot(tb, HR_order = "increase")$data
  expect_true(all(diff(up$HR[order(up$row)]) >= 0))
})

test_that("group_levels = 'block' keeps each header directly above its own levels", {
  tb <- fp_table()
  dd <- forest_plot(tb, group_levels = "block")$data
  dd <- dd[order(dd$row), , drop = FALSE]
  # every categorical row is preceded by its own covariate header or by another
  # level of the same covariate - never by another covariate's row
  for (i in seq_len(nrow(dd))) {
    if (dd$Type[i] == "cate_level") {
      above <- dd$Var1[seq_len(i - 1L)]
      expect_identical(tail(above, 1L), dd$Var1[i], info = dd$label[i])
    }
  }
  hdr <- which(dd$Type == "cate_header")
  expect_gt(length(hdr), 0L)
  # a header is immediately followed by a level of the same covariate
  for (h in hdr) expect_identical(dd$Var1[h + 1L], dd$Var1[h])
})

test_that("HR_order = 'none' reproduces the input order with its header rows", {
  tb <- fp_table()
  dd <- forest_plot(tb, HR_order = "none")$data
  dd <- dd[order(dd$row), , drop = FALSE]
  keep <- !is.na(tb$HR) | tb$Type == "cate_header"
  expect_identical(dd$Var1, tb$Var1[keep])
  expect_identical(dd$Level, tb$Level[keep])
  expect_identical(dd$Type, tb$Type[keep])
})

test_that("a table without header rows is sorted with plain labels", {
  tb <- fp_table()
  cont <- tb[tb$Type == "cont", , drop = FALSE]
  dd <- forest_plot(cont)$data
  dd <- dd[order(dd$row), , drop = FALSE]
  expect_equal(nrow(dd), nrow(cont))
  expect_true(all(diff(dd$HR) <= 0))
  expect_identical(dd$label, cont$Variates[order(cont$HR, decreasing = TRUE)])
})

# --- KM group colours: high expression red, low expression green ------------
# survminer's default hue palette assigns its first colour to the first factor
# level, and plot_km() levels the groups as c("Low", "High"), so the default
# drew LOW in the warm colour. The convention is fixed in .cpas_group_colours()
# and must survive both an explicit palette= and an options() override.

km_palette_of <- function(p) {
  sc <- p$plot$scales$get_scales("colour")
  if (is.null(sc)) return(NULL)
  sc$palette(2)
}

test_that("plot_km draws high expression red and low expression green", {
  d <- mk_cohort(150, seed = 31)
  kp <- plot_km(d, type = "OS", marker = "GAPDH")
  expect_identical(levels(kp$data.survplot$group), c("Low", "High"))
  expect_identical(unname(km_palette_of(kp)), c("#2CA02C", "#D62728"))
  ## the drawn layers must carry the same two colours
  cols <- unique(ggplot2::ggplot_build(kp$plot)$data[[1]]$colour)
  expect_setequal(cols, c("#2CA02C", "#D62728"))
})

test_that("an explicit palette still wins, and options() moves the default", {
  d <- mk_cohort(150, seed = 32)
  kp <- plot_km(d, type = "OS", marker = "GAPDH", palette = c("blue", "orange"))
  expect_identical(unname(km_palette_of(kp)), c("blue", "orange"))

  old <- options(CanPAS.km_palette = c("#004400", "#880000"))
  on.exit(options(old), add = TRUE)
  kp2 <- plot_km(d, type = "OS", marker = "GAPDH")
  expect_identical(unname(km_palette_of(kp2)), c("#004400", "#880000"))
  expect_error(options(CanPAS.km_palette = "red"),
               NA)                       # a single colour is only caught on use
  options(CanPAS.km_palette = "red")
  expect_error(plot_km(d, type = "OS", marker = "GAPDH"), "exactly two colours")
})

# --- pooled KM cut rule: 50% split or one absolute threshold ----------------
# The pooled page offers exactly two rules. "50 (top 50% high)" splits inside
# each cohort at its own median; "Custom" applies one absolute value in every
# cohort (meaningful because every mirrored matrix is log2). A per-cohort search
# for the best cut point is deliberately not offered: it would apply the
# cut-point inflation k times and then pool it.

pool_two <- function(n = 160, seed = 11) {
  set.seed(seed)
  mk <- function(pref, mu = 9) {
    data.frame(ID = paste0(pref, seq_len(n)),
               OS_time = round(rexp(n, 0.2) + 0.1, 3),
               OS_status = rbinom(n, 1, 0.6),
               FOXM1 = rnorm(n, mu, 1))
  }
  list(A = mk("a"), B = mk("b"))
}

test_that("the default cut is the within-cohort 50% split", {
  m <- pool_two()
  k <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "both")
  expect_identical(k$cut, "median")
  expect_match(k$cutpoint, "50% split within each cohort")
  expect_null(k$cut_value)
  expect_length(k$empty_cohorts, 0)
  tab <- table(k$df$dataset, k$df$group)
  expect_true(all(tab[, "High"] == tab[, "Low"]))
})

test_that("a custom threshold puts everyone above it in High and the rest in Low", {
  m <- pool_two(n = 200)
  thr <- 9.5
  k <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "ipd",
                      cut = "custom", cut_value = thr)
  expect_identical(k$cut, "custom")
  expect_identical(k$cut_value, thr)
  expect_match(k$cutpoint, "custom absolute threshold")
  expect_true(all(k$df$marker[k$df$group == "High"] > thr))
  expect_true(all(k$df$marker[k$df$group == "Low"] <= thr))
  ## the same absolute threshold is applied in both cohorts, so group sizes
  ## follow each cohort's own distribution rather than being forced equal
  per <- table(k$df$dataset, k$df$group)
  expect_identical(sum(per[, "High"]), sum(m$A$FOXM1 > thr) + sum(m$B$FOXM1 > thr))
})

test_that("a threshold that empties one side of a cohort drops that cohort, loudly", {
  m <- pool_two(n = 120)
  ## far above every value: nothing is High anywhere -> no estimable pool, and
  ## the message must say why and offer the observed medians as a guide
  err <- tryCatch(cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "ipd",
                                 cut = "custom", cut_value = 1e6),
                  error = function(e) conditionMessage(e))
  expect_match(err, "left every cohort one-sided")
  expect_match(err, "cohort medians were")
  expect_match(err, "GSE|A =|B =")
  ## a threshold above every value in one cohort but inside the other: that
  ## cohort is reported instead of being silently dropped
  m2 <- pool_two(n = 120)
  m2$B$FOXM1 <- m2$B$FOXM1 + 6            # cohort B sits far higher
  k <- cpas_km_pooled(m2, marker = "FOXM1", type = "OS", method = "ipd",
                      cut = "custom", cut_value = 10)
  expect_length(k$empty_cohorts, 1)
  expect_true(k$empty_cohorts %in% c("A", "B"))
  expect_true(nzchar(k$empty_reasons[[1]]))
  expect_identical(sort(unique(k$df$dataset)), setdiff(c("A", "B"), k$empty_cohorts))
})

test_that("a top-percent cut puts the highest x% of each cohort in High", {
  m <- pool_two(n = 200)
  k <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "ipd",
                      cut = "top_pct", top_pct = 25)
  expect_identical(k$cut, "top_pct")
  expect_identical(k$top_pct, 25)
  expect_null(k$cut_value)
  expect_match(k$cutpoint, "top 25% of each cohort by expression")
  expect_match(k$cutpoint, "75th percentile")
  expect_length(k$empty_cohorts, 0)
  ## the rule is applied inside each cohort: about a quarter of each one is High
  per <- table(k$df$dataset, k$df$group)
  expect_true(all(per[, "High"] == 50))
  expect_true(all(per[, "Low"] == 150))
  ## the threshold recorded per cohort is the (100 - x)th percentile of that cohort
  expect_equal(unname(k$cohort_thresholds[["A"]]),
               unname(stats::quantile(m$A$FOXM1, 0.75)), tolerance = 1e-12)
  for (cf in c("A", "B")) {
    hi <- k$df$dataset == cf & k$df$group == "High"
    expect_true(all(k$df$marker[hi] > k$cohort_thresholds[[cf]]))
  }
  ## and it is applied per cohort, not once on the pooled marker
  expect_false(isTRUE(all.equal(k$cohort_thresholds[["A"]],
                                k$cohort_thresholds[["B"]])))
})

test_that("top_pct = 50 reproduces the median split exactly", {
  m <- pool_two(n = 201)                    # odd n, so the two rules are comparable
  a <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "ipd", cut = "median")
  b <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "ipd",
                      cut = "top_pct", top_pct = 50)
  expect_identical(a$df$group, b$df$group)
  expect_equal(a$n_high, b$n_high)
  expect_equal(a$n_low, b$n_low)
})

test_that("a cohort that cannot be split is reported as skipped, not dropped silently", {
  m <- pool_two(n = 120)
  m$B$FOXM1 <- 7                            # every value identical in B
  k <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "ipd",
                      cut = "top_pct", top_pct = 25)
  expect_length(k$empty_cohorts, 0)
  expect_identical(k$skipped_cohorts, "B")
  expect_match(k$skipped_reasons[[1]], "single value .* no split is possible")
  expect_identical(sort(unique(k$df$dataset)), "A")
  ## the same disclosure applies to the median rule
  k2 <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "ipd")
  expect_identical(k2$skipped_cohorts, "B")
})

test_that("a top-percent cut outside 1-99 is refused", {
  m <- pool_two()
  expect_error(cpas_km_pooled(m, marker = "FOXM1", type = "OS", cut = "top_pct"),
               "needs a single finite 'top_pct' between 1 and 99")
  for (bad in list(0, 100, -5, NA_real_, c(20, 30)))
    expect_error(cpas_km_pooled(m, marker = "FOXM1", type = "OS",
                                cut = "top_pct", top_pct = bad),
                 "needs a single finite 'top_pct' between 1 and 99")
})

test_that("a custom cut without a value is refused", {
  m <- pool_two()
  expect_error(cpas_km_pooled(m, marker = "FOXM1", type = "OS", cut = "custom"),
               "needs a single finite 'cut_value'")
  expect_error(cpas_km_pooled(m, marker = "FOXM1", type = "OS", cut = "custom",
                              cut_value = NA_real_),
               "needs a single finite 'cut_value'")
})

test_that("the cut rule and the landmark table travel with the result", {
  m <- pool_two()
  k <- cpas_km_pooled(m, marker = "FOXM1", type = "OS", method = "meta",
                      cut = "custom", cut_value = 9)
  expect_true(is.character(k$cutpoint) && nzchar(k$cutpoint))
  expect_true(!is.null(k$meta_landmarks))
  expect_true(all(c("group", "time", "k", "S", "lower", "upper") %in%
                    colnames(k$meta_landmarks)))
})

test_that("plot_cpas_km can draw each route separately, and 'both' holds both", {
  m <- list(); set.seed(3)
  for (nm in c("A", "B")) {
    n <- 120
    m[[nm]] <- data.frame(ID = paste0(nm, seq_len(n)), OS_time = round(rexp(n, .2) + .1, 3),
                          OS_status = rbinom(n, 1, .6), GAPDH = rnorm(n, 12, 1))
  }
  kb <- cpas_km_pooled(m, marker = "GAPDH", type = "OS", method = "both")
  ## 'both' really does hold two routes (this is what the app now draws)
  expect_false(is.null(kb$fit))
  expect_false(is.null(kb$meta_curve))
  ## the default keeps its documented behaviour: IPD when present
  expect_s3_class(plot_cpas_km(kb), "ggsurvplot")
  ## each route can be asked for explicitly
  expect_s3_class(plot_cpas_km(kb, which = "ipd"), "ggsurvplot")
  pm <- plot_cpas_km(kb, which = "meta")
  expect_s3_class(pm, "ggplot")
  ## the meta curve carries the landmark estimates (points + CI) on the figure
  pts <- pm$layers[[which(vapply(pm$layers, function(l) inherits(l$geom, "GeomPoint"),
                                 logical(1)))[1]]]
  expect_true(!is.null(pts$data) && nrow(pts$data) == nrow(kb$meta_landmarks))
  ## asking for a route that was not run fails with a clear message
  ki <- cpas_km_pooled(m, marker = "GAPDH", type = "OS", method = "ipd")
  expect_error(plot_cpas_km(ki, which = "meta"), "no time-point meta route")
  km <- cpas_km_pooled(m, marker = "GAPDH", type = "OS", method = "meta")
  expect_error(plot_cpas_km(km, which = "ipd"), "no individual-patient-data route")
})

test_that("the pooled figure block composes both routes when both were run", {
  m <- list(); set.seed(4)
  for (nm in c("A", "B")) {
    n <- 120
    m[[nm]] <- data.frame(ID = paste0(nm, seq_len(n)), OS_time = round(rexp(n, .2) + .1, 3),
                          OS_status = rbinom(n, 1, .6), GAPDH = rnorm(n, 12, 1))
  }
  kb <- cpas_km_pooled(m, marker = "GAPDH", type = "OS", method = "both")
  pi <- plot_cpas_km(kb, which = "ipd"); pm <- plot_cpas_km(kb, which = "meta")
  skip_if_not_installed("patchwork")
  ## exactly as the app composes it, with wrap_plots and no patchwork operators
  ## (the app only loads patchwork's namespace, so operator dispatch is not
  ## guaranteed there)
  panels <- list(pi$plot + ggplot2::ggtitle("(A) IPD"),
                 pi$table + ggplot2::theme(plot.title = ggplot2::element_text(size = 9)),
                 pm + ggplot2::ggtitle("(B) meta"))
  fig <- patchwork::wrap_plots(panels, ncol = 1, heights = c(3, 1, 2))
  expect_s3_class(fig, "patchwork")
  fp <- tempfile(fileext = ".png")
  grDevices::png(fp, width = 680, height = 780)
  ok <- tryCatch({ print(fig); TRUE }, error = function(e) FALSE)
  grDevices::dev.off()
  expect_true(ok)
  expect_gt(file.info(fp)$size, 10000)      # a blank device would be far smaller
})
