# endpoint_adjudication() / endpoint_agreement() (spec B1) --------------------
# All offline: the sheet comes from the shipped endpoint-semantics table and the
# agreement calculator is fed synthesised sheets.

test_that("the rating sheet is blinded and carries the required columns", {
  s <- endpoint_adjudication(n = 25, raters = 2, seed = 7)
  expect_s3_class(s, "cpas_adjudication_sheet")
  expect_equal(nrow(s), 25L)
  expect_identical(colnames(s)[1:4], c("record_id", "cohort", "token", "evidence"))
  # the rater must not see our verdict or the semantic fields
  forbidden <- c("Family", "PoolingClass", "EventDefinition", "TimeOrigin",
                 "CensoringRule", "CompetingEvents", "SourceField", "Derived",
                 "TokenRole", "Note")
  expect_false(any(forbidden %in% colnames(s)))
  # both raters get all five empty fields
  want <- as.vector(outer(c(1, 2), .cpas_adjudication_fields,
                          function(r, f) sprintf("rater%d_%s", r, f)))
  expect_true(all(want %in% colnames(s)))
  expect_true(all(is.na(as.matrix(s[, want]))))
  expect_setequal(attr(s, "pooling_class_vocabulary"),
                  c("Exact-equivalent", "Clinically-related", "Not-poolable",
                    "Unknown", "Absent"))
  expect_identical(attr(s, "judgement_fields"), .cpas_adjudication_fields)
})

test_that("record_id is stable across n and seed, and the draw is seeded", {
  a <- endpoint_adjudication(n = 20, seed = 11)
  b <- endpoint_adjudication(n = 20, seed = 11)
  expect_identical(a, b)                                   # same seed -> same sheet
  d <- endpoint_adjudication(n = 40, seed = 99)
  both <- intersect(a$record_id, d$record_id)
  expect_true(length(both) > 0)
  # an id always names the same cohort/token/evidence
  key <- function(x) paste(x$cohort, x$token, x$evidence, sep = "|")
  expect_identical(key(a)[match(both, a$record_id)],
                   key(d)[match(both, d$record_id)])
  # and the ids are assigned before sampling, so they are not 1..n
  expect_false(identical(sort(a$record_id), sprintf("REC%04d", 1:20)))
  # id order is preserved in the sheet
  expect_identical(a$record_id, sort(a$record_id))
})

test_that("the sample is stratified by family and source and respects the filters", {
  es <- endpoint_semantics()
  # record ids are assigned over the pooled table in catalog order
  ordering <- order(as.character(es$Accession), as.character(es$Family))
  es <- es[ordering, , drop = FALSE]
  s <- endpoint_adjudication(n = 60, seed = 3)
  rows <- es[match(s$record_id, sprintf("REC%04d", seq_len(nrow(es)))), , drop = FALSE]
  # one record per (cohort, family) and the token/evidence match the table
  expect_identical(rows$Accession, s$cohort)
  expect_identical(as.character(rows$Token), s$token)
  expect_identical(as.character(rows$Evidence), s$evidence)
  # the draw is spread over all five families and more than one source
  expect_setequal(unique(as.character(rows$Family)), c("OS", "DSS", "DFS", "PFS", "MFS"))
  src <- .cpas_accession_source(s$cohort)
  expect_gt(length(unique(src)), 1L)
  # the family filter and the absent filter are honoured
  os <- endpoint_adjudication(n = 15, families = "OS", seed = 5)
  os_rows <- es[match(os$record_id, sprintf("REC%04d", seq_len(nrow(es)))), , drop = FALSE]
  expect_true(all(as.character(os_rows$Family) == "OS"))
  dd <- endpoint_adjudication(n = 20, drop_absent = TRUE, seed = 5)
  expect_false(any(dd$token == "not stated"))
  # n is capped at the number of eligible records
  expect_equal(nrow(endpoint_adjudication(n = 1e6, seed = 1)), nrow(es))
  expect_error(endpoint_adjudication(n = 0), "'n'")
  expect_error(endpoint_adjudication(raters = 0), "'raters'")
  expect_error(endpoint_adjudication(n = 5, families = "NOPE"), "No record belongs")
})

test_that("the sheet can be written to CSV with empty rating cells", {
  f <- tempfile(fileext = ".csv")
  s <- endpoint_adjudication(n = 8, raters = 3, seed = 2, path = f)
  expect_true(file.exists(f))
  expect_identical(attr(s, "csv_path"), f)
  back <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(back), 8L)
  expect_identical(colnames(back), colnames(s))
  expect_true(all(back$rater3_pooling_class == "" | is.na(back$rater3_pooling_class)))
  expect_equal(sum(grepl("^rater[0-9]+_", colnames(s))), 15L)
})

test_that("Cohen's kappa matches the hand-computed 2x2 example exactly", {
  # rater1 A x25 / B x25 ; rater2 A x30 / B x20, with 35 agreeing pairs:
  #   po = 35/50 = 0.70, pe = (25*30 + 25*20)/50^2 = 0.50, kappa = 0.40
  r1 <- c(rep("A", 25), rep("B", 25))
  r2 <- c(rep("A", 20), rep("B", 5), rep("A", 10), rep("B", 15))
  k <- CanPAS:::.cpas_cohen_kappa(r1, r2)
  expect_equal(k$po, 0.7, tolerance = 1e-12)
  expect_equal(k$pe, 0.5, tolerance = 1e-12)
  expect_equal(k$kappa, 0.4, tolerance = 1e-12)
  expect_equal(k$n, 50L)
  expect_true(is.na(k$note))
  # a second hand-computed case: 3 categories,
  #   rows [3,1,0] / [1,3,0] / [0,0,2], n = 10, po = 0.8,
  #   pe = (4*4 + 4*4 + 2*2)/100 = 0.36, kappa = 0.44/0.64 = 0.6875
  a <- c(rep("X", 4), rep("Y", 4), rep("Z", 2))
  b <- c("X", "X", "X", "Y", "X", "Y", "Y", "Y", "Z", "Z")
  k2 <- CanPAS:::.cpas_cohen_kappa(a, b)
  expect_equal(k2$po, 0.8, tolerance = 1e-12)
  expect_equal(k2$pe, 0.36, tolerance = 1e-12)
  expect_equal(k2$kappa, 0.6875, tolerance = 1e-12)
  # one category only: chance agreement 1 leaves kappa undefined, not zero
  k1 <- CanPAS:::.cpas_cohen_kappa(rep("A", 5), rep("A", 5))
  expect_equal(k1$po, 1)
  expect_true(is.na(k1$kappa))
  expect_match(k1$note, "undefined")
  # no overlap at all is chance-level agreement, systematic disagreement is negative
  k0 <- CanPAS:::.cpas_cohen_kappa(rep("A", 4), rep("B", 4))
  expect_equal(k0$po, 0)
  expect_equal(k0$kappa, 0)
  kn <- CanPAS:::.cpas_cohen_kappa(c(rep("A", 4), rep("B", 4)),
                                   c(rep("B", 4), rep("A", 4)))
  expect_equal(kn$po, 0)
  expect_equal(kn$kappa, -1)
})

test_that("kappa interpretation labels follow Landis and Koch", {
  expect_identical(CanPAS:::.cpas_kappa_label(0.4), "fair")
  expect_identical(CanPAS:::.cpas_kappa_label(0.0), "slight")
  expect_identical(CanPAS:::.cpas_kappa_label(-0.1), "poor (less than chance)")
  expect_identical(CanPAS:::.cpas_kappa_label(0.75), "substantial")
  expect_identical(CanPAS:::.cpas_kappa_label(0.95), "almost perfect")
  expect_true(is.na(CanPAS:::.cpas_kappa_label(NA_real_)))
})

test_that("the agreement calculator scores two raters end to end", {
  mk_sheet <- function(v1, v2, field = "pooling_class") {
    s <- data.frame(record_id = sprintf("REC%04d", seq_along(v1)),
                    cohort = "SYN", token = "OS", evidence = "synthetic",
                    stringsAsFactors = FALSE)
    s[[paste0("rater1_", field)]] <- v1
    s[[paste0("rater2_", field)]] <- v2
    s
  }
  r1 <- c(rep("A", 25), rep("B", 25))
  r2 <- c(rep("A", 20), rep("B", 5), rep("A", 10), rep("B", 15))
  ag <- endpoint_agreement(mk_sheet(r1, r2))
  expect_s3_class(ag, "cpas_adjudication_agreement")
  expect_equal(ag$n_records, 50L)
  expect_equal(ag$n_raters, 2L)
  expect_equal(ag$per_field$n_compared, 50L)
  expect_equal(ag$per_field$n_agree, 35L)
  expect_equal(ag$per_field$raw_agreement, 0.7, tolerance = 1e-12)
  expect_equal(ag$per_field$kappa, 0.4, tolerance = 1e-12)
  expect_identical(ag$per_field$kappa_interpretation, "fair")
  expect_equal(ag$overall_raw_agreement, 0.7, tolerance = 1e-12)
  expect_equal(ag$pooled_kappa, 0.4, tolerance = 1e-12)
  # 15 of 50 records disagree
  expect_equal(ag$adjudication_rate, 0.3, tolerance = 1e-12)
  expect_equal(nrow(ag$disagreements), 15L)
  expect_true(all(c("record_id", "cohort", "token", "evidence",
                    "fields_disagreed", "detail") %in% colnames(ag$disagreements)))
  expect_true(all(grepl("pooling_class: rater1=", ag$disagreements$detail)))
  expect_output(print(ag), "endpoint adjudication")
})

test_that("missing ratings reduce the compared count instead of biasing kappa", {
  v1 <- c("A", "A", "B", "B", NA, "", "A")
  v2 <- c("A", "B", "B", "B", "A", "B", NA)
  s <- data.frame(record_id = sprintf("REC%04d", 1:7), cohort = "SYN",
                  token = "OS", evidence = "synthetic",
                  rater1_pooling_class = v1, rater2_pooling_class = v2,
                  stringsAsFactors = FALSE)
  ag <- endpoint_agreement(s, fields = "pooling_class")
  # only rows 1-4 have both ratings present and non-empty
  expect_equal(ag$per_field$n_compared, 4L)
  expect_equal(ag$per_field$n_agree, 3L)
  expect_equal(ag$per_field$raw_agreement, 0.75, tolerance = 1e-12)
  expect_equal(ag$per_field$kappa, CanPAS:::.cpas_cohen_kappa(v1, v2)$kappa,
               tolerance = 1e-12)
  expect_equal(ag$n_comparable_records, 4L)
})

test_that("more than two raters report the pairwise mean kappa and say so", {
  set.seed(4)
  n <- 30
  s <- data.frame(record_id = sprintf("REC%04d", seq_len(n)), cohort = "SYN",
                  token = "OS", evidence = "synthetic", stringsAsFactors = FALSE)
  s$rater1_pooling_class <- sample(c("A", "B", "C"), n, TRUE)
  s$rater2_pooling_class <- s$rater1_pooling_class
  s$rater3_pooling_class <- sample(c("A", "B", "C"), n, TRUE)
  ag <- endpoint_agreement(s, fields = "pooling_class")
  expect_equal(ag$n_raters, 3L)
  expect_equal(nrow(ag$pairwise), 3L)                       # 3 rater pairs
  expect_match(ag$kappa_basis, "mean of its 3 pairwise kappas")
  expect_match(ag$per_field$kappa_basis, "mean of the 3 pairwise")
  # rater1 and rater2 agree perfectly; the mean is dragged down by rater3
  k12 <- ag$pairwise$kappa[ag$pairwise$pair == "1-2"]
  expect_equal(k12, 1, tolerance = 1e-12)
  expect_equal(ag$per_field$kappa,
               mean(ag$pairwise$kappa[is.finite(ag$pairwise$kappa)]), tolerance = 1e-12)
  expect_lt(ag$per_field$kappa, 1)
})

test_that("a single-category field reports kappa NA with a note, never NaN", {
  s <- data.frame(record_id = sprintf("REC%04d", 1:6), cohort = "SYN",
                  token = "OS", evidence = "synthetic",
                  rater1_pooling_class = rep("A", 6),
                  rater2_pooling_class = rep("A", 6),
                  stringsAsFactors = FALSE)
  ag <- endpoint_agreement(s, fields = "pooling_class")
  expect_equal(ag$per_field$raw_agreement, 1)
  expect_true(is.na(ag$per_field$kappa))
  expect_false(is.nan(ag$per_field$kappa))
  expect_match(ag$per_field$note, "undefined")
  expect_true(any(grepl("undefined", ag$notes)))
})

test_that("three judgement fields are scored together and pooled", {
  n <- 12
  s <- data.frame(record_id = sprintf("REC%04d", seq_len(n)), cohort = "SYN",
                  token = "OS", evidence = "synthetic", stringsAsFactors = FALSE)
  for (r in 1:2) for (f in c("event_definition", "time_origin")) {
    s[[sprintf("rater%d_%s", r, f)]] <- rep(c("x", "y"), length.out = n)
  }
  s$rater1_event_definition[1] <- "z"    # one disagreement
  ag <- endpoint_agreement(s, fields = c("event_definition", "time_origin"))
  expect_equal(nrow(ag$per_field), 2L)
  expect_equal(ag$per_field$n_compared, c(12L, 12L))
  expect_equal(ag$per_field$raw_agreement[1], 11 / 12, tolerance = 1e-12)
  expect_equal(ag$per_field$raw_agreement[2], 1, tolerance = 1e-12)
  expect_equal(ag$overall_raw_agreement, (11 + 12) / 24, tolerance = 1e-12)
  expect_equal(nrow(ag$disagreements), 1L)
  expect_identical(ag$disagreements$fields_disagreed, "event_definition")
})

test_that("endpoint_agreement validates its input", {
  expect_error(endpoint_agreement("nope"), "data frame")
  expect_error(endpoint_agreement(data.frame(a = 1)), "record_id")
  s <- data.frame(record_id = "REC0001", rater1_pooling_class = "A")
  expect_error(endpoint_agreement(s), "at least two raters")
  # a field that was not prepared at all
  s2 <- data.frame(record_id = "REC0001", rater1_event_definition = "a",
                   rater2_event_definition = "b")
  expect_error(endpoint_agreement(s2, fields = "time_origin"), "No rater column found")
  # a requested field missing the rating columns of one of the raters is skipped
  s3 <- data.frame(record_id = "REC0001", rater1_event_definition = "a",
                   rater2_time_origin = "b", stringsAsFactors = FALSE)
  expect_error(endpoint_agreement(s3, fields = c("event_definition", "time_origin")),
               "None of the requested")
  # a single comparable-but-unusable pair is a field with no usable pair
  s4 <- data.frame(record_id = "REC0001", rater1_time_origin = NA,
                   rater2_time_origin = "x")
  ag <- endpoint_agreement(s4, fields = "time_origin")
  expect_equal(ag$per_field$n_compared, 0L)
  expect_true(is.na(ag$per_field$raw_agreement))
  expect_true(any(grepl("no record was rated", ag$notes)))
})
