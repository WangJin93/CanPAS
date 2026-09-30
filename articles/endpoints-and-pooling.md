# Endpoints, pooling classes and shared-patient overlap

``` r
library(CanPAS)
```

## Why endpoints need a semantics layer

Different cohorts call the same clinical question by different names.
“How long until the disease comes back?” may be deposited as `DFS`,
`RFS`, `EFS` or `DFI`; “how long until progression?” as `PFS` or the
TCGA-derived `PFI`. CanPAS groups those tokens into five **pooling
families** (`OS`, `DSS`, `DFS`, `PFS`, `MFS`), resolves the token each
cohort actually provides, and records what that token means. Nothing is
pooled across families.

``` r
endpoint_family("RFS")     # a token -> its family
#> [1] "DFS"
endpoint_family("DRFS")
#> [1] "MFS"
endpoint_options("GSE31210")
#>   family token derived     label
#> 1     OS    OS   FALSE        OS
#> 2    DFS   RFS   FALSE DFS (RFS)
```

The registry behind that is the **endpoint semantics table**: one row
per cohort and family, with the token, the role of the token, what the
deposit states about the event definition / time origin / censoring /
competing events, the artefact the entry was read from, and - the column
the analysis uses - the **pooling class**.

``` r
es <- endpoint_semantics(family = "DFS")
es[es$Accession %in% c("GSE31210", "GSE13507", "TCGA-LUAD"),
   c("Accession", "Family", "Token", "TokenRole", "PoolingClass", "Evidence")]
#>     Accession Family      Token    TokenRole       PoolingClass
#> 3    GSE31210    DFS        RFS contributing Clinically-related
#> 81   GSE13507    DFS not stated   not stated             Absent
#> 130 TCGA-LUAD    DFS        DFI contributing Clinically-related
#>                                                                                                                                                         Evidence
#> 3                                                                    data/dataset_info.csv; data/processed/surv/GSE31210_surv.rds (columns RFS_status, RFS_time)
#> 81                                                                                                                                         data/dataset_info.csv
#> 130 data/dataset_info.csv; data/processed/surv/TCGA-LUAD_surv.rds (columns DFI_status, DFI_time); data/tcga/tcga_surv.rda (via pipeline/R/35_add_tcga_cohorts.R)
```

``` r
table(endpoint_semantics()$PoolingClass)
#> 
#>             Absent Clinically-related   Exact-equivalent 
#>                638                111                236
```

The vocabulary is frozen and has five values:

| class                | meaning                                                                                                                        |
|----------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `Exact-equivalent`   | the token **is** the family’s canonical definition (OS in OS, DFS in DFS, …)                                                   |
| `Clinically-related` | a different token pooled into the family by the documented rule (RFS/EFS/DFI in DFS, CSS/BCSS in DSS, PFI in PFS, DRFS in MFS) |
| `Not-poolable`       | the cohort reports an endpoint in this family that cannot be pooled                                                            |
| `Unknown`            | the deposit does not state what the endpoint means                                                                             |
| `Absent`             | the cohort has no endpoint of this family at all - the absence of an endpoint, **not** a pooling verdict                       |

``` r
endpoint_pooling_class("GSE31210", "DFS")        # RFS: related, not identical
#> [1] "Clinically-related"
endpoint_pooling_class("GSE4922_GPL96", "DFS")   # DFS itself: exact
#> [1] "Exact-equivalent"
endpoint_pooling_class("GSE13507", "DFS")        # no DFS endpoint at all -> Absent
#> [1] "Absent"
```

[`endpoint_pooling_class()`](https://wangjin93.github.io/CanPAS/reference/endpoint_pooling_class.md)
is one cell of
[`endpoint_semantics()`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md).
If the companion table is not shipped in a build, both fall back to the
package’s documented token -\> family map, which is the same rule the
table records for the exact/related split; a family the catalog does not
annotate still comes back as `Absent`.

## Two pooling modes

[`cpas_meta()`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)
has two modes and `"family"` is the documented default:

- `pooling = "family"` (default) pools `Exact-equivalent` **and**
  `Clinically-related` rows of one family. This is what makes a
  cross-token claim possible at all: the DFS pool is 98 cohorts, of
  which only about 20 report a token literally called DFS.
- `pooling = "exact"` pools only `Exact-equivalent` rows. It is
  available for a reviewer who wants the strictest reading, and it names
  every cohort it leaves out.

Both modes record the classes they pooled, per cohort and on the pooled
row.

``` r
mk_cohort <- function(n, prefix, seed, token = "OS", role = "primary") {
  set.seed(seed)
  x <- rnorm(n, 10, 1)
  hr <- exp(0.4 * (x - 10))
  d <- data.frame(ID = sprintf("%s%03d", prefix, seq_len(n)),
                  OS_time = round(rexp(n, hr / 9) + 0.05, 3),
                  OS_status = rbinom(n, 1, 0.65),
                  GAPDH = x, TP53 = rnorm(n, 9, 1.1))
  if (token != "OS") {
    d[[paste0(token, "_time")]] <- d$OS_time
    d[[paste0(token, "_status")]] <- d$OS_status
  }
  d
}
A <- mk_cohort(200, "A", 11)                 # reports OS
B <- mk_cohort(180, "B", 12, token = "RFS")  # reports RFS (DFS family)
C <- mk_cohort(150, "C", 13, token = "DFS")  # reports DFS
merged <- list(A = A, B = B, C = C)
```

A table that marks B as related and C as exact, so the two modes differ:

``` r
ct <- data.frame(
  Accession = c("A", "B", "C"),
  Family    = c("OS", "DFS", "DFS"),
  Token     = c("OS", "RFS", "DFS"),
  TokenRole = c("primary", "contributing", "primary"),
  PoolingClass = c("Exact-equivalent", "Clinically-related", "Exact-equivalent"),
  stringsAsFactors = FALSE)
ct
#>   Accession Family Token    TokenRole       PoolingClass
#> 1         A     OS    OS      primary   Exact-equivalent
#> 2         B    DFS   RFS contributing Clinically-related
#> 3         C    DFS   DFS      primary   Exact-equivalent
```

``` r
strict <- cpas_meta(c("B", "C"), marker = "GAPDH", type = "DFS", merged = merged,
                    pooling = "exact", class_table = ct)
strict$per_dataset[, c("dataset", "endpoint", "pooling_class", "n", "HR")]
#>   dataset endpoint    pooling_class   n       HR
#> 1       C      DFS Exact-equivalent 150 1.301927
strict$pooled[, c("pooling", "pooling_classes", "k", "HR", "lower", "upper")]
#>   pooling  pooling_classes k       HR    lower    upper
#> 1   exact Exact-equivalent 1 1.301927 1.051162 1.612514
```

## The pooling method: DL, REML and HK

`method` selects the estimator of the between-cohort variance tau^2 and
the interval of the pooled estimate:

- `"REML"` (the default since this version) - iterative restricted
  maximum likelihood;
- `"DL"` - DerSimonian-Laird, the estimator used before, kept for
  continuity (`"RE"` is accepted as a synonym);
- `"HK"` - the Hartung-Knapp-Sidik-Jonkman adjusted variance, on top of
  the REML tau^2, with its characteristically wider interval;
- `"FE"` - fixed effect.

Both prediction intervals are always returned: the primary construction
selected by `pi_method` and a clearly labelled alternative.
`pi_method = "t"` (the default) is the documented t(k-2) interval,
`"normal"` is the normal approximation and `"HK"` is the
Hartung-Knapp-based interval; the rule actually used is recorded in
`$pooled$pi_method` and `$pooled$pi_rule`.

``` r
rows <- lapply(c("DL", "REML", "HK", "FE"), function(mm) {
  p <- cpas_meta(c("A", "B", "C"), marker = "GAPDH", type = "OS", merged = merged,
                 method = mm)$pooled
  data.frame(method = mm, tau2 = p$tau2, HR = p$HR, lower = p$lower,
             upper = p$upper, p = p$p, pi_primary_lo = p$pi_lower,
             pi_alt_lo = p$pi_alt_lower)
})
do.call(rbind, rows)
#>   method        tau2       HR    lower    upper            p pi_primary_lo
#> 1     DL 0.002397897 1.475780 1.296641 1.679668 3.757669e-09     0.5192901
#> 2   REML 0.002201993 1.475904 1.298018 1.678169 2.839229e-09     0.5300406
#> 3     HK 0.002201993 1.475904 1.110321 1.961858 2.768388e-02     0.5266328
#> 4     FE 0.000000000 1.477599 1.314604 1.660803 5.875210e-11     0.6925971
#>   pi_alt_lo
#> 1  1.256174
#> 2  1.260240
#> 3  1.258987
#> 4  1.314604
```

Only the prediction interval changes with `pi_method`:

``` r
pis <- lapply(c("t", "normal", "HK"), function(pm) {
  p <- cpas_meta(c("A", "B", "C"), marker = "GAPDH", type = "OS", merged = merged,
                 pi_method = pm)$pooled
  data.frame(pi_method = p$pi_method, HR = p$HR, pi_lower = p$pi_lower,
             pi_upper = p$pi_upper, pi_alt_lower = p$pi_alt_lower)
})
do.call(rbind, pis)
#>   pi_method       HR  pi_lower pi_upper pi_alt_lower
#> 1         t 1.475904 0.5300406 4.109671    1.2602398
#> 2    normal 1.475904 1.2602398 1.728475    0.5300406
#> 3        HK 1.475904 1.0411323 2.092234    0.5300406
```

Leave-one-out sensitivity uses the same method as the original call:

``` r
mt <- cpas_meta(c("A", "B", "C"), marker = "GAPDH", type = "OS", merged = merged)
loo_meta(mt)[, c("left_out", "HR", "lower", "upper", "I2")]
#>   left_out       HR    lower    upper        I2
#> 1        A 1.467814 1.170212 1.841101 0.5915195
#> 2        B 1.393950 1.204829 1.612758 0.0000000
#> 3        C 1.559346 1.356242 1.792865 0.0000000
```

## Shared-patient overlap

Some cohort pairs share patients: one study split across two platforms,
a series deposited twice, or the same patients re-published. Pooling
both members counts those patients twice, which narrows the pooled
interval without adding information. The register of known pairs ships
with the package.

``` r
ov <- cohort_overlap()
dim(ov)
#> [1] 30  5
head(ov[, c("AccessionA", "AccessionB", "SharedPatients", "Basis")], 4)
#>       AccessionA     AccessionB SharedPatients         Basis
#> 1 GSE37642_GPL96 GSE37642_GPL97            422         title
#> 2  GSE9782_GPL96  GSE9782_GPL97            264         title
#> 3  GSE4922_GPL96  GSE4922_GPL97            249 platform-pair
#> 4  GSE6532_GPL96  GSE6532_GPL97            241 platform-pair
```

`cpas_meta(overlap = )` acts on the pairs that are inside the pool being
built:

- `"warn"` (default) proceeds, warns, and records the pairs in
  `$pooled$overlap_pairs` and in the manifest;
- `"refuse"` stops with an error naming the pair(s);
- `"dedupe"` keeps one member of every overlapping group - the larger
  cohort, ties broken by accession sort - and records what it dropped.

``` r
ov_small <- data.frame(AccessionA = c("A", "B"), AccessionB = c("B", "C"),
                       SharedPatients = c(40L, 12L), Basis = c("title", "platform-pair"),
                       Evidence = c("identical titles", "same ids"),
                       stringsAsFactors = FALSE)
ov_small
#>   AccessionA AccessionB SharedPatients         Basis         Evidence
#> 1          A          B             40         title identical titles
#> 2          B          C             12 platform-pair         same ids
```

``` r
w <- cpas_meta(c("A", "B", "C"), "GAPDH", "OS", merged = merged,
               overlap_table = ov_small)
w$pooled$overlap_pairs_text
#> [1] "A/B (40 shared); B/C (12 shared)"
```

``` r
cpas_meta(c("A", "B", "C"), "GAPDH", "OS", merged = merged,
          overlap = "refuse", overlap_table = ov_small)
#> Error:
#> ! Refusing to pool cohorts that share patients (overlap = "refuse"): A / B share 40 patients; B / C share 12 patients. Drop one member of each pair, or use overlap = "dedupe" to let the function keep the larger cohort of every overlapping group.
```

``` r
dd <- cpas_meta(c("A", "B", "C"), "GAPDH", "OS", merged = merged,
                overlap = "dedupe", overlap_table = ov_small)
dd$per_dataset$dataset
#> [1] "A"
dd$excluded[, c("cohort", "stage", "reason")]
#>   cohort   stage                                              reason
#> 1      B overlap shares patients with A (kept: larger N, 200 vs 180)
#> 2      C overlap shares patients with A (kept: larger N, 200 vs 150)
```

Deduplication is deterministic: the larger cohort wins, and a tie is
broken by the accession sort, so the same pool always drops the same
member.

## The analysis manifest

Everything above is recorded in one object per result, `$manifest`, also
reachable through
[`cpas_manifest()`](https://wangjin93.github.io/CanPAS/reference/cpas_manifest.md).
It is printable and convertible.

``` r
cpas_manifest(mt)
#> CanPAS analysis manifest
#> ========================
#> analysis                   : cpas_meta
#> cohorts                    : A, B, C
#> n_cohorts                  : 3
#> endpoint_family            : OS
#> resolved_token             : OS
#> resolved_tokens_per_cohort : A=OS, B=OS, C=OS
#> token_role                 : primary
#> pooling_mode               : family
#> pooling_classes_pooled     : Exact-equivalent
#> pooling_class_source       : token-map fallback
#> overlap_mode               : warn
#> overlap_source             : inst/extdata/cohort_overlap.csv
#> overlap_pairs              : none
#> overlap_dropped            : none
#> selection_rule             : the 3 cohort(s) supplied in 'datasets' were requested for family OS and the concrete token each one provides was resolved from the catalog; 3 cohort(s) entered the pool and 0 were excluded (reasons in $excluded)
#> rows_dropped               : none
#> n_rows_dropped_reported    : 0
#> covariates_dropped         : none
#> cut_rule                   : not applicable (meta-analysis of a continuous marker; no cut-point is searched)
#> cut_points_searched        : NA
#> search_adjusted_p          : NA
#> ph_test_cox_zph            : GLOBAL p = 0.319 (3 cohort x term check(s))
#> meta_method                : REML
#> auto_repair                : FALSE
#> pi_method                  : t
#> tau2                       : 0.00220199
#> I2                         : 0.183083
#> pi_primary                 : [0.530041, 4.10967]
#> pi_primary_rule            : t distribution with k-2 df on sqrt(se_pooled^2 + tau2)
#> pi_alt                     : [1.26024, 1.72847]
#> pi_alt_rule                : normal approximation (1.96) on sqrt(se_pooled^2 + tau2)
#> versions                   : R 4.5.1; CanPAS 1.0.0; survival 3.8.3; metafor 5.0.1
#> timestamp                  : 2026-09-30 20:56:53 CST
#> notes                      : auto_repair = FALSE (fail-safe): a model that would need covariates dropped, levels merged or rows excluded is reported as not estimable instead of being repaired
#> 
#> per-cohort pooling classes (3 row(s)): $manifest$pooling_table
```

``` r
as.data.frame(cpas_manifest(dd))[c("field", "value")][
  as.data.frame(cpas_manifest(dd))$field %in%
    c("cohorts", "endpoint_family", "resolved_tokens_per_cohort",
      "pooling_classes_pooled", "meta_method", "tau2", "I2", "pi_primary",
      "pi_alt", "overlap_mode", "overlap_dropped"), ]
#>                         field
#> 2                     cohorts
#> 4             endpoint_family
#> 6  resolved_tokens_per_cohort
#> 9      pooling_classes_pooled
#> 11               overlap_mode
#> 14            overlap_dropped
#> 23                meta_method
#> 26                       tau2
#> 27                         I2
#> 28                 pi_primary
#> 30                     pi_alt
#>                                                                                                              value
#> 2                                                                                                                A
#> 4                                                                                                               OS
#> 6                                                                                                             A=OS
#> 9                                                                                                 Exact-equivalent
#> 11                                                                                                          dedupe
#> 14 B: shares patients with A (kept: larger N, 200 vs 180) | C: shares patients with A (kept: larger N, 200 vs 150)
#> 23                                                                                                            REML
#> 26                                                                                                              NA
#> 27                                                                                                              NA
#> 28                                                                                                        [NA, NA]
#> 30                                                                                                        [NA, NA]
```

The structured parts stay data.frames inside the manifest, so they can
be written to CSV or compared between runs:

``` r
str(as.data.frame(cpas_manifest(mt)), max.level = 1)
#> 'data.frame':    34 obs. of  2 variables:
#>  $ field: chr  "analysis" "cohorts" "n_cohorts" "endpoint_family" ...
#>  $ value: chr  "cpas_meta" "A, B, C" "3" "OS" ...
cpas_manifest(dd)$overlap_dropped
#>   cohort kept shared_patients group
#> 1      B    A              40 A|B|C
#> 2      C    A              NA A|B|C
#>                                                reason
#> 1 shares patients with A (kept: larger N, 200 vs 180)
#> 2 shares patients with A (kept: larger N, 200 vs 150)
cpas_manifest(strict)$pooling_table
#>   accession family token token_role    pooling_class pooling_class_source
#> 1         C    DFS   DFS    primary Exact-equivalent      companion table
```

## Fail-safe: a model that would need repair is refused

The same principle applies inside a cohort. `auto_repair = FALSE` (the
default) returns an explicit `not estimable` result instead of quietly
dropping a covariate, merging a level or excluding patients;
`auto_repair = TRUE` restores the repair and records every modification.

``` r
A$flat <- 1
bad <- COX_analysis(A, type = "OS", cont_Variates = c("GAPDH", "flat"),
                    method = "multi")
bad$status
#> [1] "not estimable"
bad$offending_terms
#> [1] "flat"
```

``` r
ok <- COX_analysis(A, type = "OS", cont_Variates = c("GAPDH", "flat"),
                   method = "multi", auto_repair = TRUE)
ok$metadata$dropped_covariates[, c("variable", "reason")]
#>   variable                             reason
#> 1     flat constant across the complete cases
cpas_manifest(ok)$dropped_covariates
#>   cohort variable     detail                             reason
#> 1      A     flat continuous constant across the complete cases
```

## Checking the annotation itself: blinded adjudication

The semantics table encodes *our* reading of each deposit, so it needs
an audit that does not simply repeat it.
[`endpoint_adjudication()`](https://wangjin93.github.io/CanPAS/reference/endpoint_adjudication.md)
draws a stratified sample (by endpoint family and by source) of
cohort-endpoint records and returns a **blinded** rating sheet: the
record id, the cohort, the token and the evidence text only - the family
assignment, the pooling class and the four semantic fields are
withheld - with one empty column per judgement field and rater.

``` r
sheet <- endpoint_adjudication(n = 12, raters = 2, seed = 20260930)
sheet[, c("record_id", "cohort", "token")]
#>    record_id            cohort      token
#> 1    REC0122          GSE11121 not stated
#> 2    REC0211          GSE15459 not stated
#> 3    REC0273         GSE183088 not stated
#> 4    REC0306          GSE20685        RFS
#> 5    REC0346          GSE22219        RFS
#> 6    REC0401          GSE27020        DFS
#> 7    REC0505     GSE3494_GPL97 not stated
#> 8    REC0516    GSE37642_GPL97 not stated
#> 9    REC0541 GSE40272_GPL15973        DFS
#> 10   REC0604   GSE4716_GPL3696         OS
#> 11   REC0636           GSE5327 not stated
#> 12   REC0796           GSE9195        RFS
attr(sheet, "pooling_class_vocabulary")
#> [1] "Exact-equivalent"   "Clinically-related" "Not-poolable"      
#> [4] "Unknown"            "Absent"
```

The `record_id` is stable: it is assigned to the whole table before
filtering and sampling, so a record keeps its id between draws. Once the
raters have filled the sheet in,
[`endpoint_agreement()`](https://wangjin93.github.io/CanPAS/reference/endpoint_agreement.md)
reports the per-field raw agreement, Cohen’s kappa with its
Landis-and-Koch label, the overall pooled agreement, the adjudication
rate and the disagreeing records. The example below uses a tiny
synthesised ratings frame - **not a real study result** - whose kappa is
known by hand: 20 records where both raters say `A`, 5 where only rater
1 says `A`, 10 the other way round and 15 where both say `B`, so
`po = 0.70`, `pe = (25*30 + 25*20)/50^2 = 0.50` and `kappa = 0.40`.

``` r
synth <- data.frame(
  record_id = sprintf("REC%04d", 1:50),
  cohort = "SYNTHETIC", token = "OS", evidence = "synthetic example",
  rater1_pooling_class = c(rep("A", 25), rep("B", 25)),
  rater2_pooling_class = c(rep("A", 20), rep("B", 5), rep("A", 10), rep("B", 15)))
ag <- endpoint_agreement(synth, fields = "pooling_class")
ag$per_field
#>           field n_compared n_agree raw_agreement kappa kappa_interpretation
#> 1 pooling_class         50      35           0.7   0.4                 fair
#>                  kappa_basis note
#> 1 Cohen's kappa (two raters) <NA>
ag$overall_raw_agreement
#> [1] 0.7
ag$adjudication_rate
#> [1] 0.3
head(ag$disagreements[, c("record_id", "fields_disagreed")], 3)
#>   record_id fields_disagreed
#> 1   REC0021    pooling_class
#> 2   REC0022    pooling_class
#> 3   REC0023    pooling_class
```

With more than two raters the field’s kappa is the mean of the pairwise
Cohen’s kappas, and the result says so.

``` r
synth$rater3_pooling_class <- sample(c("A", "B"), 50, replace = TRUE)
ag3 <- endpoint_agreement(synth, fields = "pooling_class")
ag3$kappa_basis
#> [1] "Cohen's kappa per field is the mean of its 3 pairwise kappas (>2 raters)"
ag3$pairwise[, c("pair", "n_compared", "raw_agreement", "kappa")]
#>   pair n_compared raw_agreement      kappa
#> 1  1-2         50          0.70  0.4000000
#> 2  1-3         50          0.40 -0.2000000
#> 3  2-3         50          0.38 -0.2601626
```
