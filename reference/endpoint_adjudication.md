# Blinded endpoint-adjudication rating sheet

Draws a stratified sample of cohort-endpoint records from the shipped
endpoint-semantics table and returns a blinded rating sheet: one row per
record with a stable `record_id`, the cohort, the endpoint **token** and
the **evidence text** only, plus one empty column per judgement field
and rater for the raters to fill.

The sheet deliberately does **not** carry our own verdict: the `Family`
assignment, the `PoolingClass` and the four semantic fields
(`EventDefinition`, `TimeOrigin`, `CensoringRule`, `CompetingEvents`)
are withheld, so a rater judges the endpoint from the token and the
evidence text alone and cannot anchor on the annotation under review.

## Usage

``` r
endpoint_adjudication(
  n = 80,
  raters = 2,
  seed = 20260101,
  families = NULL,
  drop_absent = FALSE,
  path = NULL,
  table = NULL
)
```

## Arguments

- n:

  Number of records to draw (default 80; capped at the number of
  eligible records).

- raters:

  Number of raters the sheet is prepared for (default 2). Any number \>=
  1 is accepted; the agreement calculator reports pairwise mean kappa
  when more than two raters filled it.

- seed:

  Fixed random seed (default `20260101`) making the draw reproducible.
  The previous RNG state is restored on exit.

- families:

  Optional character vector of endpoint families to sample from
  (default: all five).

- drop_absent:

  Drop records whose token is absent for the family (`default FALSE`;
  set `TRUE` to sample only records that report an endpoint).

- path:

  Optional file path. When supplied the sheet is also written there as
  CSV with empty cells for the blank ratings (and the path is recorded
  in `attr(x, "csv_path")`).

- table:

  Optional endpoint-semantics table (the frozen schema of
  [`endpoint_semantics`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md))
  used instead of the shipped copy.

## Value

A data frame of class `cpas_adjudication_sheet` with columns
`record_id`, `cohort`, `token`, `evidence` and, for every rater `r` and
field `f`, `rater<r>_<f>` (`NA` until filled). The attributes record the
sampling design (`attr(x, "sampling")`), the pooling-class vocabulary
the raters must use (`attr(x, "pooling_class_vocabulary")`) and the
seed.

## Details

**Sampling.** Records are stratified by endpoint family and by source
(the cohort's accession prefix: GEO, TCGA, EMBL-EBI, CGGA, other), and
the `n` records are allocated across the strata in proportion to their
size (largest-remainder rounding, capped at each stratum's size). A
stratum smaller than its share is taken whole and the remainder is
redistributed, so exactly `n` records are returned whenever the table
has that many.

**Stable ids.** `record_id` is assigned to the whole cohort-endpoint
table in catalog order (`Accession`, then `Family`) as `REC0001`,
`REC0002`, ... before any filtering or sampling. A record therefore
keeps the same id when `n`, `seed`, `families`, `drop_absent` or the
rater count changes, and two sheets can be compared by id. The
consequence is that the ids of a sample are sparse (they are not
`1..n`).

**Blank fields.** Every rater column is `NA_character_` so that a
missing rating is distinguishable from a rating of `""`: the agreement
calculator compares only pairs in which both ratings are present and
non-empty.

## See also

[`endpoint_agreement`](https://wangjin93.github.io/CanPAS/reference/endpoint_agreement.md),
[`endpoint_semantics`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md),
[`endpoint_pooling_class`](https://wangjin93.github.io/CanPAS/reference/endpoint_pooling_class.md)

## Examples

``` r
## Offline: the shipped endpoint-semantics table, no network.
sheet <- endpoint_adjudication(n = 6, raters = 2, seed = 42)
sheet[, c("record_id", "cohort", "token")]
#>   record_id        cohort      token
#> 1   REC0331      GSE21653        DFS
#> 2   REC0411      GSE28735 not stated
#> 3   REC0456      GSE31519        RFS
#> 4   REC0696 GSE6532_GPL97        RFS
#> 5   REC0726      GSE71014 not stated
#> 6   REC0811       GSE9893        RFS
attr(sheet, "pooling_class_vocabulary")
#> [1] "Exact-equivalent"   "Clinically-related" "Not-poolable"      
#> [4] "Unknown"            "Absent"            

## Fill it in (here: one rater copies the second for the worked example).
sheet$rater1_pooling_class <- rep("Exact-equivalent", nrow(sheet))
sheet$rater2_pooling_class <- sheet$rater1_pooling_class
endpoint_agreement(sheet)$overall_raw_agreement
#> [1] 1
```
