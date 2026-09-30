# Agreement between endpoint raters

Takes a filled rating sheet (from
[`endpoint_adjudication`](https://wangjin93.github.io/CanPAS/reference/endpoint_adjudication.md))
and reports how far the raters agree: per-field raw agreement, **Cohen's
kappa** per field with its Landis-and-Koch interpretation, the overall
pooled agreement, the adjudication rate and the records whose raters
disagreed.

## Usage

``` r
endpoint_agreement(sheet, fields = .cpas_adjudication_fields)
```

## Arguments

- sheet:

  A filled rating sheet: the columns of
  [`endpoint_adjudication`](https://wangjin93.github.io/CanPAS/reference/endpoint_adjudication.md)
  (`record_id`, `cohort`, `token`, `evidence`) plus at least two
  `rater<r>_<field>` columns per judgement field with at least one
  non-empty pair.

- fields:

  Judgement fields to score (default: the five the sheet prepares:
  `event_definition`, `time_origin`, `censoring_rule`,
  `competing_events`, `pooling_class`).

## Value

A list of class `cpas_adjudication_agreement`:

- `n_records`, `n_raters`, `fields`::

  the size of the exercise and the fields scored

- `per_field`::

  one row per field: `field`, `n_compared` (record-pairs in which both
  ratings are present), `n_agree`, `raw_agreement`, `kappa`,
  `kappa_interpretation`, `kappa_basis`, `note`

- `overall_raw_agreement`::

  pooled over all fields and rater pairs (agreeing compared pairs / all
  compared pairs)

- `pooled_kappa`::

  the mean of the per-field kappas (the pooled agreement on the kappa
  scale)

- `adjudication_rate`::

  records with at least one disagreement divided by the records that
  carry at least one comparable rating

- `disagreements`::

  the disagreeing records, ready for adjudication: `record_id`,
  `cohort`, `token`, `evidence`, `fields_disagreed` and `detail` (which
  rater said what, field by field)

- `pairwise`::

  per field and rater pair, the raw agreement and kappa; for two raters
  this is the same number as `per_field`

- `kappa_basis`, `notes`::

  how kappa was aggregated, and any field whose kappa is undefined

## Details

**Kappa.** Cohen's kappa is computed from the nominal contingency table
of the two raters: \\\kappa = (p_o - p_e) / (1 - p_e)\\, with \\p_o\\
the observed agreement and \\p_e\\ the agreement expected from the
margins. Only rater pairs in which both ratings are present and
non-empty are used, so missing ratings reduce `n_compared` instead of
biasing the estimate, and a field in which every rating is the same
single category reports `kappa = NA` (chance agreement 1 leaves it
undefined) rather than `NaN`. With more than two raters every rater pair
is scored and the field's kappa is the **mean of the pairwise Cohen's
kappas** - kappa is defined for two raters, and this is said explicitly
in `kappa_basis` and in `notes`. For two raters there is exactly one
pair and the mean is that pair's kappa.

## See also

[`endpoint_adjudication`](https://wangjin93.github.io/CanPAS/reference/endpoint_adjudication.md)

## Examples

``` r
## A tiny synthesised sheet (no real study data): two raters, known table.
## 20 records where both say "A", 5 where rater1 says "A" and rater2 "B",
## 10 the other way round and 15 where both say "B":
##   po = 35/50 = 0.70 ; pe = (25*30 + 25*20)/50^2 = 0.50 ; kappa = 0.40
set.seed(1)
r1 <- c(rep("A", 25), rep("B", 25))
r2 <- c(rep("A", 20), rep("B", 5), rep("A", 10), rep("B", 15))
sheet <- data.frame(record_id = sprintf("REC%04d", 1:50),
                    cohort = "SYN", token = "OS", evidence = "synthetic",
                    rater1_pooling_class = r1, rater2_pooling_class = r2,
                    stringsAsFactors = FALSE)
ag <- endpoint_agreement(sheet, fields = "pooling_class")
ag$per_field[, c("field", "n_compared", "raw_agreement", "kappa",
                 "kappa_interpretation")]
#>           field n_compared raw_agreement kappa kappa_interpretation
#> 1 pooling_class         50           0.7   0.4                 fair
```
