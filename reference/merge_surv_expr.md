# Merge Survival and Expression Data of a Dataset

Pulls the survival table of a dataset from the CanPAS API and joins it
with gene expression data produced by
[`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)
on the common sample identifier (`ID`).

For datasets whose accession starts with `"TCGA-"` the survival table is
not read from the mirror: the local TCGA clinical table is used instead
([`tcga_surv_table`](https://wangjin93.github.io/CanPAS/reference/tcga_surv_table.md),
the same source the App and
[`tcga_merged`](https://wangjin93.github.io/CanPAS/reference/tcga_merged.md)
use).

## Usage

``` r
merge_surv_expr(dataset, expr_data)
```

## Arguments

- dataset:

  Dataset accession in `dataset_info` (e.g. `"GSE14814"` or
  `"TCGA-LUAD"`).

- expr_data:

  Result object of class `cpas_get_expr` returned by
  [`get_expr_data`](https://wangjin93.github.io/CanPAS/reference/get_expr_data.md)
  (an object with an `expr_data` element whose first column is `ID`) or
  a plain data.frame in that format.

## Value

List of class `cpas_merge`:

- `input_params`::

  echoed inputs and merge time

- `raw_surv_data`::

  survival table exactly as returned by the API (or the local TCGA
  clinical table)

- `raw_expr_data`::

  expression input

- `merged_data`::

  samples x (ID, `<TYPE>_time`, `<TYPE>_status`, gene columns). All
  survival times are years and all status columns are binary (0 censored
  / 1 event); rows with missing time or status are kept (see `NA`) so
  that callers may decide how to treat them, while gene columns are
  numeric

- `metadata`::

  sample counts and matched columns

## Details

Time columns are interpreted as years; status columns must be binary 0/1
(the format used by the CanPAS database). Columns of the survival table
other than `ID`, `*_status` and `*_time` (e.g. clinical covariates) are
intentionally not carried into `merged_data`: clinical data live in the
local survival RDS mirrors and can be joined by the caller.

## Examples

``` r
if (FALSE) { # \dontrun{
   ## The two prerequisite steps and the merge: expression first, survival
   ## second, giving the analysis-ready table used by every other function.
   e <- get_expr_data("GSE14814", c("GAPDH", "ACTB"), process_duplicates = "max")
   m <- merge_surv_expr("GSE14814", e)
   colnames(m$merged_data)

   ## raw_surv_data keeps the clinical columns of that cohort (age, sex, ...);
   ## cohort_merged(clin = TRUE) adds them to the merged table for you.
   colnames(m$raw_surv_data)
   head(m$merged_data)
} # }
```
