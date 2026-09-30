# Pooling class of one cohort x endpoint-family cell

Returns the `PoolingClass` of one cohort and endpoint family:

- `"Exact-equivalent"`:

  the token IS the family's canonical definition (OS in OS, DFS in DFS,
  ...)

- `"Clinically-related"`:

  a different token pooled into the family by the documented token -\>
  family map (RFS/EFS/DFI in DFS, CSS/BCSS in DSS, PFI in PFS, DRFS in
  MFS)

- `"Not-poolable"`:

  the cohort reports an endpoint in this family that cannot be pooled

- `"Unknown"`:

  the deposit does not state what the endpoint means

- `"Absent"`:

  the cohort has no endpoint of this family at all. This is the absence
  of an endpoint, **not** a pooling verdict, and such a cell is never a
  candidate for pooling. Where the companion table spells an absent cell
  `"Not-poolable"` with no token (`TokenRole = "not stated"`), this
  accessor reports `"Absent"` instead; a table that already uses
  `"Absent"` is read identically.

The value comes from the shipped companion table
([`endpoint_semantics`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md));
when that table is not available the class is derived from the package's
documented token -\> family map, which is the same rule the table
records for the exact/related split.

## Usage

``` r
endpoint_pooling_class(accession, family, class_table = NULL)
```

## Arguments

- accession:

  Cohort accession (e.g. `"GSE31210"`).

- family:

  Pooling family (`OS`, `DSS`, `DFS`, `PFS`, `MFS`) or a raw endpoint
  token.

- class_table:

  Optional companion table to use instead of the shipped one.

## Value

Character scalar: one of `"Exact-equivalent"`, `"Clinically-related"`,
`"Not-poolable"`, `"Unknown"` or `"Absent"`.

## See also

[`endpoint_semantics`](https://wangjin93.github.io/CanPAS/reference/endpoint_semantics.md),
[`cpas_meta`](https://wangjin93.github.io/CanPAS/reference/cpas_meta.md)

## Examples

``` r
  ## The token a cohort contributes decides whether it is the canonical
  ## definition of the family or a related one pooled into it by design.
  endpoint_pooling_class("GSE31210", "DFS")       # "Clinically-related" (RFS)
#> [1] "Clinically-related"
  endpoint_pooling_class("GSE4922_GPL96", "DFS")  # "Exact-equivalent"
#> [1] "Exact-equivalent"
  endpoint_pooling_class("GSE13507", "MFS")       # "Absent"
#> [1] "Absent"
```
