# Tables of a path-anomaly (HYPA) fit

Tables of a path-anomaly (HYPA) fit

## Usage

``` r
# S3 method for class 'hypergraphs_hypa'
hg_get(x, ...)

# S3 method for class 'net_hypa'
hg_get(
  x,
  what = c("anomalies", "scores", "over", "under", "pathways"),
  ...,
  type = NULL,
  order_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hypa` object from
  [`hypa()`](https://pak.dynasite.org/hypergraphs/reference/hypa.md) on
  sequences.

- ...:

  Unused.

- what:

  `"anomalies"` (default) for the anomalous paths selected by `type`,
  `"scores"` for every scored path, `"over"` or `"under"` for the
  significantly over- or under-represented paths with all their columns,
  or `"pathways"` for the anomalous paths written as pathways
  (`"a b -> c"`: the context states, then the next state), the form
  [`simplicial()`](https://pak.dynasite.org/hypergraphs/reference/simplicial.md)
  turns into simplices.

- type:

  For `"anomalies"` and `"pathways"`: `"all"`, `"over"` or `"under"`.
  Default: the `type` given to
  [`hypa()`](https://pak.dynasite.org/hypergraphs/reference/hypa.md),
  else `"all"`.

- order_by:

  `"sig"` (smallest p-value of the path's own direction first),
  `"ratio"` (largest deviation first: highest ratio for over-represented
  paths, lowest for under-represented ones, largest absolute log ratio
  when both directions are listed), `"freq"` (most observed first) or
  `"path"` (alphabetical). Default: the `order_by` given to
  [`hypa()`](https://pak.dynasite.org/hypergraphs/reference/hypa.md),
  else `"sig"`. For `"scores"`, `"sig"` sorts by the smaller of the two
  raw tails, `p_under` and `p_over`, so the paths that depart most from
  the null in either direction come first.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `order_by`. Default `NULL` returns every row.

## Value

A data.frame. For `"anomalies"`, one row per anomalous path: `path`,
`count`, `expected`, `ratio`, `p_adj` (the adjusted p-value of the
path's own direction), `direction`, and `order` when several orders were
scored. For `"scores"`, `"over"` and `"under"`, one row per path with
every column of the fit: `path`, `from`, `to`, `count`, `expected`,
`ratio`, the p-values (`p_value`, an alias of `p_under`, the lower tail
\\P(X \le f)\\; `p_over`, the strict upper tail \\P(X \> f)\\;
`p_adj_under`, `p_adj_over`), the `direction` (`"over"`, `"under"` or
`"normal"`) and the path `order`. `"scores"` holds every possible path:
the observed ones and the unobserved ones the observed transitions
permit (`count` 0). For `"pathways"`, one row per anomalous path with
the column `pathway`.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
fit <- hypa(seqs, order = 2, min_count = 1)
hg_get(fit, type = "over", order_by = "ratio", top = 5)
#>          path count expected ratio      p_adj direction
#> 1 a -> b -> c     4        2     2 0.03548389      over
#> 2 x -> b -> d     4        2     2 0.03548389      over
```
