# Edge table of a memory-network comparison

Edge table of a memory-network comparison

## Usage

``` r
# S3 method for class 'hypergraphs_comparison'
hg_get(x, what = "edges", ..., significant = FALSE, sort_by = NULL, top = NULL)
```

## Arguments

- x:

  A `hypergraphs_comparison` object from
  [`hg_compare()`](https://pak.dynasite.org/hypergraphs/reference/hg_compare.md).

- what:

  `"edges"`, the only table: one row per pooled rule edge.

- ...:

  Additional arguments (ignored).

- significant:

  Logical. `TRUE` restricts to edges whose adjusted p-value falls below
  the object's `alpha`. Default `FALSE` (all edges).

- sort_by:

  `NULL` (order/from/to, default), `"abs_diff"`, `"count"`, or
  `"p_adj"` - sort by absolute difference or count (largest first) or
  adjusted p-value (smallest first), ties broken by from/to.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per pooled rule edge (see
[`hg_compare()`](https://pak.dynasite.org/hypergraphs/reference/hg_compare.md)
for the columns).

## Examples

``` r
set.seed(1)
first_order <- replicate(6, sample(c("a", "b", "c"), 12, replace = TRUE),
                         simplify = FALSE)
second_order <- replicate(6, rep(c("a", "b", "c", "b"), 3),
                          simplify = FALSE)
by_group <- hon(c(first_order, second_order),
                group = rep(c("first", "second"), each = 6), max_order = 2)
comparison <- hg_compare(by_group, n_perm = 99, seed = 1)
hg_get(comparison, sort_by = "p_adj", top = 6)
#>   from to order count count_first count_second probability_first
#> 1    a  b     1    23           5           18         0.2500000
#> 2    a  c     1     8           8            0         0.4000000
#> 3    b  b     1    14          14            0         0.5384615
#> 4    b  c     1    24           6           18         0.2307692
#> 5    c  a     1     6           6            0         0.3000000
#> 6    c  b     1    26           8           18         0.4000000
#>   probability_second       diff p_value      p_adj significant n_perm_used
#> 1                1.0 -0.7500000    0.01 0.01333333        TRUE          99
#> 2                0.0  0.4000000    0.01 0.01333333        TRUE          99
#> 3                0.0  0.5384615    0.01 0.01333333        TRUE          99
#> 4                0.6 -0.3692308    0.01 0.01333333        TRUE          99
#> 5                0.0  0.3000000    0.01 0.01333333        TRUE          99
#> 6                1.0 -0.6000000    0.01 0.01333333        TRUE          99
```
