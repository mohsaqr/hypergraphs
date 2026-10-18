# Tidy tables of a temporal hypergraph

Tidy tables of a temporal hypergraph

## Usage

``` r
# S3 method for class 'net_temporal_hypergraph'
hg_get(x, what = c("memberships", "edges", "nodes"), ...)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/temporal_hypergraph.md).

- what:

  `"memberships"` (default) for one row per node-in-hyperedge spell
  (`member`, `edge`, `start`, `end`, `weight`), `"edges"` for one row
  per hyperedge with its clock and attributes, `"nodes"` for the node
  universe with entry times. Times are on the hypergraph's clock.

- ...:

  Ignored.

## Value

A base `data.frame`.

## Examples

``` r
seats <- data.frame(
  case = c("A", "A", "A", "B", "B", "B"),
  arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
  constituted = c(1, 1, 1, 2, 2, 2), concluded = c(4, 4, 4, 5, 5, 5)
)
thg <- temporal_hypergraph(seats, node = "arbitrator", hyperedge = "case",
                           start = "constituted", end = "concluded")
hg_get(thg)
#>   node edge start end weight
#> 1   p1    A     1   4      1
#> 2   a1    A     1   4      1
#> 3   a2    A     1   4      1
#> 4   p2    B     2   5      1
#> 5   a1    B     2   5      1
#> 6   a3    B     2   5      1
hg_get(thg, what = "edges")
#>   edge start end
#> 1    A     1   4
#> 2    B     2   5
```
