# Tables of a hypergraph

One row per hyperedge: its identifier, size (number of distinct member
states), the member states, and its weight (for hypergraphs built by
[`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md),
the number of windows collapsed into the hyperedge; `NA` for the other
constructors, whose hyperedges are unweighted).

## Usage

``` r
# S3 method for class 'net_hg'
hg_get(
  x,
  what = c("edges", "nodes", "memberships", "sets", "state_counts", "edge_data"),
  ...,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hg` object.

- what:

  `"edges"` (default) for one row per hyperedge, `"nodes"` for one row
  per node, or `"memberships"` for one row per node-in- hyperedge cell
  of the incidence matrix. A hypergraph of clustered sequences
  ([`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md)
  on a clustering of sequences) also has `"sets"` and `"state_counts"`.
  `"edge_data"` returns the per-hyperedge attribute table a constructor
  attached.

- ...:

  Additional arguments (ignored).

- sort_by:

  `NULL` (construction order, default), `"weight"`, or `"size"` - sort
  the table by that column, largest first (ties broken by hyperedge id
  so the order is deterministic).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "edges"`, one row per hyperedge with columns
`hyperedge` (character id), `size` (integer), `members` (the member
nodes, comma separated), and, for a hypergraph of windows
([`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md)),
`weight` (the number of windows with that member set). `sort_by` is
`"size"`, or `"weight"` where it exists. For `what = "nodes"`, one row
per node with columns `node` and `degree` (the number of hyperedges it
belongs to), plus `block` for a hypergraph with planted blocks
([`random_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/random_hypergraph.md)
with `type = "sbm"`); `sort_by = "degree"` orders it. For
`what = "memberships"`, one row per non-zero incidence cell with columns
`node`, `hyperedge` and `weight` (the incidence value: 1 for a binary
hypergraph, the summed weight for a weighted one), in hyperedge order;
`sort_by = "weight"` orders it. For `what = "sets"`, one row per
hyperedge with `group`, `hyperedge`, `set` (the states joined by
`" + "`), `size`, `count` (sequences of the group with exactly that set)
and `share` (`count` over the group's sequences), in group order and
decreasing count. For `what = "state_counts"`, one row per group and
state with `group`, `node`, `count` (sequences of the group containing
the state) and `share`. For `what = "edge_data"`, one row per hyperedge
with `edge` and its attributes; asking for it of a hypergraph without an
attribute table raises `hypergraphs_bad_input`.

## Examples

``` r
hg <- window_hypergraph(list(s1 = c("a", "b", "a", "c")), window = 2L)
hg_get(hg)
#>   hyperedge size members weight
#> 1        h1    2    a, b      2
#> 2        h2    2    a, c      1
hg_get(hg, sort_by = "weight")
#>   hyperedge size members weight
#> 1        h1    2    a, b      2
#> 2        h2    2    a, c      1
hg_get(hg, sort_by = "weight", top = 3)
#>   hyperedge size members weight
#> 1        h1    2    a, b      2
#> 2        h2    2    a, c      1
hg_get(hg, what = "nodes", sort_by = "degree")
#>   node degree
#> 1    a      2
#> 2    b      1
#> 3    c      1
hg_get(hg, what = "memberships")
#>   node hyperedge weight
#> 1    a        h1      2
#> 2    b        h1      2
#> 3    a        h2      1
#> 4    c        h2      1
```
