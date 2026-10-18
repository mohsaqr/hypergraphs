# Extract one snapshot from a temporal hypergraph

A snapshot is a static `net_hg` holding the hyperedges a window
measures. `mode = "active"` follows the format: an interval hyperedge is
active on its closed interval, a contact hyperedge at its instant.
`mode = "cumulative"` keeps every hyperedge begun by the end of the
window, the growing view of a contact log; with `at` unset it is the
static aggregate of everything observed.

## Usage

``` r
hg_snapshot(
  x,
  at = NULL,
  window = 0,
  mode = c("active", "cumulative"),
  multiedges = TRUE
)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/temporal_hypergraph.md).

- at:

  One time value, as a number on the hypergraph's clock or a date for a
  calendar hypergraph. `NULL` (default) is the end of observation.

- window:

  How much time the snapshot covers, in the hypergraph's time unit: `0`
  (default) evaluates the instant `at`; a positive width covers
  `[at, at + window)`; `"all"` covers the whole observation period.

- mode:

  `"active"` (default) or `"cumulative"`, see Description. The former
  `"all"` is deprecated: it is `"cumulative"` with `at` unset.

- multiedges:

  Keep distinct edge identities with identical member sets? `TRUE`
  matches a multi-hypergraph; `FALSE` collapses them to one edge and
  records their counts in `edge_multiplicity`.

## Value

A static `net_hg` usable by every hypergraphs hypergraph verb; its
`params` record `at`, `window`, `temporal_mode` and the clock.

## Examples

``` r
seats <- data.frame(
  case = c("A", "A", "A", "B", "B", "B"),
  arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
  constituted = c(1, 1, 1, 2, 2, 2), concluded = c(4, 4, 4, 5, 5, 5)
)
thg <- temporal_hypergraph(seats, node = "arbitrator", hyperedge = "case",
                           start = "constituted", end = "concluded")
hg_snapshot(thg, at = 3)
#> Hypergraph: 5 nodes, 2 hyperedges (sizes 3: 2)
#> Source: group membership (node = arbitrator, hyperedge = case)
#>  hyperedge size    members
#>          A    3 a1, a2, p1
#>          B    3 a1, a3, p2
hg_snapshot(thg, mode = "cumulative")
#> Hypergraph: 5 nodes, 2 hyperedges (sizes 3: 2)
#> Source: group membership (node = arbitrator, hyperedge = case)
#>  hyperedge size    members
#>          A    3 a1, a2, p1
#>          B    3 a1, a3, p2
```
