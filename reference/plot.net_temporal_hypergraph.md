# Plot a temporal-hypergraph snapshot

Plots the snapshot that
[`hg_snapshot()`](https://mohsaqr.github.io/hypergraphs/reference/hg_snapshot.md)
takes at `at` as a hypergraph, with
[`plot.net_hg()`](https://mohsaqr.github.io/hypergraphs/reference/plot.net_hg.md):
the hyperedges active then, as hulls by default or as the incidence
matrix with `type = "incidence"`. The snapshot keeps the hyperedge
attributes and the data's names for nodes and hyperedges, so `color_by`
can name an attribute and the legends use the data's words.

## Usage

``` r
# S3 method for class 'net_temporal_hypergraph'
plot(x, at = NULL, mode = c("active", "cumulative"), method = NULL, ...)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/temporal_hypergraph.md).

- at:

  Snapshot time; defaults to the end of observation.

- mode:

  Snapshot mode passed to
  [`hg_snapshot()`](https://mohsaqr.github.io/hypergraphs/reference/hg_snapshot.md).

- method:

  Deprecated. The snapshot used to be projected to a pairwise network;
  `plot(pairwise_network(hg_snapshot(x, at), type = method))` plots that
  network.

- ...:

  Arguments passed to
  [`plot.net_hg()`](https://mohsaqr.github.io/hypergraphs/reference/plot.net_hg.md),
  such as `type`, `color_by` or `labels`.

## Value

A ggplot object, as
[`plot.net_hg()`](https://mohsaqr.github.io/hypergraphs/reference/plot.net_hg.md)
returns.

## Conditions

`hypergraphs_bad_input` when the snapshot has no active nodes, and as
[`plot.net_hg()`](https://mohsaqr.github.io/hypergraphs/reference/plot.net_hg.md)
raises it; `hypergraphs_deprecated` (a warning) for `method`.

## Examples

``` r
seats <- data.frame(
  case = rep(c("A", "B", "C"), each = 3),
  arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
  constituted = rep(c(1, 2, 4), each = 3),
  concluded = rep(c(4, 3, 6), each = 3)
)
thg <- temporal_hypergraph(seats, node = "arbitrator", hyperedge = "case",
                           start = "constituted", end = "concluded")
plot(thg, at = 2.5)

plot(thg, at = 2.5, type = "incidence", edge_labels = TRUE)
```
