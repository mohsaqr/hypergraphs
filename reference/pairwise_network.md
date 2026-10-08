# The pairwise network of a hypergraph

Projects a hypergraph onto a network of pairs of its nodes. Two nodes
are joined when they belong to a common hyperedge, and the weight of the
edge depends on `type`. `"clique"` is the clique expansion: a hyperedge
of size \\\|e\|\\ contributes the same amount to each of its
\\\|e\|(\|e\|-1)/2\\ pairs, so the weight of a pair is the number of
hyperedges that contain both nodes (Zhou et al. 2006). `"association"`
divides the contribution of each hyperedge by \\\|e\|-1\\, following
Coupette et al. (2024):

## Usage

``` r
pairwise_network(
  hg,
  type = c("clique", "association", "citation"),
  weighted = NULL,
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL,
  directed = FALSE
)
```

## Arguments

- hg:

  A hypergraph (`net_hg`), such as one built by
  [`hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/hypergraph.md)
  or
  [`text_hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/text_hypergraph.md).

- type:

  The projection: `"clique"` (default), `"association"` or `"citation"`.

- weighted:

  `type = "clique"` only. `TRUE` (default) uses the incidence weights,
  `FALSE` their membership pattern. Setting it with another `type` is an
  error, because those projections are defined on membership.

- duplicate_edges:

  For `type = "association"` or `"citation"`: `"count"` (default) lets
  repeated hyperedges contribute repeatedly (the multi-hypergraph);
  `"collapse"` lets each distinct member set contribute once (the binary
  hypergraph).

- self_association:

  For `type = "association"`: add the source-to-member association term
  of Coupette et al. (2024). Requires one source node per hyperedge
  through `edge_source` or `hg$edge_data$source`. For source `u`, every
  membership of target `v` contributes `1 / sum_e |e|` over the
  hyperedges sourced by `u`.

- edge_source:

  The node each hyperedge comes from: the name of a column of the
  hyperedge attributes, a vector of length `n_hyperedges`, a named
  vector keyed by hyperedge, or a data frame with columns `edge` and
  `source` (a hyperedge named twice with different sources raises
  `hypergraphs_bad_input`). `NULL` uses an attribute column named
  `source`. Used with `self_association = TRUE` and `type = "citation"`.

- directed:

  For `type = "citation"`: keep the source-to-member direction. Default
  `FALSE`.

## Value

A `net_hg_pairwise`, which is also a `netobject` and a
`cograph_network`: the weighted adjacency matrix of the projection
(`weights`, zero diagonal, symmetric unless `directed = TRUE`), its
nodes and edges. `hg_get(x)` returns the edges as a data frame with one
row per pair of non-zero weight and columns `from`, `to` and `weight`,
sorted by `from` then `to`.

## Details

\$\$w(\\u,v\\) = \sum\_{e \supseteq \\u,v\\} \frac{1}{\|e\| - 1}\$\$

so the weight a hyperedge adds around each of its members is 1, and the
weighted degree of a node equals the number of hyperedges of size at
least two that contain it. A hyperedge of size one has no pairs and
contributes nothing.

`"citation"` is the graph of a hypergraph whose hyperedges have sources:
one edge from the source of every hyperedge to each of its members, so a
hypergraph of citation blocks becomes the ordinary citation graph.
`duplicate_edges = "count"` weights an edge by the number of blocks that
repeat it (the multi-graph); `"collapse"` keeps the binary graph.
`directed = TRUE` keeps the source-to-member orientation; the default
symmetrises.

The result is a network object that cograph plots and that
[`hypergraph()`](https://pak.dynasite.org/hypergraphs/reference/hypergraph.md)
reads as a network, so the cliques of the projection can be promoted
back to hyperedges. Its edges are read with
[`hg_get()`](https://pak.dynasite.org/hypergraphs/reference/hg_get.md).

## Conditions

`hypergraphs_bad_input` for an argument that does not apply to `type`,
and for `self_association` or `type = "citation"` without a source per
hyperedge.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

Zhou, D., Huang, J., & Schoelkopf, B. (2006). Learning with hypergraphs:
clustering, classification, and embedding. *NeurIPS 19*, 1601-1608.
[doi:10.7551/mitpress/7503.003.0205](https://doi.org/10.7551/mitpress/7503.003.0205)

## See also

[`hg_line_graph()`](https://pak.dynasite.org/hypergraphs/reference/hg_line_graph.md)
for the projection onto hyperedges.

## Examples

``` r
meetings <- data.frame(
  person = c("Alice", "Bob", "Carol", "Alice", "Bob", "Dave", "Eve"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m3", "m3"))
meeting_hg <- hypergraph(meetings, node = "person", hyperedge = "meeting")
meeting_network <- pairwise_network(meeting_hg)
hg_get(meeting_network)
#>    from    to weight
#> 1 Alice   Bob      2
#> 2 Alice Carol      1
#> 3   Bob Carol      1
#> 4  Dave   Eve      1
association_network <- pairwise_network(meeting_hg, type = "association")
hg_get(association_network)
#>    from    to weight
#> 1 Alice   Bob    1.5
#> 2 Alice Carol    0.5
#> 3   Bob Carol    0.5
#> 4  Dave   Eve    1.0
```
