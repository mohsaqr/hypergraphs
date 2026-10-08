# Normalized hypergraph Laplacian

Computes the normalized Laplacian of a hypergraph, either the classic
Zhou-Huang-Scholkopf form on the binary incidence pattern
(`type = "zhou"`) or the random-walk form with edge-dependent vertex
weights (`type = "random_walk"`), in which the weighted incidence cells
(e.g. the summed weights produced by
[`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md))
determine where a random walker lands inside a hyperedge, and the
resulting non-reversible walk is symmetrized through its stationary
distribution (Chung 2005). Both Laplacians are symmetric positive
semi-definite with eigenvalues in `[0, 2]`; for a binary incidence with
unit hyperedge weights they coincide.

## Usage

``` r
hg_laplacian(hg, type = c("zhou", "random_walk"), edge_weights = NULL)
```

## Arguments

- hg:

  A `net_hg` from
  [`network_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/network_hypergraph.md)
  or
  [`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md),
  dense or sparse (`sparse = TRUE`). Must be connected and have at least
  one hyperedge. Empty hyperedges (which
  [`random_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/random_hypergraph.md)
  can produce) contribute nothing.

- type:

  Character. `"zhou"` (default) for the Zhou et al. (2006) normalized
  Laplacian on the binary incidence pattern, or `"random_walk"` for the
  Hayashi et al. (2020) EDVW random-walk Laplacian on the weighted
  incidence.

- edge_weights:

  A single positive number (recycled) or a numeric vector of positive
  hyperedge weights (length `hg$n_hyperedges`), or `NULL` for the
  default. Hypergraphs built by
  [`window_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/window_hypergraph.md)
  default to their window counts (for both types). Otherwise the default
  is type-specific: unit weights for `"zhou"`; for `"random_walk"` the
  Hayashi et al. heuristic - the population standard deviation of each
  hyperedge's (non-zero) vertex weights plus one - which reduces to unit
  weights on a binary incidence.

## Value

A symmetric `n_nodes` x `n_nodes` numeric matrix (node names as
dimnames) with attributes `type` (the Laplacian type), `pi` (named
stationary distribution of the underlying random walk) and
`edge_weights` (the hyperedge weights, one per hyperedge; an empty
hyperedge's default weight is 1 and unused). For a sparse hypergraph the
matrix is a sparse symmetric `Matrix` (class `dsCMatrix`) with the same
dimnames and attributes, and the random-walk stationary distribution is
found by power iteration rather than a dense eigendecomposition. A
disconnected hypergraph raises `hypergraphs_hypergraph_disconnected`.

## References

Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *NeurIPS 19*.
[doi:10.7551/mitpress/7503.003.0205](https://doi.org/10.7551/mitpress/7503.003.0205)

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*, 495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Chung, F. (2005). Laplacians and the Cheeger inequality for directed
graphs. *Annals of Combinatorics*, 9(1), 1-19.

## Examples

``` r
events <- data.frame(
  person = c("a", "b", "c", "a", "b", "d", "c", "d", "e", "e", "a"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
              "m4", "m4"),
  hours = c(2, 1, 1, 3, 2, 1, 2, 2, 4, 1, 1)
)
hg <- group_hypergraph(events, node = "person", hyperedge = "meeting",
                       weight = "hours")
L <- hg_laplacian(hg, type = "random_walk")
range(eigen(L, symmetric = TRUE, only.values = TRUE)$values)
#> [1] -2.899699e-16  1.001279e+00
```
