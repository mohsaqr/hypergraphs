# Relations between topics (deprecated)

`hg_relations()` is the former name of
[`hg_network()`](https://mohsaqr.github.io/hypergraphs/reference/hg_network.md)
for a partition and returns the same result. It warns with the class
`hypergraphs_deprecated`.

## Usage

``` r
hg_relations(
  hg,
  clusters,
  similarity = c("none", "association", "cosine", "jaccard", "inclusion", "equivalence"),
  what = c("edges", "network")
)
```

## Arguments

- hg:

  The document hypergraph the topics were computed on.

- clusters:

  A partition: the tidy table returned by
  [`hg_cluster()`](https://mohsaqr.github.io/hypergraphs/reference/hg_cluster.md)
  (columns `node`, `cluster`), or a named vector of cluster labels.

- similarity:

  `"none"` (default, the raw count), `"association"`, `"cosine"`,
  `"jaccard"`, `"inclusion"` or `"equivalence"`. Count networks only.

- what:

  `"edges"` (default) for the edge list, `"network"` for a
  `cograph_network` built from it with
  [`cograph::as_cograph()`](https://sonsoles.me/cograph/reference/as_cograph.html),
  whose node table carries each topic's `size`.

## Value

The value of
`hg_network(hg, clusters = clusters, similarity =, what =)`.
