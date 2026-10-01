# Network of topics

Builds the topic-by-topic network of a document hypergraph, one node per
topic and one undirected weighted edge per pair of related topics. The
topics come either from a partition of the documents or from a
mixed-membership topic model.

## Usage

``` r
hg_network(
  hg,
  clusters = NULL,
  topics = NULL,
  threshold = NULL,
  cutoff = 0.01,
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

- topics:

  A mixed-membership topic model of `hg` fitted by
  [`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md).
  Give `clusters` or `topics`.

- threshold:

  Topic model only: the share at which a topic counts as present in a
  document, one number in (0, 1\]. `NULL` (default) relates topics by
  the correlation of their shares instead.

- cutoff:

  Topic model without `threshold` only: the correlation a pair must
  exceed to be kept (default `0.01`, as in `stm::topicCorr()`).

- similarity:

  `"none"` (default, the raw count), `"association"`, `"cosine"`,
  `"jaccard"`, `"inclusion"` or `"equivalence"`. Count networks only.

- what:

  `"edges"` (default) for the edge list, `"network"` for a
  `cograph_network` built from it with
  [`cograph::as_cograph()`](https://sonsoles.me/cograph/reference/as_cograph.html),
  whose node table carries each topic's `size`.

## Value

For `what = "edges"`: a base `data.frame`, one row per pair of topics
with a positive weight, columns `source`, `target`, `weight`, pairs in
the topics' natural order. For `what = "network"`: a `cograph_network`
with one node per topic (`label`, `name`, `size`) and one undirected
weighted edge per pair; `size` is the number of documents of a cluster,
the documents in which a topic is present, or (without `threshold`) a
topic's expected number of documents. Raises `hypergraphs_bad_input` for
unknown node names, for neither or both of `clusters` and `topics`, for
a topic model not fitted on `hg`, for an invalid `threshold` or
`cutoff`, and for a similarity measure on a correlation network.

## Details

For a partition (`clusters`), the network is the co-occurrence network
of the bibliometric kind: the strength between two topics is the sum,
over all words, of the product of the number of documents in each topic
that contain the word (full counting, the aggregation bibnets uses for
keyword co-occurrence).

For a topic model fitted by
[`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md)
(`topics`), two topics are related when the same documents draw on both.
Without `threshold`, the weight is the correlation of the two topics'
shares over the documents, and the pairs with a correlation above
`cutoff` are kept, the "simple" topic correlation of `stm::topicCorr()`
(Roberts et al. 2019). With `threshold`, a topic counts as present in a
document when its share is at least `threshold`, and the weight of a
pair is the number of documents in which both are present (Abuhay et al.
2017; Cassi et al. 2017).

A count network, from a partition or from a thresholded topic model, can
be normalised by the similarity measures of bibnets' `normalize()`, with
each topic's total on the diagonal (its documents with the word, or its
documents with the topic present). With `D_i` the diagonal entry of
topic i and `A_ij` the raw count: `"association"` is `A_ij / (D_i D_j)`,
`"cosine"` is `A_ij / sqrt(D_i D_j)`, `"jaccard"` is
`A_ij / (D_i + D_j - A_ij)`, `"inclusion"` is `A_ij / min(D_i, D_j)` and
`"equivalence"` is `A_ij^2 / (D_i D_j)` (van Eck & Waltman 2009).

The pairwise network of topics is the projection of the topic
combinations that
[`group_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/group_hypergraph.md)
builds from the same topic model, in which every document binds all the
topics it contains at once.

## References

van Eck, N. J., & Waltman, L. (2009). How to normalize cooccurrence
data? An analysis of some well-known similarity measures. *Journal of
the American Society for Information Science and Technology*, 60(8),
1635–1651.

Roberts, M. E., Stewart, B. M., & Tingley, D. (2019). stm: An R package
for structural topic models. *Journal of Statistical Software*, 91(2),
1-40. [doi:10.18637/jss.v091.i02](https://doi.org/10.18637/jss.v091.i02)

Abuhay, T. M., Kovalchuk, S. V., Bochenina, K., Kampis, G.,
Krzhizhanovskaya, V. V., & Lees, M. H. (2017). Analysis of computational
science papers from ICCS 2001-2016 using topic modeling and graph
theory. *Procedia Computer Science*, 108, 7-17.
[doi:10.1016/j.procs.2017.05.183](https://doi.org/10.1016/j.procs.2017.05.183)

Cassi, L., Lahatte, A., Rafols, I., Sautier, P., & de Turckheim, E.
(2017). Improving fitness: Mapping research priorities against societal
needs on obesity. *Journal of Informetrics*, 11(4), 1095-1113.
[doi:10.1016/j.joi.2017.09.010](https://doi.org/10.1016/j.joi.2017.09.010)

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
topics <- hg_cluster(hg, k = 2, seed = 1)
hg_network(hg, clusters = topics)
#>      source    target weight
#> 1 Cluster 1 Cluster 2      1
hg_network(hg, clusters = topics, similarity = "cosine")
#>      source    target     weight
#> 1 Cluster 1 Cluster 2 0.07715167

corpus <- c(
  a = "soup salt onion soup broth", b = "salt soup broth onion",
  c = "stars sky moon night", d = "sky stars night moon moon",
  e = "soup stars salt sky night broth", f = "onion salt stars broth")
corpus_hg <- text_hypergraph(corpus)
model <- hg_topics(corpus_hg, k = 2, nstart = 2)
hg_network(corpus_hg, topics = model, threshold = 0.2)
#>    source  target weight
#> 1 Topic 1 Topic 2      1
```
