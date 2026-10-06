# Agreement between two labelings of the same nodes

Compares any two tidy labelings – partitions from
[`hg_cluster()`](https://mohsaqr.github.io/hypergraphs/reference/hg_cluster.md),
predictions from
[`hg_classify()`](https://mohsaqr.github.io/hypergraphs/reference/hg_classify.md),
[`hg_neural()`](https://mohsaqr.github.io/hypergraphs/reference/hg_neural.md)
or
[`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
– joined on their shared `node` column. `agreement` is the share of
nodes with literally equal labels (meaningful when both labelings use
the same label set, e.g. a classifier scored against the clustering that
produced its seeds); `ari` is the adjusted Rand index (Hubert & Arabie
1985), which is label-permutation invariant and the right statistic when
the two label sets are arbitrary (e.g. two independent clusterings).
Adjusted mutual information (`ami`) and normalized mutual information
(`nmi`) are also available; the legal-hypergraphs workflow uses AMI to
select the medoid of repeated Infomap partitions.

## Usage

``` r
hg_agreement(
  x,
  y,
  what = c("summary", "table", "mapping"),
  method = "ari",
  node = NULL,
  label = NULL
)
```

## Arguments

- x, y:

  Tidy labelings: data.frames with a `node` (or, in the text family's
  document tables, `doc`) column and a `predicted`, `cluster`,
  `community`, `topic`, `block` or `label` column (first match in that
  order wins), one row per node. Nodes are matched by name; nodes
  present in only one labeling are dropped. A fitted model is read
  through its one-label-per-node table: a fit of
  [`hg_communities()`](https://mohsaqr.github.io/hypergraphs/reference/hg_communities.md)
  on a hypergraph through its medoid partition
  (`hg_get(fit, what = "medoid")`), a
  [`hg_mmsbm()`](https://mohsaqr.github.io/hypergraphs/reference/hg_mmsbm.md)
  fit through each node's strongest community (`what = "nodes"`), a
  [`hg_topics()`](https://mohsaqr.github.io/hypergraphs/reference/hg_topics.md)
  model through each document's dominant topic (`what = "documents"`),
  and a
  [`text_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/text_hypergraph.md)
  through its document table, so a corpus column is scored by naming it
  in `label`.

- what:

  `"summary"` (default) for the one-row comparison, `"table"` for the
  tidy contingency table of the joined labels, or `"mapping"` for one
  row per label of `x` naming the label of `y` that holds most of its
  nodes – how a partition survives a reweighting or a change of method,
  topic by topic.

- method:

  One or more label-permutation-invariant measures: `"ari"` (default),
  `"ami"`, or `"nmi"`. Ignored for `what = "table"`.

- node, label:

  Column names, one name for both labelings or two names for `x` and `y`
  in turn, that override the defaults above –
  `hg_agreement(predictions, corpus, label = c("predicted", "year"))`
  scores a classifier against a column of the corpus. `node = NULL`
  (default) uses `node`, else `doc`; `label = NULL` (default) keeps the
  `predicted` / `cluster` / `community` / `topic` / `block` / `label`
  lookup.

## Value

A base `data.frame`. For `what = "summary"`: one row with columns `n`
(nodes compared), `agreement` (share of equal labels), `aligned` (nodes
that stay with the majority of their `x` label in `y` – the sum of
`overlap` over `what = "mapping"`, so a label-name-free count of how
many nodes a re-fit keeps together) and the requested measure columns.
For `what = "table"`: one row per label pair with columns `label_x`,
`label_y` and `n`. For `what = "mapping"`: one row per label of `x` with
columns `label_x`, `n` (its nodes), `label_y` (the label of `y` holding
most of them), `overlap` (how many) and `share` (`overlap / n`), in the
natural order of `label_x`.

## References

Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
Classification*, 2, 193–218.

Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic
measures for clusterings comparison. *Journal of Machine Learning
Research*, 11, 2837–2854.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
topics <- hg_cluster(hg, k = 2, seed = 1)
fit <- hg_classify(hg, labels = c(cooking_1 = "Cluster 1",
                                  space_1 = "Cluster 2"))
hg_agreement(fit, topics)
#>   n agreement aligned ari
#> 1 4         1       4   1
hg_agreement(fit, topics, what = "table")
#>     label_x   label_y n
#> 1 Cluster 1 Cluster 1 2
#> 2 Cluster 2 Cluster 2 2
hg_agreement(fit, topics, what = "mapping")
#>     label_x n   label_y overlap share
#> 1 Cluster 1 2 Cluster 1       2     1
#> 2 Cluster 2 2 Cluster 2       2     1
# a labeling read from any table: name its label column
known <- data.frame(doc = c("cooking_1", "space_2"), theme = c("cooking", "space"))
hg_agreement(topics, known, label = c("cluster", "theme"), what = "table")
#>     label_x label_y n
#> 1 Cluster 1 cooking 1
#> 2 Cluster 2   space 1
```
