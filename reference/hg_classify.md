# Transductive label spreading on a hypergraph

Semi-supervised classification of hypergraph nodes by the regularization
framework of Zhou et al. (2006): given labels for a subset of nodes, the
scores `F = (1 - xi) * (I - xi * S)^{-1} Y` spread the labels over the
hypergraph, where `S = I - L` is the normalized similarity operator of
the chosen Laplacian and `Y` is the label indicator matrix. Each node is
assigned the class with the highest score. This is the non-neural
ancestor of hypergraph-attention text classifiers: with documents as
hyperedges over words (or vice versa) it classifies unlabeled nodes from
a handful of labeled ones.

## Usage

``` r
hg_classify(
  hg,
  labels,
  xi = 0.99,
  type = c("zhou", "random_walk"),
  normalization = c("none", "class_mass"),
  edge_weights = NULL,
  holdout = NULL,
  seed = 1L
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/text_hypergraph.md)
  (or any hypergraphs `net_hg`).

- labels:

  The known labels: a named character vector (names are node identifiers
  – documents under `nodes = "doc"` – values their class labels), a tidy
  data.frame with a `node` column and a `label`, `cluster` or
  `predicted` column, or the name of a column of the hypergraph's
  document table, which a
  [`text_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/text_hypergraph.md)
  fills with the input's other columns (`labels = "period"`).

- xi:

  Numeric in `(0, 1)`. Spreading coefficient (default `0.99`); larger
  values weight the hypergraph structure more relative to the initial
  labels.

- type, edge_weights:

  Passed to
  [`hg_laplacian()`](https://mohsaqr.github.io/hypergraphs/reference/hg_laplacian.md).

- normalization:

  Decision rule for turning spread scores into predictions: `"none"`
  (default, the raw Zhou 2006 argmax) or `"class_mass"` (class-mass
  normalization, Zhu et al. 2003). Use `"class_mass"` when the labeled
  seeds are class-imbalanced – the raw rule can collapse every
  prediction onto the majority class.

- holdout:

  `NULL` (default) uses every given label and returns the predictions. A
  share in `(0, 1)` hides that share of the labels, drawn within each
  class, predicts them from the rest, and returns an
  [hg_classification](https://mohsaqr.github.io/hypergraphs/reference/hg_get.hg_classification.md)
  that prints the held-out accuracy and balanced accuracy.

- seed:

  Seed of the held-out draw (default `1`); the caller's random-number
  stream is restored.

## Value

A base `data.frame`, one row per node, with columns `node`, `label` (the
given label or `NA`), `predicted`, `score`, and `margin`. With
`holdout`, an `hg_classification`: the same table with `label` holding
the true label, `split` (`"train"` or `"test"`) and `correct` for the
held-out documents; read its evaluation with
[`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
(`what = "accuracy"`, `"classes"`, `"confusion"`).

## Details

The result is one row per node with its given and predicted label.

## References

Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *NeurIPS 19*.

Zhu, X., Ghahramani, Z., & Lafferty, J. (2003). Semi-supervised learning
using Gaussian fields and harmonic functions. *ICML 20*.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
hg_classify(hg, labels = c(cooking_1 = "cooking", space_1 = "space"))
#>        node   label predicted     score      margin
#> 1 cooking_1 cooking   cooking 0.2665654 0.082276891
#> 2 cooking_2    <NA>   cooking 0.2538871 0.009941802
#> 3   space_1   space     space 0.3016056 0.117317071
#> 4   space_2    <NA>     space 0.2663331 0.072737125
```
