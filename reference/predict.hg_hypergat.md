# Predict documents with a fitted HyperGAT classifier

Reuses the trained network, vocabulary, preprocessing and topic
keywords. Words absent from the fitted vocabulary are ignored. A
document with no remaining words receives missing predictions and a
warning, preserving its row and id. The winning softmax probability is a
model score, not a calibrated probability of correctness.

## Usage

``` r
# S3 method for class 'hg_hypergat'
predict(
  object,
  newdata = NULL,
  column = NULL,
  id = NULL,
  labels = NULL,
  type = c("predictions", "probabilities"),
  ...
)
```

## Arguments

- object:

  A fitted `hg_hypergat`.

- newdata:

  A character vector or data frame of new documents. `NULL` returns the
  original predictions.

- column, id:

  Text and id columns of a data frame. Defaults to the columns used when
  fitting. An unnamed character vector receives sequential document ids.

- labels:

  Optional observed labels for `newdata`, in the same formats as
  [`hg_hypergat()`](https://pak.dynasite.org/hypergraphs/reference/hg_hypergat.md).
  Supplying them evaluates the predictions through
  [`hg_get.hg_classification()`](https://pak.dynasite.org/hypergraphs/reference/hg_get.hg_classification.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html). These
  labels never update the model. Unscorable labelled documents count as
  errors.

- type:

  `"predictions"` (default), a document table, or `"probabilities"`, a
  matrix with one column per fitted class.

- ...:

  Unused.

## Value

A plain data frame (`node`, `label`, `predicted`, `score`, `margin`) or
a probability matrix, in input order. No fitting occurs.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*, 4927-4936.
