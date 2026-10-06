# Read a fitted HyperGAT classifier

The classifier retains its trained network, vocabulary and topic
keywords. Reading a table or predicting new documents never trains the
network. Attention tables describe internal aggregation weights. They
are not class-specific contributions or explanations of a prediction.

## Usage

``` r
# S3 method for class 'hg_hypergat'
hg_get(
  x,
  what = NULL,
  split = NULL,
  correct = NULL,
  node = NULL,
  sort_by = NULL,
  top = NULL,
  ...
)

# S3 method for class 'hg_hypergat'
print(x, ...)

# S3 method for class 'hg_hypergat'
summary(object, ...)

# S3 method for class 'hg_hypergat'
plot(
  x,
  y,
  type = c("confusion", "history", "hyperedges", "attention"),
  node = NULL,
  top = 8L,
  ...
)
```

## Arguments

- x:

  A fitted `hg_hypergat`.

- what:

  `NULL` or `"predictions"`, `"accuracy"`, `"classes"`, `"confusion"`,
  `"history"`, `"documents"`, `"attention"`, `"hyperedges"` or
  `"hyperedge_words"`. Evaluation tables use held-out documents only.
  `"documents"` includes the text with its prediction.

- split, correct, sort_by:

  Filters and ordering for predictions and documents, as in
  [`hg_get.hg_classification()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.hg_classification.md).

- node:

  Document ids or a table with a `node` column.

- top:

  Number of rows per class (predictions and documents), per document
  (hyperedges or attention), or per hyperedge (hyperedge words).

- ...:

  Unused.

- object:

  A fitted `hg_hypergat`.

- y:

  Unused.

- type:

  `"confusion"` (default), `"history"`, `"hyperedges"` (the word
  hypergraph of one document), or `"attention"`. Attention plots show
  weights alongside their uniform normalization baseline. These weights
  do not measure a sentence's contribution to the predicted class.

## Value

A plain data frame. The attention tables include uniform normalization
baselines: `word_uniform = 1 / n_words` within an edge, and
`edge_uniform = 1 / n_edges` for a word. Hyperedge summaries include
`uniform_weight` (mean `edge_uniform`) and `exclusive_fraction` (the
fraction of words belonging only to that hyperedge). Such words have
`edge_weight = 1` regardless of the learned parameters.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*, 4927-4936.
