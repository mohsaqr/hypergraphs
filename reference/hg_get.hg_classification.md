# Read a classification result

[`hg_classify()`](https://mohsaqr.github.io/hypergraphs/reference/hg_classify.md)
and
[`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
with `holdout` return an `hg_classification`: one row per document with
its true label, its predicted label, its split and whether a held-out
prediction is correct.
[`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
reads the evaluation tables, and for
[`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
the diagnostic attention of the trained network, computed on request.

## Usage

``` r
# S3 method for class 'hg_classification'
hg_get(
  x,
  what = c("predictions", "accuracy", "classes", "confusion", "hyperedges",
    "hyperedge_words", "documents"),
  split = NULL,
  correct = NULL,
  node = NULL,
  sort_by = NULL,
  top = NULL,
  ...
)

# S3 method for class 'hg_classification'
print(x, ...)

# S3 method for class 'hg_classification'
summary(object, ...)

# S3 method for class 'hg_classification'
plot(x, y, type = c("confusion", "hyperedges"), node = NULL, top = 8L, ...)
```

## Arguments

- x:

  An `hg_classification`.

- what:

  `"predictions"` (default), the per-document table: `node`, `label`
  (the true label, `NA` for an unlabelled document), `predicted`,
  `score`, `margin`, `split` (`"train"`, `"test"` or `"unlabelled"`) and
  `correct` (for the test documents); `"accuracy"`, one row: `n_test`,
  `accuracy`, `balanced_accuracy` (the mean recall over the classes) and
  `chance` (one over the number of classes, the balanced accuracy of a
  random guess); `"classes"`, one row per class: `class`, `n_test`,
  `correct`, `recall` (the share of the class's test documents predicted
  as the class) and `precision` (the share of the test documents
  predicted as the class that belong to it); `"confusion"`, one row per
  true and predicted class: `label`, `predicted`, `n`. For a result of
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
  only: `"hyperedges"`, one row per document and hyperedge (sentence or
  topic): `node`, `label`, `predicted`, `hyperedge`, `kind`, `text`,
  `weight` (the mean edge-level attention of the hyperedge's words),
  `top_word` and `n_words`, the hyperedges of each document from the
  heaviest; or `"hyperedge_words"`, one row per document, hyperedge and
  word: `node`, `label`, `predicted`, `hyperedge`, `kind`, `text`,
  `word`, `word_weight` (the node-level attention of the word within the
  hyperedge) and `edge_weight` (the edge-level attention of the
  hyperedge for the word), the quantities of Figure 5 of Ding et al.
  (2020). These internal weights are not prediction explanations. Fitted
  HyperGAT models also support `"documents"`, the prediction and text of
  each document.

- split:

  For `"predictions"` or `"documents"`: keep only the `"train"`, the
  `"test"` or the `"unlabelled"` documents, the last being those without
  a label, whose predictions are the classifier's output.

- correct:

  For `"predictions"` or `"documents"`: `TRUE` keeps the test documents
  predicted correctly, `FALSE` those predicted wrongly.

- node:

  For `"predictions"`, `"hyperedges"` and `"hyperedge_words"`: keep only
  these documents, given as ids or as a table with a `node` column, such
  as a table
  [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  returned. For `"hyperedge_words"`, a table that also has a `hyperedge`
  column keeps only those hyperedges.

- sort_by:

  For `"predictions"`: `"margin"` or `"score"`, from the largest.

- top:

  For `"predictions"`: the first `top` rows of each true label, so
  `split = "test", correct = TRUE, sort_by = "margin", top = 1` gives
  the most confident correct prediction of each class. For
  `"hyperedges"`: the `top` heaviest hyperedges of each document. For
  `"hyperedge_words"`: the `top` words of largest `word_weight` in each
  hyperedge.

- ...:

  Unused.

- object:

  An `hg_classification`.

- y:

  Unused.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html):
  `"confusion"` (default), the counts of true against predicted labels
  on the test documents, or `"hyperedges"` (a result of
  [`hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_hypergat.md)
  only), the `top` heaviest hyperedges of each document in `node`, with
  its uniform normalization baseline. Without `node`, the first held-out
  document of each class is selected (the first training document if
  there is no holdout).

## Value

A base `data.frame` as described under `what`.
[`print()`](https://rdrr.io/r/base/print.html) shows the held-out
accuracy and the per-class table and returns `x` invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the per-class
table; [`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns
a ggplot of the confusion table (`type = "confusion"`) or of the
heaviest hyperedges of some documents (`type = "hyperedges"`).

## References

Brodersen, K. H., Ong, C. S., Stephan, K. E., & Buhmann, J. M. (2010).
The balanced accuracy and its posterior distribution. *Proceedings of
the 20th International Conference on Pattern Recognition*, 3121-3124.
[doi:10.1109/ICPR.2010.764](https://doi.org/10.1109/ICPR.2010.764)

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*, 4927-4936.

## Examples

``` r
if (FALSE) { # \dontrun{
# articles contains text and existing subject labels.
fit <- hg_hypergat(articles, column = "text", labels = "subject",
                   holdout = 0.2)
fit
hg_get(fit, what = "classes")
hg_get(fit, what = "documents", split = "test", correct = FALSE, top = 1)
} # }
```
