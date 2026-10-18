# Read a classification result

[`hg_classify()`](https://pak.dynasite.org/hypergraphs/reference/hg_classify.md)
with `holdout` returns an `hg_classification`: one row per document with
its true label, its predicted label, its split and whether a held-out
prediction is correct.
[`hg_get()`](https://pak.dynasite.org/hypergraphs/reference/hg_get.md)
reads the evaluation tables.

## Usage

``` r
# S3 method for class 'hg_classification'
hg_get(
  x,
  what = c("predictions", "accuracy", "classes", "confusion"),
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
plot(x, y, ...)
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
  true and predicted class: `label`, `predicted`, `n`.

- split:

  For `"predictions"`: keep only the `"train"`, the `"test"` or the
  `"unlabelled"` documents, the last being those without a label, whose
  predictions are the classifier's output.

- correct:

  For `"predictions"`: `TRUE` keeps the test documents predicted
  correctly, `FALSE` those predicted wrongly.

- node:

  For `"predictions"`: keep only these documents, given as ids or as a
  table with a `node` column, such as a table
  [`hg_get()`](https://pak.dynasite.org/hypergraphs/reference/hg_get.md)
  returned.

- sort_by:

  For `"predictions"`: `"margin"` or `"score"`, from the largest.

- top:

  For `"predictions"`: the first `top` rows of each true label, so
  `split = "test", correct = TRUE, sort_by = "margin", top = 1` gives
  the most confident correct prediction of each class.

- ...:

  Unused.

- object:

  An `hg_classification`.

- y:

  Unused.

## Value

A base `data.frame` as described under `what`.
[`print()`](https://rdrr.io/r/base/print.html) shows the held-out
accuracy and the per-class table and returns `x` invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the per-class
table; [`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns
a ggplot of the confusion table.

## Details

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
counts of true against predicted labels on the test documents.

## References

Brodersen, K. H., Ong, C. S., Stephan, K. E., & Buhmann, J. M. (2010).
The balanced accuracy and its posterior distribution. *Proceedings of
the 20th International Conference on Pattern Recognition*, 3121-3124.
[doi:10.1109/ICPR.2010.764](https://doi.org/10.1109/ICPR.2010.764)

## Examples

``` r
hg <- text_hypergraph(head(forum_posts, 60), column = "text")
fit <- hg_classify(hg, labels = "topic", holdout = 0.25)
fit
#> Classification (hg_classify): 60 documents, 15 held out (25% of the labelled)
#> Held-out accuracy 0.933, balanced accuracy 0.917 (chance 0.333)
#> 
#>      class n_test correct recall precision
#>  astronomy      6       6   1.00 0.8571429
#>    cooking      4       3   0.75 1.0000000
#>  gardening      5       5   1.00 1.0000000
hg_get(fit, what = "classes")
#>       class n_test correct recall precision
#> 1 astronomy      6       6   1.00 0.8571429
#> 2   cooking      4       3   0.75 1.0000000
#> 3 gardening      5       5   1.00 1.0000000
hg_get(fit, what = "confusion")
#>       label predicted n
#> 1 astronomy astronomy 6
#> 2   cooking astronomy 1
#> 3 gardening astronomy 0
#> 4 astronomy   cooking 0
#> 5   cooking   cooking 3
#> 6 gardening   cooking 0
#> 7 astronomy gardening 0
#> 8   cooking gardening 0
#> 9 gardening gardening 5
hg_get(fit, split = "test", correct = FALSE)
#>    node   label predicted     score      margin split correct
#> 1 doc_3 cooking astronomy 0.2279617 0.001093425  test   FALSE
```
