# HyperGAT document classifier

Trains the dual-attention hypergraph network of Ding et al. (2020) –
every document becomes its own hypergraph (its unique words as vertices,
its sentences as hyperedges), two attention layers aggregate words into
sentence and optional LDA-topic hyperedges and those hyperedges back
into words, and a masked mean pool feeds a linear classifier. Inductive:
only labeled documents are trained on, and every document (labeled or
not) is scored. Needs the suggested torch package. Semantics follow the
official implementation. `semantic = "none"` reproduces its
sentence-only ("w/o semantic") ablation. `semantic = "lda"` adds the
full paper path: online variational-Bayes LDA is fitted to labeled
documents only, with the topic count defaulting to the number of
classes, and each document receives one edge per topic containing the
topic's top words present in that document. The official repository's
`generate_lda.py` fits its topics on every document, test documents
included; passing its topic words as `lda_keywords` reproduces that
choice.

## Usage

``` r
hg_hypergat(
  x,
  labels,
  column = NULL,
  id = NULL,
  stop_words = stop_words_en(),
  min_count = 1L,
  lowercase = TRUE,
  semantic = c("none", "lda"),
  lda_keywords = NULL,
  lda_topics = NULL,
  lda_top_n = 10L,
  lda_max_iter = 10L,
  lda_batch_size = 128L,
  lda_offset = 50,
  lda_decay = 0.7,
  lda_seed = 0L,
  embed_dim = 300L,
  hidden = 100L,
  epochs = 10L,
  lr = 0.001,
  dropout = 0.3,
  batch_size = 8L,
  weight_decay = 1e-06,
  lr_decay = 0.1,
  lr_step = 3L,
  validation = 0.1,
  class_weights = c("balanced", "none"),
  embeddings = NULL,
  seed = 1L,
  verbose = FALSE,
  what = c("predictions", "attention", "hyperedges", "hyperedge_words"),
  holdout = NULL
)

text_hypergat(
  x,
  labels,
  column = NULL,
  id = NULL,
  stop_words = stop_words_en(),
  min_count = 1L,
  lowercase = TRUE,
  semantic = c("none", "lda"),
  lda_keywords = NULL,
  lda_topics = NULL,
  lda_top_n = 10L,
  lda_max_iter = 10L,
  lda_batch_size = 128L,
  lda_offset = 50,
  lda_decay = 0.7,
  lda_seed = 0L,
  embed_dim = 300L,
  hidden = 100L,
  epochs = 10L,
  lr = 0.001,
  dropout = 0.3,
  batch_size = 8L,
  weight_decay = 1e-06,
  lr_decay = 0.1,
  lr_step = 3L,
  validation = 0.1,
  class_weights = c("balanced", "none"),
  embeddings = NULL,
  seed = 1L,
  verbose = FALSE,
  what = c("predictions", "attention", "hyperedges", "hyperedge_words"),
  holdout = NULL
)
```

## Arguments

- x:

  A character vector of documents (names become ids) or a data.frame
  with a text column.

- labels:

  The known labels: the name of a column of `x` (a data.frame) holding
  each document's label, a named character vector (names are document
  ids, values class labels) or a tidy data.frame with a `node` column
  and a `label`, `cluster` or `predicted` column. At least two classes.

- column, id:

  When `x` is a data.frame: the text column and the optional id column,
  as in
  [`text_hypergraph()`](https://mohsaqr.github.io/hypergraphs/reference/text_hypergraph.md).

- stop_words:

  Words removed before building sentences (default
  [`stop_words_en()`](https://mohsaqr.github.io/hypergraphs/reference/stop_words_en.md)).

- min_count:

  Minimum corpus frequency for a word to become a vertex.

- lowercase:

  Lowercase the text first.

- semantic:

  `"none"` (default) for sentence hyperedges only, or `"lda"` to append
  the paper's semantic topic hyperedges.

- lda_keywords:

  Optional precomputed topic keywords: a list of character vectors, one
  per topic. With `semantic = "lda"`, `NULL` fits LDA natively;
  supplying the list bypasses fitting and reproduces a saved official
  preprocessing run.

- lda_topics:

  Number of LDA topics. `NULL` (default) uses the number of classes, as
  in the paper. Ignored when `lda_keywords` is supplied.

- lda_top_n:

  Number of highest-probability words per fitted topic (paper default
  10).

- lda_max_iter, lda_batch_size:

  Online variational-Bayes passes and minibatch size (scikit-learn
  defaults used by the official code: 10 and 128).

- lda_offset, lda_decay:

  Online learning schedule (official values 50 and 0.7).

- lda_seed:

  LDA initialization seed (official value 0), separate from the neural
  training `seed`.

- embed_dim, hidden:

  Embedding and hidden width (official defaults 300 and 100).

- epochs, lr, dropout, batch_size, weight_decay:

  Training hyperparameters (official defaults: 10, 0.001, 0.3, 8, 1e-6).

- lr_decay, lr_step:

  StepLR schedule: multiply the learning rate by `lr_decay` every
  `lr_step` epochs (official 0.1 every 3).

- validation:

  Fraction of the labeled documents held out (stratified) to pick the
  best epoch; `0` keeps the final weights.

- class_weights:

  `"balanced"` (default, inverse-frequency loss weights as in the
  official run script) or `"none"`.

- embeddings:

  Optional pretrained word-vector matrix (rownames are words,
  `embed_dim` columns); words not covered keep their random
  initialization.

- seed:

  Integer seed (R and torch); results are deterministic given a seed.

- verbose:

  Message the loss each epoch.

- what:

  Legacy extraction option. Fit once and use
  [`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
  with `what = "attention"`, `"hyperedges"` or `"hyperedge_words"`
  instead.

- holdout:

  `NULL` (default) trains on every given label. A share in `(0, 1)`
  hides that share within each class with `seed`. The fitted object
  prints its held-out accuracy and balanced accuracy.

## Value

An `hg_hypergat` fitted classifier, also an
[hg_classification](https://mohsaqr.github.io/hypergraphs/reference/hg_get.hg_classification.md)
and a data frame. [`print()`](https://rdrr.io/r/base/print.html) reports
held-out evaluation when available.
[`hg_get()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.md)
reads predictions, per-class results, confusion, document text and
training history.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows confusion
or training loss.
[`predict.hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/predict.hg_hypergat.md)
classifies new documents with the same network and frozen vocabulary.

Attention is computed on request, without training again. Word weights
sum to one within each hyperedge; edge weights sum to one over each
word's hyperedges. A word in just one edge gives it weight one
automatically. These are internal aggregation weights, not
class-specific contributions or explanations. The diagnostic tables and
their normalization baselines are described in
[`hg_get.hg_hypergat()`](https://mohsaqr.github.io/hypergraphs/reference/hg_get.hg_hypergat.md).
Softmax scores are not calibrated probabilities of correctness. The
vocabulary is estimated from known labels only; held-out and unlabelled
documents cannot change it.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*.

## Examples

``` r
if (FALSE) { # \dontrun{
# articles has text and an existing subject-label column called label.
fit <- hg_hypergat(articles, column = "text", labels = "label",
                   holdout = 0.2)
fit
plot(fit)
hg_get(fit, what = "documents", split = "test", correct = FALSE, top = 1)
predict(fit, newdata = new_articles)
} # }
```
