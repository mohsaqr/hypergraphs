# Classification results with a held-out evaluation. hg_classify() and
# hg_hypergat() with `holdout` hide a stratified share of the known labels,
# predict them, and return the per-document table with `split` and `correct`
# columns as an `hg_classification`, a data.frame whose print shows the
# evaluation and whose tables are read with hg_get().

# Labels given as the name of a column of the hypergraph's document table
# (a text hypergraph carries the input's other columns there) become a named
# vector over the documents; any other `labels` passes through unchanged.
.thg_resolve_labels <- function(hg, labels) {
  if (!(is.character(labels) && length(labels) == 1L && is.null(names(labels)))) {
    return(.thg_labels_input(labels))
  }
  documents <- if (is.list(hg$text)) hg$text$documents else NULL
  if (is.null(documents) || !labels %in% names(documents)) {
    .thg_bad_input(sprintf(paste0("`labels` = \"%s\" must name a column of the ",
                                  "hypergraph's document table"), labels))
  }
  values <- as.character(documents[[labels]])
  keep <- !is.na(values) & documents$doc %in% hg$nodes
  stats::setNames(values[keep], documents$doc[keep])
}

# Labels given as the name of a column of a data.frame `x` (the input of
# hg_hypergat()), over the document ids `ids`.
.thg_resolve_text_labels <- function(x, labels, ids) {
  if (!(is.character(labels) && length(labels) == 1L && is.null(names(labels)))) {
    return(.thg_labels_input(labels))
  }
  if (!is.data.frame(x) || !labels %in% names(x)) {
    .thg_bad_input(sprintf("`labels` = \"%s\" must name a column of `x`",
                           labels))
  }
  values <- as.character(x[[labels]])
  stats::setNames(values[!is.na(values)], ids[!is.na(values)])
}

# A stratified held-out share of the known labels: within each class,
# round(holdout * n) labels are hidden, at least one and leaving at least one
# known. The caller's random-number stream is restored.
.thg_holdout_split <- function(labels, holdout, seed) {
  stopifnot(
    "`holdout` must be a single number in (0, 1)" =
      is.numeric(holdout) && length(holdout) == 1L && holdout > 0 &&
      holdout < 1,
    "`seed` must be a single integer" = length(seed) == 1L && is.finite(seed)
  )
  if (length(unique(labels)) < 2L) {
    .thg_bad_input("`labels` must hold at least two classes to hold some out")
  }
  had_seed <- exists(".Random.seed", envir = globalenv())
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) {
    assign(".Random.seed", old_seed, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv())) {
    rm(".Random.seed", envir = globalenv())
  }, add = TRUE, after = FALSE)
  set.seed(seed)
  hidden <- unlist(lapply(split(names(labels), labels), \(nodes) {
    n_hide <- min(max(round(holdout * length(nodes)), 1L),
                  length(nodes) - 1L)
    if (n_hide < 1L) return(character())
    sample(nodes, n_hide)
  }), use.names = FALSE)
  list(known = labels[!names(labels) %in% hidden],
       hidden = labels[names(labels) %in% hidden])
}

# Assemble an hg_classification from the predictions of a fit on the known
# labels: the true label of every labelled document, its split and whether
# a held-out prediction is correct.
.thg_classification <- function(predictions, labels, hidden, method,
                                holdout) {
  out <- predictions
  out$label <- unname(labels[as.character(out$node)])
  out$split <- ifelse(out$node %in% names(hidden), "test",
                      ifelse(is.na(out$label), "unlabelled", "train"))
  out$correct <- ifelse(out$split %in% "test", !is.na(out$predicted) & out$predicted == out$label, NA)
  rownames(out) <- NULL
  structure(out, class = c("hg_classification", "data.frame"),
            method = method, holdout = holdout)
}

#' Read a classification result
#'
#' `hg_classify()` and `hg_hypergat()` with `holdout` return an
#' `hg_classification`: one row per document with its true label, its
#' predicted label, its split and whether a held-out prediction is correct.
#' `hg_get()` reads the evaluation tables, and for `hg_hypergat()` the
#' diagnostic attention of the trained network, computed on request.
#'
#' @param x An `hg_classification`.
#' @param what `"predictions"` (default), the per-document table: `node`,
#'   `label` (the true label, `NA` for an unlabelled document),
#'   `predicted`, `score`, `margin`, `split` (`"train"`, `"test"` or
#'   `"unlabelled"`) and
#'   `correct` (for the test documents); `"accuracy"`, one row: `n_test`,
#'   `accuracy`, `balanced_accuracy` (the mean recall over the classes) and
#'   `chance` (one over the number of classes, the balanced accuracy of a
#'   random guess); `"classes"`, one row per class: `class`, `n_test`,
#'   `correct`, `recall` (the share of the class's test documents predicted
#'   as the class) and `precision` (the share of the test documents
#'   predicted as the class that belong to it); `"confusion"`, one row per
#'   true and predicted class: `label`, `predicted`, `n`. For a result of
#'   `hg_hypergat()` only: `"hyperedges"`, one row per document and
#'   hyperedge (sentence or topic): `node`, `label`, `predicted`,
#'   `hyperedge`, `kind`, `text`, `weight` (the mean edge-level attention
#'   of the hyperedge's words), `top_word` and `n_words`, the hyperedges of
#'   each document from the heaviest; or `"hyperedge_words"`, one row per
#'   document, hyperedge and word: `node`, `label`, `predicted`,
#'   `hyperedge`, `kind`, `text`, `word`, `word_weight` (the node-level
#'   attention of the word within the hyperedge) and `edge_weight` (the
#'   edge-level attention of the hyperedge for the word), the quantities of
#'   Figure 5 of Ding et al. (2020). These internal weights are not
#'   prediction explanations. Fitted HyperGAT models also support
#'   `"documents"`, the prediction and text of each document.
#' @param split For `"predictions"` or `"documents"`: keep only the `"train"`, the `"test"`
#'   or the `"unlabelled"` documents, the last being those without a label,
#'   whose predictions are the classifier's output.
#' @param correct For `"predictions"` or `"documents"`: `TRUE` keeps the test documents
#'   predicted correctly, `FALSE` those predicted wrongly.
#' @param node For `"predictions"`, `"hyperedges"` and `"hyperedge_words"`:
#'   keep only these documents, given as ids or as a table with a `node`
#'   column, such as a table `hg_get()` returned. For `"hyperedge_words"`, a
#'   table that also has a `hyperedge` column keeps only those hyperedges.
#' @param sort_by For `"predictions"`: `"margin"` or `"score"`, from the
#'   largest.
#' @param top For `"predictions"`: the first `top` rows of each true label,
#'   so `split = "test", correct = TRUE, sort_by = "margin", top = 1` gives
#'   the most confident correct prediction of each class. For
#'   `"hyperedges"`: the `top` heaviest hyperedges of each document. For
#'   `"hyperedge_words"`: the `top` words of largest `word_weight` in each
#'   hyperedge.
#' @param ... Unused.
#' @return A base `data.frame` as described under `what`. `print()` shows
#'   the held-out accuracy and the per-class table and returns `x`
#'   invisibly; `summary()` returns the per-class table; `plot()` returns a
#'   ggplot of the confusion table (`type = "confusion"`) or of the heaviest
#'   hyperedges of some documents (`type = "hyperedges"`).
#' @references
#' Brodersen, K. H., Ong, C. S., Stephan, K. E., & Buhmann, J. M. (2010).
#' The balanced accuracy and its posterior distribution. \emph{Proceedings
#' of the 20th International Conference on Pattern Recognition}, 3121-3124.
#' \doi{10.1109/ICPR.2010.764}
#'
#' Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with
#' less: Hypergraph attention networks for inductive text classification.
#' \emph{EMNLP 2020}, 4927-4936.
#' @examples
#' \dontrun{
#' # articles contains text and existing subject labels.
#' fit <- hg_hypergat(articles, column = "text", labels = "subject",
#'                    holdout = 0.2)
#' fit
#' hg_get(fit, what = "classes")
#' hg_get(fit, what = "documents", split = "test", correct = FALSE, top = 1)
#' }
#' @export
hg_get.hg_classification <- function(x, what = c("predictions", "accuracy",
                                                 "classes", "confusion",
                                                 "hyperedges",
                                                 "hyperedge_words",
                                                 "documents"),
                                     split = NULL, correct = NULL,
                                     node = NULL, sort_by = NULL, top = NULL,
                                     ...) {
  what <- .ho_match_what(what)
  plain <- .ho_plain(x)
  attributes(plain) <- attributes(plain)[c("names", "row.names", "class")]
  filters <- list(split = split, correct = correct, sort_by = sort_by)
  filters <- names(filters)[!vapply(filters, is.null, logical(1L))]
  if (!what %in% c("predictions", "documents") && length(filters) > 0L) {
    .thg_bad_input(sprintf("`%s` applies only to `what = \"predictions\"`",
                           filters[[1L]]))
  }
  if (what %in% c("hyperedges", "hyperedge_words")) {
    return(.thg_classification_attention(x, plain, what, node, top))
  }
  if (what %in% c("predictions", "documents")) {
    out <- .thg_classification_predictions(plain, split, correct, node,
                                           sort_by, top)
    if (what == "documents") {
      documents <- attr(x, "documents")
      if (is.null(documents)) .thg_bad_input("this result has no document text")
      out$text <- documents$text[match(out$node, documents$node)]
    }
    return(out)
  }
  if (!is.null(node) || !is.null(top)) {
    .thg_bad_input(sprintf("`node` and `top` do not apply to `what = \"%s\"`",
                           what))
  }
  test <- plain[plain$split %in% "test", , drop = FALSE]
  classes <- attr(x, "classes") %||%
    sort(unique(stats::na.omit(plain$label)))
  if (identical(what, "confusion")) {
    predicted_classes <- c(classes, if (anyNA(test$predicted)) "(unscored)")
    predicted <- ifelse(is.na(test$predicted), "(unscored)", test$predicted)
    grid <- expand.grid(label = classes, predicted = predicted_classes,
                        stringsAsFactors = FALSE)
    counts <- table(factor(test$label, classes),
                    factor(predicted, predicted_classes))
    grid$n <- as.integer(counts[cbind(grid$label, grid$predicted)])
    return(grid)
  }
  per_class <- data.frame(
    class = classes,
    n_test = vapply(classes, \(k) sum(test$label == k), integer(1L)),
    correct = vapply(classes, \(k) sum(test$label == k & test$correct),
                     integer(1L)),
    stringsAsFactors = FALSE
  )
  per_class$recall <- ifelse(per_class$n_test > 0,
                             per_class$correct / per_class$n_test, NA_real_)
  predicted_n <- vapply(classes, \(k) sum(test$predicted == k, na.rm = TRUE),
                        integer(1L))
  per_class$precision <- ifelse(predicted_n > 0, per_class$correct /
                                  predicted_n, NA_real_)
  rownames(per_class) <- NULL
  if (identical(what, "classes")) return(per_class)
  data.frame(n_test = nrow(test),
             accuracy = if (nrow(test)) mean(test$correct) else NA_real_,
             balanced_accuracy = if (nrow(test))
               mean(per_class$recall, na.rm = TRUE) else NA_real_,
             chance = 1 / length(classes))
}

# Documents named by `node`: ids, or a table with a `node` column.
.thg_node_filter <- function(node) {
  ids <- if (is.data.frame(node)) node$node else node
  if (!is.character(ids) || length(ids) == 0L || anyNA(ids)) {
    .thg_bad_input(paste0("`node` must be document ids or a table with a ",
                          "`node` column"))
  }
  unique(ids)
}

# The first `top` rows of each group, in the current order.
.thg_top_by <- function(x, group, top) {
  if (is.null(top)) return(x)
  stopifnot(
    "`top` must be a single whole number >= 1" =
      is.numeric(top) && length(top) == 1L && is.finite(top) &&
      top == round(top) && top >= 1
  )
  key <- ifelse(is.na(group), "\r", group)
  rank <- stats::ave(seq_along(key), key, FUN = seq_along)
  out <- x[rank <= top, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.thg_classification_predictions <- function(plain, split, correct, node,
                                            sort_by, top) {
  out <- plain
  if (!is.null(split)) {
    split <- match.arg(split, c("train", "test", "unlabelled"))
    out <- out[out$split %in% split, , drop = FALSE]
  }
  if (!is.null(correct)) {
    if (!(is.logical(correct) && length(correct) == 1L && !is.na(correct))) {
      .thg_bad_input("`correct` must be TRUE or FALSE")
    }
    out <- out[out$correct %in% correct, , drop = FALSE]
  }
  if (!is.null(node)) {
    out <- out[out$node %in% .thg_node_filter(node), , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("margin", "score"))
    out <- out[order(-out[[sort_by]], out$node), , drop = FALSE]
  }
  rownames(out) <- NULL
  .thg_top_by(out, out$label, top)
}

# The attention tables of a hg_hypergat() evaluation, with each document's
# true and predicted label.
.thg_classification_attention <- function(x, plain, what, node, top) {
  words <- if (inherits(x, "hg_hypergat")) .thg_hypergat_words(x, node) else
    attr(x, "hyperedge_words")
  if (is.null(words)) {
    .thg_bad_input(paste0("`what = \"", what, "\"` needs a result of ",
                          "hg_hypergat() with `holdout`"))
  }
  out <- if (identical(what, "hyperedges")) {
    .thg_hypergat_hyperedge_summary(words)
  } else {
    words[order(words$node, words$hyperedge, -words$word_weight), ,
          drop = FALSE]
  }
  if (!is.null(node)) {
    out <- out[out$node %in% .thg_node_filter(node), , drop = FALSE]
    if (identical(what, "hyperedge_words") && is.data.frame(node) &&
        "hyperedge" %in% names(node)) {
      keep <- paste(out$node, out$hyperedge, sep = "\r") %in%
        paste(node$node, node$hyperedge, sep = "\r")
      out <- out[keep, , drop = FALSE]
    }
  }
  at <- match(out$node, plain$node)
  out <- cbind(out["node"], label = plain$label[at],
               predicted = plain$predicted[at],
               out[setdiff(names(out), "node")], stringsAsFactors = FALSE)
  rownames(out) <- NULL
  group <- if (identical(what, "hyperedges")) out$node else
    paste(out$node, out$hyperedge, sep = "\r")
  .thg_top_by(out, group, top)
}

#' @rdname hg_get.hg_classification
#' @export
print.hg_classification <- function(x, ...) {
  accuracy <- hg_get(x, what = "accuracy")
  cat(sprintf(paste0("Classification (%s): %d documents, %d held out ",
                     "(%.0f%% of the labelled)\n"),
              attr(x, "method"), nrow(x), accuracy$n_test,
              100 * attr(x, "holdout")))
  cat(sprintf(paste0("Held-out accuracy %.3f, balanced accuracy %.3f ",
                     "(chance %.3f)\n\n"),
              accuracy$accuracy, accuracy$balanced_accuracy, accuracy$chance))
  print(hg_get(x, what = "classes"), row.names = FALSE)
  invisible(x)
}

#' @rdname hg_get.hg_classification
#' @param object An `hg_classification`.
#' @export
summary.hg_classification <- function(object, ...) {
  hg_get(object, what = "classes")
}

#' @rdname hg_get.hg_classification
#' @param y Unused.
#' @param type For `plot()`: `"confusion"` (default), the counts of true
#'   against predicted labels on the test documents, or `"hyperedges"` (a
#'   result of `hg_hypergat()` only), the `top` heaviest hyperedges of each
#'   document in `node`, with its uniform normalization baseline. Without
#'   `node`, the first held-out document of each class is selected (the
#'   first training document if there is no holdout).
#' @export
plot.hg_classification <- function(x, y, type = c("confusion", "hyperedges"),
                                   node = NULL, top = 8L, ...) {
  type <- match.arg(type)
  if (identical(type, "hyperedges")) {
    return(.thg_plot_classification_hyperedges(x, node, top))
  }
  confusion <- hg_get(x, what = "confusion")
  ggplot2::ggplot(confusion, ggplot2::aes(x = .data$predicted,
                                          y = .data$label,
                                          fill = .data$n)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$n)) +
    ggplot2::scale_fill_gradient(low = "white", high = "#0072B2") +
    ggplot2::labs(x = "predicted", y = "true label", fill = "documents") +
    ggplot2::theme_minimal(base_size = 12)
}

.thg_plot_classification_hyperedges <- function(x, node, top) {
  node <- node %||% hg_get(x, split = if (any(x$split == "test"))
                            "test" else "train", top = 1)
  edges <- hg_get(x, what = "hyperedges", node = node, top = top)
  if (nrow(edges) == 0L) {
    .thg_bad_input("no hyperedge of the documents in `node`")
  }
  first_words <- ifelse(nchar(edges$text) > 90,
                        paste0(substr(edges$text, 1, 87), "..."), edges$text)
  edges$name <- vapply(paste0(edges$hyperedge, ": ", first_words),
                       \(t) paste(strwrap(t, 60), collapse = "\n"),
                       character(1L), USE.NAMES = FALSE)
  edges$panel <- sprintf("%s (true %s, predicted %s)", edges$node,
                         edges$label, edges$predicted)
  # one y position per document and hyperedge, heaviest at the top
  edges$position <- factor(paste(edges$panel, edges$name, sep = "\r"),
                           levels = rev(paste(edges$panel, edges$name,
                                              sep = "\r")))
  ggplot2::ggplot(edges, ggplot2::aes(x = .data$weight, y = .data$position,
                                      fill = .data$kind)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x = .data$uniform_weight),
                         shape = 4, size = 3, colour = "black") +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = 1,
                        scales = "free_y") +
    ggplot2::scale_y_discrete(labels = \(l) sub(".*\r", "", l)) +
    ggplot2::scale_fill_manual(values = c(sentence = "#0072B2",
                                          topic = "#E69F00")) +
    ggplot2::labs(x = "mean edge attention",
                  y = NULL, fill = "hyperedge",
                  caption = paste("Crosses: uniform normalization baseline.",
                    "A word in only one hyperedge gives it weight 1.",
                    "Attention is not a prediction explanation.")) +
    ggplot2::theme_minimal(base_size = 10)
}
