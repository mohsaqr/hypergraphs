# Classification results with a held-out evaluation. hg_classify() with
# `holdout` hides a stratified share of the known labels, predicts them, and return the per-document table with `split` and `correct`
# columns as an `hg_classification`, a data.frame whose print shows the
# evaluation and whose tables are read with hg_get().

# The one contract for a node -> label assignment (the labels of a
# classifier, the clusters of a topic summary, the groups of a prevalence
# table). Its names are node ids: present, non-empty, and naming each node
# once -- a node named twice with the same label is kept once, a node named
# twice with different labels is ambiguous and refused. A missing label
# marks the node as unlabelled, as a missing value in a document column
# does, and the entry is dropped. An unnamed vector passes through for the
# callers that take one label per node by position.
.thg_check_assignment <- function(assignment, arg = "labels") {
  ids <- names(assignment)
  if (is.null(ids)) return(assignment)
  if (anyNA(ids) || any(!nzchar(ids))) {
    .thg_bad_input(sprintf("`%s` must not name a node with a missing or empty id",
                           arg))
  }
  if (anyDuplicated(ids) > 0L) {
    value <- ifelse(is.na(assignment), "\r<NA>", as.character(assignment))
    distinct <- !duplicated(paste(ids, value, sep = "\r"))
    conflict <- unique(ids[distinct][duplicated(ids[distinct])])
    if (length(conflict) > 0L) {
      .thg_bad_input(sprintf(
        "`%s` gives node(s) more than one label: %s", arg,
        paste(utils::head(conflict, 5L), collapse = ", ")))
    }
    assignment <- assignment[!duplicated(ids)]
  }
  assignment[!is.na(assignment)]
}

# Labels given as the name of a column of the hypergraph's document table
# (a text hypergraph carries the input's other columns there) become a named
# vector over the documents; any other `labels` is read as a table or a named
# vector. Either way the assignment contract above is enforced.
.thg_resolve_labels <- function(hg, labels, arg = "labels") {
  if (!(is.character(labels) && length(labels) == 1L && is.null(names(labels)))) {
    return(.thg_check_assignment(.thg_labels_input(labels), arg))
  }
  documents <- if (is.list(hg$text)) hg$text$documents else NULL
  if (is.null(documents) || !labels %in% names(documents)) {
    .thg_bad_input(sprintf(paste0("`labels` = \"%s\" must name a column of the ",
                                  "hypergraph's document table"), labels))
  }
  values <- as.character(documents[[labels]])
  keep <- !is.na(values) & documents$doc %in% hg$nodes
  .thg_check_assignment(stats::setNames(values[keep], documents$doc[keep]),
                        arg)
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
#' `hg_classify()` with `holdout` returns an `hg_classification`: one row
#' per document with its true label, its predicted label, its split and
#' whether a held-out prediction is correct. `hg_get()` reads the
#' evaluation tables.
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
#'   true and predicted class: `label`, `predicted`, `n`.
#' @param split For `"predictions"`: keep only the `"train"`, the `"test"`
#'   or the `"unlabelled"` documents, the last being those without a label,
#'   whose predictions are the classifier's output.
#' @param correct For `"predictions"`: `TRUE` keeps the test documents
#'   predicted correctly, `FALSE` those predicted wrongly.
#' @param node For `"predictions"`: keep only these documents, given as
#'   ids or as a table with a `node` column, such as a table `hg_get()`
#'   returned.
#' @param sort_by For `"predictions"`: `"margin"` or `"score"`, from the
#'   largest.
#' @param top For `"predictions"`: the first `top` rows of each true label,
#'   so `split = "test", correct = TRUE, sort_by = "margin", top = 1` gives
#'   the most confident correct prediction of each class.
#' @param ... Unused.
#' @return A base `data.frame` as described under `what`. `print()` shows
#'   the held-out accuracy and the per-class table and returns `x`
#'   invisibly; `summary()` returns the per-class table; `plot()` returns a
#'   ggplot of the confusion table.
#' @references
#' Brodersen, K. H., Ong, C. S., Stephan, K. E., & Buhmann, J. M. (2010).
#' The balanced accuracy and its posterior distribution. \emph{Proceedings
#' of the 20th International Conference on Pattern Recognition}, 3121-3124.
#' \doi{10.1109/ICPR.2010.764}
#' @examples
#' hg <- text_hypergraph(head(forum_posts, 60), column = "text")
#' fit <- hg_classify(hg, labels = "topic", holdout = 0.25)
#' fit
#' hg_get(fit, what = "classes")
#' hg_get(fit, what = "confusion")
#' hg_get(fit, split = "test", correct = FALSE)
#' @export
hg_get.hg_classification <- function(x, what = c("predictions", "accuracy",
                                                 "classes", "confusion"),
                                     split = NULL, correct = NULL,
                                     node = NULL, sort_by = NULL, top = NULL,
                                     ...) {
  what <- .ho_match_what(what)
  plain <- .ho_plain(x)
  attributes(plain) <- attributes(plain)[c("names", "row.names", "class")]
  filters <- list(split = split, correct = correct, sort_by = sort_by)
  filters <- names(filters)[!vapply(filters, is.null, logical(1L))]
  if (!identical(what, "predictions") && length(filters) > 0L) {
    .thg_bad_input(sprintf("`%s` applies only to `what = \"predictions\"`",
                           filters[[1L]]))
  }
  if (identical(what, "predictions")) {
    return(.thg_classification_predictions(plain, split, correct, node,
                                           sort_by, top))
  }
  if (!is.null(node) || !is.null(top)) {
    .thg_bad_input(sprintf("`node` and `top` do not apply to `what = \"%s\"`",
                           what))
  }
  test <- plain[plain$split %in% "test", , drop = FALSE]
  classes <- attr(x, "classes") %||%
    sort(unique(stats::na.omit(plain$label)))
  if (identical(what, "confusion")) {
    # the column of documents without a prediction; renamed "(unscored).1"
    # if a real class already carries the name
    unscored <- utils::tail(make.unique(c(classes, "(unscored)")), 1L)
    predicted_classes <- c(classes, if (anyNA(test$predicted)) unscored)
    predicted <- ifelse(is.na(test$predicted), unscored, test$predicted)
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
#' @details `plot()` plots the counts of true against predicted labels on
#'   the test documents.
#' @export
plot.hg_classification <- function(x, y, ...) {
  .ho_no_dots(..., .for = "plot.hg_classification()")
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
