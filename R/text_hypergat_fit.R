# Readers and prediction for a single trained HyperGAT model.

#' Read a fitted HyperGAT classifier
#'
#' The classifier retains its trained network, vocabulary and topic keywords.
#' Reading a table or predicting new documents never trains the network.
#' Attention tables describe internal aggregation weights. They are not
#' class-specific contributions or explanations of a prediction.
#'
#' @param x A fitted `hg_hypergat`.
#' @param what `NULL` or `"predictions"`, `"accuracy"`, `"classes"`,
#'   `"confusion"`, `"history"`, `"documents"`, `"attention"`, `"hyperedges"`
#'   or `"hyperedge_words"`. Evaluation tables use held-out documents only.
#'   `"documents"` includes the text with its prediction.
#' @param split,correct,sort_by Filters and ordering for predictions and
#'   documents, as in [hg_get.hg_classification()].
#' @param node Document ids or a table with a `node` column.
#' @param top Number of rows per class (predictions and documents), per
#'   document (hyperedges or attention), or per hyperedge (hyperedge words).
#' @param ... Unused.
#' @return A plain data frame. The attention tables include uniform
#'   normalization baselines: `word_uniform = 1 / n_words` within an edge,
#'   and `edge_uniform = 1 / n_edges` for a word. Hyperedge summaries include
#'   `uniform_weight` (mean `edge_uniform`) and `exclusive_fraction` (the
#'   fraction of words belonging only to that hyperedge). Such words have
#'   `edge_weight = 1` regardless of the learned parameters.
#' @references
#' Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with
#' less: Hypergraph attention networks for inductive text classification.
#' \emph{EMNLP 2020}, 4927-4936.
#' @export
hg_get.hg_hypergat <- function(x, what = NULL, split = NULL, correct = NULL,
                              node = NULL, sort_by = NULL, top = NULL, ...) {
  what <- match.arg(what %||% "predictions", c("predictions", "accuracy",
    "classes", "confusion", "history", "documents", "attention",
    "hyperedges", "hyperedge_words"))
  state <- attr(x, "hypergat")
  if (what == "documents") {
    predictions <- hg_get.hg_classification(x, what = "predictions",
      split = split, correct = correct, node = node,
      sort_by = sort_by, top = top)
    predictions$text <- state$input$text[match(predictions$node,
                                               state$input$node)]
    return(predictions)
  }
  if (what %in% c("history", "attention")) {
    if (!is.null(split) || !is.null(correct) || !is.null(sort_by)) {
      .thg_bad_input("`split`, `correct` and `sort_by` apply to predictions or documents")
    }
    if (what == "history") {
      if (!is.null(node) || !is.null(top)) {
        .thg_bad_input("`node` and `top` do not apply to history")
      }
      return(attr(x, "history"))
    }
    inputs <- .thg_hypergat_diagnostic_inputs(x, node)
    state$model$eval()
    out <- .thg_hypergat_attention(state$model, inputs$documents, inputs$corpus)
    out <- out[order(out$node, -out$attention, out$word), , drop = FALSE]
    rownames(out) <- NULL
    return(.thg_top_by(out, out$node, top))
  }
  hg_get.hg_classification(x, what = what, split = split, correct = correct,
                           node = node, sort_by = sort_by, top = top)
}

.thg_hypergat_diagnostic_inputs <- function(x, node) {
  state <- attr(x, "hypergat")
  ids <- if (is.null(node)) x$node else .thg_node_filter(node)
  at <- match(ids, state$corpus$doc_id)
  at <- at[!is.na(at)]
  if (!length(at)) .thg_bad_input("no retained document in `node`")
  corpus <- state$corpus
  corpus$doc_id <- corpus$doc_id[at]
  corpus$sentences <- corpus$sentences[at]
  corpus$texts <- corpus$texts[at]
  list(corpus = corpus, documents = state$documents[at])
}

.thg_hypergat_words <- function(x, node = NULL) {
  state <- attr(x, "hypergat")
  if (!is.null(state$cache$words)) return(state$cache$words)
  inputs <- .thg_hypergat_diagnostic_inputs(x, node)
  state$model$eval()
  words <- .thg_hypergat_hyperedge_words(state$model, inputs$documents,
    inputs$corpus, attr(x, "semantic")$keywords)
  if (is.null(node)) state$cache$words <- words
  words
}

.thg_hypergat_score <- function(model, documents, classes) {
  model$eval()
  starts <- seq(1L, length(documents), by = 16L)
  probs <- lapply(starts, function(s) {
    chunk <- s:min(s + 15L, length(documents))
    batch <- .thg_hypergat_batch(documents[chunk])
    torch::with_no_grad(as.matrix(torch::nnf_softmax(
      model(batch$items, batch$adj, batch$mask), dim = 2L)))
  })
  out <- do.call(rbind, probs)
  colnames(out) <- classes
  out
}

#' Predict documents with a fitted HyperGAT classifier
#'
#' Reuses the trained network, vocabulary, preprocessing and topic keywords.
#' Words absent from the fitted vocabulary are ignored. A document with no
#' remaining words receives missing predictions and a warning, preserving
#' its row and id. The winning softmax probability is a model score, not a
#' calibrated probability of correctness.
#'
#' @param object A fitted `hg_hypergat`.
#' @param newdata A character vector or data frame of new documents. `NULL`
#'   returns the original predictions.
#' @param column,id Text and id columns of a data frame. Defaults to the
#'   columns used when fitting. An unnamed character vector receives
#'   sequential document ids.
#' @param labels Optional observed labels for `newdata`, in the same formats
#'   as [hg_hypergat()]. Supplying them evaluates the predictions through
#'   [hg_get.hg_classification()] and `plot()`. These labels never update
#'   the model. Unscorable labelled documents count as errors.
#' @param type `"predictions"` (default), a document table, or
#'   `"probabilities"`, a matrix with one column per fitted class.
#' @param ... Unused.
#' @return A plain data frame (`node`, `label`, `predicted`, `score`,
#'   `margin`) or a probability matrix, in input order. No fitting occurs.
#' @references
#' Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with
#' less: Hypergraph attention networks for inductive text classification.
#' \emph{EMNLP 2020}, 4927-4936.
#' @export
predict.hg_hypergat <- function(object, newdata = NULL, column = NULL,
                               id = NULL, labels = NULL,
                               type = c("predictions", "probabilities"), ...) {
  type <- match.arg(type)
  state <- attr(object, "hypergat")
  if (!is.null(labels) && (is.null(newdata) || type != "predictions")) {
    .thg_bad_input("`labels` requires `newdata` and `type = \"predictions\"`")
  }
  if (is.null(newdata) && type == "predictions") return(hg_get(object))
  if (is.null(newdata)) {
    at <- match(object$node, state$corpus$doc_id)
    probs <- matrix(NA_real_, nrow = nrow(object), ncol = length(state$classes),
                    dimnames = list(object$node, state$classes))
    kept <- which(!is.na(at))
    if (length(kept)) probs[kept, ] <- .thg_hypergat_score(
      state$model, state$documents[at[kept]], state$classes)
    return(probs)
  }
  if (is.data.frame(newdata)) {
    column <- column %||% state$preprocessing$column
    id <- id %||% state$preprocessing$id
    if (!is.character(column) || length(column) != 1L ||
        is.na(column) || !column %in% names(newdata)) {
      .thg_bad_input("`column` must name the text column of `newdata`")
    }
    if (!is.null(id) && (!is.character(id) || length(id) != 1L ||
                        is.na(id) || !id %in% names(newdata))) {
      .thg_bad_input("`id` must name a column of `newdata`")
    }
    text <- as.character(newdata[[column]])
    ids <- if (is.null(id)) sprintf("doc_%d", seq_len(nrow(newdata))) else
      as.character(newdata[[id]])
  } else if (is.character(newdata)) {
    text <- newdata
    ids <- names(text) %||% sprintf("doc_%d", seq_along(text))
  } else {
    .thg_bad_input("`newdata` must be a character vector or a data frame")
  }
  if (anyNA(ids) || anyDuplicated(ids)) {
    .thg_bad_input("document ids must be unique and non-missing")
  }
  probs <- matrix(NA_real_, nrow = length(text), ncol = length(state$classes),
                  dimnames = list(ids, state$classes))
  if (length(text)) {
    corpus <- .thg_hypergat_corpus(text, ids, state$preprocessing$stop_words,
      min_count = 1L, lowercase = state$preprocessing$lowercase,
      vocabulary = state$corpus$vocab)
    if (length(corpus$doc_id)) {
      documents <- corpus$sentences
      if (attr(object, "semantic")$method != "none") {
        documents <- .thg_hypergat_semantic_docs(documents, corpus$vocab,
          attr(object, "semantic")$keywords)
      }
      probs[match(corpus$doc_id, ids), ] <- .thg_hypergat_score(
        state$model, documents, state$classes)
    }
  }
  if (type == "probabilities") return(probs)
  out <- data.frame(node = ids, label = rep(NA_character_, length(ids)),
    predicted = rep(NA_character_, length(ids)), score = rep(NA_real_, length(ids)),
    margin = rep(NA_real_, length(ids)),
    stringsAsFactors = FALSE)
  kept <- which(rowSums(!is.na(probs)) > 0L)
  if (length(kept)) {
    scored <- .hl_score_predictions(probs[kept, , drop = FALSE],
                                    rep(NA_character_, length(kept)), "none")
    out[kept, c("predicted", "score", "margin")] <-
      scored[c("predicted", "score", "margin")]
  }
  if (!is.null(labels)) {
    observed <- .thg_resolve_text_labels(newdata, labels, ids)
    if (!is.character(observed) || is.null(names(observed)) ||
        anyDuplicated(names(observed)) ||
        any(!names(observed) %in% ids) ||
        any(!observed %in% state$classes)) {
      .thg_bad_input("`labels` must name new document ids and fitted classes")
    }
    out <- .thg_classification(out, observed, observed,
      method = "hg_hypergat prediction", holdout = 1)
    attr(out, "classes") <- state$classes
    attr(out, "documents") <- data.frame(node = ids, text = text,
                                          stringsAsFactors = FALSE)
  }
  out
}

#' @rdname hg_get.hg_hypergat
#' @export
print.hg_hypergat <- function(x, ...) {
  if (any(x$split == "test")) return(print.hg_classification(x, ...))
  cat(sprintf("HyperGAT classifier: %d documents, %d classes, %d labelled\n",
    nrow(x), length(attr(x, "hypergat")$classes), sum(!is.na(x$label))))
  cat("No held-out evaluation. predict() reuses the fitted model.\n")
  invisible(x)
}

#' @rdname hg_get.hg_hypergat
#' @param object A fitted `hg_hypergat`.
#' @export
summary.hg_hypergat <- function(object, ...) {
  hg_get(object, what = "classes")
}

#' @rdname hg_get.hg_hypergat
#' @param y Unused.
#' @param type `"confusion"` (default), `"history"`, `"hyperedges"` (the
#'   word hypergraph of one document), or `"attention"`.
#'   Attention plots show weights alongside their uniform normalization
#'   baseline. These weights do not measure a sentence's contribution to
#'   the predicted class.
#' @export
plot.hg_hypergat <- function(x, y,
                            type = c("confusion", "history", "hyperedges", "attention"),
                            node = NULL, top = 8L, ...) {
  type <- match.arg(type)
  if (type == "history") {
    if (!is.null(node)) .thg_bad_input("`node` does not apply to history")
    history <- hg_get(x, what = "history")
    return(ggplot2::ggplot(history, ggplot2::aes(.data$epoch, .data$loss)) +
      ggplot2::geom_line(colour = "#0072B2") +
      ggplot2::geom_point(colour = "#0072B2") +
      ggplot2::labs(x = "epoch", y = "training loss") +
      ggplot2::theme_minimal(base_size = 12))
  }
  if (type == "hyperedges") {
    node <- node %||% attr(x, "hypergat")$corpus$doc_id[1L]
    if (length(.thg_node_filter(node)) != 1L) {
      .thg_bad_input("a hypergraph plot takes one document; select it with `node`")
    }
    words <- hg_get(x, what = "hyperedge_words", node = node)
    membership <- words[c("word", "hyperedge")]
    hg <- group_hypergraph(membership, node = "word", hyperedge = "hyperedge")
    return(plot(hg, ...))
  }
  if (type == "attention") {
    return(.thg_plot_classification_hyperedges(x, node, top))
  }
  if (type == "confusion" && !any(x$split == "test")) {
    .thg_bad_input("a confusion plot needs held-out labels; fit with `holdout`")
  }
  plot.hg_classification(x, type = type, node = node, top = top, ...)
}
