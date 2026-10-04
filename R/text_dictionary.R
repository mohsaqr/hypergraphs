# Dictionary labels: a node takes the category whose terms it carries. In a
# text hypergraph with documents as nodes, a document's terms are the
# hyperedges it belongs to -- its words, or with `separator` its keywords --
# so a dictionary of keywords per category codes every paper whose authors
# named one of them.

#' Label nodes with a dictionary of terms
#'
#' Codes each node by the hyperedges it belongs to. In a text hypergraph
#' with documents as nodes, the hyperedges of a document are its words, or
#' with `text_hypergraph(separator = ";")` its keywords, so a dictionary that
#' lists the terms of each category labels every document that carries one
#' of them. This is dictionary-based content analysis (Grimmer & Stewart
#' 2013): the categories and their terms are the analyst's, and the
#' labelling is a count. A node takes the category with the most matched
#' terms. A node that matches no term, or matches two categories equally,
#' takes no label, and a message reports how many. The result is the
#' `labels` input of [hg_classify()] and [hg_hypergat()], which spread or
#' learn the labels to the nodes the dictionary cannot code.
#'
#' @param hg A `net_hg`, typically from [text_hypergraph()].
#' @param dictionary A named list of character vectors, one per category,
#'   each holding the terms of that category as they appear among the
#'   hyperedge names (lowercase for a text hypergraph built with
#'   `lowercase = TRUE`). A term may belong to one category only.
#' @return A base `data.frame`, one row per labelled node: `node`, `label`
#'   (the category), `n_terms` (how many of the category's terms the node
#'   carries) and `terms` (those terms, joined by `"; "`), ordered by
#'   `label` and `node`. Raises `hypergraphs_bad_input` for a malformed
#'   dictionary and warns with `hypergraphs_unmatched_terms` when a term of
#'   the dictionary is not a hyperedge of `hg`.
#' @references
#' Grimmer, J., & Stewart, B. M. (2013). Text as data: The promise and
#' pitfalls of automatic content analysis methods for political texts.
#' \emph{Political Analysis}, 21(3), 267-297. \doi{10.1093/pan/mps028}
#' @examples
#' papers <- data.frame(
#'   keywords = c("medical education; online learning",
#'                "teacher education; online learning",
#'                "medical students; assessment",
#'                "online learning; assessment")
#' )
#' kw <- text_hypergraph(papers, column = "keywords", separator = ";")
#' hg_dictionary(kw, dictionary = list(
#'   health = c("medical education", "medical students"),
#'   teachers = "teacher education"
#' ))
#' @export
hg_dictionary <- function(hg, dictionary) {
  .thg_check_hg(hg)
  if (!is.list(dictionary) || length(dictionary) == 0L ||
      is.null(names(dictionary)) || anyNA(names(dictionary)) ||
      !all(nzchar(names(dictionary))) || anyDuplicated(names(dictionary)) ||
      !all(vapply(dictionary, \(t) is.character(t) && length(t) > 0L &&
                    !anyNA(t), logical(1L)))) {
    .thg_bad_input(paste0("`dictionary` must be a named list of character ",
                          "vectors, one per category, with unique names"))
  }
  terms <- unlist(lapply(dictionary, unique), use.names = FALSE)
  category <- rep(names(dictionary), lengths(lapply(dictionary, unique)))
  shared <- unique(terms[duplicated(terms)])
  if (length(shared) > 0L) {
    .thg_bad_input(sprintf("terms in more than one category: %s",
                           paste(shared, collapse = ", ")))
  }
  incidence <- hg$incidence
  unmatched <- setdiff(terms, colnames(incidence))
  if (length(unmatched) > 0L) {
    warning(warningCondition(
      sprintf("%d dictionary term(s) are not hyperedges of `hg`: %s",
              length(unmatched), paste(unmatched, collapse = ", ")),
      class = "hypergraphs_unmatched_terms", call = NULL
    ))
  }
  matched <- terms %in% colnames(incidence)
  terms <- terms[matched]
  category <- category[matched]
  if (length(terms) == 0L) {
    .thg_bad_input("no term of `dictionary` is a hyperedge of `hg`")
  }
  # node x term presence, then node x category counts
  present <- as.matrix(incidence[, terms, drop = FALSE] != 0)
  counts <- vapply(names(dictionary), \(k) {
    rowSums(present[, category == k, drop = FALSE])
  }, numeric(nrow(present)))
  counts <- matrix(counts, nrow = nrow(present),
                   dimnames = list(rownames(incidence), names(dictionary)))
  best <- apply(counts, 1L, max)
  n_best <- rowSums(counts == best)
  labelled <- best > 0 & n_best == 1L
  tied <- sum(best > 0 & n_best > 1L)
  message(sprintf(paste0("%d of %d nodes labelled; %d match no term and %d ",
                         "tie between categories"),
                  sum(labelled), nrow(counts), sum(best == 0), tied))
  if (!any(labelled)) {
    .thg_bad_input("no node carries a term of `dictionary` unambiguously")
  }
  nodes <- rownames(counts)[labelled]
  label <- names(dictionary)[max.col(counts[labelled, , drop = FALSE],
                                     ties.method = "first")]
  carried <- vapply(seq_along(nodes), \(i) {
    hits <- terms[present[nodes[i], ] & category == label[i]]
    paste(sort(hits), collapse = "; ")
  }, character(1L))
  out <- data.frame(node = nodes, label = label,
                    n_terms = as.integer(best[labelled]), terms = carried,
                    stringsAsFactors = FALSE)
  out <- out[order(out$label, out$node), , drop = FALSE]
  rownames(out) <- NULL
  out
}
