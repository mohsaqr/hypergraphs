# Dual hypergraph: swap the vertex and hyperedge roles. For a bag text
# hypergraph this is exactly the nodes = "doc"/"word" flip, available on a
# fitted object without re-tokenizing.

#' Dual of a hypergraph
#'
#' Returns the dual hypergraph: every hyperedge becomes a vertex and every
#' vertex becomes a hyperedge, with the transposed weighted incidence. For a
#' bag-construction [text_hypergraph()], the dual is identical to rebuilding
#' with the opposite `nodes` orientation (tested), so document-level and
#' word-level analyses can share one constructed object. Duals of windowed
#' and kNN hypergraphs are returned as plain `net_hg` objects (their
#' corpus bookkeeping does not transpose meaningfully).
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any hypergraphs
#'   `net_hg`.
#' @return A hypergraph whose incidence is exactly the transpose of `hg`'s,
#'   sparse when `hg` is sparse, with every row and column kept: an isolated
#'   vertex of `hg` becomes an empty hyperedge of the dual and an empty
#'   hyperedge an isolated vertex, so the dual of the dual is `hg`'s
#'   incidence again. A `text_hypergraph` with flipped `nodes` for bag
#'   constructions (dense or sparse), otherwise a `net_hg`. Hyperedge
#'   metadata (`edge_data`, window counts) describes hyperedges that become
#'   vertices, so it is not carried over. Accepted by all `hg_*` verbs.
#' @references
#' Berge, C. (1989). *Hypergraphs: Combinatorics of Finite Sets*.
#' North-Holland Mathematical Library 45. North-Holland.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
#' dual <- dual_hypergraph(hg)
#' dual
#' hg_measures(dual, what = "edges")
#' @export
dual_hypergraph <- function(hg) {
  .thg_check_hg(hg)
  sparse <- .thg_is_sparse(hg)
  # The dual is the transpose of the whole incidence, so an isolated vertex
  # becomes an empty hyperedge and an empty hyperedge an isolated vertex;
  # every other field is derived from that transpose, on both storages.
  incidence <- if (sparse) Matrix::t(hg$incidence) else t(hg$incidence)
  dual <- .thg_from_incidence(
    incidence,
    params = list(source = "dual_hypergraph", sparse = sparse)
  )
  if (inherits(hg, "text_hypergraph") &&
      identical(hg$text$construction, "bag")) {
    dual$text <- hg$text
    dual$text$nodes <- if (identical(hg$text$nodes, "doc")) "word" else "doc"
    class(dual) <- c("text_hypergraph", class(dual))
  }
  dual
}
