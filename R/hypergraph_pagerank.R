# Hypergraph PageRank with edge-dependent vertex weights (Chitra & Raphael
# 2019), as a tidy verb with personalization and sparse support. The dense
# EDVW transition matrix is built here (verified at machine precision against
# the stationary distribution of the random-walk Laplacian in
# hypergraph_laplacian.R, itself HyperNetX-parity-tested) and against
# .hg_centrality_fit(type = "pagerank"); PageRank adds damping and
# personalization on top.

# EDVW transition matrix: from vertex v, pick an incident hyperedge e with
# probability proportional to omega(e), then land on w in e with
# probability proportional to the incidence weight gamma_e(w). Row-
# stochastic by construction.
.thg_transition <- function(hg, edge_weights = NULL) {
  incidence <- hg$incidence
  # Package-wide default (.hl_build, .hl_rw_transition): explicit weights,
  # else a window hypergraph's window counts, else the Hayashi et al.
  # (2020) heuristic, the HyperNetX default: population SD of the edge's
  # non-zero vertex weights, plus one (unit weights on a binary incidence).
  edge_weights <- edge_weights %||% hg$window_counts
  if (is.null(edge_weights)) {
    edge_weights <- .hl_default_edge_weights(incidence)
  } else if (!(is.numeric(edge_weights) &&
               length(edge_weights) == ncol(incidence) &&
               all(is.finite(edge_weights)) && all(edge_weights > 0))) {
    .ho_input_error("`edge_weights` must be positive and one per hyperedge")
  }
  # Empty hyperedges (allowed by the random generators) can never be picked
  # by the walk: leave them out rather than divide by their zero size.
  keep <- .hl_nonempty(incidence)
  incidence <- incidence[, keep, drop = FALSE]
  membership <- (incidence > 0) * 1
  w_keep <- edge_weights[keep]
  delta <- colSums(incidence)
  vertex_degree <- as.numeric(membership %*% w_keep)
  if (any(vertex_degree <= 0)) {
    stop(errorCondition(
      "every vertex must belong to at least one hyperedge",
      class = "hypergraphs_bad_input", call = NULL
    ))
  }
  transition <- (membership %*% (t(incidence) * (w_keep / delta))) /
    vertex_degree
  list(transition = transition, edge_weights = edge_weights)
}

#' Hypergraph PageRank with edge-dependent vertex weights
#'
#' Ranks the vertices of a weighted hypergraph by the stationary importance
#' of the Chitra & Raphael (2019) random walk: from a vertex, pick an
#' incident hyperedge with probability proportional to its edge weight,
#' then land on a member vertex with probability proportional to that
#' vertex's *edge-dependent* incidence weight (for a text hypergraph, e.g.,
#' the tf-idf of the word in the document). Damping and an optional
#' personalization vector turn the stationary distribution into PageRank
#' (Page et al. 1999): with probability `damping` the walker follows the
#' hypergraph walk, otherwise it teleports.
#'
#' Chitra & Raphael prove that with edge-*independent* vertex weights the
#' walk collapses to a random walk on a weighted clique-expansion graph —
#' so the hypergraph brings genuinely new information exactly when the
#' incidence weights differ across the hyperedges a vertex belongs to
#' (tested against the closed-form graph stationary distribution). With
#' `damping = 1` and default `edge_weights`, the result equals the
#' stationary distribution of the Hayashi et al. (2020) EDVW walk used by
#' [hg_cluster()] (tested at `1e-12` against [hg_laplacian()]), and
#' with uniform teleportation it equals
#' `hg_centrality(type = "pagerank")` (tested).
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any hypergraphs
#'   `net_hg` (connected when `damping = 1`).
#' @param damping Probability of following the hypergraph walk (default
#'   `0.85`); `1 - damping` is the teleport probability. Must be in
#'   `(0, 1]`; `damping = 1` gives the pure stationary distribution, which
#'   is unique only on a connected hypergraph: a disconnected one raises
#'   `hypergraphs_hypergraph_disconnected`.
#' @param personalized Optional restart preference: a character vector of
#'   vertex names (uniform teleport over exactly those vertices; a name
#'   given twice counts once), or a named vector of finite non-negative
#'   teleport weights over (a subset of) the vertex names, each name at
#'   most once; unnamed vertices get teleport probability 0. `NULL`
#'   (default) teleports uniformly.
#' @param edge_weights Positive hyperedge weights (one per hyperedge), or
#'   `NULL` (default): a window hypergraph's window counts when present,
#'   otherwise the Hayashi et al. heuristic used by the Laplacian engines
#'   (the population standard deviation of each edge's non-zero vertex
#'   weights plus one, which reduces to unit weights on a binary
#'   incidence).
#' @param sort_by `NULL` (default, vertex order) or `"pagerank"` to sort
#'   descending (ties broken by vertex name).
#' @param n Return only the first `n` rows after sorting: a whole number
#'   of at least 1, or `Inf` (default) for all.
#' @param max_iter,tol Power-iteration cap (a whole number of at least 1) and L1
#'   convergence tolerance (a positive number).
#'
#' @return A base `data.frame`, one row per vertex, with columns `node` and
#'   `pagerank` (non-negative, summing to 1).
#'
#' @section Conditions: Raises `hypergraphs_bad_input` for broken contracts
#'   (including a non-finite, duplicated or unknown `personalized` entry),
#'   `hypergraphs_hypergraph_disconnected` for `damping = 1` on a
#'   disconnected hypergraph, and warns with `hypergraphs_no_converge`
#'   (returning the last iterate) when `max_iter` is reached before `tol`.
#'
#' @references
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. *ICML 2019*.
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. *CIKM 2020*.
#' \doi{10.1145/3340531.3412034}
#'
#' Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank
#' citation ranking: Bringing order to the web. Stanford InfoLab.
#'
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"),
#'    weight = "tfidf")
#' hg_pagerank(hg, sort_by = "pagerank")
#' hg_pagerank(hg, personalized = c(cooking_1 = 1), sort_by = "pagerank")
#' @export
hg_pagerank <- function(hg, damping = 0.85, personalized = NULL,
                        edge_weights = NULL, sort_by = NULL, n = Inf,
                        max_iter = 1000L, tol = 1e-12) {
  .thg_check_hg(hg)
  damping <- .ho_check_number(damping, "damping", min = 0, max = 1)
  if (!(damping > 0)) .ho_input_error("`damping` must be in (0, 1]")
  n <- .ho_check_count(n, "n", allow_inf = TRUE)
  max_iter <- .ho_check_count(max_iter, "max_iter")
  tol <- .hg_check_tol(tol)
  # Undamped, the walk's stationary distribution is unique only on a
  # connected hypergraph; on a disconnected one the answer would depend on
  # where the iteration started.
  if (damping >= 1 && !.hl_connected(hg)) {
    stop(errorCondition(
      paste0("`damping = 1` needs a connected hypergraph: the walk has no ",
             "unique stationary distribution. Use `damping < 1` or analyze ",
             "components separately."),
      class = "hypergraphs_hypergraph_disconnected", call = NULL
    ))
  }

  if (.thg_is_sparse(hg)) {
    ops <- .thg_walk_operators(hg, edge_weights = edge_weights)
    step <- ops$left
    ids <- rownames(hg$incidence)
  } else {
    walk <- .thg_transition(hg, edge_weights = edge_weights)
    step <- function(v) as.numeric(v %*% walk$transition)
    ids <- rownames(walk$transition)
  }
  teleport <- .thg_teleport(personalized, ids)

  rank <- teleport
  converged <- FALSE
  iter <- 0L
  # power iteration: each step depends on the previous iterate, so the loop
  # cannot be vectorized away; bounded by max_iter with convergence surfaced
  while (iter < max_iter) {
    iter <- iter + 1L
    updated <- damping * step(rank) + (1 - damping) * teleport
    if (sum(abs(updated - rank)) < tol) {
      rank <- updated
      converged <- TRUE
      break
    }
    rank <- updated
  }
  if (!converged) {
    warning(warningCondition(
      sprintf("PageRank did not converge in %d iterations (returning the last iterate)",
              as.integer(max_iter)),
      class = "hypergraphs_no_converge"
    ))
  }

  out <- data.frame(node = ids, pagerank = rank / sum(rank),
                    row.names = NULL)
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, choices = "pagerank")
    out <- out[order(-out$pagerank, out$node), , drop = FALSE]
    rownames(out) <- NULL
  }
  if (is.finite(n) && n < nrow(out)) {
    out <- out[seq_len(n), , drop = FALSE]
  }
  out
}

# Teleport distribution of hg_pagerank() over the vertex ids. NULL is
# uniform; a character vector is a set of vertices (uniform over them); a
# named numeric vector gives finite non-negative weights, one per vertex.
# Weights are divided by their maximum before summing, so finite weights
# near the double limit cannot overflow the normalizing sum.
.thg_teleport <- function(personalized, ids) {
  if (is.null(personalized)) return(rep(1 / length(ids), length(ids)))
  if (is.character(personalized)) {
    given <- .ho_check_ids(unique(personalized), "personalized")
    personalized <- stats::setNames(rep(1, length(given)), given)
  }
  if (!is.numeric(personalized) || is.factor(personalized) ||
      is.null(names(personalized)) || anyNA(personalized) ||
      any(!is.finite(personalized)) || any(personalized < 0) ||
      !any(personalized > 0)) {
    .ho_input_error(paste0(
      "`personalized` must be a named vector of finite non-negative ",
      "weights with a positive sum, over vertex names of `hg`"
    ))
  }
  given <- .ho_check_ids(names(personalized), "personalized")
  unknown <- setdiff(given, ids)
  if (length(unknown)) {
    .ho_input_error(sprintf(
      "`personalized` names vertices not in `hg`: %s",
      paste(utils::head(unknown, 5L), collapse = ", ")
    ))
  }
  teleport <- rep(0, length(ids))
  teleport[match(given, ids)] <- personalized / max(personalized)
  teleport / sum(teleport)
}
