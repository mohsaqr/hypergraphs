# The projection tier: hypergraph -> weighted graph, and hypergraph ->
# s-line graph. Both are pure incidence algebra. Results expose graph-shaped
# matrices and edge lists at the cograph boundary: hypergraphs owns these
# hypergraph transformations; cograph owns downstream graph analysis and
# plotting.

# Membership pattern of an incidence matrix, sparse or dense, weights dropped.
.thg_binary <- function(x) (x != 0) * 1

# Scale columns of an incidence matrix by `coef`, sparse or dense. diag() with
# a length-one `coef` would build an identity of that size, hence `nrow=`.
.thg_scale_cols <- function(x, coef) {
  if (methods::is(x, "sparseMatrix")) {
    x %*% Matrix::Diagonal(x = coef)
  } else {
    x %*% diag(x = coef, nrow = length(coef))
  }
}

# Upper triangle of a symmetric weight matrix as a tidy edge list, sorted by
# (from, to) so the row order never depends on storage or input order.
.thg_tidy_pairs <- function(w, labels) {
  if (methods::is(w, "sparseMatrix")) {
    triplet <- methods::as(w, "TsparseMatrix")
    keep <- triplet@i < triplet@j & triplet@x != 0
    from <- labels[triplet@i[keep] + 1L]
    to <- labels[triplet@j[keep] + 1L]
    weight <- triplet@x[keep]
  } else {
    nz <- which(w != 0 & upper.tri(w), arr.ind = TRUE)
    from <- labels[nz[, "row"]]
    to <- labels[nz[, "col"]]
    weight <- as.numeric(w[nz])
  }
  out <- data.frame(from = from, to = to, weight = as.numeric(weight),
                    row.names = NULL)
  out <- out[order(out$from, out$to), , drop = FALSE]
  row.names(out) <- NULL
  out
}

#' The pairwise network of a hypergraph
#'
#' Projects a hypergraph onto a network of pairs of its nodes. Two nodes are
#' joined when they belong to a common hyperedge, and the weight of the
#' edge depends on `type`. `"clique"` is the clique expansion: a hyperedge
#' of size \eqn{|e|} contributes the same amount to each of its
#' \eqn{|e|(|e|-1)/2} pairs, so the weight of a pair is the number of
#' hyperedges that contain both nodes (Zhou et al. 2006). `"association"`
#' divides the contribution of each hyperedge by \eqn{|e|-1}, following
#' Coupette et al. (2024):
#'
#' \deqn{w(\{u,v\}) = \sum_{e \supseteq \{u,v\}} \frac{1}{|e| - 1}}
#'
#' so the weight a hyperedge adds around each of its members is 1, and the
#' weighted degree of a node equals the number of hyperedges of size at
#' least two that contain it. A hyperedge of size one has no pairs and
#' contributes nothing.
#'
#' `"citation"` is the graph of a hypergraph whose hyperedges have sources:
#' one edge from the source of every hyperedge to each of its members, so a
#' hypergraph of citation blocks becomes the ordinary citation graph.
#' `duplicate_edges = "count"` weights an edge by the number of blocks that
#' repeat it (the multi-graph); `"collapse"` keeps the binary graph.
#' `directed = TRUE` keeps the source-to-member orientation; the default
#' symmetrises.
#'
#' The result is a network object that cograph plots and that
#' [hypergraph()] reads as a network, so the cliques of the projection can be
#' promoted back to hyperedges. Its edges are read with `hg_get()`.
#'
#' @param hg A hypergraph (`net_hg`), such as one built by [hypergraph()] or
#'   [text_hypergraph()].
#' @param type The projection: `"clique"` (default), `"association"` or
#'   `"citation"`.
#' @param weighted `type = "clique"` only. `TRUE` (default) uses the
#'   incidence weights, `FALSE` their membership pattern. Setting it with
#'   another `type` is an error, because those projections are defined on
#'   membership.
#' @param duplicate_edges For `type = "association"` or `"citation"`:
#'   `"count"` (default) lets repeated hyperedges contribute repeatedly (the
#'   multi-hypergraph); `"collapse"` lets each distinct member set
#'   contribute once (the binary hypergraph).
#' @param self_association For `type = "association"`: add the
#'   source-to-member association term of Coupette et al. (2024). Requires one
#'   source node per hyperedge through `edge_source` or
#'   `hg$edge_data$source`. For source `u`, every membership of target `v`
#'   contributes `1 / sum_e |e|` over the hyperedges sourced by `u`.
#' @param edge_source The node each hyperedge comes from: the name of a
#'   column of the hyperedge attributes, a vector of length `n_hyperedges`, a
#'   named vector keyed by hyperedge, or a data frame with columns `edge` and
#'   `source`. `NULL` uses an attribute column named `source`. Used with
#'   `self_association = TRUE` and `type = "citation"`.
#' @param directed For `type = "citation"`: keep the source-to-member
#'   direction. Default `FALSE`.
#' @return A `net_hg_pairwise`, which is also a `netobject` and a
#'   `cograph_network`: the weighted adjacency matrix of the projection
#'   (`weights`, zero diagonal, symmetric unless `directed = TRUE`), its
#'   nodes and edges. `hg_get(x)` returns the edges as a data frame with one
#'   row per pair of non-zero weight and columns `from`, `to` and `weight`,
#'   sorted by `from` then `to`.
#' @section Conditions:
#' `hypergraphs_bad_input` for an argument that does not apply to `type`, and
#' for `self_association` or `type = "citation"` without a source per
#' hyperedge.
#' @references
#' Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
#' *Philosophical Transactions of the Royal Society A*, 382(2270), 20230141.
#' \doi{10.1098/rsta.2023.0141}
#'
#' Zhou, D., Huang, J., & Schoelkopf, B. (2006). Learning with hypergraphs:
#' clustering, classification, and embedding. *NeurIPS 19*, 1601-1608.
#' @seealso [hg_line_graph()] for the projection onto hyperedges.
#' @examples
#' meetings <- data.frame(
#'   person = c("Alice", "Bob", "Carol", "Alice", "Bob", "Dave", "Eve"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m3", "m3"))
#' meeting_hg <- hypergraph(meetings, node = "person", hyperedge = "meeting")
#' meeting_network <- pairwise_network(meeting_hg)
#' hg_get(meeting_network)
#' association_network <- pairwise_network(meeting_hg, type = "association")
#' hg_get(association_network)
#' @export
pairwise_network <- function(hg, type = c("clique", "association", "citation"),
                             weighted = NULL,
                             duplicate_edges = c("count", "collapse"),
                             self_association = FALSE, edge_source = NULL,
                             directed = FALSE) {
  .thg_check_hg(hg)
  type <- match.arg(type)
  duplicate_edges <- match.arg(duplicate_edges)
  weights <- .hg_projection(hg, type = type, weighted = weighted,
                            what = "matrix", duplicate_edges = duplicate_edges,
                            self_association = self_association,
                            edge_source = edge_source, directed = directed)
  weights <- as.matrix(weights)
  storage.mode(weights) <- "double"
  net <- .wrap_netobject(weights, method = "pairwise_network",
                         directed = isTRUE(directed))
  net$params <- list(
    source = "pairwise_network",
    type = type,
    weighted = weighted %||% TRUE,
    duplicate_edges = duplicate_edges,
    self_association = self_association,
    n_hyperedges = hg$n_hyperedges,
    hypergraph_size_distribution = hg$size_distribution
  )
  class(net) <- unique(c("net_hg_pairwise", class(net)))
  net
}

#' @rdname hg_get
#' @export
hg_get.net_hg_pairwise <- function(x, what = c("edges", "nodes"), ...,
                                  sort_by = NULL, top = NULL) {
  what <- match.arg(what)
  labels <- rownames(x$weights)
  if (identical(what, "nodes")) {
    out <- data.frame(node = labels,
                      degree = as.integer(rowSums(x$weights != 0)),
                      strength = as.numeric(rowSums(x$weights)),
                      stringsAsFactors = FALSE)
    if (!is.null(sort_by)) {
      sort_by <- match.arg(sort_by, c("degree", "strength"))
      out <- out[order(-out[[sort_by]], out$node), , drop = FALSE]
    }
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  out <- if (isTRUE(x$directed)) {
    nz <- which(x$weights != 0, arr.ind = TRUE)
    edges <- data.frame(from = labels[nz[, "row"]], to = labels[nz[, "col"]],
                        weight = as.numeric(x$weights[nz]),
                        stringsAsFactors = FALSE)
    edges[order(edges$from, edges$to), , drop = FALSE]
  } else {
    .thg_tidy_pairs(x$weights, labels)
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, "weight")
    out <- out[order(-out$weight, out$from, out$to), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Plot the pairwise network of a hypergraph
#'
#' Plots the network returned by [pairwise_network()] with [cograph::splot()]:
#' a circle layout, which keeps every label apart however strongly the
#' events are tied, edges whose width grows with the square root of their
#' weight, and nodes whose area grows with their strength, the sum of the
#' weights of their edges.
#'
#' @param x A `net_hg_pairwise` from [pairwise_network()].
#' @param ... Arguments passed to [cograph::splot()] (e.g. `layout`, `seed`,
#'   `minimum`, `node_fill`); they override the defaults set here.
#' @return `x`, invisibly. cograph plots with base graphics.
#' @examples
#' meetings <- data.frame(
#'   person = c("Alice", "Bob", "Carol", "Alice", "Bob", "Dave", "Eve"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m3", "m3"))
#' meeting_network <- pairwise_network(hypergraph(meetings, node = "person",
#'                                          hyperedge = "meeting"))
#' plot(meeting_network)
#' @export
plot.net_hg_pairwise <- function(x, ...) {
  strength <- rowSums(abs(x$weights))
  top <- max(strength)
  size <- 2 + 3 * sqrt(if (top > 0) strength / top else strength)
  defaults <- list(layout = "circle", seed = 1, tna_styling = TRUE,
                   directed = isTRUE(x$directed), edge_label_style = "none",
                   node_size = size, node_fill = "#56B4E9",
                   edge_scale_mode = "sqrt", label_size = 0.8)
  given <- list(...)
  args <- c(list(x), utils::modifyList(defaults, given))
  do.call(cograph::splot, args)
  invisible(x)
}

# The projection itself, shared by pairwise_network() and the internal callers (communities, temporal layers, representations,
# plots). `weighted = NULL` is the clique default (TRUE); a non-NULL value
# with another type is refused.
.hg_projection <- function(hg, type = c("clique", "association", "citation"),
                           weighted = NULL, what = c("matrix", "edges"),
                           duplicate_edges = c("count", "collapse"),
                           self_association = FALSE, edge_source = NULL,
                           directed = FALSE) {
  method <- match.arg(type)
  what <- match.arg(what)
  duplicate_edges <- match.arg(duplicate_edges)
  if (!is.null(weighted) && !identical(method, "clique")) {
    stop(errorCondition(
      "`weighted` applies to `type = \"clique\"` only; the association and citation weightings are defined on membership, not on incidence weights",
      class = "hypergraphs_bad_input", call = NULL
    ))
  }
  weighted <- weighted %||% TRUE
  stopifnot("`weighted` must be TRUE or FALSE" =
              length(weighted) == 1L && is.logical(weighted) &&
              !is.na(weighted),
            "`self_association` must be TRUE or FALSE" =
              length(self_association) == 1L && is.logical(self_association) &&
              !is.na(self_association),
            "`directed` must be TRUE or FALSE" =
              length(directed) == 1L && is.logical(directed) && !is.na(directed))
  if (identical(method, "clique") && !identical(duplicate_edges, "count")) {
    .thg_bad_input("`duplicate_edges` applies to `type = \"association\"` or `\"citation\"` only")
  }
  if (!identical(method, "association") && self_association) {
    .thg_bad_input("`self_association` requires `type = \"association\"`")
  }
  if (!identical(method, "citation") && directed) {
    .thg_bad_input("`directed` applies to `type = \"citation\"` only")
  }
  if (identical(method, "citation")) {
    return(.thg_citation_projection(hg, edge_source, duplicate_edges,
                                    directed, what))
  }
  incidence <- hg$incidence
  w <- if (identical(method, "clique")) {
    tcrossprod(if (weighted) incidence else .thg_binary(incidence))
  } else {
    b <- .thg_binary(incidence)
    b_projection <- b
    if (identical(duplicate_edges, "collapse") && ncol(b) > 1L) {
      signatures <- vapply(seq_len(ncol(b)), function(j) {
        paste(rownames(b)[which(b[, j] != 0)], collapse = "\r")
      }, character(1L))
      b_projection <- b[, !duplicated(signatures), drop = FALSE]
    }
    cardinality <- Matrix::colSums(b_projection)
    # A singleton hyperedge has no pairs: contribute 0 rather than divide by 0.
    coef <- ifelse(cardinality > 1, 1 / (cardinality - 1), 0)
    tcrossprod(.thg_scale_cols(b_projection, coef), b_projection)
  }
  diag(w) <- 0
  nodes <- rownames(incidence)
  dimnames(w) <- list(nodes, nodes)

  if (self_association) {
    sources <- .thg_resolve_edge_source(hg, edge_source)
    all_nodes <- union(nodes, unique(sources))
    node_index <- match(nodes, all_nodes)
    source_index <- match(sources, all_nodes)
    b_original <- .thg_binary(incidence)
    if (methods::is(b_original, "sparseMatrix")) {
      nz <- Matrix::which(b_original != 0, arr.ind = TRUE)
    } else {
      nz <- which(b_original != 0, arr.ind = TRUE)
    }
    from <- source_index[nz[, "col"]]
    to <- match(nodes[nz[, "row"]], all_nodes)
    keep <- from != to
    from <- from[keep]
    to <- to[keep]
    source_total <- rowsum(rep(1, nrow(nz)), sources[nz[, "col"]],
                           reorder = FALSE)
    denom <- stats::setNames(source_total[, 1L], rownames(source_total))
    contribution <- 1 / unname(denom[sources[nz[, "col"]][keep]])

    sparse_result <- methods::is(w, "sparseMatrix")
    full <- Matrix::Matrix(0, nrow = length(all_nodes), ncol = length(all_nodes),
                           sparse = TRUE,
                           dimnames = list(all_nodes, all_nodes))
    full[node_index, node_index] <- w
    if (length(from)) {
      self_w <- Matrix::sparseMatrix(
        i = c(from, to), j = c(to, from),
        x = c(contribution, contribution),
        dims = c(length(all_nodes), length(all_nodes)),
        dimnames = list(all_nodes, all_nodes)
      )
      full <- full + self_w
    }
    w <- if (sparse_result) full else as.matrix(full)
    nodes <- all_nodes
  }

  diag(w) <- 0
  dimnames(w) <- list(nodes, nodes)
  if (identical(what, "matrix")) return(w)
  .thg_tidy_pairs(w, nodes)
}

# Source-to-member graph of a hypergraph with sources. Rows cite columns in
# the directed matrix; the undirected matrix is its symmetrisation, so a
# pair cited in both directions sums.
.thg_citation_projection <- function(hg, edge_source, duplicate_edges,
                                     directed, what) {
  sources <- .thg_resolve_edge_source(hg, edge_source)
  nodes <- rownames(hg$incidence)
  all_nodes <- sort(union(nodes, unique(sources)))
  b <- .thg_binary(hg$incidence)
  if (methods::is(b, "sparseMatrix")) {
    triplet <- methods::as(methods::as(b, "generalMatrix"), "TsparseMatrix")
    nz <- triplet@x != 0
    row <- triplet@i[nz] + 1L
    col <- triplet@j[nz] + 1L
  } else {
    idx <- which(b != 0, arr.ind = TRUE)
    row <- idx[, 1L]
    col <- idx[, 2L]
  }
  from <- match(sources[col], all_nodes)
  to <- match(nodes[row], all_nodes)
  keep <- from != to
  from <- from[keep]
  to <- to[keep]
  n <- length(all_nodes)
  w <- Matrix::sparseMatrix(i = from, j = to, x = rep(1, length(from)),
                            dims = c(n, n), dimnames = list(all_nodes, all_nodes))
  if (identical(duplicate_edges, "collapse")) w <- (w != 0) * 1
  if (!directed) w <- w + Matrix::t(w)
  w <- Matrix::drop0(w)
  if (!methods::is(hg$incidence, "sparseMatrix")) w <- as.matrix(w)
  if (identical(what, "matrix")) return(w)
  if (directed) {
    triplet <- methods::as(methods::as(w, "sparseMatrix"), "TsparseMatrix")
    out <- data.frame(from = all_nodes[triplet@i + 1L], to = all_nodes[triplet@j + 1L],
                      weight = as.numeric(triplet@x), stringsAsFactors = FALSE)
    out <- out[order(out$from, out$to), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  .thg_tidy_pairs(w, all_nodes)
}

.thg_resolve_edge_source <- function(hg, edge_source) {
  edges <- colnames(hg$incidence) %||% paste0("h", seq_len(hg$n_hyperedges))
  source <- edge_source
  # A single string names an attribute column of the hyperedges.
  if (is.character(source) && length(source) == 1L && is.null(names(source)) &&
      !is.null(hg$edge_data) && source %in% names(hg$edge_data)) {
    source <- stats::setNames(as.character(hg$edge_data[[source]]),
                              as.character(hg$edge_data$edge))
  }
  if (is.null(source) && !is.null(hg$edge_data) &&
      all(c("edge", "source") %in% names(hg$edge_data))) {
    source <- stats::setNames(as.character(hg$edge_data$source),
                              as.character(hg$edge_data$edge))
  }
  if (is.data.frame(source)) {
    if (!all(c("edge", "source") %in% names(source))) {
      .thg_bad_input("an `edge_source` data.frame needs `edge` and `source` columns")
    }
    source <- stats::setNames(as.character(source$source),
                              as.character(source$edge))
  }
  if (is.null(source)) {
    .thg_bad_input(paste0("self-association needs one source per hyperedge; ",
                          "supply `edge_source` or construct a temporal ",
                          "hypergraph with `source =`"))
  }
  if (!is.atomic(source)) .thg_bad_input("`edge_source` must be a vector or data.frame")
  if (!is.null(names(source))) {
    missing_edges <- setdiff(edges, names(source))
    if (length(missing_edges)) {
      .thg_bad_input("the named `edge_source` vector does not cover every hyperedge")
    }
    source <- source[edges]
  } else if (length(source) != length(edges)) {
    .thg_bad_input("an unnamed `edge_source` must have one value per hyperedge")
  }
  source <- as.character(source)
  if (anyNA(source) || any(!nzchar(source))) {
    .thg_bad_input("every hyperedge must have a non-missing source")
  }
  unname(source)
}

#' The s-line graph of a hypergraph
#'
#' Turns the hyperedges into vertices: two hyperedges are adjacent when they
#' share at least `s` vertices, and the edge weight is how many they share.
#' `s = 1` is the ordinary line graph (any overlap connects); raising `s`
#' keeps only substantial overlaps, which is how walk-based measures on
#' hyperedges are parametrised (Aksoy et al. 2020).
#'
#' This is the projection of the dual: `hg_line_graph(hg, s = 1)` returns the
#' same pairs as `pairwise_network(dual_hypergraph(hg), weighted = FALSE)`, which
#' the tests assert.
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any hypergraphs
#'   `net_hg`.
#' @param s Minimum number of shared vertices for two hyperedges to be
#'   adjacent. A single integer of at least 1; `1` (default) is the ordinary
#'   line graph.
#' @param what `"edges"` (default) for the tidy edge list, or `"matrix"` for
#'   the symmetric overlap matrix to hand to a graph engine.
#' @return With `what = "edges"`, a base data.frame with one row per
#'   unordered hyperedge pair sharing at least `s` vertices, sorted by `from`
#'   then `to`, with columns `from`, `to` (hyperedge names) and `weight` (the
#'   number of shared vertices). A zero-row data.frame with those columns when
#'   no pair reaches `s`. With `what = "matrix"`, the symmetric
#'   `n_hyperedges` x `n_hyperedges` overlap matrix, zero on the diagonal and
#'   wherever the overlap is below `s`, sparse if `hg`'s incidence is sparse.
#' @references
#' Aksoy, S. G., Joslyn, C., Ortiz Marrero, C., Praggastis, B., & Purvine, E.
#' (2020). Hypernetwork science via high-order hypergraph walks.
#' *EPJ Data Science*, 9(1), 16. \doi{10.1140/epjds/s13688-020-00231-0}
#' @seealso [pairwise_network()] for the projection onto vertices,
#'   [dual_hypergraph()] for the role swap itself.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
#'                         c = "stars and salt"))
#' hg_line_graph(hg)
#' hg_line_graph(hg, s = 2)
#' @export
hg_line_graph <- function(hg, s = 1, what = c("edges", "matrix")) {
  .thg_check_hg(hg)
  what <- match.arg(what)
  stopifnot(
    "`s` must be a single whole number of at least 1" =
      length(s) == 1L && is.numeric(s) && is.finite(s) && s >= 1 &&
      isTRUE(all.equal(s, round(s)))
  )
  s <- as.integer(round(s))
  b <- .thg_binary(hg$incidence)
  overlap <- crossprod(b)
  diag(overlap) <- 0
  # Multiply by the threshold mask rather than index-assigning, which would
  # densify a sparse overlap matrix.
  overlap <- overlap * (overlap >= s)
  if (methods::is(overlap, "sparseMatrix")) overlap <- Matrix::drop0(overlap)
  edges <- colnames(hg$incidence)
  dimnames(overlap) <- list(edges, edges)
  if (identical(what, "matrix")) return(overlap)
  .thg_tidy_pairs(overlap, edges)
}
