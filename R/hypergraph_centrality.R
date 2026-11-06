# ---- Hypergraph eigenvector centralities (HON-5) -------------------------
# Three variants from Benson (2019, arXiv:1807.09644):
#   * clique-motif ("CEC") - standard eigenvector centrality on the
#                            clique-expanded pairwise graph
#   * Z-eigenvector ("Z")  - linear tensor eigenvector
#   * H-eigenvector ("H")  - H-eigenvector (power-k-1 recurrence)
#
# The clique variant is the Perron vector of the clique adjacency, solved
# directly by a symmetric eigendecomposition of each connected component;
# Z and H are solved by shifted power iteration on the hyperedge list.

#' Hypergraph eigenvector centralities
#'
#' Computes one or more eigenvector-style centralities on a
#' [net_hg][network_hypergraph]: *clique-motif* (CEC),
#' *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC). Each variant
#' captures influence differently - CEC flattens group structure via
#' clique expansion, while ZEC and HEC propagate through the
#' higher-order groups directly.
#'
#' @param hg A `net_hg` (from [network_hypergraph()],
#'   [group_hypergraph()], or [window_hypergraph()]).
#' @param type Character vector, any subset of
#'   `c("clique", "Z", "H", "pagerank", "subhypergraph")`. The default computes the three
#'   eigenvector variants; request `"pagerank"` explicitly.
#' @param max_iter Maximum number of power-iteration steps for the `"Z"`,
#'   `"H"` and `"pagerank"` variants, a whole number >= 1. Default `1000`.
#'   The `"clique"` variant is solved directly and needs no iteration.
#' @param tol Convergence tolerance on the L1 change between successive
#'   iterates, a positive number. Default `1e-8`.
#' @param normalize Logical. If `TRUE` (default), each returned
#'   centrality vector is L2-normalized to unit norm. With `FALSE` the
#'   vector of type `"clique"` is scaled to a maximum of 1, the scale of
#'   `igraph::eigen_centrality()`. Does not
#'   apply to `"pagerank"`, which always sums to 1, or
#'   `"subhypergraph"`, which is returned on its natural log scale.
#' @param damping Single numeric in (0, 1). PageRank damping factor
#'   (probability of following the walk rather than teleporting).
#'   Default `0.85`. Only used by `type = "pagerank"`.
#' @param edge_weights NULL, a single positive number (recycled to every
#'   hyperedge, e.g. `1` for unit weights), or a positive numeric vector,
#'   one per hyperedge. Hyperedge weights of the EDVW random walk behind
#'   `type = "pagerank"`; `NULL` defaults to the window counts for
#'   hypergraphs built by [window_hypergraph()], else the Hayashi et al.
#'   dispersion heuristic (unit weights on a binary incidence). Only used
#'   by `type = "pagerank"`.
#'
#' @param sort_by `NULL` (node order, default) or one of the requested
#'   `type`s - return the table ranked by that centrality, largest first,
#'   ties broken by node name so the order is deterministic.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#'
#' @return A data.frame, one row per node: `node`, followed by one column
#'   per requested `type`, in the order requested. Compare variants by
#'   reading across the columns; rank by one with `sort_by =`.
#'
#' @details
#' **Clique-motif eigenvector centrality (CEC)**: forms the
#' clique-expanded pairwise graph \eqn{W} where
#' \eqn{W_{ij} = |\{e : i, j \in e\}|} and returns the leading
#' eigenvector of \eqn{W}. With `normalize = FALSE` it equals
#' `igraph::eigen_centrality()` on [pairwise_network()] output. The vector
#' is solved directly from a symmetric eigendecomposition of each connected
#' component (Perron-Frobenius makes the leading eigenvalue of a connected
#' component simple), not by power iteration, which oscillates on bipartite
#' clique expansions. When several components share the largest eigenvalue,
#' the result is the projection of the uniform vector onto their joint
#' leading eigenspace, the limit power iteration from a uniform start
#' reaches; components with a smaller leading eigenvalue score zero.
#'
#' **No pairwise contact.** When no hyperedge has two or more members (no
#' hyperedges at all, or singletons only), the clique adjacency is zero and
#' every eigenvector of it has eigenvalue zero, so no vector is singled out.
#' `"clique"`, `"Z"` and `"H"` then return zero for every node, the same
#' convention as a hypergraph without hyperedges.
#'
#' **Z-eigenvector centrality (ZEC)**: solves the linear
#' eigen-equation on the hyperedge tensor,
#' \deqn{\lambda\, x_i \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j,}
#' via shifted power iteration (SS-HOPM). An iteration (here, for `"H"`,
#' and for `"pagerank"`) that exhausts `max_iter` warns with class
#' `hypergraphs_no_converge` and returns its last iterate.
#'
#' **H-eigenvector centrality (HEC)**: solves the power-k-1
#' eigen-equation,
#' \deqn{\lambda\, x_i^{k-1} \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j.}
#'
#' **Z and H need a uniform hypergraph.** Benson (2019) defines both for
#' \eqn{k}-uniform hypergraphs. With mixed hyperedge sizes the right-hand
#' side mixes polynomial degrees, so the equation has no scale-free
#' solution: the vector satisfies it only at the normalisation chosen inside
#' the iteration. Requesting `"Z"` or `"H"` when the hyperedges with two or
#' more members have more than one size raises
#' `hypergraphs_nonuniform_hypergraph`; singletons are ignored by the check.
#'
#' **Disconnected hypergraphs.** The Z-eigenvector is not unique when the
#' hypergraph has several components: different iterations can solve the
#' equation with different vectors, so compare `"Z"` values within one connected
#' hypergraph only.
#'
#' **Hypergraph PageRank** (`"pagerank"`): the stationary distribution of
#' the damped EDVW random walk of Chitra & Raphael (2019): from node
#' \eqn{v}, pick a hyperedge \eqn{e \ni v} with probability proportional
#' to its weight \eqn{w(e)}, then a node \eqn{u \in e} with probability
#' proportional to its edge-dependent vertex weight \eqn{\gamma_e(u)}
#' (the incidence cell, i.e. occurrence totals for
#' [window_hypergraph()]); with probability \eqn{1 - damping} teleport
#' uniformly. Their collapse theorem: with edge-*independent* vertex
#' weights (a binary incidence) the walk is equivalent to PageRank on the
#' weighted clique expansion with edge weights
#' \eqn{\sum_{e \ni u,v} w(e)/\delta(e)} -- the hypergraph adds
#' information exactly when \eqn{\gamma} is edge-dependent. Nodes left in
#' no hyperedge (possible after `min_weight`/`min_size` filtering)
#' teleport from every step and receive only teleportation mass. The
#' undamped stationary distribution of the same walk is the `pi` column
#' reported by [hg_cluster()].
#'
#' **Subhypergraph centrality** (`"subhypergraph"`): the logarithm of the
#' diagonal of the matrix exponential of the clique adjacency derived from
#' binary incidence (so entries count shared hyperedges),
#' \eqn{\log[\exp(W)]_{ii}}. It counts closed walks based at each node with a
#' factorial penalty for length and matches the implementation used by
#' HypergraphX 1.5 in the legal-hypergraphs analysis. Each diagonal entry
#' is evaluated in log space, row by row and component by component
#' (\eqn{\log \sum_j v_{ij}^2 e^{\lambda_j}} by log-sum-exp over the
#' eigenpairs of the node's own component), so a node outside the
#' component with the largest eigenvalue stays finite; a node in no
#' pairwise contact scores exactly 0.
#'
#' @seealso [network_hypergraph()], [pairwise_network()],
#'   [hg_measures()].
#'
#' @references
#' Benson, A. R. (2019). Three hypergraph eigenvector centralities.
#' \emph{SIAM Journal on Mathematics of Data Science} 1(2), 293-312.
#' arXiv:1807.09644.
#'
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. \emph{Proceedings of the 36th
#' International Conference on Machine Learning}, PMLR 97, 1172-1181.
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. \emph{Proceedings of CIKM
#' 2020}, 495-504. \doi{10.1145/3340531.3412034}
#'
#' Estrada, E., & Rodriguez-Velazquez, J. A. (2005). Complex networks as
#' hypergraphs. *arXiv preprint physics/0505137*.
#'
#' @note The `"clique"`, `"Z"`, and `"H"` variants have an independent
#'   XGI equivalence runner in `local_testing_and_equivalence/`; clique is
#'   also validated against `igraph::eigen_centrality` (cosine ~ 1).
#'   `"pagerank"` is checked against
#'   `igraph::page_rank` on the collapsed graph of the edge-independent
#'   case plus a dense linear-system solve. Clean-room list-based tensor
#'   equations remain in the unit suite for deterministic validation without
#'   Python.
#'
#' @noRd
.hg_centrality_fit <- function(hg,
                                   type     = "clique",
                                   max_iter = 1000L,
                                   tol      = 1e-8,
                                   normalize = TRUE,
                                   damping = 0.85,
                                   edge_weights = NULL,
                                   sort_by = NULL,
                                   top = NULL) {
  .thg_check_hg(hg)
  if (!(is.logical(normalize) && length(normalize) == 1L &&
        !is.na(normalize))) {
    .ho_input_error("`normalize` must be TRUE or FALSE")
  }
  damping <- .ho_check_number(damping, "damping", min = 0, max = 1)
  if (!(damping > 0 && damping < 1)) {
    .ho_input_error("`damping` must be strictly between 0 and 1")
  }
  max_iter <- .ho_check_count(max_iter, "max_iter")
  tol <- .hg_check_tol(tol)
  type <- match.arg(type, choices = c("clique", "Z", "H", "pagerank",
                                      "subhypergraph"),
                    several.ok = TRUE)

  n     <- hg$n_nodes
  nodes <- hg$nodes

  # Degenerate: no nodes -> a zero-row table with the requested columns
  if (n == 0L) {
    res <- data.frame(node = character(0L), stringsAsFactors = FALSE)
    res[type] <- rep(list(numeric(0L)), length(type))
    return(res)
  }

  # Degenerate: no hyperedges -> zero eigenvector centralities; PageRank
  # is pure teleportation, so it is uniform (it must still sum to 1)
  if (hg$n_hyperedges == 0L) {
    res <- data.frame(node = nodes, stringsAsFactors = FALSE)
    res[type] <- rep(list(rep(0, n)), length(type))
    if ("pagerank" %in% type) res$pagerank <- rep(1 / n, n)
    return(.ho_top(res, top))
  }

  # Shared initialization: uniform positive vector
  x0 <- rep(1 / sqrt(n), n)

  # Hyperedges as integer vector lists; pre-compute for Z/H
  hyperedges <- hg$hyperedges
  edge_sizes <- vapply(hyperedges, length, integer(1L))

  # Z and H are defined for uniform hypergraphs; singletons contribute nothing
  tensor_sizes <- unique(edge_sizes[edge_sizes >= 2L])
  if (any(c("Z", "H") %in% type) && length(tensor_sizes) > 1L) {
    stop(errorCondition(
      sprintf(paste0(
        "`type = \"Z\"` and `\"H\"` need a uniform hypergraph; ",
        "hyperedge sizes here are %s. Use `type = \"clique\"` or ",
        "`\"pagerank\"`, or build the hypergraph with one hyperedge size."),
        paste(sort(tensor_sizes), collapse = ", ")),
      class = "hypergraphs_nonuniform_hypergraph", call = NULL
    ))
  }
  k_max      <- max(edge_sizes)

  # No hyperedge with two or more members: the clique adjacency and the
  # tensor are zero, every vector is an eigenvector of eigenvalue 0, and
  # the convention is the zero result of the no-hyperedge branch above.
  # Guarded here, before any iteration could pick an arbitrary vector.
  pairwise <- k_max >= 2L

  # A node in no hyperedge of two or more members has f(x)_i = 0 in the
  # tensor map, so the eigen-equation lambda x_i = 0 forces x_i = 0. Its
  # start value is set to 0, where the shifted iteration keeps it exactly;
  # a positive start would only decay geometrically and stop at the
  # tolerance (about 1e-9) instead of 0.
  in_tensor <- tabulate(as.integer(unlist(hyperedges[edge_sizes >= 2L],
                                          use.names = FALSE)),
                        nbins = n) > 0
  x0_tensor <- ifelse(in_tensor, x0, 0)

  out <- list()

  # Clique adjacency (shared hyperedges per node pair, zero diagonal) and
  # its per-component eigendecomposition, shared by "clique" and
  # "subhypergraph".
  if (any(c("clique", "subhypergraph") %in% type)) {
    adjacency <- tcrossprod((hg$incidence != 0) * 1.0)
    diag(adjacency) <- 0
    spectra <- .hg_component_spectra(adjacency)
  }

  # ---- CEC: Perron vector of the clique expansion W ----
  if ("clique" %in% type) {
    x <- if (pairwise) .hg_perron_vector(spectra, x0) else rep(0, n)
    if (!normalize && any(x != 0)) x <- x / max(abs(x))
    out$clique <- stats::setNames(x, nodes)
  }

  # ---- ZEC: lambda x = sum_{e contains i} prod_{j in e,j!=i} x_j ----
  if ("Z" %in% type) {
    out$Z <- stats::setNames(
      if (pairwise) {
        .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                              exponent = 1L, max_iter = max_iter, tol = tol,
                              x0 = x0_tensor, normalize = normalize,
                              label = "Z")
      } else {
        rep(0, n)
      },
      nodes
    )
  }

  # ---- HEC: lambda x^{k-1} = sum_{e contains i} prod_{j in e,j!=i} x_j ----
  if ("H" %in% type) {
    out$H <- stats::setNames(
      if (pairwise) {
        .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                              exponent = k_max - 1L,
                              max_iter = max_iter, tol = tol,
                              x0 = x0_tensor, normalize = normalize,
                              label = "H")
      } else {
        rep(0, n)
      },
      nodes
    )
  }

  # ---- PageRank: damped EDVW random walk (Chitra & Raphael 2019) ----
  if ("pagerank" %in% type) {
    rw <- .hl_rw_transition(hg, edge_weights)
    P <- rw$P
    # Nodes in no hyperedge (d_v = 0) have undefined walk rows: they
    # teleport uniformly from every step (standard dangling-node fix)
    dangling <- rw$d_v == 0
    if (any(dangling)) P[dangling, ] <- 1 / n
    x <- rep(1 / n, n)
    converged <- FALSE
    # power iteration: an inherently sequential fixed-point loop
    for (iter in seq_len(max_iter)) {
      y <- damping * as.numeric(crossprod(P, x)) + (1 - damping) / n
      if (sum(abs(y - x)) < tol) {
        x <- y
        converged <- TRUE
        break
      }
      x <- y
    }
    if (!converged) {
      warning(warningCondition(
        sprintf(
          "PageRank did not converge in %d iterations (L1 change > %g).",
          as.integer(max_iter), tol),
        class = "hypergraphs_no_converge"
      ))
    }
    out$pagerank <- stats::setNames(x, nodes)
  }

  # ---- Log subhypergraph centrality (Estrada & Rodriguez-Velazquez) ----
  # log [exp(W)]_ii = log sum_j V_ij^2 exp(lambda_j), by a row-wise
  # log-sum-exp over the eigenpairs of the node's own component: a global
  # shift by the largest eigenvalue underflows every row with no mass on
  # the dominant eigenspace to log(0) = -Inf. A one-node component has the
  # 1 x 1 zero matrix, exp(0) = 1, so it scores exactly 0.
  if ("subhypergraph" %in% type) {
    value <- numeric(n)
    pieces <- lapply(spectra, \(piece) {
      if (length(piece$index) == 1L) return(0)
      log_terms <- sweep(2 * log(abs(piece$vectors)), 2L, piece$values, `+`)
      top <- apply(log_terms, 1L, max)
      top + log(rowSums(exp(log_terms - top)))
    })
    value[unlist(lapply(spectra, `[[`, "index"))] <- unlist(pieces)
    out$subhypergraph <- stats::setNames(value, nodes)
  }

  # One row per node, one column per requested type, in the order the user
  # asked for them.
  res <- data.frame(node = nodes, stringsAsFactors = FALSE)
  res[type] <- lapply(out[type], unname)
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, type)
    res <- res[.ho_order(-res[[sort_by]], res$node), , drop = FALSE]
    rownames(res) <- NULL
  }
  .ho_top(res, top)
}

# A convergence tolerance: one finite number > 0 (0 would never be met and
# Inf would stop after the first step). Shared by the iterative verbs of the
# hypergraph family.
.hg_check_tol <- function(tol, arg = "tol") {
  tol <- .ho_check_number(tol, arg, min = 0)
  if (!(tol > 0)) .ho_input_error(sprintf("`%s` must be > 0", arg))
  tol
}

# Connected-component labels of a symmetric adjacency by minimum-label
# propagation over its non-zero cells: every pass lets each node adopt the
# smallest label among itself and its neighbours, so the loop runs at most
# diameter + 1 times (sequential by nature). Labels are 1, 2, ... in order
# of first appearance.
.hg_component_labels <- function(adjacency) {
  n <- nrow(adjacency)
  label <- seq_len(n)
  cells <- which(adjacency != 0, arr.ind = TRUE)
  if (nrow(cells) == 0L) return(label)
  row_of <- factor(cells[, 1L], levels = seq_len(n))
  repeat {
    neighbour <- vapply(
      split(label[cells[, 2L]], row_of),
      \(v) if (length(v)) min(v) else .Machine$integer.max, integer(1L)
    )
    nxt <- pmin(label, neighbour)
    if (identical(nxt, label)) break
    label <- nxt
  }
  match(label, unique(label))
}

# Symmetric eigendecomposition of each connected component of a clique
# adjacency: a list with one element per component, `index` (its nodes),
# `values` (decreasing) and `vectors` (columns). By Perron-Frobenius the
# leading eigenvalue of a connected non-negative component is simple and its
# vector has one sign, so the leading vector is well defined component by
# component even where the whole matrix has a repeated leading eigenvalue
# (equal components) or a +/- pair (bipartite components).
.hg_component_spectra <- function(adjacency) {
  labels <- .hg_component_labels(adjacency)
  lapply(split(seq_len(nrow(adjacency)), labels), \(index) {
    eig <- eigen(adjacency[index, index, drop = FALSE], symmetric = TRUE)
    list(index = index, values = eig$values, vectors = eig$vectors)
  })
}

# Unit-L2 leading eigenvector of a clique adjacency from its component
# spectra. Components whose leading eigenvalue equals the global one (to
# rounding) share the leading eigenspace; the result is the projection of
# the start vector `x0` onto that space -- exactly the limit of power
# iteration from `x0` -- and every other node is 0.
.hg_perron_vector <- function(spectra, x0) {
  lead <- vapply(spectra, \(piece) piece$values[1L], numeric(1L))
  top_value <- max(lead)
  if (!(top_value > 0)) return(rep(0, length(x0)))
  near <- abs(lead - top_value) <= sqrt(.Machine$double.eps) * top_value
  x <- numeric(length(x0))
  parts <- lapply(spectra[near], \(piece) {
    v <- piece$vectors[, 1L]
    v <- if (sum(v) < 0) -v else v
    # a Perron vector is non-negative; clamp rounding noise below zero
    v <- pmax(v, 0)
    sum(v * x0[piece$index]) * v
  })
  x[unlist(lapply(spectra[near], `[[`, "index"))] <- unlist(parts)
  x / sqrt(sum(x^2))
}

# The tensor map f(x)_i = sum_{e contains i} prod_{j in e, j != i} x_j, with
# the (exponent)-th root taken for H. Accumulation over hyperedges with
# in-place adds into y: the per-edge leave-one-out products have no
# vectorised form across ragged edges.
#' @noRd
.hg_tensor_apply <- function(hyperedges, edge_sizes, x, n, exponent) {
  y <- numeric(n)
  for (e_idx in seq_along(hyperedges)) {
    e  <- hyperedges[[e_idx]]
    ke <- edge_sizes[e_idx]
    if (ke < 2L) next
    x_e <- x[e]
    if (any(x_e == 0)) {
      nz_idx <- which(x_e != 0)
      if (length(nz_idx) == ke - 1L) {
        zero_pos <- setdiff(seq_len(ke), nz_idx)
        y[e[zero_pos]] <- y[e[zero_pos]] + prod(x_e[nz_idx])
      }
      # With >=2 zeros in x_e, every leave-one-out product is zero.
    } else {
      y[e] <- y[e] + prod(x_e) / x_e
    }
  }
  if (exponent > 1L) {
    # k-1-root, preserving sign (for non-negative hypergraphs all >= 0)
    y <- sign(y) * abs(y)^(1 / exponent)
  }
  y
}

# Shared tensor-power-iteration kernel (Kolda & Mayo 2011 SS-HOPM shift):
# x_{t+1} ~ f(x_t) / ||f(x_t)|| + shift * x_t, whose fixed points are exactly
# the eigenvectors f(x) = c x. The update is normalised before the shift is
# added: for large hyperedges f(x) of a unit vector is a product of many
# small numbers and an unnormalised shift swamped it, so the iterate passed
# the step test without solving the equation (a 3-petal sunflower of
# 20-member hyperedges gave a hub/petal ratio of 1 instead of sqrt(3)).
# Convergence also requires the eigen-equation residual to be small;
# otherwise it warns with class `hypergraphs_no_converge` (naming `label`)
# and returns the last iterate. Same kernel as Nestimate (>= 0.9.28).
#
# exponent = 1    => Z-eigenvector (no post-root)
# exponent = k-1  => H-eigenvector (k-1-root)
#' @noRd
.hg_tensor_power_iter <- function(hyperedges, edge_sizes, n, exponent,
                                   max_iter, tol, x0, normalize,
                                   shift = 1, label = "tensor") {
  x <- x0
  step_converged <- FALSE
  # power iteration: an inherently sequential fixed-point loop
  for (iter in seq_len(max_iter)) {
    fx <- .hg_tensor_apply(hyperedges, edge_sizes, x, n, exponent)
    nf <- sqrt(sum(fx^2))
    if (!(nf > 0)) {
      x <- fx
      step_converged <- TRUE
      break
    }
    y <- fx / nf + shift * x
    y <- y / sqrt(sum(y^2))
    if (sum(abs(y - x)) < tol) {
      x <- y
      step_converged <- TRUE
      break
    }
    x <- y
  }

  # Eigen-equation residual ||f(x) - (x'f(x)) x|| / ||f(x)|| at the result
  residual_ok <- TRUE
  if (any(x != 0)) {
    fx <- .hg_tensor_apply(hyperedges, edge_sizes, x, n, exponent)
    nf <- sqrt(sum(fx^2))
    if (nf > 0) {
      residual <- sqrt(sum((fx - sum(x * fx) * x)^2)) / nf
      residual_ok <- residual <= max(sqrt(tol), 1e-6)
    }
  }
  if (!step_converged || !residual_ok) {
    warning(warningCondition(
      sprintf(paste0("%s-eigenvector centrality did not converge in %d ",
                     "iterations (tol = %g): the eigen-equation is not ",
                     "solved; returning the last iterate."),
              label, as.integer(max_iter), tol),
      class = "hypergraphs_no_converge"
    ))
  }

  if (any(x != 0) && sum(x) < 0) x <- -x
  if (!normalize && any(x != 0)) x <- x / max(abs(x))
  x
}
