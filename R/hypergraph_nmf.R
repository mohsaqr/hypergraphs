# ---- Non-negative matrix-factorization clustering ----------------------

# Internal symmetric NMF solver for Hayashi et al.'s RDC-Sym objective
#   min_{F >= 0} ||S - F F'||_F^2.
# The multiplicative update follows the positive/negative gradient split;
# a square-root step accounts for updating both tied factors simultaneously.
.hl_symnmf <- function(S, k, nstart = 1L, seed = NULL, max_iter = 500L,
                       tol = 1e-6, epsilon = 1e-12, parallel = FALSE,
                       n_cores = 2L) {
  n <- nrow(S)
  objective <- function(F) sum((S - tcrossprod(F))^2)
  # Every start's initial factor is drawn first, in start order (from its
  # own seed, or one after the other from the caller's stream); the updates
  # draw nothing, so the starts can then run in parallel with the result of
  # the serial run.
  scale <- sqrt(max(mean(S), epsilon) / k)
  initial <- lapply(seq_len(nstart), \(restart) {
    if (!is.null(seed)) set.seed(as.integer(seed) + restart - 1L)
    matrix(stats::runif(n * k, 0.5, 1.5) * scale, nrow = n)
  })
  fits <- .ho_apply(seq_len(nstart), \(restart) {
    F <- initial[[restart]]
    history <- objective(F)
    converged <- FALSE
    # Multiplicative updates: each step starts from the previous factor, so
    # the iterations cannot be vectorised; bounded by `max_iter`.
    for (iteration in seq_len(max_iter)) {
      numerator <- S %*% F
      denominator <- F %*% crossprod(F)
      candidate <- F * sqrt(pmax(numerator, 0) /
                              pmax(denominator, epsilon))
      value <- objective(candidate)
      # Numerical guard: backtrack to a monotone non-negative step.
      if (value > history[length(history)] * (1 + 1e-12)) {
        proposal <- candidate
        step <- 0.5
        repeat {
          candidate <- (1 - step) * F + step * proposal
          value <- objective(candidate)
          if (value <= history[length(history)] * (1 + 1e-12) ||
              step < 2^-20) break
          step <- step / 2
        }
      }
      history <- c(history, value)
      change <- abs(history[length(history) - 1L] - value) /
        max(history[length(history) - 1L], epsilon)
      F <- pmax(candidate, epsilon)
      if (change <= tol) {
        converged <- TRUE
        break
      }
    }
    list(factor = F, objective = utils::tail(history, 1L),
         objective_history = history, iterations = iteration,
         converged = converged, restart = restart)
  }, parallel = parallel, n_cores = n_cores)
  # the first start with the lowest objective, as the serial loop kept it
  fits[[which.min(vapply(fits, \(f) f$objective, numeric(1L)))]]
}

# General rectangular non-negative matrix factorization, used to verify the
# gamma = beta = 0 reduction of Hayashi et al.'s Eq. 18.
.hl_nmf <- function(X, k, seed = NULL, max_iter = 500L, tol = 1e-6,
                    epsilon = 1e-12) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  Z <- matrix(stats::runif(nrow(X) * k, 0.5, 1.5), nrow(X), k)
  M <- matrix(stats::runif(ncol(X) * k, 0.5, 1.5), ncol(X), k)
  objective <- function() sum((X - Z %*% t(M))^2)
  history <- objective()
  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    Z <- Z * ((X %*% M) / pmax(Z %*% crossprod(M), epsilon))
    M <- M * ((t(X) %*% Z) / pmax(M %*% crossprod(Z), epsilon))
    value <- objective()
    history <- c(history, value)
    if (abs(history[length(history) - 1L] - value) /
        max(history[length(history) - 1L], epsilon) <= tol) {
      converged <- TRUE
      break
    }
  }
  list(Z = Z, M = M, objective = utils::tail(history, 1L),
       objective_history = history, iterations = iteration,
       converged = converged)
}

#' Joint NMF clustering with an auxiliary vertex relation
#'
#' Implements the two multi-view objectives in Section 6.3.2 of Hayashi et
#' al. (2020). `method = "joint"` is Eq. 18 (J-NMF), jointly factorizing the
#' edge-dependent incidence `X` and a vertex relation `S`. `method =
#' "joint_symmetric"` is Eq. 19 (JS-NMF), jointly factorizing the
#' representative-digraph similarity `T = I - L` and `S`. The paper uses a
#' symmetric patent-citation indicator for `S`; any finite non-negative
#' node-by-node relation matrix can be supplied here.
#'
#' @param hg A connected, dense `net_hg` (a sparse one raises
#'   `hypergraphs_sparse_unsupported`: the factorization is dense).
#' @param relations A non-negative `n_nodes` by `n_nodes` numeric matrix.
#'   If it has dimnames, rows and columns are reordered to `hg$nodes`.
#' @param k Number of clusters, between 2 and `n_nodes - 1`.
#' @param method `"joint"` for J-NMF (Eq. 18) or `"joint_symmetric"` for
#'   JS-NMF (Eq. 19).
#' @param alpha Non-negative symmetry penalty in JS-NMF.
#' @param beta Non-negative agreement penalty between the content and
#'   relation vertex factors.
#' @param gamma Non-negative weight of the auxiliary-relation loss.
#' @param type,edge_weights Laplacian inputs for JS-NMF.
#' @param nstart Number of random initializations, a whole number >= 1.
#' @param seed Optional integer initialization seed.
#' @param max_iter Maximum multiplicative-update iterations, a whole
#'   number >= 1.
#' @param tol Relative objective tolerance, a positive number. A count or
#'   tolerance outside these ranges (including `Inf`) raises
#'   `hypergraphs_bad_input`.
#'
#' @return A `net_hg_cluster` object. Its `$embedding` is the fitted
#'   vertex factor `M`; `$params` contains all fitted factors, the objective
#'   trace and convergence diagnostics.
#' @references
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. \emph{CIKM 2020}, 495-504.
#' \doi{10.1145/3340531.3412034}
#' @export
hg_joint_cluster <- function(
    hg, relations, k, method = c("joint", "joint_symmetric"),
    alpha = 1, beta = 1, gamma = 1,
    type = c("random_walk", "zhou"), edge_weights = NULL,
    nstart = 10L, seed = NULL, max_iter = 500L, tol = 1e-6) {
  method <- match.arg(method)
  type <- match.arg(type)
  .hl_validate_hg(hg)
  .hl_require_dense(hg, "hg_joint_cluster()")
  n <- hg$n_nodes
  k <- .hl_check_k(k, n)
  nstart <- .ho_check_count(nstart, "nstart")
  max_iter <- .ho_check_count(max_iter, "max_iter")
  tol <- .hg_check_tol(tol)
  stopifnot(
    "`relations` must be a finite non-negative numeric matrix" =
      is.matrix(relations) && is.numeric(relations) &&
        all(is.finite(relations)) && all(relations >= 0),
    "`relations` must have n_nodes rows and columns" =
      identical(dim(relations), c(n, n)),
    "`alpha`, `beta`, and `gamma` must be finite non-negative numbers" =
      all(vapply(list(alpha, beta, gamma), function(x) {
        is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
      }, logical(1L)))
  )
  if (!is.null(rownames(relations)) || !is.null(colnames(relations))) {
    if (is.null(rownames(relations)) || is.null(colnames(relations)) ||
        !setequal(rownames(relations), hg$nodes) ||
        !setequal(colnames(relations), hg$nodes)) {
      stop("Named `relations` must contain every `hg$nodes` name in both dimensions.",
           call. = FALSE)
    }
    relations <- relations[hg$nodes, hg$nodes, drop = FALSE]
  }
  S <- relations * 1.0
  best <- NULL
  for (restart in seq_len(as.integer(nstart))) {
    restart_seed <- if (is.null(seed)) NULL else as.integer(seed) + restart - 1L
    fit <- .hl_joint_nmf_fit(
      hg, S, k, method, alpha, beta, gamma, type, edge_weights,
      restart_seed, as.integer(max_iter), tol
    )
    fit$restart <- restart
    if (is.null(best) || fit$objective < best$objective) best <- fit
  }
  assignment <- max.col(best$M, ties.method = "first")
  relabel <- match(assignment, unique(assignment))
  cluster_lab <- paste("Cluster", relabel)
  clusters <- data.frame(node = hg$nodes, cluster = cluster_lab,
                         stringsAsFactors = FALSE)
  dimnames(best$M) <- list(hg$nodes, paste0("dim", seq_len(k)))
  sizes <- as.data.frame(table(cluster = cluster_lab),
                         stringsAsFactors = FALSE)
  names(sizes) <- c("cluster", "size")
  sizes <- sizes[order(sizes$cluster), , drop = FALSE]
  rownames(sizes) <- NULL
  parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
  eigenvalues <- if (method == "joint_symmetric") {
    sort(eigen((parts$L + t(parts$L)) / 2, symmetric = TRUE,
               only.values = TRUE)$values)
  } else rep(NA_real_, n)
  structure(list(
    clusters = clusters, embedding = best$M, k = k, type = type,
    algorithm = method, eigenvalues = eigenvalues,
    eigengap = if (method == "joint_symmetric") {
      eigenvalues[k + 1L] - eigenvalues[k]
    } else NA_real_,
    sizes = sizes, pi = parts$pi, n_nodes = n,
    n_hyperedges = hg$n_hyperedges,
    params = c(list(alpha = alpha, beta = beta, gamma = gamma,
                    nstart = as.integer(nstart), seed = seed,
                    max_iter = as.integer(max_iter), tol = tol), best)
  ), class = "net_hg_cluster")
}

.hl_joint_nmf_fit <- function(hg, S, k, method, alpha, beta, gamma,
                              type, edge_weights, seed, max_iter, tol,
                              epsilon = 1e-12) {
  if (!is.null(seed)) set.seed(seed)
  n <- hg$n_nodes
  M <- matrix(stats::runif(n * k, 0.5, 1.5), n, k)
  Mtilde <- matrix(stats::runif(n * k, 0.5, 1.5), n, k)
  if (method == "joint") {
    X <- t(hg$incidence * 1.0)
    Z <- matrix(stats::runif(nrow(X) * k, 0.5, 1.5), nrow(X), k)
    objective <- function() {
      sum((X - Z %*% t(M))^2) +
        gamma * sum((S - M %*% t(Mtilde))^2) +
        beta * sum((M - Mtilde)^2)
    }
  } else {
    parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
    C <- pmax((diag(n) - parts$L + t(diag(n) - parts$L)) / 2, 0)
    Mhat <- matrix(stats::runif(n * k, 0.5, 1.5), n, k)
    objective <- function() {
      sum((C - M %*% t(Mhat))^2) + alpha * sum((M - Mhat)^2) +
        gamma * sum((S - M %*% t(Mtilde))^2) +
        beta * sum((M - Mtilde)^2)
    }
  }
  history <- objective()
  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    if (method == "joint") {
      Z <- Z * ((X %*% M) / pmax(Z %*% crossprod(M), epsilon))
      M <- M * ((t(X) %*% Z + gamma * S %*% Mtilde + beta * Mtilde) /
        pmax(M %*% crossprod(Z) + gamma * M %*% crossprod(Mtilde) +
               beta * M, epsilon))
    } else {
      M <- M * ((C %*% Mhat + alpha * Mhat + gamma * S %*% Mtilde +
                   beta * Mtilde) /
        pmax(M %*% crossprod(Mhat) + alpha * M +
               gamma * M %*% crossprod(Mtilde) + beta * M, epsilon))
      Mhat <- Mhat * ((t(C) %*% M + alpha * M) /
        pmax(Mhat %*% crossprod(M) + alpha * Mhat, epsilon))
    }
    Mtilde <- Mtilde * ((gamma * t(S) %*% M + beta * M) /
      pmax(gamma * Mtilde %*% crossprod(M) + beta * Mtilde, epsilon))
    value <- objective()
    history <- c(history, value)
    if (abs(history[length(history) - 1L] - value) /
        max(history[length(history) - 1L], epsilon) <= tol) {
      converged <- TRUE
      break
    }
  }
  out <- list(M = M, Mtilde = Mtilde, objective = utils::tail(history, 1L),
              objective_history = history, iterations = iteration,
              converged = converged)
  if (method == "joint") out$Z <- Z else out$Mhat <- Mhat
  out
}
