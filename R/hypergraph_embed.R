#' Low-Dimensional Hypergraph Embedding
#'
#' Returns the node coordinates computed by hypergraphs' existing hypergraph
#' spectral or symmetric-NMF engine without exposing the incidental k-means
#' assignments produced by [hg_cluster()].
#'
#' @param hg Any hypergraphs `net_hg`.
#' @param dimensions Number of embedding coordinates, one whole number
#'   between 2 and `n_nodes - 1`; anything else raises
#'   `hypergraphs_bad_input`.
#' @param type Laplacian type: `"zhou"` or `"random_walk"`.
#' @param method `"spectral"` (default) or Hayashi et al.'s
#'   `"symnmf"`. SymNMF requires a dense incidence matrix.
#' @param seed Optional initialization seed.
#' @param nstart Number of k-means starts for the shared engine; for spectral
#'   embeddings this affects only the discarded cluster assignment, not the
#'   returned coordinates. For SymNMF it is the number of factor restarts.
#' @param max_iter,tol SymNMF iteration controls.
#' @return A data frame with `node`, stationary mass `pi`, and
#'   `dim1` through `dim<dimensions>`.
#' @references
#' Zhou, D., Huang, J., & Schölkopf, B. (2006). Learning with hypergraphs:
#' clustering, classification, and embedding. *Advances in Neural
#' Information Processing Systems 19*. \doi{10.7551/mitpress/7503.003.0205}
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. *Proceedings of the 29th ACM
#' International Conference on Information and Knowledge Management*,
#' 495-504. \doi{10.1145/3340531.3412034}
#' @examples
#' h <- group_hypergraph(
#'   data.frame(node = c("a", "b", "c", "b", "c", "d"),
#'              edge = rep(c("e1", "e2"), each = 3)),
#'   "node", "edge"
#' )
#' hg_embed(h, dimensions = 2)
#' @export
hg_embed <- function(hg, dimensions = 2L,
                     type = c("zhou", "random_walk"),
                     method = c("spectral", "symnmf"), seed = NULL,
                     nstart = 25L, max_iter = 500L, tol = 1e-6) {
  type <- match.arg(type)
  method <- match.arg(method)
  # shared count validator: a fraction, NA, Inf or a value beyond the
  # integer range is refused with `hypergraphs_bad_input` (an as.integer()
  # comparison would turn an out-of-range value into NA and crash)
  dimensions <- .ho_check_count(dimensions, "dimensions", min = 2)
  fit <- hg_cluster(
    hg, k = dimensions, type = type, seed = seed,
    nstart = nstart, what = "embedding", algorithm = method,
    max_iter = max_iter, tol = tol
  )
  coordinates <- grep("^dim[0-9]+$", names(fit), value = TRUE)
  out <- fit[, c("node", "pi", coordinates), drop = FALSE]
  attr(out, "method") <- method
  attr(out, "laplacian") <- type
  out
}
