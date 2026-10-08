# Random hypergraph generators. The model definitions follow HyperG's four
# core samplers, but return hypergraphs' net_hg representation.

#' Random hypergraphs
#'
#' Samples a hypergraph from one of four random models, chosen with `type`.
#'
#' `"gnp"` is the Erdos--Renyi-style model: `n` nodes and `m` hyperedges, with
#' every incidence drawn independently as Bernoulli(`p`). If `m` is omitted
#' it is Poisson with mean `lambda`, or `n * p` when `lambda` is also omitted.
#' Empty and singleton hyperedges are retained, since they are outcomes of the
#' model.
#'
#' `"uniform"` makes every hyperedge contain exactly `k` distinct nodes, and
#' `"regular"` makes every node belong to exactly `k` distinct hyperedges.
#' Sampling weights `prob` can be unequal.
#'
#' `"sbm"` samples a graph stochastic block model and extends every sampled
#' pair of nodes to a hyperedge with further members from the blocks of the
#' pair; `impurity` of the added members can then be replaced by nodes outside
#' those blocks.
#'
#' Each type takes its own arguments, given by name; an argument that the
#' chosen type does not take is an error.
#'
#' @param type The model: `"uniform"` (default), `"regular"`, `"gnp"` or
#'   `"sbm"`.
#' @param ... The arguments of the model, by name:
#'   \describe{
#'     \item{`n`}{Number of nodes (all types; for `"sbm"`, `NULL` uses
#'       `sum(block_sizes)`).}
#'     \item{`m`}{Number of hyperedges (`"uniform"`, `"regular"`, `"gnp"`;
#'       for `"gnp"`, `NULL` draws it from a Poisson distribution).}
#'     \item{`k`}{Hyperedge size (`"uniform"`) or node degree
#'       (`"regular"`).}
#'     \item{`prob`}{Sampling weights, length `n` for `"uniform"` and length
#'       `m` for `"regular"`; `NULL` is equal weights.}
#'     \item{`p`}{Incidence probability in `[0, 1]` (`"gnp"`).}
#'     \item{`lambda`}{Poisson mean of `m` (`"gnp"`).}
#'     \item{`P`}{Symmetric block-to-block pair probability matrix
#'       (`"sbm"`).}
#'     \item{`block_sizes`}{Positive integer block sizes (`"sbm"`).}
#'     \item{`d`}{Hyperedge size, recycled across sampled pairs; with
#'       `variable_size = TRUE`, the Poisson mean of the size minus two
#'       (`"sbm"`).}
#'     \item{`impurity`}{Number of added members replaced by nodes outside
#'       the blocks of the pair (`"sbm"`, default 0).}
#'     \item{`variable_size`}{Draw sizes as `2 + Poisson(d)` (`"sbm"`,
#'       default `FALSE`).}
#'     \item{`absolute_purity`}{Replacements come only from outside the blocks
#'       of the pair (`"sbm"`, default `TRUE`).}
#'   }
#' @param seed Optional seed; the random number state of the caller is
#'   restored on exit.
#' @return A `net_hg` with binary incidence and the model parameters in
#'   `$params`; for `"sbm"`, `$blocks` records the planted block of each
#'   node.
#' @section Conditions:
#' `hypergraphs_bad_input` for an argument the chosen type does not take, an
#' argument given without a name, or a count (`n`, `m`, `k`, `impurity`,
#' `seed`) that is not one whole number in range -- fractions, `NA`, `Inf`
#' and values beyond the integer range are refused, never truncated.
#' @references
#' Marchette, D. J. (2021). HyperG: Hypergraphs in R. R package version
#' 1.0.0.
#' @examples
#' uniform <- random_hypergraph("uniform", n = 20, m = 8, k = 3, seed = 1)
#' hg_measures(uniform, what = "distribution", measure = "size")
#' regular <- random_hypergraph("regular", n = 20, m = 8, k = 2, seed = 1)
#' hg_measures(regular, what = "distribution", measure = "hyperdegree")
#' bernoulli <- random_hypergraph("gnp", n = 20, m = 8, p = 0.2, seed = 1)
#' bernoulli
#' blocks <- matrix(c(0.5, 0.05, 0.05, 0.5), 2, 2)
#' planted <- random_hypergraph("sbm", P = blocks, block_sizes = c(10, 10),
#'                              d = 3, seed = 1)
#' hg_get(planted, what = "nodes")
#' @export
random_hypergraph <- function(type = c("uniform", "regular", "gnp", "sbm"),
                              ..., seed = NULL) {
  type <- match.arg(type)
  sampler <- switch(type, uniform = .hgr_sample_uniform,
                    regular = .hgr_sample_regular, gnp = .hgr_sample_gnp,
                    sbm = .hgr_sample_sbm)
  args <- list(...)
  given <- names(args) %||% rep("", length(args))
  if (any(!nzchar(given))) {
    .thg_bad_input("the arguments of random_hypergraph() are given by name")
  }
  unused <- setdiff(given, setdiff(names(formals(sampler)), "seed"))
  if (length(unused)) {
    .thg_bad_input(sprintf("`type = \"%s\"` does not take %s", type,
                           paste0("`", unused, "`", collapse = ", ")))
  }
  do.call(sampler, c(args, list(seed = seed)))
}

.hgr_sample_gnp <- function(n, m = NULL, p, lambda = NULL, seed = NULL) {
  n <- .hgr_count(n, "n", minimum = 1L)
  .hgr_probability(p, "p")
  if (!is.null(lambda)) {
    if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
      .ho_input_error("`lambda` must be one non-negative number.")
    }
  }
  .hgr_seed(seed)
  if (is.null(m)) m <- stats::rpois(1L, lambda %||% (n * p))
  m <- .hgr_count(m, "m", minimum = 0L)
  # Draw in edge x node order, then transpose. This is distributionally
  # immaterial but preserves exact seeded parity with HyperG's definition.
  incidence <- if (m == 0L) {
    matrix(integer(0), nrow = n, ncol = 0L)
  } else {
    t(matrix(stats::rbinom(n * m, 1L, p), nrow = m, ncol = n))
  }
  dimnames(incidence) <- list(paste0("V", seq_len(n)),
                              if (m) paste0("h", seq_len(m)) else character())
  .hgr_from_incidence(
    incidence, "random_hypergraph",
    list(model = "gnp", n = n, m = m, p = p, lambda = lambda, seed = seed)
  )
}

.hgr_sample_sbm <- function(n = NULL, P, block_sizes, d, impurity = 0L,
                            variable_size = FALSE, absolute_purity = TRUE,
                            seed = NULL) {
  # whole-number test by rounding, never by as.integer(): a size beyond the
  # integer range would coerce to NA and crash the condition
  if (!is.numeric(block_sizes) || !length(block_sizes) ||
      any(!is.finite(block_sizes)) || any(block_sizes < 1) ||
      any(block_sizes > .Machine$integer.max) ||
      any(block_sizes != round(block_sizes))) {
    .ho_input_error("`block_sizes` must contain positive integers.")
  }
  block_sizes <- as.integer(block_sizes)
  n_blocks <- length(block_sizes)
  n <- if (is.null(n)) sum(block_sizes) else .hgr_count(n, "n", 2L)
  if (n != sum(block_sizes)) {
    .ho_input_error("`n` must equal sum(`block_sizes`).")
  }
  if (!is.matrix(P) || !is.numeric(P) ||
      !identical(dim(P), c(n_blocks, n_blocks)) || anyNA(P) ||
      any(P < 0 | P > 1) || !isTRUE(all.equal(P, t(P)))) {
    .ho_input_error(
      "`P` must be a symmetric probability matrix matching the blocks."
    )
  }
  if (!is.numeric(d) || !length(d) || any(!is.finite(d)) || any(d < 0)) {
    .ho_input_error("`d` must contain non-negative finite values.")
  }
  if (!variable_size &&
      any(d < 2 | d != round(d) | d > .Machine$integer.max)) {
    .ho_input_error("fixed `d` values must be integers >= 2.")
  }
  impurity <- .hgr_count(impurity, "impurity", 0L)
  if (!is.logical(variable_size) || length(variable_size) != 1L ||
      !is.logical(absolute_purity) || length(absolute_purity) != 1L) {
    .ho_input_error(
      "`variable_size` and `absolute_purity` must be single logicals."
    )
  }
  .hgr_seed(seed)

  blocks <- rep.int(seq_len(n_blocks), block_sizes)
  pairs <- utils::combn(n, 2L)
  probs <- P[cbind(blocks[pairs[1L, ]], blocks[pairs[2L, ]])]
  keep <- stats::rbinom(ncol(pairs), 1L, probs) == 1L
  pairs <- pairs[, keep, drop = FALSE]
  m <- ncol(pairs)
  target <- if (variable_size) {
    2L + stats::rpois(m, rep(d, length.out = m))
  } else {
    rep(as.integer(d), length.out = m)
  }
  # rpois() returns NA (with a warning) for a mean beyond the integer range
  if (anyNA(target) || any(target > n)) {
    .ho_input_error("sampled hyperedge size exceeds `n`.")
  }

  edges <- lapply(seq_len(m), function(j) {
    endpoints <- pairs[, j]
    endpoint_blocks <- unique(blocks[endpoints])
    eligible <- setdiff(which(blocks %in% endpoint_blocks), endpoints)
    need <- target[j] - 2L
    if (need > length(eligible)) {
      .ho_input_error(
        "not enough nodes in the endpoint blocks for hyperedge size `d`."
      )
    }
    # draw positions, never values: sample(x, size) on a length-one numeric
    # `x` samples from seq_len(x), which would re-draw an endpoint. For a
    # pool longer than one, x[sample.int(length(x), size)] consumes the RNG
    # exactly as sample(x, size) does, so seeded results are unchanged.
    added <- if (need) .hgr_draw(eligible, need) else integer(0)
    edge <- c(endpoints, added)
    replace_n <- min(length(added), impurity)
    # As in HyperG, cross-block dyads in a two-block model cannot receive an
    # outside-block impurity under the absolute-purity rule.
    if (replace_n && !(n_blocks == 2L && length(endpoint_blocks) == 2L)) {
      remove <- .hgr_draw(added, replace_n)
      candidates <- if (absolute_purity) {
        which(!blocks %in% endpoint_blocks)
      } else {
        setdiff(seq_len(n), edge)
      }
      if (length(candidates) < replace_n) {
        .ho_input_error(
          "not enough eligible nodes for the requested `impurity`."
        )
      }
      edge <- c(setdiff(edge, remove), .hgr_draw(candidates, replace_n))
    }
    sort(edge)
  })
  incidence <- .hgr_edges_to_incidence(edges, n)
  out <- .hgr_from_incidence(
    incidence, "random_hypergraph",
    list(model = "sbm", n = n, P = P, block_sizes = block_sizes, d = d,
         impurity = impurity, variable_size = variable_size,
         absolute_purity = absolute_purity, seed = seed)
  )
  out$blocks <- stats::setNames(blocks, out$nodes)
  out
}

.hgr_sample_uniform <- function(n, m, k, prob = NULL, seed = NULL) {
  n <- .hgr_count(n, "n", 1L)
  m <- .hgr_count(m, "m", 0L)
  k <- .hgr_count(k, "k", 0L)
  if (k > n) {
    .ho_input_error("`k` cannot exceed `n` for a uniform hypergraph.")
  }
  prob <- .hgr_sampling_prob(prob, n, k, "n")
  .hgr_seed(seed)
  incidence <- matrix(0L, n, m,
                      dimnames = list(paste0("V", seq_len(n)),
                                      if (m) paste0("h", seq_len(m)) else
                                        character()))
  for (j in seq_len(m)) incidence[sample.int(n, k, prob = prob), j] <- 1L
  .hgr_from_incidence(
    incidence, "random_hypergraph",
    list(model = "uniform", n = n, m = m, k = k, prob = prob, seed = seed)
  )
}

.hgr_sample_regular <- function(n, m, k, prob = NULL, seed = NULL) {
  n <- .hgr_count(n, "n", 1L)
  m <- .hgr_count(m, "m", 1L)
  k <- .hgr_count(k, "k", 0L)
  if (k > m) {
    .ho_input_error("`k` cannot exceed `m` for a regular hypergraph.")
  }
  prob <- .hgr_sampling_prob(prob, m, k, "m")
  .hgr_seed(seed)
  incidence <- matrix(0L, n, m,
                      dimnames = list(paste0("V", seq_len(n)),
                                      paste0("h", seq_len(m))))
  for (i in seq_len(n)) incidence[i, sample.int(m, k, prob = prob)] <- 1L
  .hgr_from_incidence(
    incidence, "random_hypergraph",
    list(model = "regular", n = n, m = m, k = k, prob = prob, seed = seed)
  )
}

# Counts go through the shared validator: fractions, NA, Inf and values
# beyond the integer range are refused with `hypergraphs_bad_input` (an
# as.integer() comparison would coerce an out-of-range value to NA).
#' @noRd
.hgr_count <- function(x, name, minimum) {
  .ho_check_count(x, name, min = minimum)
}

# `size` elements of `pool` drawn without replacement, by position, so a
# length-one pool is never read as the range seq_len(pool).
#' @noRd
.hgr_draw <- function(pool, size) {
  pool[sample.int(length(pool), size)]
}

#' @noRd
.hgr_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    .ho_input_error(sprintf("`%s` must be one probability in [0, 1].", name))
  }
  invisible(x)
}

#' @noRd
.hgr_sampling_prob <- function(prob, size, k, size_name) {
  # Explicit equal weights preserve seeded parity with HyperG's calls to
  # sample(..., prob = rep(1 / size, size)).
  if (is.null(prob)) return(rep(1 / size, size))
  if (!is.numeric(prob) || length(prob) != size || any(!is.finite(prob)) ||
      any(prob < 0) || sum(prob) <= 0 || sum(prob > 0) < k) {
    .ho_input_error(sprintf(
      paste0("`prob` must contain %s non-negative finite weights ",
             "with at least `k` positive values."),
      size_name
    ))
  }
  prob
}

# Set a local RNG scope and restore the exact caller state at function exit.
#' @noRd
.hgr_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  seed <- .hgr_count(seed, "seed", 0L)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  caller <- parent.frame()
  do.call(on.exit, list(substitute({
    if (HAD) {
      assign(".Random.seed", OLD, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, list(HAD = had_seed, OLD = old_seed)), add = TRUE), envir = caller)
  set.seed(seed)
  invisible(NULL)
}

#' @noRd
.hgr_edges_to_incidence <- function(edges, n) {
  m <- length(edges)
  out <- matrix(0L, n, m,
                dimnames = list(paste0("V", seq_len(n)),
                                if (m) paste0("h", seq_len(m)) else
                                  character()))
  if (m && sum(lengths(edges))) {
    out[cbind(unlist(edges, use.names = FALSE),
              rep.int(seq_len(m), lengths(edges)))] <- 1L
  }
  out
}

#' @noRd
.hgr_from_incidence <- function(incidence, source, params) {
  incidence <- as.matrix(incidence)
  storage.mode(incidence) <- "integer"
  n <- nrow(incidence)
  m <- ncol(incidence)
  rownames(incidence) <- rownames(incidence) %||% paste0("V", seq_len(n))
  colnames(incidence) <- colnames(incidence) %||%
    if (m) paste0("h", seq_len(m)) else character()
  edges <- lapply(seq_len(m), function(j) which(incidence[, j] > 0L))
  sizes <- lengths(edges)
  size_dist <- if (length(sizes)) {
    tab <- table(sizes)
    stats::setNames(as.integer(tab), paste0("size_", names(tab)))
  } else integer(0)
  params$source <- source
  structure(list(
    hyperedges = edges,
    incidence = incidence,
    nodes = rownames(incidence),
    n_nodes = n,
    n_hyperedges = m,
    size_distribution = size_dist,
    params = params
  ), class = "net_hg")
}
