# EDVW hypergraph PageRank: independent linear-solve reference, the
# Chitra & Raphael collapse theorem, .hg_centrality_fit parity, and
# invariants.

pr_corpus <- c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
)
pr_stop <- c("the", "with", "and", "a", "this", "at", "on", "all")

# Independent reference: build the EDVW transition entry by entry from the
# published formula, then solve the PageRank linear system directly
# (pi = (1 - d) u (I - d P)^{-1}) -- a different algorithm and code path
# from the package's power iteration.
.reference_transition <- function(hg, edge_weights) {
  R <- hg$incidence
  ids <- rownames(R)
  delta <- colSums(R)
  entry <- function(v, w) {
    e_of_v <- which(R[v, ] > 0)
    d_v <- sum(edge_weights[e_of_v])
    sum(vapply(e_of_v, \(e) edge_weights[e] / d_v * R[w, e] / delta[e],
               numeric(1)))
  }
  outer(ids, ids, Vectorize(entry))
}

test_that("pagerank matches a direct linear solve of the published formula", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop, weight = "tfidf")
  omega <- rep(1, hg$n_hyperedges)
  p_ref <- .reference_transition(hg, omega)
  d <- 0.85
  u <- rep(1 / hg$n_nodes, hg$n_nodes)
  pi_ref <- as.numeric((1 - d) * u %*% solve(diag(hg$n_nodes) - d * p_ref))
  pi_ref <- pi_ref / sum(pi_ref)

  out <- hg_pagerank(hg, damping = d, edge_weights = omega, tol = 1e-14)
  expect_identical(out$node, hg$nodes)
  expect_equal(out$pagerank, pi_ref, tolerance = 1e-10)
})

test_that("damping = 1 equals the spectral engine's stationary distribution", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop, weight = "tfidf")
  out <- hg_pagerank(hg, damping = 1, tol = 1e-15, max_iter = 10000L)
  engine <- .hg_cluster_fit(hg, k = 2, type = "random_walk",
                                          seed = 1)
  expect_equal(
    out$pagerank,
    unname(engine$pi[out$node]),
    tolerance = 1e-12
  )
})

test_that("edge-independent weights collapse to the graph walk (Chitra & Raphael)", {
  # gamma_e(v) = gamma(v) for every edge: each vertex carries one weight
  # everywhere it appears. Theorem: the hypergraph walk equals the random
  # walk on the clique-expansion graph with w(u, v) = sum_e omega_e *
  # gamma(u) gamma(v) / delta_e over shared edges (self-loops included), so
  # its stationary distribution is the normalized node strength.
  gamma <- c(a = 1, b = 2, c = 3, d = 1.5)
  long <- data.frame(
    vertex = c("a", "b", "c", "b", "c", "d", "a", "d"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3"),
    w = c(1, 2, 3, 2, 3, 1.5, 1, 1.5)
  )
  hg <- group_hypergraph(long, node = "vertex", hyperedge = "edge",
                                    weight = "w")
  omega <- c(1, 2, 0.5)

  R <- hg$incidence
  delta <- colSums(R)
  # w(u, v) = sum over shared edges of omega_e * gamma(u) * gamma(v) /
  # delta_e; with gamma_e(v) = gamma(v), R[v, e] IS gamma(v) on membership.
  strength_of <- function(u) {
    sum(vapply(seq_len(nrow(R)), \(v) {
      shared <- which(R[u, ] > 0 & R[v, ] > 0)
      sum(omega[shared] * R[u, shared] * R[v, shared] / delta[shared])
    }, numeric(1)))
  }
  strengths <- vapply(seq_len(nrow(R)), strength_of, numeric(1))
  pi_graph <- strengths / sum(strengths)

  out <- hg_pagerank(hg, damping = 1, edge_weights = omega, tol = 1e-15,
                     max_iter = 10000L)
  expect_equal(out$pagerank, pi_graph, tolerance = 1e-12)
})

test_that("pagerank is a probability vector and permutation-invariant", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop, weight = "tfidf")
  out <- hg_pagerank(hg)
  expect_equal(sum(out$pagerank), 1)
  expect_true(all(out$pagerank > 0))

  shuffled <- text_hypergraph(rev(pr_corpus), stop_words = pr_stop,
                              weight = "tfidf")
  out_shuffled <- hg_pagerank(shuffled)
  merged <- merge(out, out_shuffled, by = "node")
  expect_equal(merged$pagerank.x, merged$pagerank.y, tolerance = 1e-12)
})

test_that("personalization concentrates teleport mass", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop)
  uniform <- hg_pagerank(hg)
  focused <- hg_pagerank(hg, personalized = c(cooking_1 = 1))
  merged <- merge(uniform, focused, by = "node", suffixes = c("_u", "_p"))
  target <- subset(merged, node == "cooking_1")
  expect_gt(target$pagerank_p, target$pagerank_u)
  # a character vector of node names is uniform teleport over those nodes
  by_name <- hg_pagerank(hg, personalized = "cooking_1")
  expect_equal(by_name$pagerank, focused$pagerank, tolerance = 1e-12)
  expect_error(hg_pagerank(hg, personalized = "no_such_node"),
               class = "hypergraphs_bad_input")
})

test_that("sort_by and n select without user-side subsetting", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop)
  top <- hg_pagerank(hg, sort_by = "pagerank", n = 2)
  expect_identical(nrow(top), 2L)
  expect_true(all(diff(top$pagerank) <= 0))
})

test_that("contract violations and non-convergence are classed", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop)
  expect_error(hg_pagerank(hg, damping = 0))
  expect_error(hg_pagerank(hg, damping = 1.2))
  expect_error(hg_pagerank(hg, personalized = c(nope = 1)),
               class = "hypergraphs_bad_input")
  expect_error(hg_pagerank(hg, personalized = c(0.5)),
               class = "hypergraphs_bad_input")
  expect_error(hg_pagerank(hg, edge_weights = c(1, 2)))
  expect_error(hg_pagerank(42), class = "hypergraphs_bad_input")
  expect_warning(hg_pagerank(hg, max_iter = 1L), class = "hypergraphs_no_converge")
})

test_that("hg_pagerank equals .hg_centrality_fit(type = \"pagerank\")", {
  # same walk, same damping, uniform teleport: the verb and the engine must
  # agree to solver tolerance on binary and on weighted incidences
  long <- data.frame(
    vertex = c("a", "b", "c", "b", "c", "d", "a", "d", "e", "c", "e"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e3", "e4", "e4"),
    w = c(2, 1, 3, 1, 4, 2, 1, 1, 5, 1, 2)
  )
  binary <- group_hypergraph(long, node = "vertex", hyperedge = "edge")
  weighted <- group_hypergraph(long, node = "vertex", hyperedge = "edge",
                               weight = "w")
  for (hg in list(binary, weighted)) {
    verb <- hg_pagerank(hg, damping = 0.85, tol = 1e-14)
    engine <- .hg_centrality_fit(hg, type = "pagerank", damping = 0.85,
                                    tol = 1e-14)
    expect_identical(verb$node, engine$node)
    expect_equal(verb$pagerank, engine$pagerank, tolerance = 1e-10)
  }
})

test_that("hg_pagerank defaults to window counts like the other walk verbs", {
  # Before 0.6.0 hg_pagerank() ignored window_counts while
  # hg_centrality(type = "pagerank") and the Laplacian used them.
  hg <- window_hypergraph(ring_sequences, window = 3L)
  pr <- hg_pagerank(hg)
  cen <- hg_centrality(hg, type = "pagerank")
  # Two independent power iterations, each stopped at L1 change < 1e-8, so
  # agreement is ~1e-9, not machine precision; the old behaviour missed by
  # 0.018 on this hypergraph.
  expect_lt(max(abs(pr$pagerank - cen$pagerank)), 1e-7)
  explicit <- hg_pagerank(hg, edge_weights = as.numeric(hg$window_counts))
  expect_identical(pr, explicit)
})

# ---- Audit regressions (2026-10-06) ---------------------------------------

.pr_path <- function(sparse = FALSE) {
  group_hypergraph(
    data.frame(node = c("a", "b", "b", "c"), hyperedge = c("X", "X", "Y", "Y")),
    node = "node", hyperedge = "hyperedge", sparse = sparse
  )
}

test_that("A11: personalized refuses non-finite and duplicated weights", {
  hg <- .pr_path()
  bad <- list(c(a = Inf), c(a = NaN), c(a = -Inf, b = 1),
              c(a = 1, a = 9, c = 1), c(a = 9, a = 1, c = 1),
              stats::setNames(1, NA_character_), c(z = 1),
              NA_character_, "")
  lapply(bad, \(p) expect_error(hg_pagerank(hg, personalized = p),
                                class = "hypergraphs_bad_input"))
})

test_that("A11: a character personalization is a set; huge weights are finite", {
  hg <- .pr_path()
  expect_identical(hg_pagerank(hg, personalized = c("a", "a", "c")),
                   hg_pagerank(hg, personalized = c("a", "c")))
  expect_equal(hg_pagerank(hg, personalized = c("a", "c")),
               hg_pagerank(hg, personalized = c(a = 1, c = 1)),
               tolerance = 1e-14)
  # finite weights whose sum overflows a double
  huge <- hg_pagerank(hg, personalized = c(a = 1e308, c = 1e308))
  expect_true(all(is.finite(huge$pagerank)))
  expect_equal(huge, hg_pagerank(hg, personalized = c(a = 1, c = 1)),
               tolerance = 1e-14)
})

test_that("A12: undamped PageRank refuses a disconnected hypergraph", {
  lapply(c(FALSE, TRUE), \(sparse) {
    hg <- group_hypergraph(
      data.frame(node = c("a", "b", "c", "d"),
                 hyperedge = c("X", "X", "Y", "Y")),
      node = "node", hyperedge = "hyperedge", sparse = sparse
    )
    expect_error(hg_pagerank(hg, damping = 1, personalized = "a"),
                 class = "hypergraphs_hypergraph_disconnected")
    expect_error(hg_pagerank(hg, damping = 1),
                 class = "hypergraphs_hypergraph_disconnected")
    # damped, the disconnected graph is still valid
    expect_equal(sum(hg_pagerank(hg, damping = 0.85)$pagerank), 1)
  })
  # connected undamped input still matches the Laplacian stationary vector
  hg <- .pr_path()
  expect_equal(hg_pagerank(hg, damping = 1)$pagerank,
               unname(attr(hg_laplacian(hg, type = "random_walk"), "pi")),
               tolerance = 1e-10)
})

test_that("count and tolerance controls are validated with a classed error", {
  hg <- .pr_path()
  lapply(list(1.5, -Inf, 0, NA, "2", c(1, 2)), \(bad) {
    expect_error(hg_pagerank(hg, n = bad), class = "hypergraphs_bad_input")
  })
  lapply(list(Inf, 2.5, 0), \(bad) {
    expect_error(hg_pagerank(hg, max_iter = bad),
                 class = "hypergraphs_bad_input")
  })
  lapply(list(0, Inf, -1), \(bad) {
    expect_error(hg_pagerank(hg, tol = bad), class = "hypergraphs_bad_input")
  })
  lapply(list(0, 1.5, NA, Inf), \(bad) {
    expect_error(hg_pagerank(hg, damping = bad),
                 class = "hypergraphs_bad_input")
  })
  expect_identical(nrow(hg_pagerank(hg, n = 2)), 2L)
  expect_identical(nrow(hg_pagerank(hg, n = Inf)), 3L)
})

# ---- Audit regressions (2026-11-01) ---------------------------------------

.pr_toy <- function(sparse = FALSE, nodes = NULL) {
  d <- data.frame(node = c("a", "b", "c", "b", "c", "d", "d", "e", "a", "c", "e"),
                  edge = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3",
                           "e4", "e4", "e4"),
                  w = c(1, 2, 1, 3, 1, 1, 2, 2, 1, 1, 4))
  group_hypergraph(d, node = "node", hyperedge = "edge", weight = "w",
                   nodes = nodes, sparse = sparse)
}

test_that("a vertex in no hyperedge teleports instead of raising", {
  # hg_pagerank() used to refuse a hypergraph with an isolated vertex while
  # hg_centrality(type = "pagerank") ranked it by teleportation
  dense <- .pr_toy(nodes = letters[1:6])
  sparse <- .pr_toy(sparse = TRUE, nodes = letters[1:6])
  out <- hg_pagerank(dense)
  expect_equal(sum(out$pagerank), 1, tolerance = 1e-12)
  expect_equal(hg_pagerank(sparse), out, tolerance = 1e-10)
  expect_equal(out$pagerank,
               hg_centrality(dense, type = "pagerank")$pagerank,
               tolerance = 1e-8)
  # by hand: the isolated vertex keeps only (1 - d) / n of teleport mass
  # plus d / n of its own redistributed mass, x_f = (1 - d) / (n (1 - d / n))
  expect_equal(out$pagerank[out$node == "f"], 0.15 / (6 * (1 - 0.85 / 6)),
               tolerance = 1e-10)
  # restart on one vertex: the isolated vertex is never reached
  pers <- hg_pagerank(dense, personalized = "a")
  expect_identical(pers$pagerank[pers$node == "f"], 0)
  expect_equal(sum(pers$pagerank), 1, tolerance = 1e-12)
  # undamped, the stationary distribution is still refused
  expect_error(hg_pagerank(dense, damping = 1),
               class = "hypergraphs_hypergraph_disconnected")
})

test_that("a hypergraph without hyperedges ranks every vertex equally", {
  empty <- matrix(0, 3, 0, dimnames = list(c("a", "b", "c"), NULL))
  dense <- .thg_from_incidence(empty, params = list())
  sparse <- .thg_from_incidence(Matrix::Matrix(empty, sparse = TRUE),
                                params = list())
  expect_equal(hg_pagerank(dense)$pagerank, rep(1 / 3, 3))
  expect_equal(hg_pagerank(sparse)$pagerank, rep(1 / 3, 3))
})

test_that("edge_weights: one number recycles; bad weights are classed on both storages", {
  dense <- .pr_toy()
  sparse <- .pr_toy(sparse = TRUE)
  expect_equal(hg_pagerank(dense, edge_weights = 2),
               hg_pagerank(dense, edge_weights = rep(2, 4)))
  expect_equal(hg_pagerank(sparse, edge_weights = 2),
               hg_pagerank(dense, edge_weights = 2), tolerance = 1e-10)
  # the sparse path used an unclassed stopifnot()
  for (hg in list(dense, sparse)) {
    expect_error(hg_pagerank(hg, edge_weights = c(1, 2, 3, -1)),
                 class = "hypergraphs_bad_input")
    expect_error(hg_pagerank(hg, edge_weights = c(1, 2)),
                 class = "hypergraphs_bad_input")
  }
})
