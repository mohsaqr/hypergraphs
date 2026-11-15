# Independent mathematical checks added in the November 2026 audit.

test_that("Hayashi dispersion is stable on dense and sparse incidence", {
  b <- matrix(c(1e12, 1e12 + 1, 0, 1, 0, 2), 3,
              dimnames = list(letters[1:3], c("e1", "e2")))
  expect_equal(.hl_default_edge_weights(b), c(1.5, 1.5))
  expect_equal(.hl_default_edge_weights(Matrix::Matrix(b, sparse = TRUE)),
               c(1.5, 1.5))
  # The population SD of two numbers is half their absolute difference.
  for (scale in c(1e-200, 1, 1e200)) {
    expected <- 0.5 * scale
    expect_equal(.hl_pop_sd(c(1, 2) * scale) / expected, 1)
  }
  # Stored zeros and empty columns are not observations of vertex weights.
  sparse <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 1, 1),
                                 x = c(1e12, 1e12 + 1, 0), dims = c(3, 2))
  expect_equal(.thg_sparse_edge_weights(sparse), c(1.5, 1))
})

test_that("dispersion fixes preserve the EDVW walk on both storage paths", {
  b <- matrix(c(1e12, 1e12 + 1, 0, 1, 0, 2), 3,
              dimnames = list(letters[1:3], c("e1", "e2")))
  dense <- .thg_from_incidence(b, list())
  sparse <- .thg_from_incidence(Matrix::Matrix(b, sparse = TRUE), list())
  # Form the walk by explicitly choosing an edge and then a member.
  p <- matrix(0, 3, 3)
  for (i in 1:3) {
    incident <- which(b[i, ] > 0)
    for (e in incident) p[i, ] <- p[i, ] +
      b[, e] / sum(b[, e]) / length(incident)
  }
  expected <- as.numeric(solve(diag(3) - 0.85 * t(p), rep(0.15 / 3, 3)))
  expect_equal(hg_pagerank(dense)$pagerank, expected, tolerance = 1e-10)
  expect_equal(hg_pagerank(sparse)$pagerank, expected, tolerance = 1e-10)
  ld <- hg_laplacian(dense, type = "random_walk")
  ls <- hg_laplacian(sparse, type = "random_walk")
  expect_equal(as.numeric(ld), as.numeric(ls), tolerance = 1e-10)
  expect_equal(attr(ld, "pi"), attr(ls, "pi"), tolerance = 1e-10)
})

test_that("cosine hypergraphs are invariant to finite positive row scales", {
  x <- matrix(c(1, 2, 3, 2, 3, 4), 3,
              dimnames = list(letters[1:3], NULL))
  reference <- knn_hypergraph(x, k = 1)$incidence
  for (scale in list(rep(1e-200, 3), rep(1e200, 3), c(1e-300, 1, 1e300))) {
    expect_equal(knn_hypergraph(x * scale, k = 1)$incidence, reference,
                 tolerance = 1e-12)
  }
  x[1, ] <- 0
  expect_error(knn_hypergraph(x, k = 1), class = "hypergraphs_bad_input")
  expect_error(knn_hypergraph(matrix(numeric(), 2, 0,
                                    dimnames = list(c("a", "b"), NULL)), k = 1),
               class = "hypergraphs_bad_input")
})

test_that("duplicate association edges are identified by sets, not name separators", {
  dat <- data.frame(node = c("a", "b\rc", "a\rb", "c", "a", "b\rc"),
                    edge = rep(c("e1", "e2", "e3"), each = 2))
  for (sparse in c(FALSE, TRUE)) {
    h <- group_hypergraph(dat, node = "node", hyperedge = "edge", sparse = sparse)
    got <- pairwise_network(h, type = "association", duplicate_edges = "collapse")
    expect_equal(got$weights["a", "b\rc"], 1)
    expect_equal(got$weights["a\rb", "c"], 1)
    expect_equal(sum(got$weights), 4)
  }
})

test_that("temporal distinct-set counts are invariant to node renaming", {
  for (time in list(rep(1, 4), c(1, 2, 1, 2))) {
    dat <- data.frame(node = c("a", "b\rc", "a\rb", "c"),
                      edge = rep(c("e1", "e2"), each = 2), time = time)
    h <- temporal_hypergraph(dat, node = "node", hyperedge = "edge", time = "time")
    got <- hg_growth(h, at = 2, mode = "cumulative")
    expect_equal(got$n_edges_distinct, 2)
    expect_equal(got$n_nodes, 4)
    expect_equal(hg_snapshot(h, at = 2, mode = "cumulative",
                              multiedges = FALSE)$n_hyperedges, 2L)
  }
})

test_that("frequent-set counting accepts arbitrary distinct node names", {
  dat <- data.frame(node = c("a", "b\x01c", "a\x01b", "c", "a", "b\x01c"),
                    edge = rep(c("e1", "e2", "e3"), each = 2))
  got <- hg_get(group_hypergraph(dat, node = "node", hyperedge = "edge", top = Inf),
                 what = "sets")
  expect_equal(nrow(got), 2L)
  expect_equal(sort(got$count), c(1L, 2L))
  expect_equal(got$size, c(2L, 2L))
})

test_that("tensor centralities solve large sunflower equations without underflow", {
  k <- 300L
  members <- lapply(0:2, function(j) c("hub", paste0("p", j, "_", seq_len(k - 1L))))
  dat <- data.frame(node = unlist(members), edge = rep(paste0("e", 1:3), each = k))
  h <- group_hypergraph(dat, node = "node", hyperedge = "edge")
  expect_no_warning(got <- hg_centrality(h, type = c("Z", "H"), tol = 1e-10))
  hub <- match("hub", got$node)
  petal <- match("p0_1", got$node)
  # The tensor equations give r^2 = number of petals (Z),
  # and r^k = number of petals (H), independently of any iterative solver.
  expect_equal(got$Z[hub] / got$Z[petal], sqrt(3), tolerance = 1e-7)
  expect_equal(got$H[hub] / got$H[petal], 3^(1 / k), tolerance = 1e-7)
  expect_equal(sum(got$Z^2), 1)
  expect_equal(sum(got$H^2), 1)
})

test_that("log tensor directions respect exact zeros and H roots", {
  # A map with just one edge has a closed form: each leave-one-out
  # product is its full product divided by that coordinate.
  x <- c(0.5, rep(0.01, 299))
  for (exponent in c(1, 299)) {
    expected <- (1 / x)^(1 / exponent)
    expected <- expected / sqrt(sum(expected^2))
    expect_equal(.hg_tensor_direction(list(1:300), 300L, x, 300L, exponent),
                 expected, tolerance = 1e-12)
    xzero <- x
    xzero[1] <- 0
    expect_equal(.hg_tensor_direction(list(1:300), 300L, xzero, 300L, exponent),
                 c(1, rep(0, 299)))
    xzero[2] <- 0
    expect_equal(.hg_tensor_direction(list(1:300), 300L, xzero, 300L, exponent),
                 rep(0, 300))
  }
})

test_that("Wasserstein rejects invalid homology dimensions", {
  valid <- data.frame(dimension = 0, birth = 0, death = 1)
  for (dimension in c(-1, 0.5, Inf, -Inf, NA_real_)) {
    bad <- valid
    bad$dimension <- dimension
    expect_error(hg_wasserstein(bad, valid), class = "hypergraphs_bad_input")
    expect_error(hg_wasserstein(valid, bad), class = "hypergraphs_bad_input")
  }
})

test_that("assignment validation does not confuse missing labels with text", {
  expect_error(.thg_check_assignment(stats::setNames(c(NA, "\r<NA>"), c("a", "a"))),
               class = "hypergraphs_bad_input")
  expect_identical(.thg_check_assignment(stats::setNames(c("x", "x"), c("a", "a"))),
                   c(a = "x"))
})

test_that("a zero-length clique persistence interval is a diagonal point", {
  zero <- data.frame(dimension = 0L, birth = 0, death = 0)
  empty <- zero[FALSE, ]
  ph <- function(p) structure(list(persistence = p, mode = "clique"),
                              class = "persistent_homology")
  expect_equal(unname(hg_wasserstein(ph(zero), ph(empty))), 0)
  essential <- zero
  essential$birth <- 1
  expect_equal(unname(hg_wasserstein(ph(essential), ph(empty))), Inf)
})

test_that("unused overflowing ground costs do not make Wasserstein infinite", {
  points <- data.frame(dimension = 0L, birth = c(-1e308, 1e308),
                        death = c(-9e307, 1.1e308))
  for (p in c(1, 2, Inf)) {
    expect_equal(unname(hg_wasserstein(points, points, internal_p = p)), 0)
  }
  wide <- data.frame(dimension = 0L, birth = -1e308, death = 1e308)
  empty <- wide[FALSE, ]
  # The L-infinity distance to the diagonal is half the persistence,
  # and the L2 distance is persistence / sqrt(2), both representable here.
  expect_equal(unname(hg_wasserstein(wide, empty)) / 1e308, 1)
  expect_equal(unname(hg_wasserstein(wide, empty, internal_p = 2)) / 1e308,
               sqrt(2))
  expect_equal(unname(hg_wasserstein(wide, empty, internal_p = 1)), Inf)
})

test_that("SBM checks logical controls before using them", {
  for (bad in list(NA, NA_real_, numeric(), c(TRUE, FALSE), "yes", 1)) {
    expect_error(random_hypergraph("sbm", P = matrix(0, 1, 1),
                                   block_sizes = 3, d = 2, variable_size = bad),
                 class = "hypergraphs_bad_input")
    expect_error(random_hypergraph("sbm", P = matrix(0, 1, 1),
                                   block_sizes = 3, d = 2, absolute_purity = bad),
                 class = "hypergraphs_bad_input")
  }
  expect_error(random_hypergraph("sbm", P = matrix(0, 1, 1),
                                 block_sizes = 1, d = 2),
               class = "hypergraphs_bad_input")
})
