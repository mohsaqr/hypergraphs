test_that("G(n,p) generator follows Bernoulli incidence model", {
  h <- random_hypergraph("gnp", n = 12, m = 7, p = 0.3, seed = 17)
  expect_s3_class(h, "net_hg")
  expect_equal(dim(h$incidence), c(12, 7))
  expect_true(all(h$incidence %in% 0:1))
  expect_identical(h, random_hypergraph("gnp", n = 12, m = 7, p = 0.3, seed = 17))
})

test_that("G(n,p) keeps valid limiting and Poisson cases", {
  empty <- random_hypergraph("gnp", n = 5, m = 3, p = 0, seed = 1)
  expect_true(all(empty$incidence == 0))
  full <- random_hypergraph("gnp", n = 5, m = 3, p = 1, seed = 1)
  expect_true(all(full$incidence == 1))
  h <- random_hypergraph("gnp", n = 5, p = 0.2, lambda = 0, seed = 1)
  expect_equal(h$n_hyperedges, 0L)
})

test_that("uniform generator fixes every hyperedge size", {
  h <- random_hypergraph("uniform", n = 20, m = 30, k = 4, seed = 4)
  expect_equal(unname(colSums(h$incidence)), rep(4, 30))
  expect_equal(h$size_distribution, c(size_4 = 30L))
})

test_that("regular generator fixes every node degree", {
  h <- random_hypergraph("regular", n = 20, m = 9, k = 3, seed = 4)
  expect_equal(unname(rowSums(h$incidence)), rep(3, 20))
})

test_that("SBM generator records planted blocks and fixed sizes", {
  P <- matrix(c(0.8, 0.05, 0.05, 0.8), 2, 2)
  h <- random_hypergraph("sbm", P = P, block_sizes = c(8, 8), d = 3, seed = 8)
  expect_s3_class(h, "net_hg")
  expect_equal(as.integer(table(h$blocks)), c(8L, 8L))
  expect_true(all(colSums(h$incidence) == 3))
  expect_identical(h, random_hypergraph("sbm", P = P, block_sizes = c(8, 8),
                                       d = 3, seed = 8))
})

test_that("SBM variable sizes are at least dyadic", {
  P <- matrix(1, 2, 2)
  h <- random_hypergraph("sbm", P = P, block_sizes = c(10, 10), d = 1,
                         variable_size = TRUE, seed = 2)
  expect_true(all(colSums(h$incidence) >= 2))
})

test_that("random generators restore caller RNG state", {
  set.seed(101)
  before <- .Random.seed
  random_hypergraph("gnp", n = 10, m = 4, p = .2, seed = 9)
  expect_identical(.Random.seed, before)
  random_hypergraph("uniform", n = 10, m = 4, k = 2, seed = 9)
  expect_identical(.Random.seed, before)
  random_hypergraph("regular", n = 10, m = 4, k = 2, seed = 9)
  expect_identical(.Random.seed, before)
})

test_that("random generators reject impossible contracts", {
  expect_error(random_hypergraph("gnp", n = 4, m = 2, p = 2), "probability")
  expect_error(random_hypergraph("uniform", n = 4, m = 2, k = 5), "cannot exceed")
  expect_error(random_hypergraph("regular", n = 4, m = 2, k = 3), "cannot exceed")
  expect_error(random_hypergraph("sbm", P = diag(2), block_sizes = c(2, 2), d = 3),
               "not enough nodes")
})

test_that("random_hypergraph() refuses arguments its type does not take", {
  expect_error(random_hypergraph("uniform", n = 10, m = 4, k = 2, p = 0.5),
               class = "hypergraphs_bad_input")
  expect_error(random_hypergraph("gnp", 10, m = 4, p = 0.2),
               class = "hypergraphs_bad_input")
  expect_error(random_hypergraph("dense", n = 10), "should be one of")
})

test_that("random_hypergraph() records its model and prints it", {
  h <- random_hypergraph("regular", n = 6, m = 4, k = 2, seed = 1)
  expect_identical(h$params$source, "random_hypergraph")
  expect_identical(h$params$model, "regular")
  expect_output(print(h), "random regular model")
})
