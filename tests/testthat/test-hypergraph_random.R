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

# Members of hyperedge j (as node indices) of a generated hypergraph.
.sbm_members <- function(h) lapply(seq_len(h$n_hyperedges),
                                   \(j) which(h$incidence[, j] > 0))

test_that("REGRESSION A02: a single eligible added node keeps the size exact", {
  # one block of three, d = 3: every pair has exactly one eligible third node
  for (seed in 1:15) {
    h <- random_hypergraph("sbm", P = matrix(1, 1, 1), block_sizes = 3,
                           d = 3, seed = seed)
    expect_identical(as.integer(colSums(h$incidence)), rep(3L, 3L),
                     info = seed)
  }
})

test_that("REGRESSION A02: a single removable node keeps size and purity", {
  # three blocks, within-block pairs only; d = 3 adds one node, which the
  # impurity of 1 replaces by a node from another block
  P <- diag(3)
  for (seed in 1:15) {
    h <- random_hypergraph("sbm", P = P, block_sizes = c(3, 3, 3), d = 3,
                           impurity = 1L, seed = seed)
    members <- .sbm_members(h)
    expect_true(all(lengths(members) == 3L), info = seed)
    outside <- vapply(members, \(v) {
      block_count <- table(h$blocks[v])
      as.integer(sum(block_count) - max(block_count))
    }, integer(1))
    expect_true(all(outside == 1L), info = seed)
  }
})

test_that("REGRESSION A02: a single outside-block candidate is the one drawn", {
  # pairs only inside block 1 (size 3); block 2 holds one node, V4, the only
  # candidate for the impurity replacement
  P <- matrix(c(1, 0, 0, 0), 2, 2)
  for (seed in 1:15) {
    h <- random_hypergraph("sbm", P = P, block_sizes = c(3, 1), d = 3,
                           impurity = 1L, seed = seed)
    members <- .sbm_members(h)
    expect_true(all(lengths(members) == 3L), info = seed)
    expect_true(all(vapply(members, \(v) 4L %in% v, logical(1))), info = seed)
    expect_true(all(vapply(members, \(v) sum(v <= 3L) == 2L, logical(1))),
                info = seed)
  }
})

test_that("random generator counts refuse overflow and fractions by class", {
  bad <- list(
    list("uniform", n = 3e9, m = 1, k = 1),
    list("uniform", n = 2.5, m = 1, k = 1),
    list("uniform", n = 4, m = NA_real_, k = 1),
    list("regular", n = 4, m = 2, k = Inf),
    list("gnp", n = 4, m = 3e9, p = 0.5),
    list("uniform", n = 4, m = 1, k = 1, seed = 3e9),
    list("sbm", P = diag(1), block_sizes = 3e9, d = 2),
    list("sbm", P = diag(1), block_sizes = 4, d = 3e9),
    list("sbm", P = diag(1), block_sizes = 4, d = 2, impurity = 0.5)
  )
  for (args in bad) {
    expect_error(do.call(random_hypergraph, args),
                 class = "hypergraphs_bad_input",
                 info = paste(names(args), args, collapse = " "))
  }
  expect_error(random_hypergraph("sbm", P = diag(2), block_sizes = c(2, 2),
                                 d = 3), class = "hypergraphs_bad_input")
})
