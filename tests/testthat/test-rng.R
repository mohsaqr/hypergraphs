skip_on_cran()

# A verb that seeds itself gives the caller's random stream back unchanged
# (CRAN policy; the cran-comments say so).

rng_corpus <- c(
  a = "soup salt onion soup broth night", b = "salt soup broth onion night",
  c = "stars sky moon night", d = "sky stars night moon moon",
  e = "soup stars salt sky night broth", f = "onion salt stars broth night")

leaves_rng <- function(expr) {
  set.seed(42)
  before <- .Random.seed
  force(expr)
  identical(before, .Random.seed)
}

test_that("seeded verbs restore the caller's random stream", {
  thg <- text_hypergraph(rng_corpus)
  sparse_thg <- text_hypergraph(rng_corpus, sparse = TRUE)
  seqs <- replicate(6, c("a", "b", "c", "a", "b"), simplify = FALSE)
  model <- hon(seqs, group = rep(c("x", "y"), each = 3), max_order = 2)
  clusters <- hg_cluster(thg, k = 2, seed = 1)
  adjacency <- matrix(1, 5, 5, dimnames = list(letters[1:5], letters[1:5]))
  diag(adjacency) <- 0
  expect_true(leaves_rng(hg_cluster(thg, k = 2, seed = 1)))
  expect_true(leaves_rng(hg_cluster(thg, k = 2, seed = 1,
                                    algorithm = "symnmf")))
  expect_true(leaves_rng(hg_cluster(sparse_thg, k = 2, seed = 1)))
  expect_true(leaves_rng(network_hypergraph(adjacency, p = 0.5, seed = 1)))
  expect_true(leaves_rng(hg_bootstrap(seqs, n_boot = 5, max_order = 2,
                                      seed = 1)))
  expect_true(leaves_rng(hg_compare(model, n_perm = 5, seed = 1)))
  # verbs that seed internally with no seed argument of their own
  expect_true(leaves_rng(hg_membership(thg, clusters)))
  expect_true(leaves_rng(hg_stability(thg, k = 2, resample = "seeds")))
})

test_that("a NULL seed draws from the caller's stream", {
  seqs <- replicate(6, c("a", "b", "c", "a", "b"), simplify = FALSE)
  expect_false(leaves_rng(hg_bootstrap(seqs, n_boot = 5, max_order = 2)))
})
