testthat::skip_on_cran()

# ---- One verb per idea: dispatch on the input class ----
#
# hg_centrality(), hg_communities() and hypa() pick their estimator from
# the class of `x`. Each method keeps its own arguments and defaults; an
# argument only the other method takes raises hypergraphs_bad_input.

.dp_seqs <- function() {
  list(c("a", "h", "b", "a", "h", "b", "a"),
       c("c", "h", "d", "c", "h", "d", "c"),
       c("a", "h", "b", "a", "h", "b"),
       c("c", "h", "d", "c", "h", "d"))
}

.dp_hg <- function() {
  group_hypergraph(data.frame(
    member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z",
               "a", "x"),
    edge = c(rep(paste0("e", 1:4), each = 3), "e5", "e5")
  ), "member", "edge")
}

test_that("hg_centrality dispatches on memory networks and hypergraphs", {
  h <- hon(.dp_seqs(), max_order = 2L)
  mem <- hg_centrality(h)
  expect_named(mem, c("state", "pagerank", "betweenness", "closeness"))
  expect_identical(hg_centrality(h, type = "pagerank", project = FALSE),
                   hg_centrality.net_hon(h, type = "pagerank",
                                         project = FALSE))
  hg <- .dp_hg()
  expect_named(hg_centrality(hg, type = "clique"), c("node", "clique"))
})

test_that("hg_centrality refuses an argument of the other method", {
  h <- hon(.dp_seqs(), max_order = 2L)
  expect_error(hg_centrality(h, alpha = 0.1), class = "hypergraphs_bad_input")
  expect_error(hg_centrality(h, n = 3), class = "hypergraphs_bad_input")
  expect_error(hg_centrality(.dp_hg(), project = FALSE),
               class = "hypergraphs_bad_input")
  expect_error(hg_centrality(.dp_hg(), top = 2),
               class = "hypergraphs_bad_input")
  expect_error(hg_centrality(data.frame(a = 1)), class = "hypergraphs_bad_input")
})

test_that("hg_communities dispatches on memory networks and hypergraphs", {
  h <- hon(.dp_seqs(), max_order = 2L)
  comm <- hg_communities(h, trials = 2L)
  expect_s3_class(comm, "hypergraphs_memory_communities")
  expect_identical(comm, hg_communities.net_hon(h, trials = 2L))
  skip_if_not_installed("igraph")
  fit <- hg_communities(.dp_hg(), n_runs = 2L, trials = 2L, seeds = 1:2)
  expect_s3_class(fit, "hg_communities")
})

test_that("hg_communities refuses an argument of the other method", {
  h <- hon(.dp_seqs(), max_order = 2L)
  expect_error(hg_communities(h, n_runs = 2L), class = "hypergraphs_bad_input")
  expect_error(hg_communities(h, type = "irmm"), class = "hypergraphs_bad_input")
  expect_error(hg_communities(.dp_hg(), teleportation = 0.1),
               class = "hypergraphs_bad_input")
  expect_error(hg_communities(.dp_hg(), partition = NULL),
               class = "hypergraphs_bad_input")
  expect_error(hg_communities(matrix(1)), class = "hypergraphs_bad_input")
})

test_that("hypa dispatches on hypergraphs and on sequences", {
  pairs <- hypa(.dp_hg(), min_count = 1L)
  expect_s3_class(pairs, "data.frame")
  expect_true(all(c("from", "to", "observed", "expected") %in% names(pairs)))
  paths <- hypa(.dp_seqs(), min_count = 1L)
  expect_s3_class(paths, "net_hypa")
  # the expected counts redistribute the observed ones
  expect_equal(sum(hg_get(paths, what = "scores")$expected),
               sum(hg_get(paths, what = "scores")$count))
  # worked by hand: the 3-grams are a h b 4, h b a 3, b a h 2 and the same on
  # the c/d side. a h and c h each send 4 into h b and h d, which each take 4,
  # so all four routes expect 4 * 4 / 8 = 2; every other route is the only way
  # out of its node and expects what it saw
  scores <- hg_get(paths, what = "scores")
  hand <- c("a -> h -> b" = 2, "c -> h -> d" = 2, "a -> h -> d" = 2,
            "c -> h -> b" = 2, "b -> a -> h" = 2, "d -> c -> h" = 2,
            "h -> b -> a" = 3, "h -> d -> c" = 3)
  expect_identical(sort(scores$path), sort(names(hand)))
  expect_equal(scores$expected, unname(hand[scores$path]))
  expect_equal(scores$count, c(4, 4, 0, 0, 2, 2, 3, 3)[match(scores$path, names(hand))])
})

test_that("hypa refuses an argument of the other method", {
  expect_error(hypa(.dp_hg(), order = 2L), class = "hypergraphs_bad_input")
  expect_error(hypa(.dp_hg(), action = "a"), class = "hypergraphs_bad_input")
  expect_error(hypa(.dp_seqs(), top = 2L), class = "hypergraphs_bad_input")
})

test_that("the dispatching verbs keep each method's defaults", {
  expect_identical(formals(hg_centrality.net_hg)$tol, 1e-8)
  expect_identical(formals(hg_centrality.net_hon)$tol, 1e-12)
  expect_identical(formals(hg_communities.net_hg)$trials, 100L)
  expect_identical(formals(hg_communities.net_hon)$trials, 10L)
  expect_identical(formals(hypa.net_hg)$min_count, 5L)
  expect_identical(formals(hypa.default)$min_count, 5L)
  # both HYPA methods follow LaRock et al. (2020): no correction by default
  expect_identical(formals(hypa.net_hg)$p_adjust, "none")
  expect_identical(formals(hypa.default)$p_adjust, "none")
})

test_that("hypa() on sequences checks its arguments with classed errors", {
  seqs <- .dp_seqs()
  bad <- list(list(alpha = 0.5), list(alpha = 0), list(alpha = NA),
              list(alpha = "0.05"), list(min_count = 0), list(min_count = 1.5),
              list(order = 0), list(order = 1.5), list(order = "2"),
              list(order = integer()), list(p_adjust = "bogus"),
              list(p_adjust = c("BH", "none")), list(n = 0))
  lapply(bad, \(args) {
    expect_error(do.call(hypa, c(list(seqs), args)),
                 class = "hypergraphs_bad_input")
  })
  # the values reach the estimator
  stricter <- hypa(seqs, min_count = 1, alpha = 0.1, p_adjust = "BH")
  expect_identical(stricter$alpha, 0.1)
  expect_identical(stricter$p_adjust, "BH")
  expect_equal(sum(hg_get(stricter, what = "scores")$expected),
               sum(hg_get(stricter, what = "scores")$count))
  # the expected counts do not depend on alpha or the correction: 2 on the four
  # routes that share h (4 * 4 / 8) and the observed count elsewhere
  scores <- hg_get(stricter, what = "scores")
  hand <- c("a -> h -> b" = 2, "c -> h -> d" = 2, "a -> h -> d" = 2,
            "c -> h -> b" = 2, "b -> a -> h" = 2, "d -> c -> h" = 2,
            "h -> b -> a" = 3, "h -> d -> c" = 3)
  expect_equal(scores$expected, unname(hand[scores$path]))
  # BH over 8 tests lifts the two smallest over-tail p-values (tied at ranks 1
  # and 2) to 8 / 2 = 4 times their size, above alpha = 0.1, so nothing is
  # anomalous; uncorrected, both are
  expect_equal(scores$p_adj_over[scores$path == "a -> h -> b"],
               4 * scores$p_over[scores$path == "a -> h -> b"])
  expect_identical(stricter$n_anomalous, 0L)
  expect_identical(hypa(seqs, min_count = 1, alpha = 0.1)$n_anomalous, 2L)
  # a network without its sequences cannot be scored
  net <- Nestimate::build_network(.ho_wide_sequences(seqs), method = "relative")
  net$data <- NULL
  expect_error(hypa(net), class = "hypergraphs_bad_input")
})
