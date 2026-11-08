# ---- hg_centrality() tests ----------------------------------------------

.hc_cycle <- function(n = 4L) {
  # a -> b -> c -> a, first order only
  hon(replicate(n, rep(c("a", "b", "c"), 5), simplify = FALSE),
            max_order = 1L)
}

.hc_planted <- function(n = 6L) {
  hon(c(replicate(n, rep(c("a", "b", "c"), 4), simplify = FALSE),
              replicate(n, rep(c("x", "b", "d"), 4), simplify = FALSE)),
            max_order = 2L)
}

test_that("directed 3-cycle matches the hand-computed centralities", {
  # Nodes a, b, c with edges a->b, b->c, c->a.
  #   PageRank: symmetric cycle -> uniform 1/3.
  #   Betweenness: each ordered pair has one shortest path; the only
  #     two-hop paths are a->b->c, b->c->a, c->a->b, so every node is an
  #     intermediate exactly once.
  #   Harmonic closeness (incoming): distances into a are 1 (from c) and
  #     2 (from b) -> 1 + 1/2 = 1.5.
  cen <- hg_centrality(.hc_cycle())
  expect_identical(cen$state, c("a", "b", "c"))
  expect_equal(cen$pagerank, rep(1 / 3, 3), tolerance = 1e-10)
  expect_equal(cen$betweenness, rep(1, 3), tolerance = 1e-12)
  expect_equal(cen$closeness, rep(1.5, 3), tolerance = 1e-12)
})

test_that("a first-order network projects to itself", {
  hon <- .hc_cycle()
  projected <- hg_centrality(hon, project = TRUE)
  raw <- hg_centrality(hon, project = FALSE)
  expect_identical(projected$state, raw$node)
  expect_identical(raw$order, rep(1L, 3L))
  expect_equal(projected$pagerank, raw$pagerank, tolerance = 1e-12)
  expect_equal(projected$betweenness, raw$betweenness, tolerance = 1e-12)
  expect_equal(projected$closeness, raw$closeness, tolerance = 1e-12)
})

test_that("projected PageRank conserves mass except under 'all'", {
  hon <- .hc_planted()
  results <- lapply(c("scaled", "last", "first"), function(pj) {
    cen <- hg_centrality(hon, type = "pagerank", projection = pj)
    expect_equal(sum(cen$pagerank), 1, tolerance = 1e-9)
    cen
  })
  expect_length(results, 3L)
  # 'all' hands the full value to every state on the path, so a network
  # with higher-order nodes must exceed 1
  all_pr <- hg_centrality(hon, type = "pagerank", projection = "all")
  expect_gt(sum(all_pr$pagerank), 1)
  # the three mass-conserving projections are not all identical here
  expect_false(isTRUE(all.equal(results[[1L]]$pagerank,
                                results[[2L]]$pagerank)))
})

test_that("projections place higher-order mass on the intended states", {
  hon <- .hc_planted()
  raw <- hg_centrality(hon, type = "pagerank", project = FALSE)
  ho <- subset(raw, order == 2L)
  expect_identical(sort(ho$node), c("a -> b", "x -> b"))
  last_pr <- hg_centrality(hon, type = "pagerank", projection = "last")
  first_pr <- hg_centrality(hon, type = "pagerank", projection = "first")
  # "a -> b" and "x -> b" both end in b, so 'last' concentrates their
  # mass on b while 'first' splits it between a and x
  expect_gt(subset(last_pr, state == "b")$pagerank,
            subset(first_pr, state == "b")$pagerank)
  expect_gt(subset(first_pr, state == "a")$pagerank,
            subset(last_pr, state == "a")$pagerank)
})

test_that("memory changes the ranking relative to the first-order network", {
  seqs <- c(replicate(8, rep(c("a", "b", "c"), 4), simplify = FALSE),
            replicate(8, rep(c("x", "b", "d"), 4), simplify = FALSE))
  hon_ho <- hon(seqs, max_order = 2L)
  ho <- hg_centrality(hon_ho)
  hon_fo <- hon(seqs, max_order = 1L)
  fo <- hg_centrality(hon_fo)
  expect_identical(ho$state, fo$state)
  # b is the shared hub: with memory its flow splits across two contexts
  expect_false(isTRUE(all.equal(ho$pagerank, fo$pagerank)))
  expect_lt(subset(ho, state == "b")$pagerank,
            subset(fo, state == "b")$pagerank)
})

test_that("centralities are invariant under state relabeling", {
  relabel <- c(a = "p", b = "q", c = "r")
  seqs <- replicate(6, rep(c("a", "b", "c"), 5), simplify = FALSE)
  hon1 <- hon(seqs, max_order = 2L)
  cen1 <- hg_centrality(hon1)
  hon2 <- hon(lapply(seqs, function(s) unname(relabel[s])),
                    max_order = 2L)
  cen2 <- hg_centrality(hon2)
  expect_identical(unname(relabel[cen1$state]), cen2$state)
  expect_equal(cen1$pagerank, cen2$pagerank, tolerance = 1e-12)
  expect_equal(cen1$betweenness, cen2$betweenness, tolerance = 1e-12)
  expect_equal(cen1$closeness, cen2$closeness, tolerance = 1e-12)
})

test_that("results are non-negative, finite and deterministically ordered", {
  set.seed(31)
  results <- lapply(1:8, function(i) {
    seqs <- replicate(10, sample(letters[1:5], sample(10:25, 1),
                                 replace = TRUE), simplify = FALSE)
    cen <- hg_centrality(hon(seqs, max_order = 3L))
    expect_identical(cen$state, sort(cen$state))
    expect_true(all(is.finite(cen$pagerank)))
    expect_true(all(cen$betweenness >= 0))
    expect_true(all(cen$closeness >= 0))
    expect_equal(sum(cen$pagerank), 1, tolerance = 1e-9)
    cen
  })
  expect_length(results, 8L)
})

test_that("type selection returns exactly the requested columns", {
  hon <- .hc_planted()
  pr_only <- hg_centrality(hon, type = "pagerank")
  expect_identical(names(pr_only), c("state", "pagerank"))
  cl_bt <- hg_centrality(hon, type = c("closeness", "betweenness"))
  expect_identical(names(cl_bt), c("state", "betweenness", "closeness"))
  pr_raw <- hg_centrality(hon, project = FALSE, type = "pagerank")
  expect_identical(names(pr_raw), c("node", "order", "pagerank"))
})

test_that("damping and weighting change PageRank as expected", {
  hon <- .hc_planted()
  d_low <- hg_centrality(hon, type = "pagerank", damping = 0.2)
  d_high <- hg_centrality(hon, type = "pagerank", damping = 0.95)
  expect_false(isTRUE(all.equal(d_low$pagerank, d_high$pagerank)))
  # as damping -> 0 the walk is pure teleportation, so the higher-order
  # PageRank tends to the uniform distribution over higher-order nodes
  # (the projection then splits it across each node's states)
  raw_tiny <- hg_centrality(hon, type = "pagerank", project = FALSE,
                             damping = 1e-6)
  expect_equal(raw_tiny$pagerank, rep(1 / nrow(raw_tiny), nrow(raw_tiny)),
               tolerance = 1e-5)
  wt <- hg_centrality(hon, type = "pagerank", weighted = TRUE)
  expect_equal(sum(wt$pagerank), 1, tolerance = 1e-9)
})

test_that("error paths and the path-enumeration guard", {
  hon <- .hc_planted()
  expect_error(hg_centrality(list(a = 1)), "net_hon")
  expect_error(hg_centrality(hon, damping = 0), "damping")
  expect_error(hg_centrality(hon, damping = 1), "damping")
  expect_error(hg_centrality(hon, tol = 0), "tol")
  expect_error(hg_centrality(hon, project = NA), "project")
  expect_error(hg_centrality(hon, type = "nope"), "arg")
  expect_error(hg_centrality(hon, type = "betweenness", max_paths = 1),
               class = "hypergraphs_too_many_paths")
  expect_warning(
    hg_centrality(hon, type = "pagerank", max_iter = 1L),
    class = "hypergraphs_no_converge")
})

test_that("first-order distances take the minimum over realizations (hand)", {
  # Nodes v1 = "a", v2 = "b -> c" (order 2), v3 = "c"; hop distances
  # v1 -> v2 = 1, v2 -> v3 = 1, v1 -> v3 = 2, diagonal 0, the rest Inf.
  # A hop v -> w realizes first(v) -> last(w) at dist + order(v) - 1:
  #   (a, c): v1 -> v2 gives 1 + 0 = 1, v1 -> v3 gives 2 + 0 = 2 -> min 1
  #   (b, c): v2 itself gives 0 + 1 = 1, v2 -> v3 gives 1 + 1 = 2 -> min 1
  # Every other off-diagonal cell has no realization -> Inf; diagonal 0.
  # The larger candidate of each cell is enumerated second, so keeping the
  # last value written would give 2 in both cells.
  dist_ho <- matrix(c(0, 1, 2,
                      Inf, 0, 1,
                      Inf, Inf, 0), 3L, byrow = TRUE)
  got <- hypergraphs:::.hoc_dist_first(
    dist_ho, hypergraphs:::.hoc_node_paths(c("a", "b -> c", "c")),
    c("a", "b", "c"))
  expect_identical(got, matrix(c(0, Inf, 1,
                                 Inf, 0, 1,
                                 Inf, Inf, 0), 3L, byrow = TRUE,
                               dimnames = list(c("a", "b", "c"),
                                               c("a", "b", "c"))))
})

test_that("sort_by orders by a computed centrality, largest first", {
  hon <- .hc_planted()
  by_pr <- hg_centrality(hon, type = "pagerank", project = FALSE,
                          sort_by = "pagerank")
  expect_true(all(diff(by_pr$pagerank) <= 1e-12))
  by_cl <- hg_centrality(hon, type = c("pagerank", "closeness"),
                          sort_by = "closeness")
  expect_true(all(diff(by_cl$closeness) <= 1e-12))
  # only a requested type may be used as the sort key
  expect_error(hg_centrality(hon, type = "pagerank", sort_by = "closeness"),
               "arg")
})

# ---- State labels holding the arrow (M07) ----------------------------------

test_that("hg_centrality() refuses a network whose states hold the arrow", {
  # a net_hon built outside hon(), which refuses such labels itself
  net <- Nestimate::build_hon(rep(list(c("a -> b", "c", "a -> b", "c")), 2L),
                              max_order = 1L)
  expect_error(hg_centrality(net, type = "pagerank", projection = "last"),
               class = "hypergraphs_bad_input")
})

test_that("projected PageRank conserves its mass on ordinary labels", {
  net <- hon(rep(list(c("a", "b", "c", "a", "b", "d")), 3L), max_order = 2L)
  pr <- hg_centrality(net, type = "pagerank", projection = "last")
  expect_equal(sum(pr$pagerank), 1)
})

# ---- Iteration and path controls (M13) -------------------------------------

test_that("hg_centrality() refuses fractional or infinite iteration counts", {
  net <- hon(rep(list(c("a", "b", "c", "a")), 3L), max_order = 1L)
  expect_error(hg_centrality(net, type = "pagerank", max_iter = 2.5),
               class = "hypergraphs_bad_input")
  expect_error(hg_centrality(net, type = "pagerank", max_iter = Inf),
               class = "hypergraphs_bad_input")
  expect_error(hg_centrality(net, type = "pagerank", tol = Inf),
               class = "hypergraphs_bad_input")
  expect_error(hg_centrality(net, type = "betweenness", max_paths = 0.5),
               class = "hypergraphs_bad_input")
  # Inf lifts the path cap and is accepted
  expect_identical(
    hg_centrality(net, type = "betweenness", max_paths = Inf),
    hg_centrality(net, type = "betweenness"))
})
