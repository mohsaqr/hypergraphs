test_that("s-edge centralities equal graph centrality on the s-line graph", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "b", "c", "c", "d"),
               event = rep(c("e1", "e2", "e3"), each = 2)),
    "member", "event"
  )
  out <- hg_edge_centrality(h, s = 1)
  between <- out$value[out$measure == "betweenness"]
  names(between) <- out$edge[out$measure == "betweenness"]
  close <- out$value[out$measure == "closeness"]
  names(close) <- out$edge[out$measure == "closeness"]
  expect_equal(unname(between[c("e1", "e2", "e3")]), c(0, 1, 0))
  expect_equal(unname(close[c("e1", "e2", "e3")]), c(2 / 3, 1, 2 / 3))
})

test_that("edge centrality applies over temporal snapshots", {
  dat <- data.frame(
    member = c("a", "b", "b", "c", "c", "d"),
    event = rep(c("e1", "e2", "e3"), each = 2),
    time = rep(1:3, each = 2)
  )
  thg <- temporal_hypergraph(dat, node = "member", hyperedge = "event", time = "time")
  out <- hg_edge_centrality(thg, measure = "closeness", snapshot_mode = "cumulative")
  # a numeric clock reports the number on the clock
  expect_identical(unique(out$time), c(1, 2, 3))
  expect_equal(nrow(out), 1 + 2 + 3)
})

test_that("a calendar temporal hypergraph reports snapshot times as dates", {
  seats <- data.frame(
    case = rep(c("A", "B"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
    constituted = rep(as.Date(c("2003-01-15", "2004-06-01")), each = 3),
    concluded = rep(as.Date(c("2008-01-01", "2009-01-01")), each = 3)
  )
  thg <- temporal_hypergraph(seats, node = "arbitrator", hyperedge = "case",
                             start = "constituted", end = "concluded")
  at <- as.Date(c("2005-01-01", "2007-01-01"))
  centrality <- hg_edge_centrality(thg, measure = "closeness", at = at)
  expect_s3_class(centrality$time, "Date")
  expect_identical(unique(centrality$time), at)
  motifs <- hg_motifs(thg, n = 5, seed = 1, at = at)
  expect_s3_class(motifs$time, "Date")
  expect_identical(unique(motifs$time), at)
})

test_that("edge centrality uses the paper's disconnected-graph normalization", {
  h <- group_hypergraph(
    data.frame(
      member = c("a", "b", "b", "c", "x", "y"),
      event = rep(c("e1", "e2", "isolated"), each = 2)
    ),
    "member", "event"
  )
  out <- hg_edge_centrality(h, s = 1, normalized = TRUE)
  close <- out[out$measure == "closeness", c("edge", "value")]
  close <- stats::setNames(close$value, close$edge)

  # NetworkX/HypergraphX use the Wasserman--Faust correction: the two
  # connected vertices have raw closeness 1, multiplied by (2 - 1)/(3 - 1).
  expect_equal(unname(close[c("e1", "e2", "isolated")]), c(0.5, 0.5, 0))
})

test_that("subhypergraph centrality has the closed form for one dyad", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b"), event = "e1"), "member", "event"
  )
  out <- .hg_centrality_fit(h, type = "subhypergraph")
  expect_equal(out$subhypergraph, rep(log(cosh(1)), 2), tolerance = 1e-12)
})

test_that("edge centrality counts refuse fractions and overflow by class", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "b", "c", "c", "d"),
               event = rep(c("e1", "e2", "e3"), each = 2)),
    "member", "event"
  )
  for (bad in list(0, 1.5, 3e9, NA_real_)) {
    expect_error(hg_edge_centrality(h, s = bad),
                 class = "hypergraphs_bad_input")
    expect_error(hg_edge_centrality(h, s = 1, top = bad),
                 class = "hypergraphs_bad_input")
  }
  # top applies per measure: two measures, two hyperedges each
  expect_identical(nrow(hg_edge_centrality(h, s = 1, top = 2)), 4L)
})

test_that("a sparse incidence gives the dense s-centralities, isolates at 0", {
  long <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "d", "e", "a", "e", "f", "x", "y"),
    event = c("X", "X", "X", "Y", "Y", "Y", "Z", "Z", "W", "W", "W", "I", "I")
  )
  dense <- group_hypergraph(long, "member", "event")
  sparse <- group_hypergraph(long, "member", "event", sparse = TRUE)
  lapply(1:2, \(s) {
    expect_equal(hg_edge_centrality(sparse, s = s),
                 hg_edge_centrality(dense, s = s), tolerance = 1e-12)
  })
  # hyperedge I shares no member with any other: both measures are 0
  isolated <- subset(hg_edge_centrality(dense, s = 1), edge == "I")
  expect_identical(isolated$value, c(0, 0))
  # s = 3: no pair of hyperedges shares three members, every value is 0
  expect_true(all(hg_edge_centrality(sparse, s = 3)$value == 0))
})

test_that("edge-list components equal the adjacency components", {
  set.seed(3)
  lapply(1:5, \(i) {
    n <- 30L
    pairs <- unique(t(apply(matrix(sample(n, 30, replace = TRUE), ncol = 2), 1,
                            sort)))
    pairs <- pairs[pairs[, 1] != pairs[, 2], , drop = FALSE]
    adjacency <- matrix(0, n, n)
    adjacency[pairs] <- 1
    adjacency <- adjacency + t(adjacency)
    fast <- .thg_edge_list_components(pairs[, 1], pairs[, 2], n)
    slow <- .thg_components(adjacency)
    # same partition: the labels pair up one to one
    expect_identical(nrow(unique(cbind(fast, slow))), length(unique(slow)))
    expect_identical(length(unique(fast)), length(unique(slow)))
  })
  expect_identical(.thg_edge_list_components(integer(), integer(), 4L), 1:4)
})

test_that("many isolated hyperedges stay cheap and score 0", {
  skip_on_cran()
  long <- data.frame(member = paste0("v", 1:1500), event = paste0("e", 1:1500))
  hg <- group_hypergraph(long, "member", "event")
  elapsed <- system.time(out <- hg_edge_centrality(hg, s = 1))[["elapsed"]]
  expect_true(all(out$value == 0))
  expect_lt(elapsed, 2)
})

test_that("a temporal hypergraph with an empty snapshot gives zero rows for it", {
  events <- data.frame(from = c("a", "b"), to = c("b", "c"), time = c(1, 3))
  thg <- temporal_hypergraph(events, from = "from", to = "to", time = "time")
  out <- hg_edge_centrality(thg, at = 1:3)
  expect_identical(sort(unique(out$time)), c(1, 3))
  expect_identical(names(out), c("time", "edge", "s", "measure", "value"))
})
