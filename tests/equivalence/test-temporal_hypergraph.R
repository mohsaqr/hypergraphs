# temporal_hypergraph(), hg_snapshot(), hg_snapshots(), hg_growth(), hg_subset(): hyperedges that begin at a
# time, and the hypergraph seen in a window. A positive window from t holds the hyperedges that begin in
# [t, t + window); the cumulative view at t holds every hyperedge begun before t + window, the window defaulting
# to the step (both rules as documented in ?hg_snapshots). Reference: the membership table filtered by those rules
# in base R.

events <- data.frame(
  node = c("a", "b", "c", "a", "d", "b", "c", "e", "a", "e", "c", "d", "f", "b", "f"),
  hyperedge = c("h1", "h1", "h1", "h2", "h2", "h3", "h3", "h3", "h4", "h4", "h5", "h5", "h5", "h6", "h6"),
  time = c(1, 1, 1, 3, 3, 5, 5, 5, 8, 8, 8, 8, 8, 11, 11))
temporal <- temporal_hypergraph(events, node = "node", hyperedge = "hyperedge", time = "time")
grid <- seq(0, 10, by = 2)

#' The memberships of the hyperedges that begin in [from, to)
in_window <- function(from, to) subset(events, time >= from & time < to)

#' Counts of a membership table
counts <- function(members) {
  c(n_nodes = length(unique(members$node)), n_edges = length(unique(members$hyperedge)), n_memberships = nrow(members))
}

test_that("the stored memberships and start times equal the event table", {
  expect_reference(hg_get(temporal, what = "memberships"),
                   data.frame(node = events$node, edge = events$hyperedge, start = events$time, weight = 1),
                   reference = "independent: the event table", covers = c("temporal_hypergraph()", "hg_get.net_temporal_hypergraph()"),
                   quantity = "memberships", tolerance = 0, key = c("node", "edge"))
  first <- stats::aggregate(list(start = events$time), list(edge = events$hyperedge), min)
  expect_reference(hg_get(temporal, what = "edges"), first, reference = "independent: the event table",
                   covers = c("temporal_hypergraph()", "hg_get.net_temporal_hypergraph()"), quantity = "hyperedge start",
                   tolerance = 0, key = "edge")
})

test_that("a snapshot holds the hyperedges begun in [t, t + window)", {
  invisible(lapply(c(1, 5, 8), \(at) {
    members <- in_window(at, at + 2)
    sizes <- table(members$hyperedge)
    expected <- data.frame(hyperedge = names(sizes), size = as.vector(sizes))
    expect_reference(hg_get(hg_snapshot(temporal, at = at, window = 2), what = "edges"), expected,
                     reference = "independent: the event table filtered to the window", covers = "hg_snapshot()",
                     quantity = sprintf("snapshot at %g, window 2: hyperedge sizes", at), tolerance = 0, key = "hyperedge")
  }))
})

test_that("snapshots on a grid hold the hyperedges begun in each window", {
  snapshots <- hg_snapshots(temporal, start = 0, end = 10, step = 2, window = 2)
  # the window from 6 holds no hyperedge: its empty table is left out (the growth test below counts it)
  shown <- grid[grid != 6]
  invisible(Map(\(snapshot, at) {
    members <- in_window(at, at + 2)
    sizes <- table(members$hyperedge)
    expect_reference(hg_get(snapshot, what = "memberships"), data.frame(node = members$node, hyperedge = members$hyperedge, weight = 1),
                     reference = "independent: the event table filtered to the window", covers = "hg_snapshots()",
                     quantity = sprintf("snapshot at %g, window 2: memberships", at), tolerance = 0, key = c("node", "hyperedge"))
  }, snapshots[as.character(shown)], shown))
})

test_that("growth counts per window and cumulatively equal the filtered event table", {
  active <- do.call(rbind, lapply(grid, \(at) {
    members <- in_window(at, at + 2)
    so_far <- subset(events, time < at + 2)
    data.frame(time = at, as.list(counts(members)),
               n_nodes_cumulative = length(unique(so_far$node)), n_edges_cumulative = length(unique(so_far$hyperedge)))
  }))
  expect_reference(hg_growth(temporal, start = 0, end = 10, step = 2, window = 2), active,
                   reference = "independent: the event table filtered to each window", covers = "hg_growth()",
                   quantity = "window 2 growth", tolerance = 0, key = "time")
  # the window defaults to the step (documented), and the cumulative view holds what began before t + window
  cumulative <- do.call(rbind, lapply(grid, \(at) data.frame(time = at, as.list(counts(subset(events, time < at + 2))))))
  expect_reference(hg_growth(temporal, start = 0, end = 10, step = 2, mode = "cumulative"), cumulative,
                   reference = "independent: the event table up to the end of each window", covers = "hg_growth()",
                   quantity = "cumulative growth", tolerance = 0, key = "time")
})

test_that("a subset by size keeps the hyperedges of that size", {
  hg <- group_hypergraph(events, node = "node", hyperedge = "hyperedge")
  sizes <- table(events$hyperedge)
  kept <- names(sizes)[sizes == 3L]
  expected <- subset(events, hyperedge %in% kept, select = c(node, hyperedge))
  expect_reference(transform(hg_get(hg_subset(hg, size = 3), what = "memberships"), member = TRUE),
                   transform(expected, member = TRUE), reference = "independent: the membership table filtered by size",
                   covers = c("hg_subset()", "hg_get.net_hg()"), quantity = "hyperedges of size 3", tolerance = 0,
                   key = c("node", "hyperedge"))
})
