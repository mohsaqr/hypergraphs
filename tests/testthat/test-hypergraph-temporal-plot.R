skip_on_cran()

# plot() of a temporal hypergraph plots its snapshot as a hypergraph, and the
# snapshot keeps the data's names for nodes and hyperedges.

.seats <- function() {
  temporal_hypergraph(
    data.frame(
      case = rep(c("A", "B", "C"), each = 3),
      arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
      sector = rep(c("energy", "mining", "energy"), each = 3),
      constituted = rep(c(1, 2, 4), each = 3),
      concluded = rep(c(4, 3, 6), each = 3)
    ),
    node = "arbitrator", hyperedge = "case",
    start = "constituted", end = "concluded"
  )
}

test_that("a snapshot carries the data's node and hyperedge names", {
  snapshot <- hg_snapshot(.seats(), at = 2.5)
  expect_identical(snapshot$params$node, "arbitrator")
  expect_identical(snapshot$params$hyperedge, "case")
  # the incidence itself is unchanged: cases A and B are active at 2.5
  expect_identical(sort(colnames(snapshot$incidence)), c("A", "B"))
  expect_identical(sort(snapshot$nodes), c("a1", "a2", "a3", "p1", "p2"))
})

test_that("an edge-list temporal hypergraph keeps the default names", {
  contacts <- temporal_hypergraph(
    data.frame(from = c("x", "y"), to = c("y", "z"), time = c(1, 2)),
    from = "from", to = "to", time = "time"
  )
  snapshot <- hg_snapshot(contacts, at = 1)
  expect_identical(snapshot$params$node, "node")
  expect_identical(snapshot$params$hyperedge, "edge")
})

test_that("plot() of a temporal hypergraph is the hypergraph plot of its snapshot", {
  thg <- .seats()
  p <- plot(thg, at = 2.5)
  expect_s3_class(p, "ggplot")
  expect_identical(
    vapply(p$layers, \(l) class(l$geom)[1L], character(1)),
    vapply(plot(hg_snapshot(thg, at = 2.5))$layers,
           \(l) class(l$geom)[1L], character(1))
  )
  # arguments reach plot.net_hg(): the incidence view and an attribute colour
  incidence <- plot(thg, at = 2.5, type = "incidence", color_by = "sector")
  expect_identical(incidence$scales$get_scales("shape")$name, "sector")
})

test_that("the legends use the data's words", {
  # the snapshot that read "Edges" / "Edges with the node" before the fix
  cases <- subset(icsid_tribunals, respondent == "Argentine Republic")
  tribunals <- temporal_hypergraph(cases, node = "arbitrator",
                                   hyperedge = "case", start = "constituted",
                                   end = "concluded")
  p <- plot(tribunals, at = as.Date("2006-01-01"))
  titles <- unlist(lapply(p$scales$scales, \(s) s$name))
  expect_identical(titles, c("Cases", "Cases with the arbitrator"))
  expect_match(p$labels$title, "^All cases")
})

test_that("a snapshot's numeric clock column sorts the incidence columns in time order", {
  thg <- .seats()
  snapshot <- hg_snapshot(thg, at = 2.5)
  expect_true(is.numeric(snapshot$edge_data$start))
  p <- plot(thg, at = 2.5, type = "incidence", sort_by = "start")
  members <- Filter(\(l) inherits(l$geom, "GeomPoint") &&
                      "node" %in% names(l$data), p$layers)[[1L]]$data
  cols <- unique(members[c("edge", "col")])
  # A starts at 1, B at 2
  expect_identical(cols$edge[order(cols$col)], c("A", "B"))
})

test_that("method is deprecated and an empty snapshot is refused", {
  thg <- .seats()
  expect_warning(plot(thg, at = 2.5, method = "clique"),
                 class = "hypergraphs_deprecated")
  expect_error(plot(thg, at = 100), class = "hypergraphs_bad_input")
})
