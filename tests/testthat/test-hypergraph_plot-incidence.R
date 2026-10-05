skip_on_cran()

# plot(hg, type = "incidence"): rows are nodes by hyperdegree, columns are
# hyperedges, a bar joins the members of each hyperedge.

.inc_data <- function() {
  data.frame(member = c("a", "b", "c", "b", "c", "d", "d", "e", "f", "a"),
             event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3",
                       "e4", "e4"))
}

.inc_hg <- function(data = .inc_data(), ...) {
  group_hypergraph(data, node = "member", hyperedge = "event", ...)
}

# the member points: the GeomPoint layer whose data names the nodes
.inc_members <- function(p) {
  layers <- Filter(\(l) inherits(l$geom, "GeomPoint") &&
                     "node" %in% names(l$data), p$layers)
  members <- layers[[1L]]$data
  members[order(members$edge, members$node), c("edge", "node", "row", "col")]
}

.inc_geoms <- function(p, geom) {
  sum(vapply(p$layers, \(l) inherits(l$geom, geom), logical(1)))
}

.argentina <- function() {
  h <- group_hypergraph(icsid_tribunals, node = "arbitrator",
                        hyperedge = "case")
  hg_subset(h, where = list(respondent = "Argentine Republic"),
            component = "largest")
}

test_that("one point per membership, rows by hyperdegree, columns grouped", {
  hg <- .inc_hg()
  p <- plot(hg, type = "incidence")
  expect_s3_class(p, "ggplot")
  members <- .inc_members(p)
  expect_identical(nrow(members), 10L)
  # hyperdegree a, b, c, d = 2; e, f = 1; ties by name
  rows <- unique(members[order(members$row), c("node", "row")])
  expect_identical(rows$node, c("a", "b", "c", "d", "e", "f"))
  # columns by highest member row, then name: e1 and e4 hold a (row 1),
  # e2 holds b (row 2), e3 holds d (row 4)
  cols <- unique(members[order(members$col), c("edge", "col")])
  expect_identical(cols$edge, c("e1", "e4", "e2", "e3"))
  # one bar per hyperedge with two or more members
  expect_identical(.inc_geoms(p, "GeomSegment"), 1L)
  spans <- Filter(\(l) inherits(l$geom, "GeomSegment"), p$layers)[[1L]]$data
  expect_identical(nrow(spans), 4L)
})

test_that("the layout does not depend on the order of the input rows", {
  shuffled <- .inc_data()[c(10, 3, 7, 1, 9, 5, 2, 8, 4, 6), ]
  expect_identical(.inc_members(plot(.inc_hg(), type = "incidence")),
                   .inc_members(plot(.inc_hg(shuffled), type = "incidence")))
})

test_that("a sparse hypergraph plots the same incidence view", {
  expect_identical(
    .inc_members(plot(.inc_hg(), type = "incidence")),
    .inc_members(plot(.inc_hg(sparse = TRUE), type = "incidence"))
  )
})

test_that("sort_by orders the columns by an edge-metadata column", {
  argentina <- .argentina()
  p <- plot(argentina, type = "incidence", sort_by = "constituted")
  cols <- unique(.inc_members(p)[c("edge", "col")])
  cols <- cols[order(cols$col), ]
  dates <- argentina$edge_data$constituted[
    match(cols$edge, argentina$edge_data$edge)]
  expect_false(is.unsorted(dates))
  expect_identical(nrow(cols), argentina$n_hyperedges)
})

test_that("a numeric sort_by puts the largest first", {
  p <- plot(.inc_hg(), type = "incidence", sort_by = "size")
  cols <- unique(.inc_members(p)[c("edge", "col")])
  # sizes e1 = 3, e2 = 3, e3 = 2, e4 = 2; ties by highest member row
  expect_identical(cols$edge[order(cols$col)], c("e1", "e2", "e4", "e3"))
})

test_that("discrete color_by pairs every colour with a shape", {
  p <- plot(.argentina(), type = "incidence", color_by = "economic_sector")
  shape <- p$scales$get_scales("shape")
  colour <- p$scales$get_scales("colour")
  expect_false(is.null(shape))
  expect_identical(shape$name, "economic_sector")
  sectors <- sort(unique(.argentina()$edge_data$economic_sector))
  expect_identical(anyDuplicated(shape$palette(length(sectors))), 0L)
  expect_identical(unname(colour$palette(length(sectors))),
                   unname(.thg_okabe_ito[seq_along(sectors)]))
  numeric_colour <- plot(.inc_hg(), type = "incidence", color_by = "size")
  expect_null(numeric_colour$scales$get_scales("shape"))
  expect_s3_class(numeric_colour$scales$get_scales("colour"),
                  "ScaleContinuous")
})

test_that("labels and edge_labels set the axis text", {
  hg <- .inc_hg()
  plain <- plot(hg, type = "incidence")
  expect_identical(plain$scales$get_scales("x")$labels, rep("", 4L))
  named <- plot(hg, type = "incidence", edge_labels = TRUE,
                labels = c(a = "Alpha"))
  expect_identical(named$scales$get_scales("x")$labels,
                   c("e1", "e4", "e2", "e3"))
  expect_identical(named$scales$get_scales("y")$labels,
                   c("Alpha", "b", "c", "d", "e", "f"))
})

test_that("size bars appear only when hyperedge sizes differ", {
  # sizes 3, 3, 2, 2 here; every ICSID tribunal has three members
  expect_identical(.inc_geoms(plot(.inc_hg(), type = "incidence"),
                              "GeomRect"), 2L)
  expect_identical(.inc_geoms(plot(.argentina(), type = "incidence"),
                              "GeomRect"), 1L)
})

test_that("the grid of empty cells is left out of a large figure", {
  small <- plot(.inc_hg(), type = "incidence")
  expect_identical(.inc_geoms(small, "GeomPoint"), 2L)
  everything <- group_hypergraph(icsid_tribunals, node = "arbitrator",
                                 hyperedge = "case")
  expect_gt(everything$n_nodes * everything$n_hyperedges, 20000)
  large <- plot(everything, type = "incidence")
  expect_identical(.inc_geoms(large, "GeomPoint"), 1L)
})

test_that("hull-only arguments and bad input raise hypergraphs_bad_input", {
  hg <- .inc_hg()
  expect_error(plot(hg, type = "incidence", layout = "spring"),
               class = "hypergraphs_bad_input")
  expect_error(plot(hg, type = "incidence", alpha = 0.3),
               class = "hypergraphs_bad_input")
  expect_error(plot(hg, type = "incidence", dismantled = TRUE),
               class = "hypergraphs_bad_input")
  expect_error(plot(hg, sort_by = "size"), class = "hypergraphs_bad_input")
  expect_error(plot(hg, type = "incidence", sort_by = "nope"),
               class = "hypergraphs_bad_input")
  expect_error(plot(hg, type = "incidence", label_size = -1),
               class = "hypergraphs_bad_input")
  expect_error(plot(hg, type = "incidence", labels = c("x", "y")),
               class = "hypergraphs_bad_input")
  expect_error(plot(hg, type = "nope"))
})

test_that("the hull plot is unchanged by the new arguments", {
  expect_s3_class(plot(.inc_hg()), "ggplot")
  expect_identical(.inc_geoms(plot(.inc_hg()), "GeomSegment"), 0L)
})
