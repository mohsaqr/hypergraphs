skip_on_cran()

# plot(temporal_hypergraph, type = "storyline"): one line per node from its
# first to its last hyperedge, hyperedges as columns in order of start.

.story_data <- function() {
  data.frame(
    case = rep(c("A", "B", "C", "D"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3",
                   "p1", "a4", "a2", "p2", "a3", "a4"),
    constituted = rep(c(1, 2, 4, 5), each = 3),
    concluded = rep(c(4, 3, 6, 8), each = 3)
  )
}

.story_tg <- function(data = .story_data()) {
  temporal_hypergraph(data, node = "arbitrator", hyperedge = "case",
                      start = "constituted", end = "concluded")
}

.story_layer <- function(p, geom) {
  Filter(\(l) inherits(l$geom, geom), p$layers)
}

# the rows of the lines: one per node and column it is alive at
.story_rows <- function(p) {
  path <- .story_layer(p, "GeomPath")[[1L]]$data
  rows <- unique(transform(path, col = round(x))[c("node", "col", "y")])
  rows[order(rows$node, rows$col), , drop = FALSE]
}

test_that("columns are hyperedges in order of start and lines span first to last", {
  p <- plot(.story_tg(), type = "storyline", top = NULL)
  expect_s3_class(p, "ggplot")
  expect_identical(p$scales$get_scales("x")$labels,
                   c("A  1", "B  2", "C  4", "D  5"))
  rows <- .story_rows(p)
  # a1 sits on A and B only, so its line spans columns 1 and 2
  expect_identical(subset(rows, node == "a1")$col, c(1, 2))
  # p1 sits on A and C, so it is alive at columns 1 to 3
  expect_identical(subset(rows, node == "p1")$col, c(1, 2, 3))
  expect_setequal(unique(rows$node), c("a1", "a2", "a3", "a4", "p1", "p2"))
})

test_that("INVARIANT: the members of every hyperedge sit in adjacent rows", {
  argentina <- subset(icsid_tribunals, respondent == "Argentine Republic")
  tg <- temporal_hypergraph(argentina, node = "arbitrator", hyperedge = "case",
                            start = "constituted", end = "concluded")
  p <- plot(tg, type = "storyline", top = 12)
  hits <- .story_layer(p, "GeomPoint")[[1L]]$data
  spans <- aggregate(y ~ col, data = hits,
                     FUN = \(y) max(y) - min(y) + 1 - length(y))
  expect_true(all(spans$y == 0))
  expect_identical(length(unique(.story_rows(p)$node)), 12L)
})

test_that("the layout does not depend on the order of the input rows", {
  shuffled <- .story_data()[c(7, 2, 12, 4, 9, 1, 11, 5, 3, 10, 6, 8), ]
  expect_identical(.story_rows(plot(.story_tg(), type = "storyline")),
                   .story_rows(plot(.story_tg(shuffled), type = "storyline")))
})

test_that("the kept layout has no more crossings than the first sweep", {
  argentina <- subset(icsid_tribunals, respondent == "Argentine Republic")
  m <- unique(subset(argentina, select = c(arbitrator, case, constituted)))
  busiest <- names(sort(table(m$arbitrator), decreasing = TRUE))[1:12]
  m <- subset(m, arbitrator %in% busiest)
  edge_order <- unique(m[order(m$constituted, m$case), "case"])
  col <- match(m$case, edge_order)
  best <- .thg_storyline_layout(m$arbitrator, col, length(edge_order))
  first_only <- .thg_storyline_layout(m$arbitrator, col, length(edge_order),
                                      sweeps = 0L)
  expect_lte(attr(best, "crossings"), attr(first_only, "crossings"))
})

test_that("top keeps the nodes with most hyperedges and start/end set the period", {
  tg <- .story_tg()
  # four nodes sit on two hyperedges each; ties break by name
  p <- plot(tg, type = "storyline", top = 2)
  expect_setequal(unique(.story_rows(p)$node), c("a1", "a2"))
  windowed <- plot(tg, type = "storyline", top = NULL, start = 2, end = 4)
  expect_identical(windowed$scales$get_scales("x")$labels, c("B  2", "C  4"))
})

test_that("a calendar hypergraph labels the columns with dates", {
  seats <- transform(.story_data(),
                     constituted = as.Date("2003-01-01") + constituted * 100,
                     concluded = as.Date("2003-01-01") + concluded * 100)
  p <- plot(.story_tg(seats), type = "storyline",
            start = as.Date("2003-06-01"))
  expect_identical(p$scales$get_scales("x")$labels,
                   c("B  2003-07-20", "C  2004-02-05", "D  2004-05-15"))
})

test_that("arguments of other views and bad input raise hypergraphs_bad_input", {
  tg <- .story_tg()
  expect_error(plot(tg, type = "storyline", at = 2),
               class = "hypergraphs_bad_input")
  expect_error(plot(tg, type = "storyline", color_by = "size"),
               class = "hypergraphs_bad_input")
  expect_error(plot(tg, type = "storyline", top = 0),
               class = "hypergraphs_bad_input")
  expect_error(plot(tg, type = "storyline", start = 50),
               class = "hypergraphs_bad_input")
  expect_error(plot(tg, at = 2, top = 3), class = "hypergraphs_bad_input")
  expect_error(plot(tg, type = "storyline", start = as.Date("2003-01-01")),
               class = "hypergraphs_bad_input")
})

test_that("strength spacing scales each gap linearly with shared hyperedges", {
  tg <- .story_tg()
  even <- .story_rows(plot(tg, type = "storyline", top = NULL))
  close <- .story_rows(plot(tg, type = "storyline", top = NULL,
                            spacing = "strength"))
  # INVARIANT: the order of the lines in every column is unchanged
  rank_in_col <- function(r) ave(r$y, r$col, FUN = rank)
  expect_identical(rank_in_col(even), rank_in_col(close))
  # shared hyperedges, by hand: A = {p1, a1, a2}, B = {p2, a1, a3},
  # C = {p1, a4, a2}, D = {p2, a3, a4}
  shared <- function(u, v) {
    sets <- list(c("p1", "a1", "a2"), c("p2", "a1", "a3"),
                 c("p1", "a4", "a2"), c("p2", "a3", "a4"))
    sum(vapply(sets, \(s) all(c(u, v) %in% s), logical(1)))
  }
  # the most any pair shares here is 2 (p1-a2 and p2-a3; checked below)
  gaps <- do.call(rbind, lapply(split(close, close$col), \(r) {
    r <- r[order(r$y), ]
    s <- mapply(shared, head(r$node, -1), tail(r$node, -1))
    data.frame(gap = diff(r$y), expected = 1 - 0.65 * s / 2)
  }))
  expect_equal(gaps$gap, gaps$expected)
  # p1 and a2 share the most (two) and sit 0.35 apart; p1 and a1 share one,
  # half the strongest tie, and sit 1 - 0.65 / 2 = 0.675 apart
  first <- subset(close, col == 1)
  first <- first[order(first$y), ]
  expect_identical(first$node, c("a2", "p1", "a1"))
  expect_equal(diff(first$y), c(0.35, 0.675))
  pairs <- utils::combn(sort(unique(close$node)), 2)
  expect_identical(max(apply(pairs, 2, \(p) shared(p[1], p[2]))), 2L)
})

test_that("an invalid spacing raises an error", {
  expect_error(plot(.story_tg(), type = "storyline", spacing = "tight"))
})

test_that("INVARIANT: no two of the first 72 lines share colour and shape", {
  pairs <- paste(.thg_storyline_colours(72), .thg_storyline_shapes(72))
  expect_identical(anyDuplicated(pairs), 0L)
  expect_gt(anyDuplicated(paste(.thg_storyline_colours(73),
                                .thg_storyline_shapes(73))), 0L)
})

test_that("width_by = \"degree\" widens lines with the node's hyperedges", {
  argentina <- subset(icsid_tribunals, respondent == "Argentine Republic")
  tg <- temporal_hypergraph(argentina, node = "arbitrator", hyperedge = "case",
                            start = "constituted", end = "concluded")
  p <- plot(tg, type = "storyline", width_by = "degree")
  path_layer <- .story_layer(p, "GeomPath")[[1L]]
  path <- unique(path_layer$data[c("node", "degree")])
  # degree is the number of cases of each arbitrator, counted by hand
  by_hand <- table(unique(argentina[c("arbitrator", "case")])$arbitrator)
  expect_identical(path$degree, as.integer(by_hand[path$node]))
  # INVARIANT: built widths never decrease as the degree grows
  built <- ggplot2::ggplot_build(p)$data[[match(TRUE, vapply(
    p$layers, identical, logical(1), path_layer))]]
  widths <- unique(data.frame(group = built$group, linewidth = built$linewidth))
  degree_of_group <- path$degree[order(path$node)]
  widths <- widths[order(widths$group), ]
  ordered <- widths$linewidth[order(degree_of_group)]
  expect_true(all(diff(ordered) >= 0))
  expect_gt(max(widths$linewidth), min(widths$linewidth))
  expect_identical(p$scales$get_scales("linewidth")$name, "cases")
  # the default draws every line alike
  plain <- plot(tg, type = "storyline")
  expect_null(plain$scales$get_scales("linewidth"))
  expect_error(plot(tg, type = "storyline", width_by = "size"),
               class = "hypergraphs_bad_input")
})

test_that("window_hypergraph(collapse = FALSE) keeps every window in order", {
  steps <- list(c("plan", "monitor", "discuss", "plan", "adapt", "monitor"))
  windows <- window_hypergraph(steps, window = 3, collapse = FALSE)
  expect_identical(colnames(windows$incidence), c("1-3", "2-4", "3-5", "4-6"))
  expect_identical(windows$edge_data$start, 1:4)
  expect_identical(windows$edge_data$end, 3:6)
  expect_identical(windows$window_counts, rep(1L, 4))
  # window 1-3 holds plan, monitor and discuss once each
  expect_identical(unname(windows$incidence[c("discuss", "monitor", "plan"), "1-3"]),
                   c(1, 1, 1))
  # plan occurs at steps 1 and 4; window 2-4 holds only the second
  expect_identical(unname(windows$incidence["plan", "2-4"]), 1)
  two <- window_hypergraph(list(a = c("x", "y", "z"), b = c("z", "y", "x")),
                           window = 2, collapse = FALSE)
  expect_identical(colnames(two$incidence), c("a:1-2", "a:2-3", "b:1-2", "b:2-3"))
  expect_error(window_hypergraph(steps, window = 3, collapse = FALSE,
                                 min_weight = 2),
               class = "hypergraphs_bad_input")
  # the default still collapses identical sets
  expect_identical(window_hypergraph(steps, window = 3)$params$collapse, NULL)
})

test_that("a hypergraph of ordered windows plots as a storyline", {
  steps <- list(c("plan", "monitor", "discuss", "plan", "adapt", "monitor"))
  windows <- window_hypergraph(steps, window = 3, collapse = FALSE)
  p <- plot(windows, type = "storyline", top = NULL)
  expect_s3_class(p, "ggplot")
  # columns in stored window order, named by default
  expect_identical(p$scales$get_scales("x")$labels, c("1-3", "2-4", "3-5", "4-6"))
  expect_identical(p$scales$get_scales("colour")$name, "state")
  rows <- .story_rows(p)
  # adapt occurs only at step 5, so it is in windows 3-5 and 4-6
  expect_identical(subset(rows, node == "adapt")$col, c(3, 4))
  # sort_by reorders the columns ascending
  backwards <- plot(windows, type = "storyline", top = NULL,
                    sort_by = c("1-3" = 4, "2-4" = 3, "3-5" = 2, "4-6" = 1))
  expect_identical(backwards$scales$get_scales("x")$labels,
                   c("4-6", "3-5", "2-4", "1-3"))
  expect_identical(plot(windows, type = "storyline", edge_labels = FALSE)$scales$get_scales("x")$labels,
                   rep("", 4))
  expect_error(plot(windows, type = "storyline", alpha = 0.3),
               class = "hypergraphs_bad_input")
  expect_error(plot(windows, type = "storyline", colour_me = TRUE),
               class = "hypergraphs_bad_input")
})
