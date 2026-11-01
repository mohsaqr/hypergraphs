# hg_edge_centrality(): s-betweenness and s-closeness of hyperedges on the
# s-line graph, checked against values computed by hand from the NetworkX
# conventions the paper uses (Coupette et al. 2024).

.ec_toy <- function(sparse = FALSE) {
  d <- data.frame(node = c("a", "b", "c", "b", "c", "d", "d", "e",
                           "a", "c", "e"),
                  edge = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3",
                           "e4", "e4", "e4"))
  group_hypergraph(d, node = "node", hyperedge = "edge", sparse = sparse)
}

.ec_value <- function(out, s, measure) {
  hit <- out[out$s == s & out$measure == measure, , drop = FALSE]
  stats::setNames(hit$value, hit$edge)
}

test_that("s = 1 and s = 2 match the hand-computed NetworkX values", {
  out <- hg_edge_centrality(.ec_toy(), s = 1:2,
                            measure = c("betweenness", "closeness"))
  expect_named(out, c("edge", "s", "measure", "value"))
  # s = 1 line graph: e1-e2, e1-e4, e2-e3, e2-e4, e3-e4 (connected)
  between_1 <- .ec_value(out, 1L, "betweenness")
  expect_equal(between_1[c("e1", "e2", "e3", "e4")], c(e1 = 0, e2 = 1 / 6,
                                                       e3 = 0, e4 = 1 / 6))
  close_1 <- .ec_value(out, 1L, "closeness")
  expect_equal(close_1[c("e1", "e2", "e3", "e4")],
               c(e1 = 3 / 4, e2 = 1, e3 = 3 / 4, e4 = 1))
  # s = 2 line graph: e2 - e1 - e4, e3 isolated; Wasserman-Faust closeness
  between_2 <- .ec_value(out, 2L, "betweenness")
  expect_equal(between_2[["e1"]], 1 / 3)
  close_2 <- .ec_value(out, 2L, "closeness")
  expect_equal(close_2[c("e1", "e2", "e3", "e4")],
               c(e1 = 2 / 3, e2 = 4 / 9, e3 = 0, e4 = 4 / 9))
})

test_that("dense and sparse incidences agree; rows are ranked per block", {
  dense <- hg_edge_centrality(.ec_toy(), s = 1:2,
                              measure = c("betweenness", "closeness"))
  sparse <- hg_edge_centrality(.ec_toy(sparse = TRUE), s = 1:2,
                               measure = c("betweenness", "closeness"))
  expect_identical(sparse, dense)
  top <- hg_edge_centrality(.ec_toy(), s = 1:2, top = 1)
  # one row per (s, measure) block: betweenness then closeness, s = 1 then 2
  expect_identical(nrow(top), 4L)
  expect_identical(top$edge, c("e2", "e2", "e1", "e1"))
})

test_that("broken arguments are classed", {
  hg <- .ec_toy()
  expect_error(hg_edge_centrality(hg, s = 0), class = "hypergraphs_bad_input")
  expect_error(hg_edge_centrality(hg, normalized = NA),
               class = "hypergraphs_bad_input")
  expect_error(hg_edge_centrality(hg, top = 0), class = "hypergraphs_bad_input")
  expect_error(hg_edge_centrality(list()), class = "hypergraphs_bad_input")
})
