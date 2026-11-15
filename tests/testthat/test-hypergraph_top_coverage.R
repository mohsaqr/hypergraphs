# Hand-counted transactions: ab occurs five times, ac three, and bc twice.
.top_coverage_data <- function(counts = c(5L, 3L, 2L),
                               sets = list(c("a", "b"), c("a", "c"),
                                           c("b", "c"))) {
  transactions <- rep(seq_along(sets), counts)
  do.call(rbind, lapply(seq_along(transactions), function(i) {
    data.frame(item = sets[[transactions[i]]], basket = i)
  }))
}

test_that("fractional top retains the smallest prefix reaching coverage", {
  d <- .top_coverage_data()
  read_sets <- function(top) hg_get(group_hypergraph(
    d, node = "item", hyperedge = "basket", top = top), what = "sets")
  expect_identical(read_sets(0.5)$count, 5L)
  expect_identical(read_sets(0.6)$count, c(5L, 3L))
  expect_identical(read_sets(0.8)$count, c(5L, 3L))
  expect_identical(read_sets(0.95)$count, c(5L, 3L, 2L))
  expect_identical(read_sets(1)$count, 5L)
  expect_identical(read_sets(2)$count, c(5L, 3L))
  expect_identical(read_sets(Inf)$count, c(5L, 3L, 2L))
  expect_equal(read_sets(0.6)$share, c(0.5, 0.3))
  expect_identical(hypergraph(d, action = "item", session = "basket",
                              top = 0.8),
                   group_hypergraph(d, node = "item", hyperedge = "basket",
                                    top = 0.8))
})

test_that("coverage is calculated separately within each group", {
  first <- transform(.top_coverage_data(), cohort = "first")
  second <- transform(.top_coverage_data(c(8L, 1L, 1L)), cohort = "second",
                      basket = basket + 10L)
  sets <- hg_get(group_hypergraph(rbind(first, second), node = "item",
                                  hyperedge = "basket", group = "cohort",
                                  top = 0.8), what = "sets")
  expect_identical(sets$group, c("first", "first", "second"))
  expect_identical(sets$count, c(5L, 3L, 8L))
  expect_equal(sets$share, c(0.5, 0.3, 0.8))
})

test_that("ties use set names and decimal boundaries do not add a set", {
  tied <- .top_coverage_data(c(4L, 3L, 3L))
  sets <- hg_get(group_hypergraph(tied, node = "item", hyperedge = "basket",
                                  top = 0.7), what = "sets")
  expect_identical(sets$set, c("a + b", "a + c"))
  shuffled <- tied[rev(seq_len(nrow(tied))), ]
  expect_identical(hg_get(group_hypergraph(shuffled, node = "item",
                                          hyperedge = "basket", top = 0.7),
                          what = "sets"), sets)
  boundary <- .top_coverage_data(c(28L, 27L, 25L, 20L),
    list(c("a", "b"), c("a", "c"), c("b", "c"), c("a", "d")))
  exact <- group_hypergraph(boundary, node = "item", hyperedge = "basket",
                            top = 0.28)
  above <- group_hypergraph(boundary, node = "item", hyperedge = "basket",
                            top = 0.2801)
  expect_identical(hg_get(exact, what = "sets")$count, 28L)
  expect_identical(hg_get(above, what = "sets")$count, c(28L, 27L))
})

test_that("coverage uses eligible transactions and shares keep the full denominator", {
  d <- .top_coverage_data(c(7L, 2L, 1L),
                         list("a", c("a", "b"), c("b", "c")))
  sizes <- hg_get(group_hypergraph(d, node = "item", hyperedge = "basket",
                                   top = 0.6, min_size = 2), what = "sets")
  expect_identical(sizes$set, "a + b")
  expect_equal(sizes$share, 0.2)
  support <- hg_get(group_hypergraph(.top_coverage_data(), node = "item",
                                     hyperedge = "basket", top = 0.6,
                                     min_share = 0.25), what = "sets")
  expect_identical(support$count, 5L)
  states <- hg_get(group_hypergraph(.top_coverage_data(), node = "item",
                                    hyperedge = "basket", top = 0.8,
                                    states = c("a", "b")), what = "sets")
  expect_identical(states$set, c("a + b", "a"))
  expect_identical(states$count, c(5L, 3L))
})

test_that("fractional top works for clustered sequences and topic combinations", {
  rows <- rep(1:3, c(5L, 3L, 2L))
  sequences <- data.frame(V1 = c("a", "a", "b")[rows],
                          V2 = c("b", "c", "c")[rows])
  clustered <- structure(list(data = sequences, assignments = rep(1L, 10L),
                              k = 1L), class = "net_mmm")
  hg <- group_hypergraph(clustered, top = 0.8)
  expect_identical(hg_get(hg, what = "sets")$count, c(5L, 3L))
  expect_identical(hg$params$top, 0.8)
  expect_output(print(hg), "cumulative 80% coverage")
  expect_no_warning(group_hypergraph(clustered, top = Inf))
  d <- .top_coverage_data()
  topics <- structure(list(shares = data.frame(
    node = as.character(d$basket), topic = paste("Topic", match(d$item, letters)),
    share = 0.5)), class = "net_hg_topics")
  topic_hg <- group_hypergraph(topics, threshold = 0.2, top = 0.8)
  expect_identical(hg_get(topic_hg, what = "sets")$count, c(5L, 3L))
  expect_output(print(topic_hg), "cumulative 80% coverage")
})

test_that("invalid top values raise classed errors", {
  d <- .top_coverage_data()
  for (top in list(0, -0.5, -Inf, NA_real_, NaN, 1.5, 2 - 1e-10,
                   0.95 + 0i, c(0.5, 0.9), numeric(), "0.95", TRUE)) {
    expect_error(group_hypergraph(d, node = "item", hyperedge = "basket",
                                  top = top), class = "hypergraphs_bad_input")
  }
})
