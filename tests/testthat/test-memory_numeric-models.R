# Numeric state labels are values in Nestimate models, not alphabet positions.
testthat::skip_on_cran()

.nm_sequences <- function() {
  as.data.frame(rbind(c(1, 3, 5, 1, 3, 5), c(5, 3, 1, 5, 3, 1),
                     c(1, 5, 3, 1, 5, 3), c(3, 5, 1, 3, 5, 1)))
}

.nm_tna <- function(data, labels = c("5", "1", "3")) {
  codes <- vapply(data, function(col) match(as.character(col), labels),
                   integer(nrow(data)))
  attr(codes, "labels") <- labels
  structure(list(data = codes, labels = labels), class = "tna")
}

test_that("numeric model states reach all memory constructors unchanged", {
  raw <- .nm_sequences()
  net <- Nestimate::build_network(raw, method = "relative")
  expect_identical(.coerce_sequence_input(net), raw)
  expect_identical(.hon_sequence_values(net), lapply(seq_len(nrow(raw)),
                   function(i) as.character(unlist(raw[i, ], use.names = FALSE))))

  raw_hypa <- hypa(raw, order = 2L, min_count = 1L)
  model_hypa <- hypa(net, order = 2L, min_count = 1L)
  expect_gt(nrow(hg_get(raw_hypa, what = "scores")), 0L)
  expect_equal(hg_get(model_hypa, what = "scores"),
                 hg_get(raw_hypa, what = "scores"))
  expect_equal(hg_get(hon(net, max_order = 2L), what = "nodes"),
                 hg_get(hon(raw, max_order = 2L), what = "nodes"))
  expect_equal(hg_get(mogen(net, max_order = 2L), what = "paths", k = 2L),
                 hg_get(mogen(raw, max_order = 2L), what = "paths", k = 2L))
  expect_identical(.hi_parse(net), .hi_parse(raw))
  expect_equal(hg_get(hg_bootstrap(net, n_boot = 5L, max_order = 2L, seed = 8L)),
                 hg_get(hg_bootstrap(raw, n_boot = 5L, max_order = 2L, seed = 8L)))
})

test_that("tna encoding is decoded by provenance even for numeric labels", {
  raw <- .nm_sequences()
  literal <- as.data.frame(lapply(raw, as.character))
  coded <- .nm_tna(raw)
  expect_equal(.coerce_sequence_input(coded), literal)
  expect_equal(hg_get(hypa(coded, order = 2L, min_count = 1L), what = "scores"),
                 hg_get(hypa(raw, order = 2L, min_count = 1L), what = "scores"))

  # Every code is itself a label here; membership cannot identify encoding.
  digits <- as.data.frame(matrix(as.character(1:12), 2L, 12L, byrow = TRUE))
  coded_digits <- .nm_tna(digits, labels = sort(as.character(1:12)))
  expect_equal(.coerce_sequence_input(coded_digits), digits)
  expect_equal(hg_get(mogen(coded_digits, max_order = 2L), what = "paths", k = 2L),
                 hg_get(mogen(digits, max_order = 2L), what = "paths", k = 2L))
})

test_that("generic cograph sequence values are literal unless marked encoded", {
  raw <- .nm_sequences()
  net <- Nestimate::build_network(raw, method = "relative")
  plain <- net
  class(plain) <- "cograph_network"
  expect_identical(.coerce_sequence_input(plain), raw)
  expect_equal(hg_get(hypa(plain, order = 2L, min_count = 1L), what = "scores"),
                 hg_get(hypa(raw, order = 2L, min_count = 1L), what = "scores"))

  encoded <- plain
  labels <- rownames(net$weights)
  encoded$data <- .nm_tna(raw, labels)$data
  encoded$meta$source <- "tna"
  expect_equal(.coerce_sequence_input(encoded), as.data.frame(lapply(raw, as.character)))
  expect_equal(hg_get(hypa(encoded, order = 2L, min_count = 1L), what = "scores"),
                 hg_get(hypa(raw, order = 2L, min_count = 1L), what = "scores"))
})

test_that("numeric mixture and grouped model states remain literal", {
  raw <- .nm_sequences()
  # Match the fields stored by build_mmm(), with deliberately nonsequential
  # numeric states; component dimnames are not an encoding declaration.
  net <- Nestimate::build_network(raw, method = "relative")
  fit <- structure(list(data = raw, assignments = c(1L, 1L, 2L, 2L), k = 2L,
                         models = list(net, net)), class = "net_mmm")
  grouped <- .coerce_grouped_sequences(fit)
  expected <- lapply(seq_len(nrow(raw)),
                      function(i) as.character(unlist(raw[i, ], use.names = FALSE)))
  expect_identical(grouped, list(`Cluster 1` = expected[1:2],
                                  `Cluster 2` = expected[3:4]))
  expect_setequal(group_hypergraph(fit)$nodes, c("1", "3", "5"))

  # Real fit: stored numeric data obey the same representation contract.
  actual <- suppressWarnings(Nestimate::build_mmm(raw[rep(1:4, 4L), ],
                              k = 2L, n_starts = 1L, max_iter = 20L, seed = 4L))
  sequences <- .coerce_grouped_sequences(actual)
  expect_setequal(unique(unlist(sequences, use.names = FALSE)), c("1", "3", "5"))
  expect_true(all(lengths(unlist(sequences, recursive = FALSE)) == ncol(raw)))
})
