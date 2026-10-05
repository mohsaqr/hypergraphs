skip_on_cran()

# The hull plot's node-label size: 4.2 mm up to a dozen labels, then
# 4.2 * sqrt(12 / n) to a floor of 2.2 mm.

.label_text_size <- function(p) {
  texts <- Filter(\(layer) inherits(layer$geom, "GeomText"), p$layers)
  utils::tail(texts, 1L)[[1L]]$aes_params$size
}

.argentina <- function() {
  h <- group_hypergraph(icsid_tribunals, node = "arbitrator",
                        hyperedge = "case")
  hg_subset(h, where = list(respondent = "Argentine Republic"),
            component = "largest")
}

test_that(".thg_label_size matches its formula at hand-computed points", {
  expect_identical(.thg_label_size(1), 4.2)
  expect_identical(.thg_label_size(12), 4.2)
  expect_equal(.thg_label_size(20), 4.2 * sqrt(12 / 20))
  expect_equal(.thg_label_size(37), 2.391878, tolerance = 1e-6)
  expect_identical(.thg_label_size(44), 2.2)
  expect_identical(.thg_label_size(500), 2.2)
  expect_identical(.thg_label_size(0), 4.2)
})

test_that("INVARIANT: label size never grows with the number of labels", {
  sizes <- vapply(seq_len(300), .thg_label_size, numeric(1))
  expect_true(all(diff(sizes) <= 0))
  expect_true(all(sizes >= 2.2 & sizes <= 4.2))
})

test_that("the hull plot scales the default to the labelled nodes", {
  small <- group_hypergraph(
    data.frame(node = c("a", "b", "c", "b", "c", "d"),
               hyperedge = rep(c("e1", "e2"), each = 3)),
    node = "node", hyperedge = "hyperedge")
  expect_identical(.label_text_size(plot(small)), 4.2)
  argentina <- .argentina()
  # 44 labels or more reach the floor
  expect_gte(argentina$n_nodes, 44L)
  expect_identical(.label_text_size(plot(argentina)), 2.2)
  # only the labels written count: ten names and the rest blank size as ten
  ten <- stats::setNames(
    c(head(argentina$nodes, 10L), rep("", argentina$n_nodes - 10L)),
    argentina$nodes
  )
  expect_identical(.label_text_size(plot(argentina, labels = ten)), 4.2)
  # an explicit size is used as given
  expect_identical(.label_text_size(plot(argentina, label_size = 3.5)), 3.5)
})

test_that("an invalid label_size raises hypergraphs_bad_input", {
  hg <- .argentina()
  expect_error(plot(hg, label_size = -1), class = "hypergraphs_bad_input")
  expect_error(plot(hg, label_size = c(2, 3)), class = "hypergraphs_bad_input")
  expect_error(plot(hg, label_size = "big"), class = "hypergraphs_bad_input")
})
