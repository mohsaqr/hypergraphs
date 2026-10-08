# Dual hypergraph: transpose identity, involution, and equality with the
# opposite-orientation construction.

dual_corpus <- c(a = "salt and soup", b = "soup and stars",
                 c = "stars and salt soup")

test_that("the dual incidence is exactly the transpose", {
  hg <- text_hypergraph(dual_corpus, weight = "tfidf")
  dual <- dual_hypergraph(hg)
  expect_identical(dual$incidence, t(hg$incidence))
  expect_identical(dual$n_nodes, hg$n_hyperedges)
  expect_identical(dual$n_hyperedges, hg$n_nodes)
})

test_that("the dual of the dual restores the original incidence", {
  hg <- text_hypergraph(dual_corpus, weight = "tfidf")
  dual <- dual_hypergraph(hg)
  back <- dual_hypergraph(dual)
  expect_identical(back$incidence, hg$incidence)
})

test_that("the dual of a bag hypergraph equals the opposite orientation", {
  doc_hg <- text_hypergraph(dual_corpus, nodes = "doc", weight = "tfidf")
  word_hg <- text_hypergraph(dual_corpus, nodes = "word", weight = "tfidf")
  dual <- dual_hypergraph(doc_hg)
  expect_identical(dual$incidence, word_hg$incidence)
  expect_identical(dual$text$nodes, "word")
  expect_s3_class(dual, "text_hypergraph")
  # and every verb sees the same hypergraph
  dual_measures <- hg_measures(dual)
  word_measures <- hg_measures(word_hg)
  expect_identical(dual_measures, word_measures)
})

test_that("duals of non-bag constructions drop the text layer", {
  win <- text_hypergraph(c(d = "a b c a b"), construction = "window",
                         window = 2)
  dual <- dual_hypergraph(win)
  expect_false(inherits(dual, "text_hypergraph"))
  expect_identical(dual$incidence, t(win$incidence))
})

test_that("dual rejects non-hypergraphs with a classed error", {
  expect_error(dual_hypergraph("x"), class = "hypergraphs_bad_input")
})

# A general (dgC) sparse copy of a dense incidence; Matrix::Matrix() would
# pick symmetric storage for a square symmetric matrix.
.dual_sparse <- function(m) {
  methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix")
}

# Isolated vertex z, empty hyperedge E, plus the degenerate shapes.
.dual_fixtures <- function() {
  full <- matrix(c(1, 1, 1, 0,
                   1, 0, 1, 0,
                   0, 0, 0, 0), nrow = 4,
                 dimnames = list(c("a", "b", "c", "z"), c("X", "Y", "E")))
  list(
    isolated_and_empty = full,
    zero_edges = full[, 0L, drop = FALSE],
    zero_nodes = full[0L, , drop = FALSE]
  )
}

test_that("REGRESSION R08: the dual keeps isolated vertices and empty edges", {
  for (nm in names(.dual_fixtures())) {
    m <- .dual_fixtures()[[nm]]
    for (sparse in c(FALSE, TRUE)) {
      inc <- if (sparse) .dual_sparse(m) else m
      hg <- hypergraphs:::.thg_from_incidence(inc, params = list())
      dual <- dual_hypergraph(hg)
      expected <- if (sparse) Matrix::t(inc) else t(inc)
      expect_identical(dual$incidence, expected, info = nm)
      expect_identical(dual$n_nodes, ncol(m), info = nm)
      expect_identical(dual$n_hyperedges, nrow(m), info = nm)
      # a dimension of length zero carries no names in base R
      expect_identical(dual$nodes, colnames(m) %||% character(0), info = nm)
      # double dual restores the incidence exactly
      expect_identical(dual_hypergraph(dual)$incidence, inc, info = nm)
      # readers run on every shape
      expect_identical(nrow(hg_get(dual, what = "nodes")), ncol(m), info = nm)
      expect_identical(nrow(hg_get(dual)), nrow(m), info = nm)
    }
  }
  # the isolated vertex z becomes an empty hyperedge of the dual
  hg <- hypergraphs:::.thg_from_incidence(
    .dual_fixtures()$isolated_and_empty, params = list()
  )
  edges <- hg_get(dual_hypergraph(hg))
  expect_identical(edges$size, c(2L, 1L, 2L, 0L))
  expect_identical(edges$hyperedge, c("a", "b", "c", "z"))
})

test_that("REGRESSION R07: a sparse dual satisfies the full net_hg contract", {
  m <- .dual_fixtures()$isolated_and_empty
  dense <- dual_hypergraph(hypergraphs:::.thg_from_incidence(m, list()))
  sparse <- dual_hypergraph(
    hypergraphs:::.thg_from_incidence(.dual_sparse(m), list())
  )
  expect_identical(names(sparse), names(dense))
  expect_identical(sparse$hyperedges, dense$hyperedges)
  expect_identical(sparse$size_distribution, dense$size_distribution)
  expect_identical(hg_get(sparse), hg_get(dense))
  expect_identical(hg_get(sparse, what = "nodes"), hg_get(dense, what = "nodes"))
  for (what in c("nodes", "edges", "overlap", "summary")) {
    expect_equal(hg_measures(sparse, what = what),
                 hg_measures(dense, what = what), info = what)
  }
  # a group hypergraph built sparse, as in the audit
  d <- data.frame(node = c("a", "b", "c", "a", "c"),
                  edge = c("X", "X", "X", "Y", "Y"))
  hs <- group_hypergraph(d, node = "node", hyperedge = "edge", sparse = TRUE)
  hd <- group_hypergraph(d, node = "node", hyperedge = "edge")
  expect_identical(hg_get(dual_hypergraph(hs)), hg_get(dual_hypergraph(hd)))
})

test_that("REGRESSION R07: a sparse bag dual keeps the text layer and flips it", {
  for (orientation in c("doc", "word")) {
    dense <- text_hypergraph(dual_corpus, nodes = orientation, weight = "tfidf")
    sparse <- text_hypergraph(dual_corpus, nodes = orientation,
                              weight = "tfidf", sparse = TRUE)
    dual_sparse <- dual_hypergraph(sparse)
    dual_dense <- dual_hypergraph(dense)
    expect_s3_class(dual_sparse, "text_hypergraph")
    expect_s4_class(dual_sparse$incidence, "sparseMatrix")
    expect_identical(dual_sparse$text$nodes, dual_dense$text$nodes)
    expect_false(identical(dual_sparse$text$nodes, orientation))
    expect_identical(dual_sparse$text$documents, sparse$text$documents)
    expect_identical(dual_sparse$text$vocabulary, sparse$text$vocabulary)
    expect_identical(as.matrix(dual_sparse$incidence), dual_dense$incidence)
    expect_identical(dual_sparse$hyperedges, dual_dense$hyperedges)
  }
})
