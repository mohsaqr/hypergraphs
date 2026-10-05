testthat::skip_on_cran()

# e1 = {a, b, c}, e2 = {b, c, d}, e3 = {d, e}; sources s1, s1, s2.
.subset_fixture <- function(sparse = FALSE) {
  dat <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "d", "e"),
    event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3")
  )
  hg <- group_hypergraph(dat, "member", "event", sparse = sparse)
  hg$edge_data <- data.frame(edge = c("e1", "e2", "e3"),
                             source = c("s1", "s1", "s2"),
                             stringsAsFactors = FALSE)
  hg$edge_multiplicity <- c(1L, 2L, 1L)
  hg
}

test_that("hg_subset keeps named hyperedges and their nodes", {
  out <- hg_subset(.subset_fixture(), edges = c("e1", "e2"))
  expect_s3_class(out, "net_hg")
  expect_identical(colnames(out$incidence), c("e1", "e2"))
  expect_identical(out$nodes, c("a", "b", "c", "d"))
  expect_identical(out$n_hyperedges, 2L)
  expect_identical(out$n_nodes, 4L)
  expect_identical(out$hyperedges, list(1:3, 2:4))
  expect_identical(out$size_distribution, c(size_3 = 2L))
  expect_identical(out$edge_data$edge, c("e1", "e2"))
  expect_identical(out$edge_multiplicity, c(1L, 2L))
})

test_that("hg_subset by nodes keeps the induced sub-hypergraph", {
  out <- hg_subset(.subset_fixture(), nodes = c("b", "c", "d", "e"))
  expect_identical(colnames(out$incidence), c("e2", "e3"))
  expect_identical(out$nodes, c("b", "c", "d", "e"))
  # a node set with no hyperedge inside gives an edgeless hypergraph
  empty <- hg_subset(.subset_fixture(), nodes = c("a", "e"))
  expect_identical(empty$n_hyperedges, 0L)
  expect_identical(empty$n_nodes, 2L)
})

test_that("hg_subset by hyperedge attribute uses the edge metadata", {
  out <- hg_subset(.subset_fixture(), where = c(source = "s2"))
  expect_identical(colnames(out$incidence), "e3")
  expect_identical(out$nodes, c("d", "e"))
  kept <- hg_subset(.subset_fixture(), where = list(source = "s2"),
                    drop_isolated = FALSE)
  expect_identical(kept$nodes, c("a", "b", "c", "d", "e"))
  both <- hg_subset(.subset_fixture(), where = list(source = c("s1", "s2")))
  expect_identical(both$n_hyperedges, 3L)
  no_source <- group_hypergraph(
    data.frame(member = c("a", "b"), event = "e1"), "member", "event"
  )
  expect_error(hg_subset(no_source, where = c(source = "s1")), class = "hypergraphs_bad_input")
  expect_error(hg_subset(.subset_fixture(), where = "s1"), class = "hypergraphs_bad_input")
})

test_that("a where value that no hyperedge takes is an error, not an empty result", {
  expect_error(hg_subset(.subset_fixture(), where = c(source = "s3")),
               class = "hypergraphs_bad_input")
  # one absent value among present ones is still an error
  expect_error(hg_subset(.subset_fixture(), where = list(source = c("s1", "s9"))),
               class = "hypergraphs_bad_input")
  # the message offers the close values of the attribute
  tribunals <- group_hypergraph(icsid_tribunals, node = "arbitrator",
                                hyperedge = "case")
  err <- tryCatch(hg_subset(tribunals, where = list(respondent = "Argentina")),
                  hypergraphs_bad_input = identity)
  expect_match(conditionMessage(err), "Argentine Republic", fixed = TRUE)
  # filters that are each satisfiable may still leave nothing: that is a
  # result, not an error
  none <- hg_subset(.subset_fixture(), edges = "e1", where = c(source = "s2"))
  expect_identical(none$n_hyperedges, 0L)
})

test_that("INVARIANT: sparse and dense subsets agree", {
  dense <- hg_subset(.subset_fixture(), edges = c("e2", "e3"))
  sparse <- hg_subset(.subset_fixture(sparse = TRUE), edges = c("e2", "e3"))
  expect_true(methods::is(sparse$incidence, "sparseMatrix"))
  expect_identical(as.matrix(sparse$incidence) * 1, dense$incidence * 1)
  expect_identical(sparse$hyperedges, dense$hyperedges)
  expect_identical(sparse$nodes, dense$nodes)
})

test_that("hg_subset takes a ranking table through its edge column", {
  hg <- .subset_fixture()
  ranking <- hg_edge_centrality(hg, s = 1, measure = "betweenness", top = 2)
  from_table <- hg_subset(hg, edges = ranking)
  from_names <- hg_subset(hg, edges = ranking$edge)
  expect_identical(from_table, from_names)
  expect_error(hg_subset(hg, edges = data.frame(x = 1)), class = "hypergraphs_bad_input")
})

test_that("hg_subset rejects unknown names and empty selectors", {
  hg <- .subset_fixture()
  expect_error(hg_subset(hg), class = "hypergraphs_bad_input")
  expect_error(hg_subset(hg, edges = "nope"), class = "hypergraphs_bad_input")
  expect_error(hg_subset(hg, nodes = "nope"), class = "hypergraphs_bad_input")
  expect_error(hg_subset(42, edges = "e1"), class = "hypergraphs_bad_input")
})

test_that("hg_subset() keeps the window counts aligned with the hyperedges", {
  sessions <- list(c("a", "a", "b", "c", "a", "b"), c("a", "b", "c", "d"))
  hg <- window_hypergraph(sessions, window = 3)
  kept <- hg_subset(hg, size = 3)
  expect_length(kept$window_counts, kept$n_hyperedges)
  full <- hg_get(hg)
  three <- subset(full, size == 3)
  expect_identical(hg_get(kept)$weight, three$weight)
  expect_identical(hg_get(kept)$members, three$members)
})

test_that("hg_subset(component = \"largest\") keeps the largest connected component", {
  dat <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "x", "y"),
    event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3")
  )
  hg <- group_hypergraph(dat, node = "member", hyperedge = "event")
  big <- hg_subset(hg, component = "largest")
  expect_setequal(big$nodes, c("a", "b", "c", "d"))
  expect_setequal(big$params$subset$removed, c("x", "y"))
  # a connected hypergraph comes back unchanged
  expect_identical(hg_subset(big, component = "largest"), big)
  # after another filter
  expect_setequal(hg_subset(hg, edges = c("e1", "e3"), component = "largest")$nodes,
                  c("a", "b", "c"))
  expect_error(hg_subset(hg, component = "biggest"),
               class = "hypergraphs_bad_input")
  expect_error(hg_subset(hg), class = "hypergraphs_bad_input")
})

test_that("labels of nodes removed by hg_subset() are set aside with a warning", {
  # "night" bridges the two themes; f shares no word with the rest
  corpus <- data.frame(
    node = c("a", "b", "c", "d", "e", "f"),
    text = c("soup salt night", "salt soup broth", "broth soup salt",
             "stars sky night", "sky stars moon", "zebra quilt"),
    label = c("food", "food", "food", "sky", "sky", "sky")
  )
  hg <- text_hypergraph(corpus, column = "text", id = "node")
  connected <- hg_subset(hg, component = "largest")
  expect_setequal(connected$nodes, c("a", "b", "c", "d", "e"))
  expect_warning(spread <- hg_classify(connected, corpus),
                 class = "hypergraphs_dropped_documents")
  expect_false("f" %in% spread$node)
  expect_error(suppressWarnings(hg_classify(connected,
                                            rbind(corpus, data.frame(node = "zz", text = "q", label = "sky")))),
               class = "hypergraphs_bad_input")
})

test_that("hg_subset() cuts the text layer of a text hypergraph", {
  corpus <- data.frame(
    node = c("a", "b", "c", "d", "e", "f"),
    text = c("soup salt night", "salt soup broth", "broth soup salt",
             "stars sky night", "sky stars moon", "zebra quilt"),
    label = c("food", "food", "food", "sky", "sky", "sky")
  )
  hg <- text_hypergraph(corpus, column = "text", id = "node")
  connected <- hg_subset(hg, component = "largest")
  weights <- hg_get(connected, what = "weights")
  expect_false(any(weights$doc == "f"))
  expect_false(any(c("zebra", "quilt") %in% weights$word))
  expect_false("f" %in% hg_get(connected, what = "documents")$doc)
  vocabulary <- hg_get(connected, what = "vocabulary")
  expect_setequal(vocabulary$word, colnames(connected$incidence))
  # recounted over the kept documents
  expect_identical(subset(vocabulary, word == "salt")$count, 3L)
  # keyword verbs that read the token layer work on the subset
  expect_warning(kw <- hg_keywords(connected, corpus, n = 2, type = "ctfidf"),
                 class = "hypergraphs_dropped_documents")
  expect_setequal(unique(kw$cluster), c("food", "sky"))
})
