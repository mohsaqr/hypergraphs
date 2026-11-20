# dual_hypergraph(): nodes and hyperedges exchange roles; a node of the dual is a hyperedge of the original and
# a hyperedge of the dual holds the hyperedges a node belongs to. References: the transposed incidence matrix in
# base R and HyperG::dual_hypergraph().

memberships <- ref_memberships()
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")

test_that("the dual's memberships equal the transposed incidence matrix", {
  incidence <- ref_incidence(memberships)
  expected <- subset(ref_long(t(incidence), name = "weight", row = "node", column = "hyperedge"), weight > 0)
  expect_reference(hg_get(dual_hypergraph(hg), what = "memberships"), expected,
                   reference = "independent: the transposed incidence matrix, base R",
                   covers = c("dual_hypergraph()", "hg_get.net_hg()"), quantity = "dual memberships", tolerance = 0,
                   key = c("node", "hyperedge"))
})

test_that("the dual equals HyperG::dual_hypergraph()", {
  skip_without_reference("HyperG")
  incidence <- ref_incidence(memberships) > 0
  # HyperG reads hyperedges in rows and vertices in columns
  reference <- HyperG::dual_hypergraph(HyperG::hypergraph_from_incidence_matrix(t(incidence * 1)))
  dual_incidence <- as.matrix(HyperG::incidence_matrix(reference))
  dimnames(dual_incidence) <- list(rownames(incidence), colnames(incidence))
  expected <- subset(ref_long(dual_incidence > 0, name = "member", row = "hyperedge", column = "node"), member)
  ours <- transform(hg_get(dual_hypergraph(hg), what = "memberships"), member = TRUE)
  expect_reference(ours, expected, reference = "HyperG", covers = c("dual_hypergraph()", "hg_get.net_hg()"),
                   quantity = "dual, binary memberships", tolerance = 0, key = c("node", "hyperedge"))
})
