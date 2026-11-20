# pairwise_network(), hg_get.net_hg_pairwise(): the clique (two-section) expansion of a hypergraph. Two nodes are
# joined when they share a hyperedge, weighted by the number of hyperedges they share; node degree and strength
# are read on that network. References: igraph's bipartite projection of the incidence graph with
# multiplicities, and the off-diagonal of H H' in base R.

memberships <- ref_memberships()
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
network <- pairwise_network(hg, type = "clique")

test_that("the clique expansion's edge weights equal the off-diagonal of H H'", {
  shared <- tcrossprod(ref_incidence(memberships))
  # Counts of shared hyperedges: integers on both sides.
  expect_reference(hg_get(network, what = "edges"), ref_pairs_upper(shared, name = "weight"),
                   reference = "independent: H H' of the incidence matrix, base R",
                   covers = c("pairwise_network()", "hg_get.net_hg_pairwise()", "hg_get()"),
                   quantity = "clique expansion weight", tolerance = 0, key = c("from", "to"))
})

test_that("the clique expansion's edge weights and node strengths equal igraph's bipartite projection", {
  skip_without_reference("igraph")
  bipartite <- igraph::graph_from_data_frame(memberships, directed = FALSE)
  igraph::V(bipartite)$type <- igraph::V(bipartite)$name %in% memberships$edge
  projection <- igraph::bipartite_projection(bipartite, multiplicity = TRUE, which = "false")
  pairs <- igraph::as_data_frame(projection, what = "edges")
  expected <- cbind(ref_ordered_pairs(pairs$from, pairs$to), weight = pairs$weight)
  expect_reference(hg_get(network, what = "edges"), expected, reference = "igraph",
                   covers = c("pairwise_network()", "hg_get.net_hg_pairwise()", "hg_get()"),
                   quantity = "clique expansion weight", tolerance = 0, key = c("from", "to"))
  nodes <- data.frame(node = igraph::V(projection)$name, degree = unname(igraph::degree(projection)),
                      strength = unname(igraph::strength(projection, weights = igraph::E(projection)$weight)))
  expect_reference(hg_get(network, what = "nodes"), nodes, reference = "igraph",
                   covers = c("pairwise_network()", "hg_get.net_hg_pairwise()"),
                   quantity = "clique expansion degree and strength", tolerance = 0, key = "node")
})

test_that("the clique expansion of the ICSID tribunals equals igraph's projection", {
  skip_without_reference("igraph")
  tribunals <- unique(subset(icsid_tribunals, select = c(arbitrator, case)))
  bipartite <- igraph::graph_from_data_frame(tribunals, directed = FALSE)
  igraph::V(bipartite)$type <- igraph::V(bipartite)$name %in% tribunals$case
  projection <- igraph::bipartite_projection(bipartite, multiplicity = TRUE, which = "false")
  pairs <- igraph::as_data_frame(projection, what = "edges")
  expected <- cbind(ref_ordered_pairs(pairs$from, pairs$to), weight = pairs$weight)
  ours <- pairwise_network(group_hypergraph(tribunals, node = "arbitrator", hyperedge = "case"), type = "clique")
  observed <- hg_get(ours, what = "edges")
  # the package writes a pair in the order of its node list; the reference with the smaller name first
  observed <- cbind(ref_ordered_pairs(observed$from, observed$to), weight = observed$weight)
  expect_reference(observed, expected, reference = "igraph", covers = c("pairwise_network()", "hg_get.net_hg_pairwise()"),
                   quantity = "ICSID tribunals, clique expansion weight", tolerance = 0, key = c("from", "to"))
})
