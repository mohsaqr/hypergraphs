# network_hypergraph(type = "clique"): a network's cliques as hyperedges, of sizes 2 (its edges, when
# include_pairwise = TRUE) up to max_size. Reference: igraph::cliques() on the same network, and every node
# triple tested in base R.

#' A hyperedge set as one sorted string per hyperedge
as_sets <- function(sets) data.frame(members = vapply(sets, \(set) paste(sort(set), collapse = ", "), character(1)), present = TRUE)

#' A random undirected network of 12 nodes with a planted 4-clique
random_network <- function(seed) {
  set.seed(seed)
  nodes <- sprintf("n%02d", seq_len(12L))
  adjacency <- matrix(0, 12L, 12L, dimnames = list(nodes, nodes))
  adjacency[upper.tri(adjacency)] <- stats::rbinom(66L, 1L, 0.3)
  adjacency[1:4, 1:4] <- 1
  adjacency[lower.tri(adjacency)] <- t(adjacency)[lower.tri(adjacency)]
  diag(adjacency) <- 0
  adjacency
}

test_that("hyperedges of sizes 2 to 4 equal igraph's cliques", {
  skip_without_reference("igraph")
  invisible(lapply(c(11L, 12L, 13L), \(seed) {
    adjacency <- random_network(seed)
    graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected")
    cliques <- lapply(igraph::cliques(graph, min = 2, max = 4), \(clique) igraph::V(graph)$name[clique])
    ours <- hg_get(network_hypergraph(adjacency, type = "clique", include_pairwise = TRUE, max_size = 4L), what = "edges")
    expect_reference(transform(ours, present = TRUE), as_sets(cliques), reference = "igraph",
                     covers = "network_hypergraph()", quantity = sprintf("network %d, cliques of 2 to 4 nodes", seed),
                     tolerance = 0, key = "members")
  }))
})

test_that("hyperedges of three nodes equal the triangles found by testing every triple", {
  adjacency <- random_network(21L)
  nodes <- rownames(adjacency)
  triples <- utils::combn(nodes, 3L, simplify = FALSE)
  triangles <- Filter(\(triple) all(adjacency[triple, triple][upper.tri(diag(3L))] > 0), triples)
  ours <- hg_get(network_hypergraph(adjacency, type = "clique", include_pairwise = FALSE, max_size = 3L), what = "edges")
  expect_reference(transform(ours, present = TRUE), as_sets(triangles),
                   reference = "independent: every node triple tested, base R", covers = "network_hypergraph()",
                   quantity = "triangles", tolerance = 0, key = "members")
})
