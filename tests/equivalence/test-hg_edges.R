# hg_edges(): one row per hyperedge (size, weight, number of hyperedges it meets, number of outside nodes adjacent
# to its members), the distribution of a measure and its summary. References: igraph on the bipartite incidence
# graph and the incidence matrix in base R, on the small seven-edge hypergraph and on the ICSID tribunals.

memberships <- ref_memberships()
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")

#' The four hyperedge measures from an incidence matrix
edge_measures <- function(incidence) {
  together <- tcrossprod(incidence) > 0
  diag(together) <- FALSE
  data.frame(edge = colnames(incidence), size = colSums(incidence), weight = colSums(incidence),
             n_incident_edges = rowSums(crossprod(incidence) > 0) - 1,
             n_neighbors = colSums(together %*% incidence > 0 & incidence == 0))
}

test_that("hyperedge size, weight, incident hyperedges and outside neighbours equal the incidence matrix", {
  expect_reference(hg_edges(hg, what = "edges"), edge_measures(ref_incidence(memberships)),
                   reference = "independent: the incidence matrix, base R", covers = "hg_edges()",
                   quantity = "hyperedge measures", tolerance = 0, key = "edge")
  tribunals <- unique(subset(icsid_tribunals, select = c(arbitrator, case)))
  expect_reference(hg_edges(group_hypergraph(tribunals, node = "arbitrator", hyperedge = "case"), what = "edges"),
                   edge_measures(ref_incidence(tribunals, node = "arbitrator", edge = "case")),
                   reference = "independent: the incidence matrix, base R", covers = c("hg_edges()", "group_hypergraph()"),
                   quantity = "ICSID tribunals, hyperedge measures", tolerance = 0, key = "edge")
})

test_that("hyperedge size, weight, incident hyperedges and outside neighbours equal igraph on the incidence graph", {
  skip_without_reference("igraph")
  bipartite <- igraph::graph_from_data_frame(memberships, directed = FALSE)
  igraph::E(bipartite)$weight <- 1
  edges <- unique(memberships$edge)
  expected <- data.frame(edge = edges, size = unname(igraph::degree(bipartite, edges)),
                         weight = unname(igraph::strength(bipartite, edges)),
                         n_incident_edges = igraph::ego_size(bipartite, order = 2, nodes = edges, mindist = 2),
                         n_neighbors = igraph::ego_size(bipartite, order = 3, nodes = edges, mindist = 3))
  expect_reference(hg_edges(hg, what = "edges"), expected, reference = "igraph", covers = "hg_edges()",
                   quantity = "hyperedge measures", tolerance = 0, key = "edge")
})

test_that("the distribution and summary of hyperedge sizes equal their definitions", {
  tribunals <- unique(subset(icsid_tribunals, select = c(arbitrator, case)))
  sizes <- colSums(ref_incidence(tribunals, node = "arbitrator", edge = "case"))
  counts <- table(sizes)
  values <- as.numeric(names(counts))
  ours <- group_hypergraph(tribunals, node = "arbitrator", hyperedge = "case")
  expect_reference(hg_edges(ours, what = "distribution", measure = "size"),
                   data.frame(value = values, n = as.vector(counts), proportion = as.vector(counts) / length(sizes),
                              ccdf = vapply(values, \(v) mean(sizes >= v), numeric(1))),
                   reference = "independent: the hyperedge sizes, base R", covers = "hg_edges()",
                   quantity = "ICSID tribunals, size distribution", tolerance = 1e-12, key = "value")
  quartiles <- stats::quantile(sizes, c(0.25, 0.5, 0.75), names = FALSE)
  # the summary is one row: compared as a named vector, one number per column
  expect_reference(unlist(hg_edges(ours, what = "summary", measure = "size")),
                   c(n_edges = length(sizes), mean = mean(sizes), sd = stats::sd(sizes), min = min(sizes),
                     q25 = quartiles[[1L]], median = quartiles[[2L]], q75 = quartiles[[3L]], max = max(sizes)),
                   reference = "stats", covers = "hg_edges()", quantity = "ICSID tribunals, size summary",
                   tolerance = 1e-12, note = "stats::sd() and stats::quantile(type = 7).")
})
