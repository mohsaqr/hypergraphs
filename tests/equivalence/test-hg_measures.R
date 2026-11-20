# hg_measures(), group_hypergraph(): the structural measures of a hypergraph built from a membership table.
# Node degrees on the bundled ICSID tribunals (441 arbitrators, 742 cases as hyperedges) against igraph on
# the bipartite incidence graph and against the incidence matrix in base R; the hyperedge overlaps, the size
# distribution and the summary against their definitions on the small seven-edge hypergraph.

tribunals <- unique(subset(icsid_tribunals, select = c(arbitrator, case)))
memberships <- ref_memberships()

test_that("node hyperdegree, strength, largest hyperedge and neighbours equal the incidence matrix", {
  incidence <- ref_incidence(tribunals, node = "arbitrator", edge = "case")
  sizes <- colSums(incidence)
  together <- tcrossprod(incidence) > 0
  diag(together) <- FALSE
  expected <- data.frame(
    node = rownames(incidence),
    hyperdegree = rowSums(incidence),
    strength = as.vector(incidence %*% sizes),
    max_edge_size = apply(incidence, 1L, \(row) max(sizes[row > 0])),
    n_neighbors = rowSums(together))
  hg <- group_hypergraph(tribunals, node = "arbitrator", hyperedge = "case")
  # Counts: integers on both sides.
  expect_reference(hg_measures(hg, what = "nodes"), expected, reference = "independent: the incidence matrix, base R",
                   covers = c("group_hypergraph()", "hg_measures()"), quantity = "ICSID tribunals, node measures",
                   tolerance = 0, key = "node")
})

test_that("node hyperdegree, strength, largest hyperedge and neighbours equal igraph on the incidence graph", {
  skip_without_reference("igraph")
  bipartite <- igraph::graph_from_data_frame(tribunals, directed = FALSE)
  arbitrators <- unique(tribunals$arbitrator)
  cases <- igraph::neighborhood(bipartite, order = 1, nodes = arbitrators, mindist = 1)
  expected <- data.frame(
    node = arbitrators,
    hyperdegree = unname(igraph::degree(bipartite, arbitrators)),
    strength = vapply(cases, \(own) sum(igraph::degree(bipartite, own)), numeric(1)),
    max_edge_size = vapply(cases, \(own) max(igraph::degree(bipartite, own)), numeric(1)),
    n_neighbors = igraph::ego_size(bipartite, order = 2, nodes = arbitrators, mindist = 2))
  hg <- group_hypergraph(tribunals, node = "arbitrator", hyperedge = "case")
  expect_reference(hg_measures(hg, what = "nodes"), expected, reference = "igraph",
                   covers = c("group_hypergraph()", "hg_measures()"), quantity = "ICSID tribunals, node measures",
                   tolerance = 0, key = "node")
})

test_that("pairwise overlaps of hyperedges equal the intersections of their member sets", {
  incidence <- ref_incidence(memberships)
  shared <- crossprod(incidence)
  sizes <- colSums(incidence)
  # every pair of hyperedges, those that share no node included
  overlap <- ref_pairs_upper(shared, keep = \(m) m >= 0, name = "overlap")
  smaller <- pmin(sizes[overlap$from], sizes[overlap$to])
  union <- sizes[overlap$from] + sizes[overlap$to] - overlap$overlap
  expected <- data.frame(edge_1 = overlap$from, edge_2 = overlap$to, overlap = overlap$overlap,
                         overlap_coefficient = overlap$overlap / smaller, jaccard = overlap$overlap / union)
  hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
  expect_reference(hg_measures(hg, what = "overlap"), expected, reference = "independent: the incidence matrix, base R",
                   covers = "hg_measures()", quantity = "hyperedge overlaps", tolerance = 1e-12,
                   key = c("edge_1", "edge_2"))
})

test_that("the size distribution and the summary equal their definitions", {
  incidence <- ref_incidence(memberships)
  sizes <- colSums(incidence)
  counts <- table(sizes)
  values <- as.numeric(names(counts))
  distribution <- data.frame(value = values, n = as.vector(counts), proportion = as.vector(counts) / length(sizes),
                             ccdf = vapply(values, \(v) mean(sizes >= v), numeric(1)))
  hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
  expect_reference(hg_measures(hg, what = "distribution"), distribution,
                   reference = "independent: the hyperedge sizes, base R", covers = "hg_measures()",
                   quantity = "hyperedge size distribution", tolerance = 1e-12, key = "value")
  summary_table <- data.frame(
    measure = c("n_nodes", "n_hyperedges", "density", "avg_edge_size"),
    value = c(nrow(incidence), ncol(incidence), sum(incidence) / length(incidence), mean(sizes)))
  expect_reference(subset(hg_measures(hg, what = "summary"), measure %in% summary_table$measure), summary_table,
                   reference = "independent: the incidence matrix, base R", covers = "hg_measures()",
                   quantity = "summary", tolerance = 1e-12, key = "measure")
})
