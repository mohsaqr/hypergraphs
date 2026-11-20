# hg_transitivity(): local clustering of the nodes of a hypergraph. "projection" is the local clustering
# coefficient of the simple clique expansion; "extra_overlap" is Klimm et al. (2021)'s local coefficient;
# "two_node_<kind>" are the two-node coefficients of Gallagher and Goldberg (2013). References: igraph and
# diag(A^3) / (k (k - 1)) in base R for the projection, and XGI for every type on six random hypergraphs.

memberships <- ref_memberships()

#' Six random membership tables, each with one hyperedge repeated and some of size 1
transitivity_configs <- function() {
  one <- function(seed, n_nodes, n_edges, max_size) {
    set.seed(seed)
    sizes <- sample(seq_len(max_size), n_edges, replace = TRUE, prob = c(0.1, rep(0.9 / (max_size - 1L), max_size - 1L)))
    edges <- lapply(sizes, \(s) sample(n_nodes, s))
    edges <- c(edges, edges[1L])
    data.frame(member = paste0("n", unlist(edges)), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges))))
  }
  list(one(21L, 8L, 10L, 4L), one(22L, 12L, 18L, 5L), one(23L, 15L, 25L, 3L),
       one(24L, 6L, 5L, 4L), one(25L, 20L, 30L, 6L), one(26L, 10L, 40L, 2L))
}

test_that("projection transitivity equals diag(A^3) / (k (k - 1)) and igraph's local transitivity", {
  incidence <- ref_incidence(memberships)
  adjacency <- (tcrossprod(incidence) > 0) * 1
  diag(adjacency) <- 0
  neighbours <- rowSums(adjacency)
  expected <- data.frame(node = rownames(adjacency),
                         projection = diag(adjacency %*% adjacency %*% adjacency) / (neighbours * (neighbours - 1)))
  ours <- hg_transitivity(group_hypergraph(memberships, node = "node", hyperedge = "edge"), type = "projection", what = "nodes")
  # A ratio of two small integers.
  expect_reference(ours, expected, reference = "independent: diag(A^3) / (k (k - 1)) on the clique expansion, base R",
                   covers = "hg_transitivity()", quantity = "projection", tolerance = 1e-12, key = "node")
  skip_without_reference("igraph")
  bipartite <- igraph::graph_from_data_frame(memberships, directed = FALSE)
  igraph::V(bipartite)$type <- igraph::V(bipartite)$name %in% memberships$edge
  expansion <- igraph::bipartite_projection(bipartite, multiplicity = FALSE, which = "false")
  expect_reference(ours, data.frame(node = igraph::V(expansion)$name, projection = igraph::transitivity(expansion, type = "local")),
                   reference = "igraph", covers = "hg_transitivity()", quantity = "projection", tolerance = 1e-12, key = "node")
})

test_that("every type of local clustering equals XGI on six random hypergraphs", {
  skip_without_python(c("xgi", "pandas"))
  script <- c(
    "import sys, json, xgi",
    "spec = json.load(open(sys.argv[1]))",
    "groups = {}",
    "for n, e in zip(spec['member'], spec['edge']): groups.setdefault(e, []).append(n)",
    "H = xgi.Hypergraph([groups[e] for e in sorted(groups)])",
    "cols = {'projection': xgi.clustering_coefficient(H), 'extra_overlap': xgi.local_clustering_coefficient(H),",
    "        'two_node_union': xgi.two_node_clustering_coefficient(H, kind='union'),",
    "        'two_node_min': xgi.two_node_clustering_coefficient(H, kind='min'),",
    "        'two_node_max': xgi.two_node_clustering_coefficient(H, kind='max')}",
    "nodes = sorted(H.nodes)",
    "out = {'node': nodes}",
    "out.update({k: [float(v[n]) for n in nodes] for k, v in cols.items()})",
    "json.dump(out, open(sys.argv[2], 'w'))")
  types <- c("projection", "extra_overlap", "two_node_union", "two_node_min", "two_node_max")
  version <- ref_python_version("xgi")
  # every hypergraph, node and type as one row, the six hypergraphs stacked
  long <- \(table, index) do.call(rbind, lapply(types, \(type) data.frame(hypergraph = index, node = table$node, type = type, value = table[[type]])))
  cells <- do.call(rbind, Map(\(config, index) {
    theirs <- as.data.frame(ref_python_run(script, as.list(config)))
    ours <- hg_transitivity(group_hypergraph(config, node = "member", hyperedge = "edge"), type = types, what = "nodes")
    merge(long(ours, index), long(theirs, index), by = c("hypergraph", "node", "type"), suffixes = c("_ours", "_xgi"), all = TRUE)
  }, transitivity_configs(), seq_along(transitivity_configs())))
  defined <- subset(cells, !is.na(value_ours))
  undefined <- subset(cells, is.na(value_ours))
  key <- c("hypergraph", "node", "type")
  # A handful of rational operations per value.
  expect_reference(setNames(defined[c(key, "value_ours")], c(key, "value")), setNames(defined[c(key, "value_xgi")], c(key, "value")),
                   reference = "xgi", covers = "hg_transitivity()", quantity = "six random hypergraphs, every type",
                   tolerance = 1e-12, key = key, note = version)
  expect_known_difference(setNames(undefined[c(key, "value_ours")], c(key, "value")), setNames(undefined[c(key, "value_xgi")], c(key, "value")),
                          reference = "xgi", covers = "hg_transitivity()", quantity = "six random hypergraphs, where undefined",
                          tolerance = 1e-12, key = key, note = version,
                          reason = "Where the coefficient is undefined (fewer than two neighbours in the projection, no neighbour for the two-node coefficients) hypergraphs returns NA and XGI returns 0.")
})
