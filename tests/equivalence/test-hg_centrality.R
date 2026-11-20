# hg_centrality() on a hypergraph, hg_pagerank(): node centralities. Clique eigenvector centrality is the leading
# eigenvector of the weighted clique expansion at unit length; PageRank (damping 0.85) is the hypergraph random
# walk, which for a binary incidence is the walk on W[u, v] = sum over shared hyperedges of 1 / size, self-loops
# included (Chitra and Raphael 2019); Katz is (I - alpha A)^-1 1 - 1 on the weighted clique expansion. Z- and
# H-eigenvector centralities of a uniform hypergraph (Benson 2019) against XGI; the stationary distribution of
# the weighted walk (damping 1) against HyperNetX.

memberships <- ref_memberships()
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")

test_that("clique eigenvector, PageRank and Katz centralities equal their closed forms", {
  incidence <- ref_incidence(memberships)
  nodes <- nrow(incidence)
  adjacency <- tcrossprod(incidence)
  diag(adjacency) <- 0
  leading <- abs(eigen(adjacency, symmetric = TRUE)$vectors[, 1L])
  walk <- incidence %*% (t(incidence) / colSums(incidence))
  transition <- walk / rowSums(walk)
  expected <- data.frame(node = rownames(incidence), clique = leading / sqrt(sum(leading^2)),
                         pagerank = solve(diag(nodes) - 0.85 * t(transition), rep((1 - 0.85) / nodes, nodes)),
                         katz = solve(diag(nodes) - 0.05 * adjacency, rep(1, nodes)) - 1)
  # The package iterates (power iteration to tol = 1e-12); the error of an iteration stopped at 1e-12 is that
  # divided by the spectral gap, so 1e-8 is the agreement that can be asked of it.
  expect_reference(hg_centrality(hg, type = c("clique", "pagerank", "katz"), damping = 0.85, alpha = 0.05,
                                 tol = 1e-12, max_iter = 10000L), expected,
                   reference = "independent: eigen() and solve() on the incidence matrix, base R",
                   covers = c("hg_centrality()", "hg_centrality.net_hg()"), quantity = "clique, PageRank, Katz",
                   tolerance = 1e-8, key = "node")
  expect_reference(hg_pagerank(hg, damping = 0.85, tol = 1e-12), subset(expected, select = c(node, pagerank)),
                   reference = "independent: eigen() and solve() on the incidence matrix, base R", covers = "hg_pagerank()",
                   quantity = "PageRank", tolerance = 1e-8, key = "node")
})

test_that("clique eigenvector, PageRank and Katz centralities equal igraph", {
  skip_without_reference("igraph")
  bipartite <- igraph::graph_from_data_frame(memberships, directed = FALSE)
  igraph::V(bipartite)$type <- igraph::V(bipartite)$name %in% memberships$edge
  expansion <- igraph::bipartite_projection(bipartite, multiplicity = TRUE, which = "false")
  eigenvector <- igraph::eigen_centrality(expansion, weights = igraph::E(expansion)$weight)$vector
  incidence <- igraph::as_biadjacency_matrix(bipartite)
  walk <- incidence %*% (t(incidence) / colSums(incidence))
  walk_graph <- igraph::graph_from_adjacency_matrix(walk, mode = "directed", weighted = TRUE, diag = TRUE)
  names <- igraph::V(expansion)$name
  expected <- data.frame(node = names, clique = unname(eigenvector / sqrt(sum(eigenvector^2))),
                         pagerank = unname(igraph::page_rank(walk_graph, damping = 0.85)$vector[names]),
                         katz = unname(igraph::alpha_centrality(expansion, alpha = 0.05, exo = 1,
                                                                weights = igraph::E(expansion)$weight)) - 1)
  # Iterative on both sides (power iteration here, ARPACK and PRPACK there): 1e-8, as above.
  expect_reference(hg_centrality(hg, type = c("clique", "pagerank", "katz"), damping = 0.85, alpha = 0.05,
                                 tol = 1e-12, max_iter = 10000L), expected, reference = "igraph",
                   covers = c("hg_centrality()", "hg_centrality.net_hg()"), quantity = "clique, PageRank, Katz",
                   tolerance = 1e-8, key = "node")
  expect_reference(hg_pagerank(hg, damping = 0.85, tol = 1e-12), subset(expected, select = c(node, pagerank)),
                   reference = "igraph", covers = "hg_pagerank()", quantity = "PageRank", tolerance = 1e-8, key = "node")
})

test_that("clique, Z- and H-eigenvector centralities of a 3-uniform hypergraph equal XGI", {
  skip_without_python(c("xgi", "numpy"))
  uniform <- data.frame(node = c("a", "b", "c", "a", "b", "d", "a", "c", "e", "b", "d", "e"),
                        edge = rep(paste0("e", seq_len(4L)), each = 3L))
  script <- c(
    "import sys, json, xgi",
    "spec = json.load(open(sys.argv[1]))",
    "edges = {}",
    "for n, e in zip(spec['node'], spec['edge']): edges.setdefault(e, []).append(n)",
    "H = xgi.Hypergraph(edges)",
    "clique = xgi.clique_eigenvector_centrality(H, tol=1e-12)",
    "z = xgi.z_eigenvector_centrality(H, max_iter=100000, tol=1e-12)",
    "h = xgi.h_eigenvector_centrality(H, max_iter=100000, tol=1e-12)",
    "nodes = sorted(H.nodes)",
    "json.dump({'node': nodes, 'clique': [float(clique[n]) for n in nodes], 'Z': [float(z[n]) for n in nodes],",
    "           'H': [float(h[n]) for n in nodes]}, open(sys.argv[2], 'w'))")
  theirs <- as.data.frame(ref_python_run(script, as.list(uniform)))
  ours <- hg_centrality(group_hypergraph(uniform, node = "node", hyperedge = "edge"), type = c("clique", "Z", "H"),
                        tol = 1e-12, max_iter = 10000L)
  # Eigenvectors are defined up to scale: both sides are divided by their sum before the comparison.
  to_shares <- \(table) data.frame(node = table$node, clique = table$clique / sum(table$clique),
                                   Z = table$Z / sum(table$Z), H = table$H / sum(table$H))
  # Iterative on both sides, each stopped at 1e-12: 1e-8 as above.
  expect_reference(to_shares(ours), to_shares(theirs), reference = "xgi",
                   covers = c("hg_centrality()", "hg_centrality.net_hg()"), quantity = "3-uniform: clique, Z, H (shares)",
                   tolerance = 1e-8, key = "node", note = ref_python_version("xgi"))
})

test_that("the stationary distribution of the weighted walk equals HyperNetX's get_pi", {
  skip_without_python(c("hypernetx", "pandas", "numpy"))
  set.seed(7)
  weighted <- data.frame(vertex = c("a", "b", "c", "b", "c", "d", "a", "d", "e", "c", "e"),
                         edge = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e3", "e4", "e4"),
                         w = round(stats::runif(11L, 0.5, 3), 3))
  script <- c(
    "import sys, json, numpy as np, pandas as pd, hypernetx as hnx",
    "from hypernetx.algorithms.clustering import laplacians_clustering as lc",
    "spec = json.load(open(sys.argv[1]))",
    "df = pd.DataFrame({'edge': spec['edge'], 'node': spec['vertex'], 'weight': spec['w']})",
    "H = hnx.Hypergraph(df, edge_col='edge', node_col='node', cell_weight_col='weight')",
    "P, idx = lc.prob_trans(H, weights=True)",
    "pi = [float(v) for v in np.asarray(lc.get_pi(P)).ravel()]",
    "json.dump({'node': [str(i) for i in idx], 'pagerank': pi}, open(sys.argv[2], 'w'))")
  theirs <- as.data.frame(ref_python_run(script, as.list(weighted)))
  ours <- hg_pagerank(group_hypergraph(weighted, node = "vertex", hyperedge = "edge", weight = "w"), damping = 1,
                      tol = 1e-15, max_iter = 100000L)
  # A power iteration stopped at 1e-15 against an eigen-solver: the tolerance of the standing battery.
  expect_reference(ours, theirs, reference = "hypernetx", covers = "hg_pagerank()",
                   quantity = "edge-dependent vertex weights, damping 1", tolerance = 1e-10, key = "node",
                   note = ref_python_version("hypernetx"))
})
