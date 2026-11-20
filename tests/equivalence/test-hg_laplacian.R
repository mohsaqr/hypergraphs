# hg_laplacian(): the normalised hypergraph Laplacian of Zhou, Huang and Scholkopf (2007),
# I - Dv^-1/2 H W De^-1 H' Dv^-1/2, and the random-walk Laplacian of Hayashi et al. (2020) with edge-dependent
# vertex weights. References: the formula in base R; XGI's normalized_hypergraph_laplacian(); HyperNetX's
# norm_lap() (laplacians_clustering, the authors' group at PNNL).

#' Eight random connected hypergraphs of 5 to 12 nodes, the first with one hyperedge repeated, and edge weights
laplacian_configs <- function() {
  draw <- function(seed) {
    set.seed(seed)
    n <- sample(c(5L, 8L, 12L), 1L)
    members <- lapply(seq_len(sample(c(4L, 7L, 11L), 1L)), \(j) sort(sample(n, sample(2:min(5L, n), 1L))))
    list(n = n, members = members, w = round(stats::runif(length(members), 0.25, 4), 3))
  }
  connected <- function(config) {
    incidence <- vapply(config$members, \(edge) seq_len(config$n) %in% edge, logical(config$n))
    graph <- igraph::graph_from_adjacency_matrix((tcrossprod(incidence * 1) > 0) * 1, mode = "undirected")
    all(rowSums(incidence) > 0) && igraph::is_connected(graph)
  }
  configs <- utils::head(Filter(connected, lapply(seq_len(200L) * 17L, draw)), 8L)
  configs[[1L]]$members <- c(configs[[1L]]$members, configs[[1L]]$members[1L])
  configs[[1L]]$w <- c(configs[[1L]]$w, 1.5)
  lapply(configs, \(config) {
    data.frame(node = paste0("v", unlist(config$members)),
               edge = sprintf("h%03d", rep(seq_along(config$members), lengths(config$members))),
               w = rep(config$w, lengths(config$members)))
  })
}

test_that("the Zhou Laplacian equals its formula on the incidence matrix", {
  memberships <- ref_memberships()
  incidence <- ref_incidence(memberships)
  degree_root <- diag(1 / sqrt(rowSums(incidence)))
  size_inverse <- diag(1 / colSums(incidence))
  expected <- diag(nrow(incidence)) - degree_root %*% incidence %*% size_inverse %*% t(incidence) %*% degree_root
  dimnames(expected) <- list(rownames(incidence), rownames(incidence))
  # A few products of numbers of order 1.
  expect_reference(hg_laplacian(group_hypergraph(memberships, node = "node", hyperedge = "edge"), type = "zhou"),
                   expected, reference = "independent: the formula of Zhou, Huang and Scholkopf (2007), base R",
                   covers = "hg_laplacian()", quantity = "Zhou Laplacian", tolerance = 1e-12)
})

test_that("the Zhou Laplacian equals XGI and HyperNetX with unit weights, and XGI with edge weights", {
  skip_without_reference("igraph")
  skip_without_python(c("xgi", "hypernetx", "numpy", "pandas"))
  script <- c(
    "import sys, json, numpy as np, pandas as pd, xgi, hypernetx as hnx",
    "from hypernetx.algorithms.clustering.laplacians_clustering import norm_lap",
    "spec = json.load(open(sys.argv[1]))",
    "edges = {}",
    "for n, e in zip(spec['node'], spec['edge']): edges.setdefault(e, []).append(n)",
    "weights = dict(zip(spec['edge'], spec['w']))",
    "H = xgi.Hypergraph(edges)",
    "L0, idx = xgi.normalized_hypergraph_laplacian(H, weighted=False, sparse=False, index=True)",
    "H.set_edge_attributes(weights, name='weight')",
    "L1, idx1 = xgi.normalized_hypergraph_laplacian(H, weighted=True, sparse=False, index=True)",
    "HX = hnx.Hypergraph(pd.DataFrame({'edges': spec['edge'], 'nodes': spec['node']}), edge_col='edges', node_col='nodes')",
    "L2, idx2 = norm_lap(HX, weights=False, index=True)",
    "L2 = L2.toarray()",
    "labels = [str(idx[i]) for i in range(len(idx))]",
    "at2 = {str(v): i for i, v in enumerate(idx2)}",
    "out = {'from': [], 'to': [], 'xgi_unit': [], 'xgi_weighted': [], 'hnx_unit': []}",
    "for i, a in enumerate(labels):",
    "    for j, b in enumerate(labels):",
    "        out['from'].append(a); out['to'].append(b)",
    "        out['xgi_unit'].append(float(L0[i, j])); out['xgi_weighted'].append(float(L1[i, j]))",
    "        out['hnx_unit'].append(float(L2[at2[a], at2[b]]))",
    "json.dump(out, open(sys.argv[2], 'w'))")
  xgi_version <- ref_python_version("xgi")
  hnx_version <- ref_python_version("hypernetx")
  invisible(Map(\(config, index) {
    theirs <- as.data.frame(ref_python_run(script, as.list(config)))
    hg <- group_hypergraph(subset(config, select = c(node, edge)), node = "node", hyperedge = "edge")
    edge_weights <- unique(subset(config, select = c(edge, w)))
    ours <- ref_long(hg_laplacian(hg, type = "zhou"), name = "value")
    label <- sprintf("hypergraph %d", index)
    # A few products of numbers of order 1 on each side.
    expect_reference(ours, data.frame(from = theirs$from, to = theirs$to, value = theirs$xgi_unit), reference = "xgi",
                     covers = "hg_laplacian()", quantity = paste(label, "Zhou Laplacian, unit weights"), tolerance = 1e-12,
                     key = c("from", "to"), note = xgi_version)
    expect_reference(ours, data.frame(from = theirs$from, to = theirs$to, value = theirs$hnx_unit), reference = "hypernetx",
                     covers = "hg_laplacian()", quantity = paste(label, "Zhou Laplacian, unit weights"), tolerance = 1e-12,
                     key = c("from", "to"), note = hnx_version)
    weighted <- ref_long(hg_laplacian(hg, type = "zhou", edge_weights = edge_weights$w[match(hg_get(hg, what = "edges")$hyperedge, edge_weights$edge)]),
                         name = "value")
    expect_known_difference(weighted, data.frame(from = theirs$from, to = theirs$to, value = theirs$xgi_weighted),
                            reference = "xgi", covers = "hg_laplacian()", quantity = paste(label, "Zhou Laplacian, edge weights"),
                            tolerance = 1e-12, key = c("from", "to"), note = xgi_version,
                            reason = "With edge weights, hypergraphs takes the node degree Dv as the sum of the weights of the hyperedges a node is in, d(v) = sum w(e) as Zhou et al. (2007) define it; XGI's normalized_hypergraph_laplacian(weighted = True) weights the hyperedges but keeps the unweighted degree (the standing battery found the two equal after rescaling by sqrt(count / weighted degree)).")
  }, laplacian_configs(), seq_along(laplacian_configs())))
})

test_that("the random-walk Laplacian and its stationary distribution equal HyperNetX's norm_lap and get_pi", {
  skip_without_reference("igraph")
  skip_without_python(c("hypernetx", "numpy", "pandas"))
  script <- c(
    "import sys, json, numpy as np, pandas as pd, hypernetx as hnx",
    "from hypernetx.algorithms.clustering import laplacians_clustering as lc",
    "spec = json.load(open(sys.argv[1]))",
    "out = {}",
    "for weighted in (False, True):",
    "    df = pd.DataFrame({'edges': spec['edge'], 'nodes': spec['node'], 'cell_weights': spec['w'] if weighted else [1.0] * len(spec['w'])})",
    "    H = hnx.Hypergraph(df, edge_col='edges', node_col='nodes', cell_weight_col='cell_weights')",
    "    L, idx = lc.norm_lap(H, weights=weighted, index=True)",
    "    L = np.asarray(L.todense()) if hasattr(L, 'todense') else np.asarray(L)",
    "    labels = [str(v) for v in idx]",
    "    key = 'weighted' if weighted else 'unit'",
    "    out[key] = {'from': [a for a in labels for b in labels], 'to': [b for a in labels for b in labels],",
    "                'value': [float(L[i, j]) for i in range(len(labels)) for j in range(len(labels))]}",
    "json.dump(out, open(sys.argv[2], 'w'))")
  version <- ref_python_version("hypernetx")
  invisible(Map(\(config, index) {
    theirs <- ref_python_run(script, as.list(config))
    unit <- group_hypergraph(subset(config, select = c(node, edge)), node = "node", hyperedge = "edge")
    weighted <- group_hypergraph(config, node = "node", hyperedge = "edge", weight = "w")
    # A few products of numbers of order 1 on each side.
    expect_reference(ref_long(hg_laplacian(unit, type = "random_walk"), name = "value"), as.data.frame(theirs$unit),
                     reference = "hypernetx", covers = "hg_laplacian()",
                     quantity = sprintf("hypergraph %d, random-walk Laplacian, unit weights", index), tolerance = 1e-12,
                     key = c("from", "to"), note = version)
    expect_reference(ref_long(hg_laplacian(weighted, type = "random_walk"), name = "value"), as.data.frame(theirs$weighted),
                     reference = "hypernetx", covers = "hg_laplacian()",
                     quantity = sprintf("hypergraph %d, random-walk Laplacian, vertex weights", index), tolerance = 1e-12,
                     key = c("from", "to"), note = version)
  }, laplacian_configs(), seq_along(laplacian_configs())))
})
