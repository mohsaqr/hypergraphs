# hg_cluster(), hg_membership(): spectral clustering of a hypergraph on the Zhou Laplacian (Zhou, Huang and
# Scholkopf 2007) and the soft membership of every node in every cluster (fuzzy c-means weights with fuzziness 2
# around the cluster centres in the row-normalised embedding, Bezdek 1981). References: the spectrum of the
# Laplacian by eigen() in base R; HyperNetX's spec_clus() on a planted hypergraph; XGI's Laplacian with numpy's
# eigh() and scikit-fuzzy's cmeans_predict() for the memberships.

#' A planted hypergraph: k blocks with a few bridging hyperedges, so the clusters are real and it is connected
planted <- function(seed, n_nodes, n_edges, k) {
  set.seed(seed)
  block <- rep(seq_len(k), length.out = n_nodes)
  edges <- lapply(seq_len(n_edges), \(i) {
    pool <- if (i %% 5L == 0L) seq_len(n_nodes) else which(block == sample(k, 1L))
    sample(pool, min(sample(2:4, 1L), length(pool)))
  })
  list(memberships = data.frame(member = sprintf("n%02d", unlist(edges)), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges)))),
       k = k)
}

configs <- list(planted(41L, 12L, 18L, 2L), planted(42L, 18L, 30L, 3L), planted(43L, 24L, 40L, 4L), planted(44L, 30L, 45L, 3L))

#' Each node labelled by the first node of its cluster, so two partitions compare label by label
canonical <- function(node, cluster) {
  first <- tapply(node, cluster, min)
  data.frame(node = node, cluster = unname(first[as.character(cluster)]))
}

test_that("the spectrum of the Zhou Laplacian equals eigen() of its formula", {
  invisible(Map(\(config, index) {
    incidence <- ref_incidence(config$memberships, node = "member", edge = "edge") > 0
    degree_root <- diag(1 / sqrt(rowSums(incidence)))
    laplacian <- diag(nrow(incidence)) - degree_root %*% incidence %*% diag(1 / colSums(incidence)) %*% t(incidence) %*% degree_root
    values <- sort(eigen(laplacian, symmetric = TRUE, only.values = TRUE)$values)
    hg <- group_hypergraph(config$memberships, node = "member", hyperedge = "edge")
    # Eigenvalues of a symmetric matrix of order 30 by LAPACK on both sides: rounding near eps times the norm.
    expect_reference(hg_cluster(hg, k = config$k, seed = 1L, what = "eigenvalues"),
                     data.frame(index = seq_along(values), value = values),
                     reference = "independent: eigen() of the Zhou Laplacian, base R", covers = "hg_cluster()",
                     quantity = sprintf("hypergraph %d, Laplacian spectrum", index), tolerance = 1e-10, key = "index")
  }, configs, seq_along(configs)))
})

test_that("the clusters equal HyperNetX's spectral clustering on planted hypergraphs", {
  skip_without_python(c("hypernetx", "pandas", "scikit-learn" = "sklearn"))
  script <- c(
    "import sys, json, pandas as pd, hypernetx as hnx",
    "from hypernetx.algorithms.clustering import laplacians_clustering as lc",
    "spec = json.load(open(sys.argv[1]))",
    "H = hnx.Hypergraph(pd.DataFrame({'edges': spec['edge'], 'nodes': spec['member']}), edge_col='edges', node_col='nodes')",
    "clusters = lc.spec_clus(H, spec['k'], weights=False)",
    "node, cluster = [], []",
    "for c, members in clusters.items():",
    "    for n in members: node.append(str(n)); cluster.append(int(c))",
    "json.dump({'node': node, 'cluster': cluster}, open(sys.argv[2], 'w'))")
  version <- ref_python_version("hypernetx")
  invisible(Map(\(config, index) {
    theirs <- ref_python_run(script, c(as.list(config$memberships), list(k = config$k)))
    ours <- hg_cluster(group_hypergraph(config$memberships, node = "member", hyperedge = "edge"), k = config$k, seed = 1L)
    expect_reference(canonical(ours$node, ours$cluster), canonical(theirs$node, theirs$cluster), reference = "hypernetx",
                     covers = "hg_cluster()", quantity = sprintf("planted hypergraph %d, cluster of each node", index),
                     tolerance = 0, key = "node",
                     note = paste(version, "laplacians_clustering.spec_clus (k-means on the spectral embedding, scikit-learn)"))
  }, configs, seq_along(configs)))
})

test_that("soft memberships equal scikit-fuzzy's c-means memberships in XGI's spectral embedding", {
  skip_without_python(c("xgi", "numpy", "scikit-fuzzy" = "skfuzzy"))
  script <- c(
    "import sys, json, numpy as np, xgi, skfuzzy",
    "spec = json.load(open(sys.argv[1]))",
    "groups = {}",
    "for n, e in zip(spec['member'], spec['edge']): groups.setdefault(e, []).append(n)",
    "H = xgi.Hypergraph([groups[e] for e in sorted(groups)])",
    "L, idx = xgi.normalized_hypergraph_laplacian(H, weighted=False, sparse=False, index=True)",
    "nodes = [idx[i] for i in range(len(idx))]",
    "k = spec['k']",
    "vals, vecs = np.linalg.eigh((L + L.T) / 2)",
    "U = vecs[:, :k]",
    "norm = np.linalg.norm(U, axis=1)",
    "U = U / np.where(norm > 0, norm, 1)[:, None]",
    "label = dict(zip(spec['node'], spec['cluster']))",
    "labels = np.array([label[n] for n in nodes])",
    "clusters = sorted(set(labels))",
    "centres = np.vstack([U[labels == c].mean(axis=0) for c in clusters])",
    "u, u0, dist, jm, p, fpc = skfuzzy.cmeans_predict(U.T, centres, 2.0, error=0.0, maxiter=5)",
    "json.dump({'node': [n for n in nodes for c in clusters], 'topic': [c for n in nodes for c in clusters],",
    "           'membership': [float(u[j, i]) for i in range(len(nodes)) for j in range(len(clusters))],",
    "           'gap': float(vals[k] - vals[k - 1])}, open(sys.argv[2], 'w'))")
  version <- paste(ref_python_version("scikit-fuzzy"), "with", ref_python_version("xgi"))
  invisible(Map(\(config, index) {
    hg <- group_hypergraph(config$memberships, node = "member", hyperedge = "edge")
    clusters <- hg_cluster(hg, k = config$k, seed = 1L)
    theirs <- ref_python_run(script, c(as.list(config$memberships), list(k = config$k, node = clusters$node, cluster = clusters$cluster)))
    # Both sides solve a dense symmetric eigenproblem (LAPACK through R and through numpy); memberships are smooth
    # in the embedding and the eigenvector error is of order eps / gap: the tolerance of the standing battery.
    expect_reference(hg_membership(hg, clusters), data.frame(node = theirs$node, topic = theirs$topic, membership = theirs$membership),
                     reference = "scikit-fuzzy", covers = "hg_membership()",
                     quantity = sprintf("planted hypergraph %d, membership", index), tolerance = 1e-8,
                     key = c("node", "topic"), note = sprintf("%s; eigengap after k: %.3g", version, theirs$gap))
  }, configs, seq_along(configs)))
})

test_that("the stationary distribution of the embedding and the seeds equal their definitions", {
  invisible(Map(\(config, index) {
    hg <- group_hypergraph(config$memberships, node = "member", hyperedge = "edge")
    embedding <- hg_cluster(hg, k = config$k, seed = 1L, what = "embedding")
    incidence <- ref_incidence(config$memberships, node = "member", edge = "edge") > 0
    # Zhou's walk with unit hyperedge weights: pi(v) = d(v) / sum d
    stationary <- rowSums(incidence) / sum(incidence)
    # hg_cluster(what = "embedding") returns the table hg_get.net_hg_cluster() builds from the fit (it adds the `pi`
    # column), so this comparison is evidence for that method too: the only route by which the exported API reaches it.
    expect_reference(subset(embedding, select = c(node, pi)), data.frame(node = names(stationary), pi = unname(stationary)),
                     reference = "independent: degree over total degree, the stationary distribution of Zhou et al. (2007)",
                     covers = c("hg_cluster()", "hg_get.net_hg_cluster()"), quantity = sprintf("hypergraph %d, stationary distribution", index), tolerance = 1e-12,
                     key = "node")
    # the seeds of a cluster: its n nodes of largest pi, ties broken by node name
    ranked <- embedding[order(embedding$cluster, -round(stationary[embedding$node], 12), embedding$node), ]
    chosen <- do.call(rbind, lapply(split(ranked, ranked$cluster), utils::head, 3L))
    expect_reference(data.frame(node = names(hg_seeds(embedding, n = 3L)), cluster = unname(hg_seeds(embedding, n = 3L))),
                     data.frame(node = chosen$node, cluster = chosen$cluster),
                     reference = "independent: the n nodes of largest stationary probability per cluster, base R",
                     covers = "hg_seeds()", quantity = sprintf("hypergraph %d, three seeds per cluster", index), tolerance = 0,
                     key = "node")
  }, configs, seq_along(configs)))
})
