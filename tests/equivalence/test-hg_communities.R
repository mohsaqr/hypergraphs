# hg_communities() on a hypergraph (Infomap on the association graph, where each hyperedge of m nodes adds
# 1 / (m - 1) to each of its pairs), hg_community_quality() (coverage, performance, modularity and conductance of a
# partition of that graph) and hg_compare_communities() (agreement of two partitions). References: Python infomap
# on the association graph built from the definition, NetworkX's partition_quality(), modularity() and
# conductance(), and scikit-learn's AMI, ARI and NMI.

#' Four hypergraphs with communities, three hyperedges repeated, and a random labelling of their nodes
quality_configs <- function() {
  one <- function(seed, n_nodes, n_edges, k) {
    set.seed(seed)
    block <- rep(seq_len(k), length.out = n_nodes)
    edges <- lapply(seq_len(n_edges), \(i) {
      pool <- if (i %% 4L == 0L) seq_len(n_nodes) else which(block == sample(k, 1L))
      sample(pool, min(sample(2:4, 1L), length(pool)))
    })
    edges <- c(edges, edges[1:3])
    list(memberships = data.frame(member = sprintf("n%02d", unlist(edges)),
                                  edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges)))),
         labels = stats::setNames(sprintf("c%d", sample(k, n_nodes, replace = TRUE)), sprintf("n%02d", seq_len(n_nodes))))
  }
  list(one(61L, 15L, 20L, 3L), one(62L, 24L, 35L, 4L), one(63L, 30L, 40L, 2L), one(64L, 20L, 30L, 5L))
}

#' A planted hypergraph: three blocks of six nodes, ten hyperedges of three nodes inside each block, two across
planted <- function() {
  set.seed(61)
  blocks <- rep(1:3, each = 6L)
  nodes <- sprintf("n%02d", seq_len(18L))
  edges <- c(lapply(seq_len(30L), \(i) sample(nodes[blocks == ((i - 1L) %% 3L + 1L)], 3L)),
             list(c("n01", "n07"), c("n08", "n13")))
  data.frame(member = unlist(edges), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges))))
}

#' Each node labelled by the first node of its community, so two partitions compare label by label
canonical <- function(node, community) {
  first <- tapply(node, community, min)
  data.frame(node = node, community = unname(first[as.character(community)]))
}

association_graph <- c(
  "def association(spec, collapse=False):",
  "    groups = {}",
  "    for n, e in zip(spec['member'], spec['edge']): groups.setdefault(e, set()).add(n)",
  "    sets = [tuple(sorted(g)) for _, g in sorted(groups.items())]",
  "    if collapse: sets = sorted(set(sets))",
  "    weights = {}",
  "    for g in sets:",
  "        for i in range(len(g)):",
  "            for j in range(i + 1, len(g)):",
  "                weights[(g[i], g[j])] = weights.get((g[i], g[j]), 0.0) + 1.0 / (len(g) - 1)",
  "    return sorted(set(spec['member'])), weights")

test_that("coverage, performance, modularity and conductance equal NetworkX on both association graphs", {
  skip_without_python(c("networkx", "pandas"))
  script <- c(
    "import sys, json, networkx as nx",
    association_graph,
    "spec = json.load(open(sys.argv[1]))",
    "out = {}",
    "for mode in ('count', 'collapse'):",
    "    nodes, weights = association(spec, collapse=(mode == 'collapse'))",
    "    G = nx.Graph(); G.add_nodes_from(nodes)",
    "    for (u, v), w in weights.items(): G.add_edge(u, v, weight=w)",
    "    labels = dict(zip(spec['node'], spec['label']))",
    "    parts = [set(n for n in nodes if labels[n] == c) for c in sorted(set(labels.values()))]",
    "    B = nx.Graph(); B.add_nodes_from(G.nodes); B.add_edges_from(G.edges())",
    "    coverage, performance = nx.community.partition_quality(B, parts)",
    "    total = sum(w for _, _, w in G.edges(data='weight'))",
    "    inside = sum(w for u, v, w in G.edges(data='weight') if labels[u] == labels[v])",
    "    vol = dict(G.degree(weight='weight')); total_vol = sum(vol.values())",
    "    cond = [nx.conductance(G, S, weight='weight') for S in parts if 0 < sum(vol[n] for n in S) < total_vol]",
    "    out[mode] = {'coverage': coverage, 'weighted_coverage': inside / total, 'performance': performance,",
    "                 'modularity': nx.community.modularity(G, parts, weight='weight'),",
    "                 'conductance': max(cond), 'n_communities': len(parts)}",
    "json.dump(out, open(sys.argv[2], 'w'))")
  version <- ref_python_version("networkx")
  invisible(Map(\(config, index) {
    hg <- group_hypergraph(config$memberships, node = "member", hyperedge = "edge")
    nodes <- hg_get(hg, what = "nodes")$node
    labels <- config$labels[nodes]
    theirs <- ref_python_run(script, c(as.list(config$memberships), list(node = names(labels), label = unname(labels))))
    invisible(lapply(c("count", "collapse"), \(mode) {
      # Sums of rational weights over a few hundred pairs.
      expect_reference(unlist(hg_community_quality(hg, labels, duplicate_edges = mode)), unlist(theirs[[mode]]),
                       reference = "networkx", covers = "hg_community_quality()",
                       quantity = sprintf("hypergraph %d, duplicate_edges = %s", index, mode), tolerance = 1e-12,
                       note = version)
    }))
  }, quality_configs(), seq_along(quality_configs())))
})

test_that("Infomap's partition and codelength on a planted hypergraph equal Python infomap on the association graph", {
  skip_without_python("infomap")
  script <- c(
    "import sys, json, infomap",
    association_graph,
    "spec = json.load(open(sys.argv[1]))",
    "nodes, weights = association(spec)",
    "index = {n: i for i, n in enumerate(nodes)}",
    "im = infomap.Infomap('--two-level --silent --num-trials 3 --seed 1')",
    "for (u, v), w in weights.items(): im.add_link(index[u], index[v], w)",
    "im.run()",
    "modules = im.get_modules()",
    "json.dump({'node': nodes, 'community': [modules[index[n]] for n in nodes], 'codelength': im.codelength},",
    "          open(sys.argv[2], 'w'))")
  memberships <- planted()
  theirs <- ref_python_run(script, as.list(memberships))
  found <- hg_communities(group_hypergraph(memberships, node = "member", hyperedge = "edge"), n_runs = 3, trials = 3,
                          seeds = 1:3)
  version <- ref_python_version("infomap")
  medoid <- hg_get(found, what = "medoid")
  expect_reference(canonical(medoid$node, medoid$community), canonical(theirs$node, theirs$community), reference = "infomap",
                   covers = c("hg_communities()", "hg_communities.net_hg()", "hg_get.hg_communities()"),
                   quantity = "planted hypergraph, community of each node", tolerance = 0, key = "node", note = version)
  summarised <- hg_get(summary(found), what = "medoid")
  expect_reference(canonical(summarised$node, summarised$community), canonical(theirs$node, theirs$community), reference = "infomap",
                   covers = c("summary.hg_communities()", "hg_get.hypergraphs_summary()"),
                   quantity = "planted hypergraph, community of each node through summary()", tolerance = 0, key = "node",
                   note = version)
  # The codelength is the two-level map equation (Rosvall and Bergstrom 2008) of each run's partition, with current
  # Infomap's flow model: on this undirected graph, no teleportation. Python Infomap found the same partition, so its
  # own optimum is the same quantity.
  expect_reference(hg_get(found, what = "runs"), data.frame(run = 1:3, codelength = theirs$codelength),
                   reference = "infomap", covers = c("hg_communities()", "hg_communities.net_hg()", "hg_get.hg_communities()"),
                   quantity = "planted hypergraph, codelength", tolerance = 1e-10, key = "run", note = version)
  expect_reference(hg_get(summary(found), what = "runs"), data.frame(run = 1:3, codelength = theirs$codelength),
                   reference = "infomap", covers = c("summary.hg_communities()", "hg_get.hypergraphs_summary()"),
                   quantity = "planted hypergraph, codelength through summary()", tolerance = 1e-10, key = "run",
                   note = version)
  # the two-level map equation (Rosvall and Bergstrom 2008) of Infomap's partition, in bits: node flow p = strength /
  # total strength, module exit flow q_i, L = q log q - 2 sum q_i log q_i - sum p log p + sum (q_i + p_i) log (q_i + p_i)
  incidence <- ref_incidence(memberships, node = "member", edge = "edge")
  association <- incidence %*% diag(1 / (colSums(incidence) - 1)) %*% t(incidence)
  diag(association) <- 0
  module <- theirs$community[match(rownames(association), theirs$node)]
  total <- sum(association)
  flow <- rowSums(association) / total
  exits <- vapply(sort(unique(module)), \(m) sum(association[module == m, module != m]) / total, numeric(1))
  inside <- vapply(sort(unique(module)), \(m) sum(flow[module == m]), numeric(1))
  plogp <- \(x) ifelse(x > 0, x * log2(x), 0)
  map_equation <- plogp(sum(exits)) - 2 * sum(plogp(exits)) - sum(plogp(flow)) + sum(plogp(exits + inside))
  expect_reference(hg_get(found, what = "runs"), data.frame(run = 1:3, codelength = map_equation),
                   reference = "independent: two-level map equation, base R",
                   covers = c("hg_communities()", "hg_communities.net_hg()", "hg_get.hg_communities()"),
                   quantity = "planted hypergraph, codelength", tolerance = 1e-10, key = "run")
})

test_that("the codelength of every run's partition equals Python infomap's on the same partition, directed and undirected", {
  skip_without_python("infomap")
  # Infomap with --no-infomap scores a given partition without optimising: the codelength of exactly our partition.
  script <- c(
    "import sys, json, infomap",
    "s = json.load(open(sys.argv[1]))",
    "out = []",
    "for m in s['modules']:",
    "    im = infomap.Infomap('--two-level --silent --no-infomap' + (' --directed' if s['directed'] else ''))",
    "    for k in range(s['n']): im.add_node(k)",
    "    for a, b, w in zip(s['from'], s['to'], s['weight']): im.add_link(a, b, w)",
    "    im.run(initial_partition={k: c for k, c in enumerate(m)})",
    "    out.append(im.codelength)",
    "json.dump({'codelength': out}, open(sys.argv[2], 'w'))")
  memberships <- planted()
  hg <- group_hypergraph(memberships, node = "member", hyperedge = "edge")
  # each hyperedge cites from its first member (a column constant within a hyperedge is a hyperedge attribute)
  first_member <- memberships$member[!duplicated(memberships$edge)]
  cited <- group_hypergraph(within(memberships, source <- first_member[match(edge, unique(edge))]),
                            node = "member", hyperedge = "edge")
  version <- ref_python_version("infomap")
  settings <- list(
    list(hg = hg, args = list(), label = "association graph"),
    list(hg = hg, args = list(duplicate_edges = "collapse"), label = "association graph, duplicates collapsed"),
    list(hg = cited, args = list(method = "citation"), label = "undirected citation graph"),
    list(hg = cited, args = list(method = "citation", directed = TRUE), label = "directed citation graph"))
  invisible(lapply(settings, \(setting) {
    found <- do.call(hg_communities, c(list(setting$hg, n_runs = 3, trials = 3, seeds = 1:3), setting$args))
    graph <- found$projection
    directed <- isTRUE(setting$args$directed)
    links <- which(graph != 0 & (directed | upper.tri(graph, diag = TRUE)), arr.ind = TRUE)
    partitions <- hg_get(found, what = "partitions")
    modules <- lapply(split(partitions, partitions$run), \(p) {
      as.integer(factor(p$community[match(rownames(graph), p$node)])) - 1L
    })
    theirs <- ref_python_run(script, list(n = nrow(graph), from = unname(links[, 1L]) - 1L, to = unname(links[, 2L]) - 1L,
                                          weight = graph[links], modules = unname(modules), directed = directed))
    expect_reference(hg_get(found, what = "runs"), data.frame(run = 1:3, codelength = unlist(theirs$codelength)),
                     reference = "infomap", covers = c("hg_communities()", "hg_communities.net_hg()", "hg_get.hg_communities()"),
                     quantity = sprintf("planted hypergraph, %s, codelength of each run's partition", setting$label),
                     tolerance = 1e-10, key = "run", note = version)
  }))
})

test_that("the agreement of two partitions equals scikit-learn's AMI, ARI and NMI", {
  skip_without_python(c("scikit-learn" = "sklearn"))
  script <- c(
    "import sys, json",
    "from sklearn.metrics import adjusted_mutual_info_score as ami, adjusted_rand_score as ari, normalized_mutual_info_score as nmi",
    "spec = json.load(open(sys.argv[1]))",
    "json.dump({'ami': ami(spec['a'], spec['b']), 'ari': ari(spec['a'], spec['b']), 'nmi': nmi(spec['a'], spec['b'])},",
    "          open(sys.argv[2], 'w'))")
  version <- ref_python_version("scikit-learn")
  invisible(Map(\(config, index) {
    hg <- group_hypergraph(config$memberships, node = "member", hyperedge = "edge")
    counted <- hg_communities(hg, n_runs = 3, trials = 3, seeds = 1:3)
    collapsed <- hg_communities(hg, n_runs = 3, trials = 3, seeds = 1:3, duplicate_edges = "collapse")
    comparison <- hg_compare_communities(mh = counted, bh = collapsed, hg = hg)
    shared <- merge(hg_get(counted, what = "medoid"), hg_get(collapsed, what = "medoid"), by = "node", suffixes = c("_a", "_b"))
    theirs <- ref_python_run(script, list(a = shared$community_a, b = shared$community_b))
    # Entropies and expected mutual information of small contingency tables.
    expect_reference(unlist(subset(hg_get(comparison, what = "similarity"), select = c(ami, ari, nmi))), unlist(theirs),
                     reference = "scikit-learn",
                     covers = c("hg_compare_communities()", "hg_get.hypergraphs_community_comparison()"),
                     quantity = sprintf("hypergraph %d, counted against collapsed duplicates", index), tolerance = 1e-12,
                     note = version)
    expect_reference(unlist(subset(hg_get(summary(comparison), what = "similarity"), select = c(ami, ari, nmi))), unlist(theirs),
                     reference = "scikit-learn", covers = c("summary.hypergraphs_community_comparison()", "hg_get.hypergraphs_summary()"),
                     quantity = sprintf("hypergraph %d, through summary()", index), tolerance = 1e-12, note = version)
  }, quality_configs(), seq_along(quality_configs())))
})
