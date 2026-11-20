# hg_centrality() on a higher-order network: PageRank, betweenness and closeness of the higher-order nodes,
# projected to the first-order states (Scholtes, Wider and Garas 2016), or unprojected. References: pathpy, the
# authors' implementation, on pathpy's own k-th order topology (wrapped as a net_hon so that both compute on one
# graph), and igraph on the unprojected network of hon().

pathpy_script <- c(
  "import sys, json, scipy, numpy",
  "# pathpy 2.2.0 imports names that SciPy has since removed; they are taken from NumPy, pathpy is not changed",
  "for name in dir(numpy):",
  "    if not hasattr(scipy, name):",
  "        try: setattr(scipy, name, getattr(numpy, name))",
  "        except Exception: pass",
  "import pathpy as pp",
  "import pathpy.algorithms.centralities as centralities",
  "spec = json.load(open(sys.argv[1]))",
  "paths = pp.Paths()",
  "for s in spec['sequences']: paths.add_path(','.join(s), frequency=1)",
  "network = pp.HigherOrderNetwork(paths, k=int(spec['k']))",
  "# the power iteration is run to 1e-14, where pathpy's default stops at 1e-6",
  "pagerank = dict(centralities.pagerank(network, max_iter=100000, tol=1e-14))",
  "betweenness = dict(centralities.betweenness(network))",
  "closeness = dict(centralities.closeness(network))",
  "states = sorted(paths.nodes)",
  "json.dump({'nodes': sorted(network.nodes), 'from': [a for a, b in network.edges], 'to': [b for a, b in network.edges],",
  "           'state': states, 'pagerank': [float(pagerank[s]) for s in states],",
  "           'betweenness': [float(betweenness[s]) for s in states], 'closeness': [float(closeness[s]) for s in states]},",
  "          open(sys.argv[2], 'w'))")

#' Sequences with two second-order patterns and random ones
corpus <- function(seed) {
  set.seed(seed)
  c(replicate(10L, rep(c("a", "b", "c"), 4L), simplify = FALSE), replicate(10L, rep(c("x", "b", "d"), 4L), simplify = FALSE),
    replicate(10L, sample(c("a", "b", "c", "d", "x"), 12L, replace = TRUE), simplify = FALSE))
}

#' pathpy's topology as a net_hon: its higher-order nodes and edges, unit weights
as_hon <- function(theirs) {
  nodes <- gsub(",", " -> ", theirs$nodes, fixed = TRUE)
  adjacency <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  adjacency[cbind(gsub(",", " -> ", theirs$from, fixed = TRUE), gsub(",", " -> ", theirs$to, fixed = TRUE))] <- 1
  structure(list(matrix = adjacency, first_order_states = sort(theirs$state)), class = c("net_hon", "cograph_network"))
}

test_that("projected PageRank, betweenness and closeness equal pathpy on its k-th order topologies", {
  skip_without_python("pathpy")
  version <- ref_python_version("pathpy")
  inputs <- c(lapply(1:3, \(k) list(sequences = corpus(42L), k = k)),
              lapply(1:4, \(i) {
                set.seed(700L + i)
                states <- letters[seq_len(sample(3:6, 1L))]
                list(sequences = replicate(sample(8:20, 1L), sample(states, sample(8:25, 1L), replace = TRUE), simplify = FALSE),
                     k = sample(1:3, 1L))
              }))
  invisible(Map(\(input, index) {
    theirs <- ref_python_run(pathpy_script, input)
    ours <- hg_centrality(as_hon(theirs), tol = 1e-14, max_iter = 100000L)
    # Shortest paths on both sides; PageRank by power iteration on both sides, each run to 1e-14: 1e-8, as for
    # the other iterative centralities.
    expect_reference(subset(ours, select = c(state, pagerank, betweenness, closeness)),
                     data.frame(state = theirs$state, pagerank = theirs$pagerank, betweenness = theirs$betweenness,
                                closeness = theirs$closeness),
                     reference = "pathpy", covers = c("hg_centrality()", "hg_centrality.net_hon()"),
                     quantity = sprintf("corpus %d, order %d, projected centralities", index, input$k), tolerance = 1e-8,
                     key = "state", note = version)
  }, inputs, seq_along(inputs)))
})

test_that("unprojected PageRank, betweenness and harmonic closeness equal igraph on the network of hon()", {
  skip_without_reference("igraph")
  invisible(lapply(1:4, \(i) {
    set.seed(800L + i)
    sequences <- replicate(12L, sample(letters[1:5], sample(10:30, 1L), replace = TRUE), simplify = FALSE)
    network <- hon(sequences, max_order = 2L, min_freq = 1L)
    rules <- hg_get(network, what = "rules")
    states <- strsplit(rules$path, " -> ", fixed = TRUE)
    target <- mapply(\(path, order) paste(utils::tail(path, order), collapse = " -> "), states, rules$to_order)
    graph <- igraph::graph_from_data_frame(data.frame(from = rules$from, to = target), directed = TRUE,
                                           vertices = hg_get(network, what = "nodes")$node)
    names <- igraph::V(graph)$name
    expected <- data.frame(node = names, pagerank = unname(igraph::page_rank(graph, damping = 0.85)$vector),
                           betweenness = unname(igraph::betweenness(graph, directed = TRUE)),
                           closeness = unname(igraph::harmonic_centrality(graph, mode = "in", normalized = FALSE)))
    # PageRank by PRPACK against a power iteration to 1e-12: 1e-8.
    expect_reference(hg_centrality(network, project = FALSE, tol = 1e-12), expected, reference = "igraph",
                     covers = c("hg_centrality()", "hg_centrality.net_hon()"),
                     quantity = sprintf("network %d, unprojected centralities", i), tolerance = 1e-8, key = "node")
  }))
})
