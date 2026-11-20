# hg_line_graph(), hg_edge_centrality(): the s-line graph of a hypergraph (hyperedges joined when they share at
# least s nodes, weighted by the number they share) and the shortest-path centralities of the hyperedges on it
# (betweenness scaled by 2 / ((n - 1)(n - 2)), closeness as NetworkX defines it). References: H'H in base R,
# igraph's bipartite projection, and HyperNetX on six random hypergraphs for s = 1, 2, 3.

memberships <- ref_memberships()
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")

#' Six random membership tables of 8 to 30 hyperedges, each with one hyperedge repeated
line_configs <- function() {
  one <- function(seed, n_nodes, n_edges, max_size) {
    set.seed(seed)
    sizes <- sample(2:max_size, n_edges, replace = TRUE)
    edges <- lapply(sizes, \(s) sample(n_nodes, s))
    edges <- c(edges, edges[1L])
    data.frame(member = paste0("n", unlist(edges)), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges))),
               stringsAsFactors = FALSE)
  }
  list(one(31L, 10L, 8L, 4L), one(32L, 15L, 14L, 5L), one(33L, 25L, 20L, 6L),
       one(34L, 40L, 12L, 3L), one(35L, 12L, 25L, 6L), one(36L, 30L, 30L, 4L))
}

test_that("the 1- and 2-line graphs equal the off-diagonal of H'H", {
  shared <- crossprod(ref_incidence(memberships))
  invisible(lapply(1:2, \(s) {
    expect_reference(hg_line_graph(hg, s = s, what = "edges"), ref_pairs_upper(shared, keep = \(m) m >= s),
                     reference = "independent: H'H of the incidence matrix, base R", covers = "hg_line_graph()",
                     quantity = sprintf("s = %d line graph weight", s), tolerance = 0, key = c("from", "to"))
  }))
})

test_that("the 1- and 2-line graphs equal igraph's projection on the hyperedges", {
  skip_without_reference("igraph")
  bipartite <- igraph::graph_from_data_frame(memberships, directed = FALSE)
  igraph::V(bipartite)$type <- igraph::V(bipartite)$name %in% memberships$edge
  projection <- igraph::bipartite_projection(bipartite, multiplicity = TRUE, which = "true")
  pairs <- igraph::as_data_frame(projection, what = "edges")
  expected <- cbind(ref_ordered_pairs(pairs$from, pairs$to), weight = pairs$weight)
  expect_reference(hg_line_graph(hg, s = 1, what = "edges"), expected, reference = "igraph", covers = "hg_line_graph()",
                   quantity = "s = 1 line graph weight", tolerance = 0, key = c("from", "to"))
  expect_reference(hg_line_graph(hg, s = 2, what = "edges"), subset(expected, weight >= 2), reference = "igraph",
                   covers = "hg_line_graph()", quantity = "s = 2 line graph weight", tolerance = 0, key = c("from", "to"))
})

test_that("s-betweenness and s-closeness of the hyperedges equal igraph on the 1-line graph", {
  skip_without_reference("igraph")
  bipartite <- igraph::graph_from_data_frame(memberships, directed = FALSE)
  igraph::V(bipartite)$type <- igraph::V(bipartite)$name %in% memberships$edge
  line <- igraph::bipartite_projection(bipartite, multiplicity = FALSE, which = "true")
  # the line graph of this hypergraph is connected, so igraph's closeness and NetworkX's agree here
  expected <- rbind(
    data.frame(edge = igraph::V(line)$name, measure = "betweenness", value = unname(igraph::betweenness(line, normalized = TRUE))),
    data.frame(edge = igraph::V(line)$name, measure = "closeness", value = unname(igraph::closeness(line, normalized = TRUE))))
  # Shortest-path counts on a seven-node graph: rounding of a few divisions.
  expect_reference(hg_edge_centrality(hg, s = 1, measure = c("betweenness", "closeness"), normalized = TRUE), expected,
                   reference = "igraph", covers = "hg_edge_centrality()", quantity = "s = 1 edge centrality",
                   tolerance = 1e-12, key = c("edge", "measure"))
})

test_that("line graphs and edge centralities equal HyperNetX for s = 1, 2, 3 on six random hypergraphs", {
  skip_without_python(c("hypernetx", "networkx", "pandas"))
  script <- c(
    "import sys, json, networkx as nx, pandas as pd, hypernetx as hnx",
    "spec = json.load(open(sys.argv[1]))",
    "d = pd.DataFrame({'edges': spec['edge'], 'nodes': spec['member']})",
    "H = hnx.Hypergraph(d, edge_col='edges', node_col='nodes')",
    "members = {e: set(v) for e, v in H.incidence_dict.items()}",
    "names = sorted(members)",
    "out = {}",
    "for s in (1, 2, 3):",
    "    L = H.get_linegraph(s=s, edges=True)",
    "    L.add_nodes_from(names)",
    "    pairs = sorted(tuple(sorted((str(a), str(b)))) for a, b in L.edges())",
    "    btw = hnx.s_betweenness_centrality(H, s=s, edges=True, normalized=True)",
    "    clo = nx.closeness_centrality(nx.Graph(L))",
    "    out[str(s)] = {'from': [p[0] for p in pairs], 'to': [p[1] for p in pairs],",
    "                   'weight': [len(members[p[0]] & members[p[1]]) for p in pairs],",
    "                   'edge': names, 'betweenness': [float(btw.get(e, 0.0)) for e in names],",
    "                   'closeness': [float(clo[e]) for e in names]}",
    "json.dump(out, open(sys.argv[2], 'w'))")
  version <- ref_python_version("hypernetx")
  # every hypergraph and s stacked in one table, so that an s-line graph without pairs adds no row
  both <- Map(\(config, index) {
    theirs <- ref_python_run(script, as.list(config))
    ours <- group_hypergraph(config, node = "member", hyperedge = "edge")
    lapply(1:3, \(s) {
      reference <- theirs[[as.character(s)]]
      line <- hg_line_graph(ours, s = s, what = "edges")
      list(pairs = list(ours = cbind(hypergraph = rep(index, nrow(line)), s = rep(s, nrow(line)), line),
                        theirs = data.frame(hypergraph = rep(index, length(reference$from)), s = rep(s, length(reference$from)),
                                            from = as.character(reference$from), to = as.character(reference$to),
                                            weight = as.numeric(reference$weight))),
           centrality = list(ours = cbind(hypergraph = index, hg_edge_centrality(ours, s = s, measure = c("betweenness", "closeness"),
                                                                                 normalized = TRUE)),
                             theirs = rbind(data.frame(hypergraph = index, s = s, edge = reference$edge, measure = "betweenness",
                                                       value = reference$betweenness),
                                            data.frame(hypergraph = index, s = s, edge = reference$edge, measure = "closeness",
                                                       value = reference$closeness))))
    })
  }, line_configs(), seq_along(line_configs()))
  stacked <- \(part, side) do.call(rbind, lapply(unlist(both, recursive = FALSE), \(one) one[[part]][[side]]))
  expect_reference(stacked("pairs", "ours"), stacked("pairs", "theirs"), reference = "hypernetx", covers = "hg_line_graph()",
                   quantity = "six hypergraphs, s = 1 to 3, line graph weight", tolerance = 0,
                   key = c("hypergraph", "s", "from", "to"), note = version)
  centrality_ours <- stacked("centrality", "ours")
  centrality_theirs <- stacked("centrality", "theirs")
  key <- c("hypergraph", "s", "edge", "measure")
  # Shortest-path counts on graphs of up to 31 vertices: rounding of a few divisions.
  expect_reference(subset(centrality_ours, measure == "betweenness"), subset(centrality_theirs, measure == "betweenness"),
                   reference = "hypernetx", covers = "hg_edge_centrality()", quantity = "six hypergraphs, s = 1 to 3, s-betweenness",
                   tolerance = 1e-12, key = key, note = version)
  expect_reference(subset(centrality_ours, measure == "closeness"), subset(centrality_theirs, measure == "closeness"),
                   reference = "networkx", covers = "hg_edge_centrality()", quantity = "six hypergraphs, s = 1 to 3, s-closeness",
                   tolerance = 1e-12, key = key,
                   note = paste(ref_python_version("networkx"), "closeness_centrality on HyperNetX's s-line graph"))
})
