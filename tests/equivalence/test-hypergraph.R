# hypergraph() and its methods: the main constructor routes each input to one constructor. Each route is checked
# on the numbers it produces: a membership table (incidence matrix in base R), a list of sequences (windows
# enumerated in base R), a network as a matrix, sparse Matrix, cograph network, Nestimate network or tna model
# (igraph::cliques()), and a clustering of sequences as a mixture Markov fit, a distance clustering or per-cluster
# networks (the sets of distinct states counted per cluster in base R: frequent-itemset counting restricted to
# exact sets, Agrawal and Srikant 1994).

#' The network every network route reads: 8 nodes, two triangles and a 4-clique sharing nodes
network_matrix <- function() {
  nodes <- letters[1:8]
  adjacency <- matrix(0, 8L, 8L, dimnames = list(nodes, nodes))
  links <- rbind(c("a", "b"), c("a", "c"), c("b", "c"), c("c", "d"), c("d", "e"), c("d", "f"), c("e", "f"),
                 c("e", "g"), c("f", "g"), c("d", "g"), c("g", "h"))
  adjacency[links] <- 1
  adjacency[links[, 2:1]] <- 1
  adjacency
}

#' Every clique of 2 and 3 nodes of an adjacency matrix, by igraph
igraph_cliques <- function(adjacency) {
  graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected")
  sets <- lapply(igraph::cliques(graph, min = 2, max = 3), \(clique) sort(igraph::V(graph)$name[clique]))
  data.frame(members = vapply(sets, paste, character(1), collapse = ", "), present = TRUE)
}

test_that("a membership table gives the incidence matrix's node degrees", {
  memberships <- ref_memberships()
  incidence <- ref_incidence(memberships)
  expect_reference(hg_get(hypergraph(memberships, node = "node", hyperedge = "edge"), what = "nodes"),
                   data.frame(node = rownames(incidence), degree = rowSums(incidence)),
                   reference = "independent: the incidence matrix, base R",
                   covers = c("hypergraph()", "hypergraph.data.frame()", "hg_get.net_hg()"),
                   quantity = "membership table, node degree", tolerance = 0, key = "node")
  # summary() returns the same per-node table, invisibly, after printing
  expect_reference(summary(group_hypergraph(memberships, node = "node", hyperedge = "edge")),
                   data.frame(node = rownames(incidence), degree = rowSums(incidence)),
                   reference = "independent: the incidence matrix, base R", covers = "summary.net_hg()",
                   quantity = "summary, node degree", tolerance = 0, key = "node")
})

test_that("a list of sequences gives the hyperedges of its windows", {
  sequences <- list(c("a", "b", "a", "c", "c", "b"), c("b", "c", "d", "a"))
  starts <- list(1:4, 1:2)
  windows <- unlist(Map(\(events, first) lapply(first, \(start) sort(unique(events[start:(start + 2L)]))), sequences, starts),
                    recursive = FALSE)
  sets <- table(vapply(windows, paste, character(1), collapse = ", "))
  expect_reference(hg_get(hypergraph(sequences, window = 3L, step = 1L), what = "edges"),
                   data.frame(members = names(sets), weight = as.vector(sets)),
                   reference = "independent: windows enumerated over the sequences, base R",
                   covers = c("hypergraph()", "hypergraph.list()"), quantity = "windows of 3, hyperedge weight",
                   tolerance = 0, key = "members")
})

test_that("a network as a matrix, Matrix, cograph network, Nestimate network or tna model gives igraph's cliques", {
  skip_without_reference("igraph")
  adjacency <- network_matrix()
  expected <- igraph_cliques(adjacency)
  routes <- list(
    "hypergraph.matrix()" = adjacency,
    "hypergraph.Matrix()" = Matrix::Matrix(adjacency, sparse = TRUE),
    "hypergraph.cograph_network()" = cograph::as_cograph(adjacency),
    "hypergraph.netobject()" = Nestimate::build_network(adjacency, method = "relative"))
  invisible(Map(\(input, method) {
    ours <- transform(hg_get(hypergraph(input), what = "edges"), present = TRUE)
    expect_reference(ours, expected, reference = "igraph", covers = c("hypergraph()", method),
                     quantity = sprintf("%s, cliques of 2 and 3 nodes", method), tolerance = 0, key = "members")
  }, routes, names(routes)))
})

test_that("a tna model gives igraph's cliques of its network", {
  skip_without_reference(c("igraph", "tna"))
  adjacency <- network_matrix()
  ours <- transform(hg_get(hypergraph(tna::tna(adjacency)), what = "edges"), present = TRUE)
  expect_reference(ours, igraph_cliques(adjacency), reference = "igraph", covers = c("hypergraph()", "hypergraph.tna()"),
                   quantity = "tna model, cliques of 2 and 3 nodes", tolerance = 0, key = "members")
})

test_that("a clustering of sequences gives the exact sets of distinct states counted per cluster", {
  set.seed(3)
  states <- matrix(sample(c("a", "b", "c", "d"), 40L * 5L, replace = TRUE), 40L, 5L)
  states[sample(200L, 20L)] <- NA
  sequences <- as.data.frame(states)
  assignments <- rep(1:2, 20L)
  sets <- vapply(seq_len(nrow(states)), \(row) paste(sort(unique(stats::na.omit(states[row, ]))), collapse = " + "), character(1))
  counted <- as.data.frame(table(cluster = assignments, set = sets), stringsAsFactors = FALSE)
  counted <- subset(counted, Freq > 0)
  totals <- table(assignments)
  expected <- data.frame(group = paste("Cluster", counted$cluster), set = counted$set, count = counted$Freq,
                         share = counted$Freq / as.vector(totals[counted$cluster]),
                         size = lengths(strsplit(counted$set, " + ", fixed = TRUE)))
  fits <- list(
    "hypergraph.net_mmm()" = structure(list(data = sequences, assignments = assignments, k = 2L), class = "net_mmm"),
    "hypergraph.net_clustering()" = structure(list(data = sequences, assignments = assignments, k = 2L), class = "net_clustering"))
  invisible(Map(\(fit, method) {
    expect_reference(hg_get(hypergraph(fit, top = Inf), what = "sets"), expected,
                     reference = "independent: exact sets of distinct states counted per cluster, base R",
                     covers = c("hypergraph()", method, "group_hypergraph()", "hg_get.net_hg()"),
                     quantity = sprintf("%s, set counts", method), tolerance = 1e-12, key = c("group", "set"))
  }, fits, names(fits)))
  grouped <- Nestimate::build_network(transform(sequences, cluster = paste("Cluster", assignments)), method = "relative",
                                      group = "cluster")
  expect_reference(hg_get(hypergraph(grouped, top = Inf), what = "sets"), expected,
                   reference = "independent: exact sets of distinct states counted per cluster, base R",
                   covers = c("hypergraph()", "hypergraph.netobject_group()"), quantity = "netobject_group, set counts",
                   tolerance = 1e-12, key = c("group", "set"))
})

test_that("a membership table gives xgi's node degrees", {
  skip_without_python("xgi")
  memberships <- ref_memberships()
  script <- c(
    "import sys, json, xgi",
    "data = json.load(open(sys.argv[1]))",
    "edges = {}",
    "for n, e in zip(data['node'], data['edge']): edges.setdefault(e, []).append(n)",
    "H = xgi.Hypergraph([edges[e] for e in sorted(edges)])",
    "d = H.nodes.degree.asdict()",
    "json.dump({'node': sorted(d), 'degree': [d[n] for n in sorted(d)]}, open(sys.argv[2], 'w'))")
  theirs <- as.data.frame(ref_python_run(script, as.list(memberships)))
  expect_reference(hg_get(hypergraph(memberships, node = "node", hyperedge = "edge"), what = "nodes"), theirs,
                   reference = "xgi", covers = c("hypergraph()", "hypergraph.data.frame()", "hg_get.net_hg()"),
                   quantity = "membership table, node degree, xgi", tolerance = 0, key = "node", note = ref_python_version("xgi"))
  expect_reference(summary(group_hypergraph(memberships, node = "node", hyperedge = "edge")), theirs,
                   reference = "xgi", covers = "summary.net_hg()", quantity = "summary, node degree, xgi", tolerance = 0,
                   key = "node")
})
