# hg_modularity(): the hypergraph modularity of Kaminski et al. (2019) of a node partition, with the hyperedge
# weighting "linear", "majority" or "strict". Reference: HyperNetX's modularity() (adapted from Theberge's
# reference code) on three planted hypergraphs, five partitions each, unit and integer hyperedge weights.

#' A planted hypergraph: three blocks of eight nodes, ten hyperedges inside each block and four across
planted <- function(seed) {
  set.seed(seed)
  nodes <- sprintf("v%02d", seq_len(24L))
  block <- rep(1:3, each = 8L)
  within <- unlist(lapply(1:3, \(b) lapply(seq_len(10L), \(i) sample(nodes[block == b], sample(2:5, 1L)))), recursive = FALSE)
  members <- c(within, lapply(seq_len(4L), \(i) sample(nodes, sample(2:5, 1L))))
  edges <- sprintf("e%03d", seq_along(members))
  set.seed(seed + 1L)
  partitions <- list(truth = block, random3 = sample(3L, 24L, replace = TRUE), random5 = sample(5L, 24L, replace = TRUE),
                     one = rep(1L, 24L), merged = pmin(block, 2L))
  list(memberships = data.frame(edge = rep(edges, lengths(members)), node = unlist(members)),
       weights = data.frame(edge = edges, unit = 1, integer = sample(1:4, length(edges), replace = TRUE)),
       nodes = nodes, partitions = partitions)
}

script <- c(
  "import sys, json",
  "import hypernetx as hnx",
  "import hypernetx.algorithms as A",
  "spec = json.load(open(sys.argv[1]))",
  "members = {}",
  "for n, e in zip(spec['node'], spec['edge']): members.setdefault(e, []).append(n)",
  "rows = {'weights': [], 'partition': [], 'type': [], 'modularity': []}",
  "for scheme in ('unit', 'integer'):",
  "    HG = hnx.Hypergraph(dict(members))",
  "    for e, x in zip(spec['weight_edge'], spec[scheme]): HG.edges[e].weight = int(x)",
  "    for name, labels in spec['partitions'].items():",
  "        groups = {}",
  "        for n, c in zip(spec['nodes'], labels): groups.setdefault(c, set()).add(n)",
  "        for kind, f in (('linear', A.linear), ('majority', A.majority), ('strict', A.strict)):",
  "            rows['weights'].append(scheme); rows['partition'].append(name); rows['type'].append(kind)",
  "            rows['modularity'].append(float(A.modularity(HG, list(groups.values()), wdc=f)))",
  "json.dump(rows, open(sys.argv[2], 'w'))")

test_that("the modularity of five partitions equals HyperNetX on three planted hypergraphs", {
  skip_without_python("hypernetx")
  version <- ref_python_version("hypernetx")
  invisible(lapply(c(11L, 23L, 37L), \(seed) {
    fixture <- planted(seed)
    theirs <- as.data.frame(ref_python_run(script, list(node = fixture$memberships$node, edge = fixture$memberships$edge,
                                                        weight_edge = fixture$weights$edge, unit = fixture$weights$unit,
                                                        integer = fixture$weights$integer, nodes = fixture$nodes,
                                                        partitions = fixture$partitions)))
    hg <- group_hypergraph(fixture$memberships, node = "node", hyperedge = "edge")
    edge_order <- hg_get(hg, what = "edges")$hyperedge
    grid <- expand.grid(weights = c("unit", "integer"), partition = names(fixture$partitions),
                        type = c("linear", "majority", "strict"), stringsAsFactors = FALSE)
    grid$modularity <- mapply(\(weights, partition, type) {
      score <- hg_modularity(hg, partition = stats::setNames(fixture$partitions[[partition]], fixture$nodes), type = type,
                             edge_weights = fixture$weights[[weights]][match(edge_order, fixture$weights$edge)])
      score$modularity
    }, grid$weights, grid$partition, grid$type)
    # Sums of a few hundred rational terms.
    expect_reference(grid, theirs, reference = "hypernetx", covers = "hg_modularity()",
                     quantity = sprintf("planted hypergraph %d", seed), tolerance = 1e-12,
                     key = c("weights", "partition", "type"), note = version)
  }))
})
