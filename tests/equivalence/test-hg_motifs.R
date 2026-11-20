# hg_motifs(): the census of the Y, T and O motifs of three 3-node hyperedges on four nodes (hyperedges sharing
# two nodes: two, three or four of the triples on four nodes). The observed counts against hypergraphx's
# compute_motifs(order = 4) on six random 3-uniform hypergraphs (each with one hyperedge repeated). The null model
# (configuration-model draws inside the function) is not compared.

#' Six random 3-uniform hypergraphs, each with its first hyperedge repeated
motif_configs <- function() {
  one <- function(seed, n_nodes, n_edges) {
    set.seed(seed)
    triples <- lapply(seq_len(n_edges), \(j) sample(n_nodes, 3L))
    triples <- c(triples, triples[1L])
    data.frame(member = unlist(triples), edge = rep(seq_along(triples), each = 3L))
  }
  list(one(11L, 7L, 14L), one(12L, 8L, 25L), one(13L, 10L, 30L), one(14L, 6L, 12L), one(15L, 12L, 20L), one(16L, 9L, 40L))
}

test_that("observed Y, T and O counts equal hypergraphx's motif census", {
  skip_without_python("hypergraphx")
  script <- c(
    "import sys, json",
    "import hypergraphx as hgx",
    "from hypergraphx.motifs import compute_motifs",
    "spec = json.load(open(sys.argv[1]))",
    "groups = {}",
    "for n, e in zip(spec['member'], spec['edge']): groups.setdefault(e, []).append(int(n))",
    "edges = [tuple(sorted(groups[e])) for e in sorted(groups)]",
    "H = hgx.Hypergraph(edges)",
    "observed = dict(compute_motifs(H, order=4, runs_config_model=0)['observed'])",
    "keys = {'Y': ((1, 2, 3), (1, 2, 4)), 'T': ((1, 2, 3), (1, 2, 4), (1, 3, 4)),",
    "        'O': ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4))}",
    "other = sum(v for k, v in observed.items() if k not in keys.values())",
    "json.dump({'motif': list(keys) + ['other'], 'count': [float(observed.get(v, 0)) for v in keys.values()] + [float(other)]},",
    "          open(sys.argv[2], 'w'))")
  version <- ref_python_version("hypergraphx")
  invisible(Map(\(config, index) {
    theirs <- as.data.frame(ref_python_run(script, as.list(config)))
    found <- hg_motifs(group_hypergraph(transform(config, member = paste0("n", member), edge = sprintf("e%03d", edge)),
                                        node = "member", hyperedge = "edge"), n = 20L, seed = 1L)
    # a census of three-edge motifs on four nodes counts only Y, T and O: hypergraphx's count of any other
    # motif of order 4 made of 3-node hyperedges is compared with zero
    counts <- rbind(subset(hg_get(found, what = "test"), select = c(motif, count)), data.frame(motif = "other", count = 0))
    expect_reference(counts, theirs, reference = "hypergraphx", covers = c("hg_motifs()", "hg_get.hypergraphs_motifs()"),
                     quantity = sprintf("hypergraph %d, observed motifs", index), tolerance = 0, key = "motif", note = version)
    expect_reference(subset(hg_get(summary(found), what = "test"), select = c(motif, count)), subset(theirs, motif != "other"),
                     reference = "hypergraphx", covers = c("summary.hypergraphs_motifs()", "hg_get.hypergraphs_summary()"),
                     quantity = sprintf("hypergraph %d, observed motifs through summary()", index), tolerance = 0, key = "motif",
                     note = version)
  }, motif_configs(), seq_along(motif_configs())))
})
