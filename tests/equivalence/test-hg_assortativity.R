# hg_assortativity(): degree assortativity of a hypergraph (Chodrow 2020): the weighted Pearson correlation of the
# degrees (or degree ranks) of node pairs within hyperedges, every ordered pair weighted 1 / (m (m - 1))
# ("uniform"), the two largest ("top_2") or the largest and the smallest ("top_bottom"). hg_degree_correlation():
# the correlation of node degrees between hyperedges of two sizes. References: the definitions in base R, igraph's
# Newman assortativity on a graph, XGI's degree_assortativity() and hypergraphx's degree_correlation().

#' A random hypergraph of hyperedges of the given sizes, exact duplicates dropped
assortativity_memberships <- function(seed, n_nodes, n_edges, sizes) {
  set.seed(seed)
  size <- if (length(sizes) == 1L) rep(sizes, n_edges) else sample(sizes, n_edges, replace = TRUE)
  edges <- lapply(size, \(s) sample(n_nodes, s))
  edges <- edges[!duplicated(vapply(edges, \(e) paste(sort(e), collapse = "-"), character(1)))]
  data.frame(member = paste0("n", unlist(edges)), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges))))
}

#' Chodrow's assortativity from its definition: the weighted correlation of the pair scores
by_definition <- function(memberships, scale) {
  degree <- table(memberships$member)
  score <- if (scale == "degree") as.numeric(degree) else rank(as.numeric(degree))
  names(score) <- names(degree)
  groups <- split(memberships$member, memberships$edge)
  pairs <- function(type) do.call(rbind, lapply(Filter(\(g) length(g) >= 2L, groups), \(group) {
    s <- sort(score[group])
    m <- length(group)
    switch(type,
      uniform = {
        ordered <- subset(expand.grid(i = seq_len(m), j = seq_len(m)), i != j)
        data.frame(x = score[group][ordered$i], y = score[group][ordered$j], w = 1 / (m * (m - 1)))
      },
      top_2 = data.frame(x = s[c(m, m - 1L)], y = s[c(m - 1L, m)], w = 0.5),
      top_bottom = data.frame(x = s[[m]], y = s[[1L]], w = 1))
  }))
  weighted_correlation <- function(p) {
    w <- p$w / sum(p$w)
    mx <- sum(w * p$x)
    my <- sum(w * p$y)
    sum(w * (p$x - mx) * (p$y - my)) / sqrt(sum(w * (p$x - mx)^2) * sum(w * (p$y - my)^2))
  }
  types <- c("uniform", "top_2", "top_bottom")
  data.frame(type = types, scale = scale, assortativity = vapply(types, \(type) weighted_correlation(pairs(type)), numeric(1)))
}

configs <- list(mixed_a = assortativity_memberships(31L, 30L, 60L, 2:5), mixed_b = assortativity_memberships(32L, 50L, 120L, 2:6),
                three_uniform = assortativity_memberships(33L, 25L, 70L, 3L), four_uniform = assortativity_memberships(34L, 40L, 90L, 4L))

test_that("every type and scale equals the weighted correlation of Chodrow's definition", {
  invisible(Map(\(memberships, name) {
    hg <- group_hypergraph(memberships, node = "member", hyperedge = "edge")
    expected <- rbind(by_definition(memberships, "degree"), by_definition(memberships, "rank"))
    ours <- rbind(hg_assortativity(hg, type = c("uniform", "top_2", "top_bottom"), scale = "degree"),
                  hg_assortativity(hg, type = c("uniform", "top_2", "top_bottom"), scale = "rank"))
    # Weighted moments of a few hundred small integers.
    expect_reference(ours, expected, reference = "independent: Chodrow (2020) weighted pair correlation, base R",
                     covers = "hg_assortativity()", quantity = paste(name, "assortativity"), tolerance = 1e-12,
                     key = c("type", "scale"))
  }, configs, names(configs)))
})

test_that("on a graph (every hyperedge of size 2) uniform and top-2 assortativity equal igraph's Newman assortativity", {
  skip_without_reference("igraph")
  set.seed(35)
  graph <- igraph::sample_gnm(40, 90)
  ends <- igraph::as_edgelist(graph)
  memberships <- data.frame(member = paste0("v", c(t(ends))), edge = sprintf("e%03d", rep(seq_len(nrow(ends)), each = 2L)))
  newman <- igraph::assortativity_degree(graph, directed = FALSE)
  expect_reference(hg_assortativity(group_hypergraph(memberships, node = "member", hyperedge = "edge"),
                                    type = c("uniform", "top_2"), scale = "degree"),
                   data.frame(type = c("uniform", "top_2"), assortativity = newman), reference = "igraph",
                   covers = "hg_assortativity()", quantity = "graph, degree assortativity", tolerance = 1e-12, key = "type")
})

test_that("top-2 assortativity equals XGI on every hypergraph, uniform assortativity on uniform hypergraphs", {
  skip_without_python("xgi")
  script <- c(
    "import sys, json, xgi",
    "spec = json.load(open(sys.argv[1]))",
    "groups = {}",
    "for n, e in zip(spec['member'], spec['edge']): groups.setdefault(e, []).append(n)",
    "H = xgi.Hypergraph([groups[e] for e in sorted(groups)])",
    "json.dump({'type': ['top_2', 'uniform'],",
    "           'assortativity': [xgi.degree_assortativity(H, kind='top-2', exact=True),",
    "                             xgi.degree_assortativity(H, kind='uniform', exact=True)]}, open(sys.argv[2], 'w'))")
  version <- ref_python_version("xgi")
  invisible(Map(\(memberships, name) {
    theirs <- as.data.frame(ref_python_run(script, as.list(memberships)))
    ours <- hg_assortativity(group_hypergraph(memberships, node = "member", hyperedge = "edge"), type = c("top_2", "uniform"),
                             scale = "degree")
    expect_reference(subset(ours, type == "top_2"), subset(theirs, type == "top_2"), reference = "xgi",
                     covers = "hg_assortativity()", quantity = paste(name, "top-2"), tolerance = 1e-12, key = "type", note = version)
  }, configs, names(configs)))
  uniform <- configs[c("three_uniform", "four_uniform")]
  invisible(Map(\(memberships, name) {
    theirs <- as.data.frame(ref_python_run(script, as.list(memberships)))
    ours <- hg_assortativity(group_hypergraph(memberships, node = "member", hyperedge = "edge"), type = "uniform", scale = "degree")
    expect_reference(ours, subset(theirs, type == "uniform"), reference = "xgi", covers = "hg_assortativity()",
                     quantity = paste(name, "uniform"), tolerance = 1e-12, key = "type", note = version)
  }, uniform, names(uniform)))
  mixed <- configs[c("mixed_a", "mixed_b")]
  invisible(Map(\(memberships, name) {
    theirs <- as.data.frame(ref_python_run(script, as.list(memberships)))
    ours <- hg_assortativity(group_hypergraph(memberships, node = "member", hyperedge = "edge"), type = "uniform", scale = "degree")
    expect_known_difference(ours, subset(theirs, type == "uniform"), reference = "xgi", covers = "hg_assortativity()",
                            quantity = paste(name, "uniform"), tolerance = 1e-12, key = "type", note = version,
                            reason = "On hyperedges of mixed sizes hypergraphs weights each ordered pair 1 / (m (m - 1)), so that every hyperedge counts once (Chodrow 2020); XGI's exact uniform kind weights every pair alike. The two coincide on a uniform hypergraph (compared above).")
  }, mixed, names(mixed)))
})

test_that("degree correlations between hyperedge sizes equal hypergraphx", {
  skip_without_python("hypergraphx")
  script <- c(
    "import sys, json",
    "from hypergraphx import Hypergraph",
    "from hypergraphx.measures.degree import degree_correlation",
    "spec = json.load(open(sys.argv[1]))",
    "groups = {}",
    "for n, e in zip(spec['member'], spec['edge']): groups.setdefault(e, []).append(n)",
    "H = Hypergraph([tuple(sorted(groups[e])) for e in sorted(groups)])",
    "M = degree_correlation(H)",
    "rows = [(i + 2, j + 2, float(M[i, j])) for i in range(M.shape[0]) for j in range(i + 1, M.shape[0])]",
    "json.dump({'size_1': [r[0] for r in rows], 'size_2': [r[1] for r in rows], 'correlation': [r[2] for r in rows]},",
    "          open(sys.argv[2], 'w'))")
  version <- ref_python_version("hypergraphx")
  inputs <- list(assortativity_memberships(41L, 30L, 80L, 2:4), assortativity_memberships(42L, 60L, 150L, 2:5),
                 assortativity_memberships(43L, 20L, 40L, c(2L, 4L)))
  invisible(Map(\(memberships, index) {
    theirs <- as.data.frame(ref_python_run(script, as.list(memberships)))
    ours <- hg_degree_correlation(group_hypergraph(memberships, node = "member", hyperedge = "edge"))
    # hypergraphx returns a square table over the sizes 2 to the largest; the comparison is on the pairs of sizes
    # hypergraphs reports (sizes present in the hypergraph)
    expect_reference(ours, subset(theirs, paste(size_1, size_2) %in% paste(ours$size_1, ours$size_2)), reference = "hypergraphx",
                     covers = "hg_degree_correlation()", quantity = sprintf("hypergraph %d, degree correlation", index),
                     tolerance = 1e-12, key = c("size_1", "size_2"), note = version)
  }, inputs, seq_along(inputs)))
})
