# random_hypergraph(): samples from four random models with draws made inside the function (a seed only). HyperG's
# generators share no stream, so the comparison is statistical and by properties, declared as such:
# - properties of the definition that hold for any draw (?random_hypergraph): a "uniform" hyperedge has exactly k
#   distinct members, a "regular" node lies in exactly k hyperedges, an "sbm" hyperedge with a block-diagonal P and
#   no impurity lies inside one block;
# - Monte Carlo summaries against the model's expectation: the share of incidences of "gnp" is p, each node of a
#   "uniform" model with equal weights lies in a share k / n of the hyperedges, the number of "sbm" hyperedges has
#   mean sum_b choose(n_b, 2) P_bb. Tolerances are 4 binomial standard errors, fixed before the run; a sweep of 20
#   seeds before fixing them passed 20 of 20 for every check.

rh_seeds <- c(3L, 4L, 5L)

test_that("every uniform hyperedge has k members and every regular node k hyperedges, for any draw", {
  invisible(lapply(rh_seeds, \(seed) {
    uniform <- hg_get(random_hypergraph("uniform", n = 25, m = 40, k = 4, seed = seed), what = "edges")
    expect_reference(data.frame(hyperedge = uniform$hyperedge, size = uniform$size),
                     data.frame(hyperedge = uniform$hyperedge, size = 4), reference = "independent: the model's definition",
                     covers = "random_hypergraph()", quantity = sprintf("uniform hyperedge size, seed %d", seed),
                     tolerance = 0, key = "hyperedge")
    regular <- hg_get(random_hypergraph("regular", n = 30, m = 12, k = 3, seed = seed), what = "nodes")
    expect_reference(data.frame(node = regular$node, degree = regular$degree), data.frame(node = regular$node, degree = 3),
                     reference = "independent: the model's definition", covers = "random_hypergraph()",
                     quantity = sprintf("regular node degree, seed %d", seed), tolerance = 0, key = "node")
  }))
})

test_that("an SBM with a block-diagonal P keeps every hyperedge inside one block, for any draw", {
  invisible(lapply(rh_seeds, \(seed) {
    hg <- random_hypergraph("sbm", P = diag(c(0.3, 0.3)), block_sizes = c(20L, 20L), d = 3, seed = seed)
    block <- hg_get(hg, what = "nodes")
    memberships <- hg_get(hg, what = "memberships")
    blocks_per_edge <- tapply(block$block[match(memberships$node, block$node)], memberships$hyperedge,
                              \(b) length(unique(b)))
    expect_reference(data.frame(hyperedge = names(blocks_per_edge), blocks = as.vector(blocks_per_edge)),
                     data.frame(hyperedge = names(blocks_per_edge), blocks = 1), reference = "independent: the model's definition",
                     covers = "random_hypergraph()", quantity = sprintf("blocks spanned by a hyperedge, seed %d", seed),
                     tolerance = 0, key = "hyperedge")
  }))
})

test_that("incidence shares and hyperedge counts follow the models' expectations (Monte Carlo summary)", {
  invisible(lapply(rh_seeds, \(seed) {
    gnp <- hg_get(random_hypergraph("gnp", n = 200, m = 300, p = 0.1, seed = seed), what = "memberships")
    expect_reference(nrow(gnp) / (200 * 300), 0.1, reference = "independent: Bernoulli(p) incidences",
                     covers = "random_hypergraph()", quantity = sprintf("Monte Carlo summary: gnp incidence share, seed %d", seed),
                     tolerance = 4 * sqrt(0.1 * 0.9 / (200 * 300)))
    uniform <- hg_get(random_hypergraph("uniform", n = 10, m = 2000, k = 3, seed = seed), what = "nodes")
    expect_reference(data.frame(node = uniform$node, share = uniform$degree / 2000),
                     data.frame(node = uniform$node, share = 3 / 10), reference = "independent: k of n nodes drawn equally",
                     covers = "random_hypergraph()", quantity = sprintf("Monte Carlo summary: uniform node share, seed %d", seed),
                     tolerance = 4 * sqrt(0.3 * 0.7 / 2000), key = "node")
    sbm <- hg_get(random_hypergraph("sbm", P = diag(c(0.3, 0.3)), block_sizes = c(20L, 20L), d = 3, seed = seed),
                  what = "edges")
    pairs <- 2 * choose(20, 2)
    expect_reference(nrow(sbm), pairs * 0.3, reference = "independent: Bernoulli(P) pairs of the block model",
                     covers = "random_hypergraph()", quantity = sprintf("Monte Carlo summary: SBM hyperedge count, seed %d", seed),
                     tolerance = 4 * sqrt(pairs * 0.3 * 0.7))
  }))
})
