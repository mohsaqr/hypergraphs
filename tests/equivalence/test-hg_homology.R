# hg_homology(): persistent homology of the Vietoris-Rips filtration of a point cloud's distances; hg_landscape():
# its persistence landscapes (Bubenik 2015), lambda_k(t) the k-th largest of max(0, min(t - b, d - t)) over the
# pairs; hg_bottleneck() and hg_wasserstein(): distances between diagrams, a point matched to the diagonal at
# (d - b) / 2 in the sup norm. References: TDAstats (ripser) and TDA (Dionysus: ripsDiag, bottleneck,
# wasserstein, landscape), SciPy's linear_sum_assignment for the order-2 Wasserstein matching, the landscape's
# definition in base R. TDAstats reports no essential class (death beyond the threshold): the essential classes
# hypergraphs reports with death = Inf are left out of that comparison.

#' Finite pairs of a diagram up to `top`, ranked within each dimension by birth then death
ranked <- function(diagram, top = 1L) {
  diagram <- subset(diagram, is.finite(death) & dimension <= top, select = c(dimension, birth, death))
  diagram <- diagram[order(diagram$dimension, diagram$birth, diagram$death), ]
  diagram$rank <- stats::ave(diagram$birth, diagram$dimension, FUN = seq_along)
  diagram
}

clouds <- lapply(1:4, \(seed) {
  set.seed(seed)
  matrix(stats::rnorm(c(8L, 12L, 15L, 12L)[[seed]] * 2L), ncol = 2L)
})

test_that("finite persistence pairs in dimensions 0 and 1 equal TDAstats (ripser)", {
  skip_without_reference("TDAstats")
  invisible(Map(\(points, index) {
    distances <- as.matrix(stats::dist(points))
    threshold <- max(distances) * 1.05
    theirs <- as.data.frame(TDAstats::calculate_homology(points, dim = 1L, threshold = threshold))
    ours <- hg_homology(distances, n_steps = 5L, max_dim = 2L, type = "vr", max_scale = threshold)
    # Births and deaths are pairwise distances, read on both sides: rounding only.
    expect_reference(ranked(hg_get(ours, what = "persistence")), ranked(theirs), reference = "TDAstats",
                     covers = c("hg_homology()", "hg_get.persistent_homology()", "hg_get.hypergraphs_result()"),
                     quantity = sprintf("cloud %d, finite pairs", index), tolerance = 1e-10, key = c("dimension", "rank"))
  }, clouds, seq_along(clouds)))
})

test_that("finite persistence pairs equal TDA's Rips diagram (Dionysus)", {
  skip_without_reference("TDA")
  invisible(Map(\(points, index) {
    distances <- as.matrix(stats::dist(points))
    threshold <- max(distances) * 1.05
    diagram <- TDA::ripsDiag(points, maxdimension = 1L, maxscale = threshold, library = "Dionysus")$diagram
    theirs <- data.frame(dimension = diagram[, "dimension"], birth = diagram[, "Birth"], death = diagram[, "Death"])
    # TDA writes an essential class with death = maxscale: those are left out like the infinite ones
    theirs <- subset(theirs, death < threshold)
    ours <- hg_homology(distances, n_steps = 5L, max_dim = 2L, type = "vr", max_scale = threshold)
    expect_reference(ranked(hg_get(ours, what = "persistence")), ranked(theirs), reference = "TDA",
                     covers = c("hg_homology()", "hg_get.persistent_homology()"),
                     quantity = sprintf("cloud %d, finite pairs", index), tolerance = 1e-10, key = c("dimension", "rank"))
  }, clouds, seq_along(clouds)))
})

test_that("persistence landscapes equal their definition and TDA's landscape()", {
  points <- clouds[[3L]]
  distances <- as.matrix(stats::dist(points))
  homology <- hg_homology(distances, max_dim = 2L, type = "vr", max_scale = max(distances) * 1.05)
  landscape <- hg_landscape(homology, k_max = 3L, dimension = 1L)
  values <- hg_get(landscape)
  pairs <- subset(hg_get(homology, what = "persistence"), dimension == 1L & is.finite(death))
  grid <- sort(unique(values$t))
  # Rows are keyed by their place on the grid, not by t: t is a computed double, and a key label printed from it
  # changes with its last bits on another machine (a re-run then cannot match the row). t is compared as a value.
  values$point <- match(values$t, grid)
  tents <- vapply(grid, \(t) sort(pmax(0, pmin(t - pairs$birth, pairs$death - t)), decreasing = TRUE)[1:3], numeric(3))
  expected <- data.frame(k = rep(1:3, times = length(grid)), t = rep(grid, each = 3L), value = as.vector(tents),
                         point = rep(seq_along(grid), each = 3L))
  expected$value[is.na(expected$value)] <- 0
  expect_reference(values, expected, reference = "independent: lambda_k(t) of Bubenik (2015), base R",
                   covers = c("hg_landscape()", "hg_get.persistence_landscape()"), quantity = "dimension 1, k = 1 to 3",
                   tolerance = 1e-12, key = c("k", "point"))
  skip_without_reference("TDA")
  diagram <- cbind(dimension = pairs$dimension, Birth = pairs$birth, Death = pairs$death)
  theirs <- TDA::landscape(diagram, dimension = 1L, KK = 1:3, tseq = grid)
  expect_reference(values, data.frame(k = rep(1:3, each = length(grid)), t = rep(grid, times = 3L), value = as.vector(theirs),
                                      point = rep(seq_along(grid), times = 3L)),
                   reference = "TDA", covers = c("hg_landscape()", "hg_get.persistence_landscape()"),
                   quantity = "dimension 1, k = 1 to 3", tolerance = 1e-12, key = c("k", "point"))
})

#' Two random diagrams of one dimension, all pairs finite
diagrams <- function(seed) {
  set.seed(seed)
  list(data.frame(dimension = 0L, birth = stats::runif(12, 0, 2), death = stats::runif(12, 3, 6)),
       data.frame(dimension = 0L, birth = stats::runif(9, 0, 2), death = stats::runif(9, 3, 6)))
}

test_that("bottleneck and order-1 Wasserstein distances equal TDA's (Dionysus)", {
  skip_without_reference("TDA")
  invisible(lapply(c(20260902L, 11L, 12L), \(seed) {
    pair <- diagrams(seed)
    as_tda <- \(d) cbind(dimension = d$dimension, Birth = d$birth, Death = d$death)
    # Sup-norm costs and one optimal matching on each side.
    expect_reference(unname(hg_bottleneck(pair[[1L]], pair[[2L]], dimension = 0L)),
                     TDA::bottleneck(as_tda(pair[[1L]]), as_tda(pair[[2L]]), dimension = 0L), reference = "TDA",
                     covers = "hg_bottleneck()", quantity = sprintf("diagrams %d, bottleneck", seed), tolerance = 1e-10)
    expect_reference(unname(hg_wasserstein(pair[[1L]], pair[[2L]], dimension = 0L, order = 1)),
                     TDA::wasserstein(as_tda(pair[[1L]]), as_tda(pair[[2L]]), p = 1, dimension = 0L), reference = "TDA",
                     covers = "hg_wasserstein()", quantity = sprintf("diagrams %d, Wasserstein order 1", seed), tolerance = 1e-10)
  }))
})

test_that("the order-2 Wasserstein distance equals SciPy's optimal assignment", {
  skip_without_python(c("scipy", "numpy"))
  script <- c(
    "import sys, json, numpy as np",
    "from scipy.optimize import linear_sum_assignment",
    "spec = json.load(open(sys.argv[1]))",
    "a = np.column_stack([spec['b1'], spec['d1']]); b = np.column_stack([spec['b2'], spec['d2']])",
    "n, m = len(a), len(b)",
    "pair = np.max(np.abs(a[:, None, :] - b[None, :, :]), axis=2)",
    "da = np.abs(a[:, 1] - a[:, 0]) / 2; db = np.abs(b[:, 1] - b[:, 0]) / 2",
    "cost = np.zeros((n + m, n + m))",
    "cost[:n, :m] = pair ** 2",
    "cost[:n, m:] = np.repeat((da ** 2)[:, None], n, axis=1)",
    "cost[n:, :m] = np.repeat((db ** 2)[None, :], m, axis=0)",
    "rows, cols = linear_sum_assignment(cost)",
    "json.dump({'distance': float(np.sqrt(cost[rows, cols].sum()))}, open(sys.argv[2], 'w'))")
  version <- ref_python_version("scipy")
  invisible(lapply(c(20260902L, 11L, 12L), \(seed) {
    pair <- diagrams(seed)
    theirs <- ref_python_run(script, list(b1 = pair[[1L]]$birth, d1 = pair[[1L]]$death, b2 = pair[[2L]]$birth, d2 = pair[[2L]]$death))
    expect_reference(unname(hg_wasserstein(pair[[1L]], pair[[2L]], order = 2)), theirs$distance, reference = "scipy",
                     covers = "hg_wasserstein()", quantity = sprintf("diagrams %d, Wasserstein order 2", seed), tolerance = 1e-12,
                     note = paste(version, "scipy.optimize.linear_sum_assignment on the augmented cost matrix"))
  }))
})
