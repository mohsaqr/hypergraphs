# Degree-preserving null-model test: margin preservation, invariances,
# calibration, and detection of blatant structure.

blocky <- c(
  a = "soup salt onion", b = "soup salt pepper", c = "salt onion pepper",
  d = "stars sky moon", e = "stars sky comet", f = "sky moon comet"
)

test_that("checkerboard swaps preserve both margins exactly", {
  hg <- text_hypergraph(blocky)
  m <- (hg$incidence > 0) * 1L
  set.seed(1)
  swapped <- hypergraphs:::.thg_swap_chain(m, attempts = 5000L)
  expect_identical(rowSums(swapped), rowSums(m))
  expect_identical(colSums(swapped), colSums(m))
  expect_true(all(swapped %in% c(0L, 1L)))
  # and the chain actually moves
  expect_gt(sum(swapped != m), 0)
})

test_that("avg_edge_size is invariant under the null, by construction", {
  hg <- text_hypergraph(blocky)
  out <- hg_null_test(hg, statistic = "avg_edge_size", n = 29, seed = 1)
  expect_equal(out$observed, out$null_mean)
  expect_identical(out$p_value, 1)
  expect_true(is.na(out$z))
})

test_that("the test is deterministic given a seed and restores the RNG", {
  hg <- text_hypergraph(blocky)
  set.seed(999)
  before <- get(".Random.seed", envir = globalenv())
  a <- hg_null_test(hg, statistic = "pairwise_participation", n = 29,
                    seed = 42)
  after <- get(".Random.seed", envir = globalenv())
  expect_identical(before, after)
  b <- hg_null_test(hg, statistic = "pairwise_participation", n = 29,
                    seed = 42)
  expect_identical(a, b)
})

test_that("observed values match hg_measures and p stays in bounds", {
  hg <- text_hypergraph(blocky)
  out <- hg_null_test(hg,
                      statistic = c("pairwise_participation", "density"),
                      n = 29, seed = 7)
  summary_tab <- hg_measures(hg, what = "summary")
  expect_equal(
    out$observed[out$statistic == "pairwise_participation"],
    summary_tab$value[summary_tab$measure == "pairwise_participation"]
  )
  expect_true(all(out$p_value >= 1 / 30 & out$p_value <= 1))
  expect_identical(names(out),
                   c("statistic", "observed", "null_mean", "null_lo",
                     "null_hi", "z", "p_value", "n", "method"))
})

test_that("blatant block structure is detected against the null", {
  # two vocabulary blocks that never mix: observed pairwise participation
  # (0.4) sits far below what degree-preserving rewiring produces
  hg <- text_hypergraph(blocky)
  out <- hg_null_test(hg, statistic = "pairwise_participation", n = 99,
                      seed = 11, alternative = "less")
  expect_lt(out$observed, out$null_mean)
  expect_lt(out$p_value, 0.05)
})

test_that("contract violations raise classed or plain errors", {
  hg <- text_hypergraph(blocky)
  expect_error(hg_null_test(42), class = "hypergraphs_bad_input")
  expect_error(hg_null_test(hg, n = 5))
})

test_that("direct null statistics equal the delegated measures path", {
  set.seed(11)
  m <- matrix(rbinom(12 * 30, 1, 0.3), nrow = 12,
              dimnames = list(paste0("d", seq_len(12)),
                              paste0("w", seq_len(30))))
  # drop empty rows/columns so the rebuilt hypergraph matches exactly
  m <- m[rowSums(m) > 0, colSums(m) > 0]
  stats <- c("density", "avg_edge_size", "pairwise_participation",
             "avg_jaccard")
  fast <- .thg_null_statistics(m, stats)
  nz <- which(m > 0, arr.ind = TRUE)
  long <- data.frame(vertex = rownames(m)[nz[, "row"]],
                     edge = colnames(m)[nz[, "col"]], w = 1)
  ref_hg <- group_hypergraph(long, node = "vertex",
                                        hyperedge = "edge", weight = "w")
  s_tab <- hg_measures(ref_hg, what = "summary")
  expect_equal(unname(fast["density"]),
               subset(s_tab, measure == "density")$value,
               tolerance = 1e-12)
  expect_equal(unname(fast["avg_edge_size"]),
               subset(s_tab, measure == "avg_edge_size")$value,
               tolerance = 1e-12)
  expect_equal(unname(fast["pairwise_participation"]),
               subset(s_tab, measure == "pairwise_participation")$value,
               tolerance = 1e-12)
  overlap_tab <- hg_measures(ref_hg, what = "overlap")
  expect_equal(unname(fast["avg_jaccard"]),
               mean(overlap_tab$jaccard),
               tolerance = 1e-12)
})

test_that("the configuration null runs and reports itself", {
  hg <- text_hypergraph(c(
    a = "salt and soup and onions", b = "soup and salt",
    c = "stars and sky and salt", d = "stars and sky"
  ))
  out <- suppressWarnings(hg_null_test(hg,
                      statistic = "pairwise_participation",
                      method = "configuration", n = 49, seed = 1))
  expect_identical(out$method, "configuration")
  expect_identical(out$n, 49L)
  expect_true(out$p_value >= 1 / 50 && out$p_value <= 1)
  # The observed statistic does not depend on the null.
  swap <- hg_null_test(hg, statistic = "pairwise_participation",
                       method = "swap", n = 49, seed = 1)
  expect_equal(out$observed, swap$observed)
})

test_that("the configuration null is reproducible under a seed", {
  hg <- text_hypergraph(c(
    a = "salt and soup and onions", b = "soup and salt",
    c = "stars and sky and salt", d = "stars and sky"
  ))
  args <- list(hg, statistic = "avg_jaccard", method = "configuration",
               n = 29, seed = 42)
  expect_equal(suppressWarnings(do.call(hg_null_test, args)),
               suppressWarnings(do.call(hg_null_test, args)))
})

test_that("INVARIANT: configuration draws never exceed the observed margins", {
  # Stub matching preserves each margin up to collapse, so a draw's degrees
  # and sizes can fall but must never rise, and the total membership count
  # can only fall.
  hg <- text_hypergraph(c(
    a = "salt and soup and onions and night", b = "soup and salt and night",
    c = "stars and sky and salt", d = "stars and sky and night"
  ))
  m <- (hg$incidence > 0) * 1L
  set.seed(7)
  draws <- replicate(30, .thg_configuration_draw(m), simplify = FALSE)
  expect_true(all(vapply(draws, \(d) all(rowSums(d) <= rowSums(m)),
                         logical(1))))
  expect_true(all(vapply(draws, \(d) all(colSums(d) <= colSums(m)),
                         logical(1))))
  expect_true(all(vapply(draws, \(d) sum(d) <= sum(m), logical(1))))
  # and the margins are preserved in expectation, so collapse is rare
})

test_that("collapse in the configuration null is warned about, not hidden", {
  # Few vertices, hyperedges nearly as large as the vertex set: stub matching
  # collides often, so the margins are only approximately preserved.
  hg <- text_hypergraph(c(
    a = "salt and soup and onions and night", b = "soup and salt and night",
    c = "stars and sky and salt", d = "stars and sky and night"
  ))
  expect_warning(
    hg_null_test(hg, statistic = "density", method = "configuration",
                 n = 29, seed = 3),
    class = "hypergraphs_configuration_collapse"
  )
  # The swap null preserves both margins exactly and never warns.
  expect_no_warning(
    hg_null_test(hg, statistic = "density", method = "swap",
                 n = 29, seed = 3)
  )
})

test_that("the configuration null actually randomises", {
  # A draw that is not shuffled would reproduce the observed membership
  # exactly, making every p-value 1 and the null vacuous.
  hg <- text_hypergraph(blocky)
  m <- (hg$incidence > 0) * 1L
  set.seed(5)
  draws <- replicate(20, .thg_configuration_draw(m), simplify = FALSE)
  expect_true(any(vapply(draws, \(d) !identical(d, m), logical(1))))
  spread <- stats::sd(vapply(draws, \(d) sum(tcrossprod(d) > 0), numeric(1)))
  expect_gt(spread, 0)
})

test_that("the configuration null sees block structure, with less power", {
  # The swap null reaches p < 0.05 on this fixture; the configuration null
  # does not. Measured over seeds 11-16 at n = 199: p in [0.090, 0.115],
  # z in [-1.69, -1.49]. Collapse and repeated incidences blur the blocks,
  # so the direction is reliable here but the p-value is not significant --
  # asserting p < 0.05 would be asserting something false.
  hg <- text_hypergraph(blocky)
  out <- suppressWarnings(hg_null_test(
    hg, statistic = "pairwise_participation", method = "configuration",
    n = 199, seed = 11, alternative = "less"
  ))
  expect_lt(out$observed, out$null_mean)
  expect_lt(out$z, -1)
})

# ---- regression tests from the 2026-10-06 audit -------------------------

.null_hg <- function(members) {
  # members: named list of hyperedges -> dense net_hg via group_hypergraph()
  group_hypergraph(
    data.frame(node = unlist(members, use.names = FALSE),
               edge = rep(names(members), lengths(members))),
    node = "node", hyperedge = "edge"
  )
}

test_that("REGRESSION A03: a sole stub is permuted, never multiplied", {
  for (seed in 1:10) {
    set.seed(seed)
    draw <- .thg_configuration_draw(matrix(c(0L, 0L, 1L), 3, 1))
    expect_identical(draw, matrix(c(0L, 0L, 1L), 3, 1), info = seed)
  }
  # zero stubs: nothing to place
  set.seed(1)
  expect_identical(.thg_configuration_draw(matrix(0L, 3, 2)),
                   matrix(0L, 3, 2))
})

test_that("INVARIANT A03: configuration margins never exceed the input", {
  # singleton hyperedges and isolated vertices included
  fixtures <- list(
    matrix(c(0L, 0L, 1L), 3, 1),
    matrix(c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L), 4, 3),
    matrix(c(1L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L), 4, 3),
    matrix(1L, 1, 1)
  )
  for (i in seq_along(fixtures)) {
    m <- fixtures[[i]]
    for (seed in 1:20) {
      set.seed(seed)
      draw <- .thg_configuration_draw(m)
      expect_true(all(rowSums(draw) <= rowSums(m)), info = paste(i, seed))
      expect_true(all(colSums(draw) <= colSums(m)), info = paste(i, seed))
      expect_true(all(draw %in% c(0L, 1L)), info = paste(i, seed))
    }
  }
})

test_that("REGRESSION A08: the sample count is one whole number >= 19", {
  path <- .null_hg(list(e1 = c("a", "b"), e2 = c("b", "c")))
  for (bad in list(19.5, 5, NA_real_, Inf, -Inf, 3e9, c(19, 20), "19")) {
    expect_error(hg_null_test(path, statistic = "avg_edge_size", n = bad,
                              seed = 1),
                 class = "hypergraphs_bad_input", info = format(bad))
  }
})

test_that("INVARIANT A08: an invariant statistic has p = 1 under every null", {
  path <- .null_hg(list(e1 = c("a", "b"), e2 = c("b", "c"),
                        e3 = c("c", "d", "a")))
  for (method in c("swap", "assignment")) {
    out <- hg_null_test(path, statistic = c("avg_edge_size", "density"),
                        method = method, n = 19, seed = 1)
    expect_identical(out$p_value, c(1, 1), info = method)
    expect_identical(out$n, c(19L, 19L), info = method)
    expect_equal(out$observed, out$null_mean, info = method)
  }
})

test_that("REGRESSION A09: null density is the hg_measures() density", {
  fixtures <- list(
    triangle = list(e1 = c("a", "b"), e2 = c("a", "c"), e3 = c("b", "c")),
    uniform_3 = list(e1 = c("a", "b", "c"), e2 = c("b", "c", "d"),
                     e3 = c("a", "c", "d"), e4 = c("a", "b", "d")),
    mixed = list(e1 = c("a", "b", "c"), e2 = c("c", "d"), e3 = c("d", "e")),
    duplicates = list(e1 = c("a", "b"), e2 = c("a", "b"), e3 = c("b", "c"))
  )
  for (nm in names(fixtures)) {
    hg <- .null_hg(fixtures[[nm]])
    out <- hg_null_test(hg, statistic = c("density", "avg_edge_size",
                                          "pairwise_participation"),
                        n = 19, seed = 2)
    summary_tab <- hg_measures(hg, what = "summary")
    expect_equal(out$observed,
                 summary_tab$value[match(out$statistic,
                                         summary_tab$measure)],
                 tolerance = 1e-12, info = nm)
  }
  # the audit's case: a simple triangle is complete, density 1
  triangle <- .null_hg(fixtures$triangle)
  expect_identical(hg_null_test(triangle, statistic = "density", n = 19,
                                seed = 1)$observed, 1)
  # an empty hypergraph has density 0 on both paths
  empty <- hg_subset(triangle, edges = character(0), drop_isolated = FALSE)
  expect_identical(
    unname(.thg_null_statistics((empty$incidence > 0) * 1L, "density")),
    hg_measures(empty, what = "summary")$value[3L]
  )
})

test_that("REGRESSION A10: undefined statistics are NA with a classed warning", {
  one <- .null_hg(list(e1 = c("a", "b")))
  expect_warning(
    out <- hg_null_test(one, n = 19, seed = 1),
    class = "hypergraphs_undefined_statistic"
  )
  jac <- out[out$statistic == "avg_jaccard", ]
  expect_true(is.na(jac$observed) && is.na(jac$p_value) && is.na(jac$null_lo))
  expect_identical(jac$n, 19L)
  # the defined statistics are still tested
  expect_false(anyNA(out$p_value[out$statistic != "avg_jaccard"]))
  # a well-defined statistic requested on its own raises nothing
  expect_no_warning(hg_null_test(one, statistic = "density", n = 19,
                                 seed = 1))
})

test_that("REGRESSION A10: degenerate hypergraphs never reach quantile()", {
  base <- .null_hg(list(e1 = c("a", "b"), e2 = c("b", "c")))
  zero_edges <- hg_subset(base, edges = character(0), drop_isolated = FALSE)
  one_node <- .null_hg(list(e1 = "a"))
  all_empty <- hypergraphs:::.thg_from_incidence(
    matrix(0, 2, 2, dimnames = list(c("a", "b"), c("e1", "e2"))), list()
  )
  for (hg in list(zero_edges, one_node, all_empty)) {
    for (method in c("swap", "configuration", "assignment")) {
      out <- withCallingHandlers(
        hg_null_test(hg, method = method, n = 19, seed = 1),
        hypergraphs_undefined_statistic = \(w) invokeRestart("muffleWarning")
      )
      expect_identical(nrow(out), 6L)
      defined <- !is.na(out$observed)
      expect_true(all(out$p_value[defined] >= 1 / 20 &
                        out$p_value[defined] <= 1))
    }
  }
  # zero/one-node pairwise participation is 0, as hg_measures() reports
  out <- hg_null_test(one_node, statistic = "pairwise_participation",
                      n = 19, seed = 1)
  expect_identical(out$observed, 0)
})
