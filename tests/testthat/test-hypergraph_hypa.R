skip_on_cran()

# ---- fixtures ---------------------------------------------------------------

planted_hg <- function(seed = 1) {
  set.seed(seed)
  bg <- data.frame(actor = sample(letters[1:12], 400, TRUE),
                   group = paste0("g", sample(60, 400, TRUE)))
  planted <- data.frame(actor = rep(c("x", "y"), 25),
                        group = rep(paste0("p", 1:25), each = 2))
  group_hypergraph(rbind(bg, planted), node = "actor", hyperedge = "group")
}

# Two hubs that each join 60 background hyperedges but never the same one:
# their strengths imply about a dozen shared hyperedges, they share none.
avoiding_hubs_hg <- function(seed = 1) {
  set.seed(seed)
  background <- sprintf("b%02d", 1:20)
  edges <- lapply(seq_len(300), \(e) sample(background, 3))
  edges[1:60] <- lapply(edges[1:60], \(e) c(e, "h1"))
  edges[61:120] <- lapply(edges[61:120], \(e) c(e, "h2"))
  group_hypergraph(
    data.frame(node = unlist(edges),
               hyperedge = rep(sprintf("e%03d", seq_along(edges)),
                               lengths(edges))),
    node = "node", hyperedge = "hyperedge")
}

# Random incidence: node i joins hyperedge e independently with probability
# p_i q_e, so E[c_ij] is proportional to p_i p_j, the product form of the
# null, with no pair-level structure.
random_incidence_hg <- function(seed, n = 40L, n_edges = 1000L) {
  set.seed(seed)
  p <- stats::runif(n, 0.2, 1)
  q <- stats::runif(n_edges, 0.03, 0.12)
  inc <- matrix(stats::rbinom(n * n_edges, 1L, outer(p, q)), n, n_edges)
  member <- which(inc == 1L, arr.ind = TRUE)
  group_hypergraph(
    data.frame(node = sprintf("n%02d", member[, "row"]),
               hyperedge = sprintf("e%04d", member[, "col"])),
    node = "node", hyperedge = "hyperedge")
}

star_hg <- function(leaves = 8L) {
  group_hypergraph(
    data.frame(node = c(rbind("hub", paste0("l", seq_len(leaves)))),
               hyperedge = rep(paste0("e", seq_len(leaves)), each = 2L)),
    node = "node", hyperedge = "hyperedge")
}

pair_key <- function(d) {
  paste(pmin(d$from, d$to), pmax(d$from, d$to), sep = "|")
}

# ---- the null model ------------------------------------------------------

test_that("node strengths equal the row sums of the co-occurrence matrix", {
  hg <- random_incidence_hg(1L)
  co <- .hg_hypa_cooccurrence(hg$incidence)
  b <- (as.matrix(hg$incidence) != 0) * 1
  dense <- tcrossprod(b)
  diag(dense) <- 0
  expect_equal(co$strength, unname(rowSums(dense)))
  expect_equal(sum(co$pairs$co), sum(dense) / 2)
})

test_that("expected counts preserve every node's strength", {
  strengths <- lapply(1:5, \(seed) {
    co <- .hg_hypa_cooccurrence(random_incidence_hg(seed, n = 30L,
                                                    n_edges = 300L)$incidence)
    co$strength[co$strength > 0]
  })
  # a planted fixture with a heavy, isolated pair as well
  co <- .hg_hypa_cooccurrence(planted_hg()$incidence)
  strengths <- c(strengths, list(co$strength[co$strength > 0]))
  invisible(lapply(strengths, \(s) {
    fit <- .hg_hypa_fit(s)
    expect_true(fit$converged)
    expected <- outer(fit$x, fit$x)
    diag(expected) <- 0
    expect_equal(rowSums(expected), s, tolerance = 1e-8)
    expect_equal(sum(expected) / 2, fit$m, tolerance = 1e-8)
  }))
})

test_that("the fitted null equals symmetric matrix scaling (independent oracle)", {
  co <- .hg_hypa_cooccurrence(random_incidence_hg(2L, n = 25L,
                                                  n_edges = 300L)$incidence)
  s <- co$strength[co$strength > 0]
  fit <- .hg_hypa_fit(s)
  # Symmetric Sinkhorn on the all-ones matrix with a zero diagonal: the
  # strength-preserving symmetric scaling is unique, so both must agree.
  start <- matrix(1, length(s), length(s))
  diag(start) <- 0
  oracle <- Reduce(\(e, i) {
    d <- sqrt(s / rowSums(e))
    e * outer(d, d)
  }, seq_len(3000L), start)
  fitted <- outer(fit$x, fit$x)
  diag(fitted) <- 0
  expect_equal(fitted, oracle, tolerance = 1e-6)
})

test_that("a star has no fitted null and warns by class", {
  expect_warning(out <- hypa(star_hg(), min_count = 1L),
                 class = "hypergraphs_no_converge")
  fit <- attr(out, "hypa_fit")
  expect_false(fit$converged)
  expect_identical(fit$iterations, 10000L)
  expect_gt(fit$max_error, 1e-10 * fit$n_cooccurrences)
})

test_that("a converged fit is recorded on the result", {
  out <- hypa(planted_hg(), min_count = 3L)
  fit <- attr(out, "hypa_fit")
  expect_true(fit$converged)
  expect_lt(fit$max_error, 1e-10 * fit$n_cooccurrences)
  expect_lt(fit$iterations, 10000L)
})

# ---- planted anomalies ---------------------------------------------------

test_that("hypa recovers a planted over-represented pair", {
  out <- hypa(planted_hg(), min_count = 3L)
  hit <- subset(out, from == "x" & to == "y")
  expect_identical(nrow(hit), 1L)
  expect_identical(hit$anomaly, "over")
  expect_gt(hit$ratio, 2)
  # the most over-represented pair comes first
  expect_identical(pair_key(out)[1L], "x|y")
  # and it does not cry wolf on the random background
  background <- subset(out, !(from == "x" & to == "y"))
  expect_gt(nrow(background), 0L)
  expect_true(all(background$anomaly == "normal"))
})

test_that("two high-strength nodes that never co-occur are under-represented", {
  out <- hypa(avoiding_hubs_hg(), min_count = 5L)
  hubs <- subset(out, from %in% c("h1", "h2") & to %in% c("h1", "h2"))
  expect_identical(nrow(hubs), 1L)
  expect_identical(hubs$observed, 0L)
  expect_gt(hubs$expected, 5)
  expect_lt(hubs$p_under, 1e-4)
  expect_identical(hubs$anomaly, "under")
})

# ---- tails, classification, scored set ------------------------------------

test_that("the tails are the paper's: P(X <= f) and the strict P(X > f)", {
  out <- hypa(planted_hg(), min_count = 2L)
  expect_equal(out$p_under + out$p_over, rep(1, nrow(out)), tolerance = 1e-12)
  m <- attr(out, "hypa_fit")$n_cooccurrences
  k <- round(m * out$expected)
  expect_equal(out$p_over,
               stats::phyper(out$observed, k, m^2 - k, m, lower.tail = FALSE))
  expect_equal(out$p_under, stats::phyper(out$observed, k, m^2 - k, m))
})

test_that("classification applies the under rule first, each with its floor", {
  out <- hypa(avoiding_hubs_hg(), min_count = 5L)
  under <- out$expected >= 5 & out$p_adjusted_under < 0.05
  over <- !under & out$observed >= 5 & out$p_adjusted_over < 0.05
  expect_identical(out$anomaly,
                   ifelse(under, "under", ifelse(over, "over", "normal")))
  expect_true(any(out$anomaly == "under"))
})

test_that("the scored set is exactly the pairs that can be classified", {
  hg <- avoiding_hubs_hg()
  out <- hypa(hg, min_count = 4L)
  # brute force over every possible pair of nodes with positive strength
  co <- .hg_hypa_cooccurrence(hg$incidence)
  active <- which(co$strength > 0)
  fit <- .hg_hypa_fit(co$strength[active])
  b <- (as.matrix(hg$incidence) != 0) * 1
  dense <- tcrossprod(b)[active, active]
  idx <- which(upper.tri(dense), arr.ind = TRUE)
  all_pairs <- data.frame(from = hg$nodes[active[idx[, 1]]],
                          to = hg$nodes[active[idx[, 2]]],
                          observed = dense[idx],
                          expected = fit$x[idx[, 1]] * fit$x[idx[, 2]])
  eligible <- subset(all_pairs, observed >= 4 | expected >= 4)
  expect_gt(sum(eligible$observed < 4), 0L)  # some are under-eligible only
  expect_setequal(pair_key(out), pair_key(eligible))
  matched <- eligible[match(pair_key(out), pair_key(eligible)), ]
  expect_equal(out$observed, as.integer(matched$observed))
  expect_equal(out$expected, matched$expected)
})

test_that("the null is anchored to all pairs, not to the survivors", {
  hg <- planted_hg()
  low <- hypa(hg, min_count = 2L)
  high <- hypa(hg, min_count = 8L)
  a <- subset(low, from == "x" & to == "y")$expected
  b <- subset(high, from == "x" & to == "y")$expected
  # raising the floor changes which pairs are tested, never the null itself
  expect_equal(a, b)
  expect_gt(nrow(low), nrow(high))
})

test_that("hypa returns a tidy frame with ratio = observed / expected", {
  out <- hypa(planted_hg(), min_count = 3L)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out),
    c("from", "to", "observed", "expected", "ratio", "p_under", "p_over",
      "p_adjusted_under", "p_adjusted_over", "anomaly"))
  expect_equal(out$ratio, out$observed / out$expected)
  # ordered by significance, not by ratio (ratio favours pairs at the floor)
  expect_false(is.unsorted(out$p_adjusted_over))
  expect_true(all(out$observed >= 3L | out$expected >= 3))
  expect_identical(anyDuplicated(pair_key(out)), 0L)
})

test_that("p_adjust defaults to none; a correction runs over every possible pair", {
  hg <- planted_hg()
  expect_identical(formals(hypa.net_hg)$p_adjust, "none")
  raw <- hypa(hg, min_count = 3L)
  expect_identical(raw$p_adjusted_over, raw$p_over)
  expect_identical(raw$p_adjusted_under, raw$p_under)
  bh <- hypa(hg, min_count = 3L, p_adjust = "BH")
  co <- .hg_hypa_cooccurrence(hg$incidence)
  n_possible <- choose(sum(co$strength > 0), 2)
  expect_gt(n_possible, nrow(raw))
  bh <- bh[match(pair_key(raw), pair_key(bh)), ]
  expect_equal(bh$p_adjusted_over,
               stats::p.adjust(raw$p_over, "BH", n = n_possible))
  expect_equal(bh$p_adjusted_under,
               stats::p.adjust(raw$p_under, "BH", n = n_possible))
  expect_identical(nrow(hypa(hg, min_count = 3L, top = 4L)), 4L)
})

# ---- calibration ---------------------------------------------------------

test_that("random incidence gives few anomalies across 20 seeds", {
  counts <- vapply(seq_len(20L), \(seed) {
    hg <- random_incidence_hg(seed)
    raw <- hypa(hg, min_count = 5L)
    bh <- hypa(hg, min_count = 5L, p_adjust = "BH")
    n_possible <- choose(sum(.hg_hypa_cooccurrence(hg$incidence)$strength > 0),
                         2)
    c(possible = n_possible, over = sum(raw$anomaly == "over"),
      under = sum(raw$anomaly == "under"), bh = sum(bh$anomaly != "normal"))
  }, numeric(4L))
  # See LEARNINGS.md (2026-10-12) for the measured counts. The strict upper
  # tail of the paper is anti-conservative with small counts, so the
  # uncorrected over-rate sits a little under alpha per possible pair.
  rate <- sum(counts[c("over", "under"), ]) / sum(counts["possible", ])
  expect_lt(rate, 0.06)
  expect_lt(sum(counts["under", ]) / sum(counts["possible", ]), 0.01)
  expect_lte(sum(counts["bh", ]), 3)
})

# ---- invariance ------------------------------------------------------------

test_that("hypa is invariant to relabelling and reordering the nodes", {
  hg <- avoiding_hubs_hg()
  out <- hypa(hg, min_count = 5L)
  # build the same hypergraph from shuffled memberships with new names
  inc <- as.matrix(hg$incidence)
  member <- which(inc != 0, arr.ind = TRUE)
  set.seed(7)
  shuffle <- sample(nrow(member))
  map <- stats::setNames(paste0("N", sample(nrow(inc))), rownames(inc))
  relabelled <- group_hypergraph(
    data.frame(node = unname(map[rownames(inc)[member[shuffle, "row"]]]),
               hyperedge = colnames(inc)[member[shuffle, "col"]]),
    node = "node", hyperedge = "hyperedge")
  expect_false(identical(relabelled$nodes, unname(map[hg$nodes])))
  out2 <- hypa(relabelled, min_count = 5L)
  renamed <- transform(out, from = unname(map[from]), to = unname(map[to]))
  expect_setequal(pair_key(out2), pair_key(renamed))
  out2 <- out2[match(pair_key(renamed), pair_key(out2)), ]
  expect_equal(out2$expected, renamed$expected, tolerance = 1e-8)
  expect_equal(out2$p_under, renamed$p_under, tolerance = 1e-8)
  expect_identical(out2$observed, renamed$observed)
  expect_identical(out2$anomaly, renamed$anomaly)
})

# ---- validation ------------------------------------------------------------

test_that("hypa rejects bad input by class", {
  hg <- planted_hg()
  bad <- "hypergraphs_bad_input"
  expect_error(hypa(hg, min_count = 0), class = bad)
  expect_error(hypa(hg, min_count = 1.5), class = bad)
  expect_error(hypa(hg, min_count = c(2, 3)), class = bad)
  expect_error(hypa(hg, min_count = NA), class = bad)
  expect_error(hypa(hg, alpha = 0.5), class = bad)
  expect_error(hypa(hg, alpha = 0), class = bad)
  expect_error(hypa(hg, alpha = c(0.01, 0.05)), class = bad)
  expect_error(hypa(hg, alpha = NA_real_), class = bad)
  expect_error(hypa(hg, p_adjust = "nope"), class = bad)
  expect_error(hypa(hg, top = 0), class = bad)
  expect_error(hypa(hg, order = 2L), class = bad)
  expect_error(hypa(hg, min_count = 10000L),
               class = "hypergraphs_empty_result")
  singletons <- group_hypergraph(
    data.frame(node = c("a", "b"), hyperedge = c("e1", "e2")),
    node = "node", hyperedge = "hyperedge")
  expect_error(hypa(singletons, min_count = 1L),
               class = "hypergraphs_empty_result")
})
