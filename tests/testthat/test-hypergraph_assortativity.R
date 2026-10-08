# hg_assortativity() and hg_degree_correlation(): an enumeration reference
# written from Chodrow's (2020) definition (a different code path from the
# closed-form moments in the package), the graph reduction to Newman's
# coefficient, invariants and error classes. External oracles (XGI,
# hypergraphx) live in local_testing_and_equivalence/.

.as_hg <- function(edges, sparse = FALSE) {
  long <- data.frame(
    member = unlist(edges),
    group = rep(names(edges), lengths(edges)),
    stringsAsFactors = FALSE
  )
  group_hypergraph(long, node = "member", hyperedge = "group", sparse = sparse)
}

# Enumerate the chosen pairs with their probability weights, then take the
# weighted Pearson correlation of the two coordinates.
.as_reference <- function(edges, type, scale) {
  nodes <- sort(unique(unlist(edges)))
  degree <- vapply(nodes, \(v) sum(vapply(edges, \(e) v %in% e, logical(1))),
                   numeric(1))
  score <- if (scale == "rank") rank(degree) else degree
  rows <- lapply(Filter(\(e) length(e) >= 2L, edges), \(e) {
    s <- score[e]
    m <- length(s)
    if (type == "uniform") {
      grid <- expand.grid(i = seq_len(m), j = seq_len(m))
      grid <- subset(grid, i != j)
      data.frame(x = s[grid$i], y = s[grid$j], w = 1 / (m * (m - 1)))
    } else if (type == "top_2") {
      top <- sort(s, decreasing = TRUE)[1:2]
      data.frame(x = top, y = rev(top), w = 0.5)
    } else {
      data.frame(x = max(s), y = min(s), w = 1)
    }
  })
  pairs <- do.call(rbind, rows)
  stats::cov.wt(pairs[c("x", "y")], wt = pairs$w, cor = TRUE)$cor[1, 2]
}

set.seed(19)
as_edges <- lapply(seq_len(40), \(i) paste0("n", sample(20, sample(2:5, 1))))
names(as_edges) <- sprintf("e%02d", seq_along(as_edges))

test_that("every type and scale matches the enumerated definition", {
  hg <- .as_hg(as_edges)
  types <- c("uniform", "top_2", "top_bottom")
  vapply(c("degree", "rank"), \(scale) {
    out <- hg_assortativity(hg, type = types, scale = scale)
    expect_identical(out$type, types)
    expect_identical(out$n_edges, rep(40L, 3L))
    reference <- vapply(types, \(t) .as_reference(as_edges, t, scale),
                        numeric(1))
    expect_equal(out$assortativity, unname(reference), tolerance = 1e-12)
    TRUE
  }, logical(1))
})

test_that("a star graph is perfectly disassortative", {
  edges <- lapply(paste0("leaf", 1:6), \(l) c("hub", l))
  names(edges) <- paste0("e", 1:6)
  out <- hg_assortativity(.as_hg(edges), type = c("uniform", "top_bottom"),
                          scale = "degree")
  expect_equal(out$assortativity[1], -1, tolerance = 1e-12)
  # top-bottom always picks (hub, leaf): both coordinates are constant
  expect_true(is.na(out$assortativity[2]))
})

test_that("on a graph, uniform equals Newman's degree assortativity", {
  skip_if_not_installed("igraph")
  set.seed(5)
  g <- igraph::sample_gnm(30, 60)
  ends <- igraph::as_edgelist(g)
  edges <- lapply(seq_len(nrow(ends)), \(i) paste0("v", ends[i, ]))
  names(edges) <- sprintf("e%03d", seq_along(edges))
  out <- hg_assortativity(.as_hg(edges), type = c("uniform", "top_2"),
                          scale = "degree")
  newman <- igraph::assortativity_degree(g, directed = FALSE)
  expect_equal(out$assortativity, c(newman, newman), tolerance = 1e-12)
})

test_that("undefined coefficients are NA (regular hypergraph)", {
  # every node has hyperdegree 2
  edges <- list(e1 = c("a", "b", "c"), e2 = c("d", "e", "f"),
                e3 = c("a", "d"), e4 = c("b", "e"), e5 = c("c", "f"))
  out <- hg_assortativity(.as_hg(edges), type = c("uniform", "top_bottom"))
  expect_true(all(is.na(out$assortativity)))
})

test_that("dense and sparse agree; relabelling does not matter", {
  types <- c("uniform", "top_2", "top_bottom")
  dense <- hg_assortativity(.as_hg(as_edges), type = types)
  sparse <- hg_assortativity(.as_hg(as_edges, sparse = TRUE), type = types)
  expect_equal(sparse, dense, tolerance = 1e-14)

  relabel <- setNames(paste0("z", sample(20)), paste0("n", seq_len(20)))
  moved <- lapply(rev(as_edges), \(e) unname(relabel[e]))
  names(moved) <- paste0("f", seq_along(moved))
  expect_equal(hg_assortativity(.as_hg(moved), type = types), dense,
               tolerance = 1e-12)
  values <- dense$assortativity
  expect_true(all(values >= -1 & values <= 1))
})

test_that("singleton hyperedges are ignored", {
  with_singletons <- c(as_edges, list(s1 = "n1", s2 = "n2"))
  a <- hg_assortativity(.as_hg(as_edges))
  b <- hg_assortativity(.as_hg(with_singletons), scale = "rank")
  expect_identical(b$n_edges, 40L)
  # degrees change (n1, n2 gain a hyperedge), so only the edge count is fixed
  expect_true(is.finite(b$assortativity))
  expect_true(is.finite(a$assortativity))
})

test_that("hg_degree_correlation: hand-computed Pearson across orders", {
  edges <- list(p1 = c("a", "b"), p2 = c("a", "c"), p3 = c("b", "c"),
                p4 = c("a", "d"), t1 = c("a", "b", "c"),
                t2 = c("b", "c", "d"), q1 = c("a", "b", "c", "d"))
  out <- hg_degree_correlation(.as_hg(edges))
  expect_identical(out$size_1, c(2L, 2L, 3L))
  expect_identical(out$size_2, c(3L, 4L, 4L))
  pair_degree <- c(a = 3, b = 2, c = 2, d = 1)
  triple_degree <- c(a = 1, b = 2, c = 2, d = 1)
  expect_equal(out$correlation[1], stats::cor(pair_degree, triple_degree),
               tolerance = 1e-14)
  # size 4: every node has degree 1 -> constant -> NA
  expect_true(all(is.na(out$correlation[2:3])))
  expect_identical(out$n_nodes, rep(4L, 3L))
})

test_that("broken contracts raise classed errors", {
  expect_error(hg_assortativity(list()), class = "hypergraphs_bad_input")
  expect_error(hg_assortativity(.as_hg(list(s1 = "a", s2 = "b"))),
               class = "hypergraphs_bad_input")
  expect_error(hg_degree_correlation(.as_hg(list(e1 = c("a", "b"),
                                                 e2 = c("b", "c")))),
               class = "hypergraphs_bad_input")
  expect_error(hg_degree_correlation(list()), class = "hypergraphs_bad_input")
})

# ---- Audit regressions (2026-10-06) ---------------------------------------

# Independent oracle: enumerate the chosen ordered pairs of every hyperedge
# with their probabilities (each hyperedge carries total weight 1) and take
# the weighted Pearson correlation with stats::cov.wt().
.as_enumerated <- function(hg, type, scale = "degree") {
  inc <- (as.matrix(hg$incidence) > 0) * 1
  degree <- rowSums(inc)
  x <- if (identical(scale, "rank")) rank(degree) else degree
  edges <- lapply(seq_len(ncol(inc)), \(j) which(inc[, j] > 0))
  edges <- edges[lengths(edges) >= 2L]
  pairs <- do.call(rbind, lapply(edges, \(v) {
    s <- x[v]
    if (identical(type, "uniform")) {
      g <- expand.grid(i = seq_along(v), j = seq_along(v))
      g <- g[g$i != g$j, ]
      data.frame(a = s[g$i], b = s[g$j], w = 1 / nrow(g))
    } else if (identical(type, "top_2")) {
      top <- sort(s, decreasing = TRUE)[1:2]
      data.frame(a = top, b = rev(top), w = 0.5)
    } else {
      data.frame(a = max(s), b = min(s), w = 1)
    }
  }))
  cw <- stats::cov.wt(as.matrix(pairs[c("a", "b")]), wt = pairs$w / sum(pairs$w),
                      cor = TRUE, method = "ML")
  cw$cor[1L, 2L]
}

test_that("A07: near-regular large degrees keep a defined coefficient", {
  n <- 10000L
  d <- data.frame(
    node = c(rep(c("a", "b", "c"), n), "a", "b"),
    hyperedge = c(rep(sprintf("e%05d", seq_len(n)), each = 3), "z", "z")
  )
  hg <- group_hypergraph(d, node = "node", hyperedge = "hyperedge")
  got <- hg_assortativity(hg, type = "uniform", scale = "degree")
  expect_false(is.na(got$assortativity))
  # degrees (10001, 10001, 10000) shift exactly to the 0/1 scores (1, 1, 0).
  # Closed form over the n triangles (ordered pairs (1,1) x 2, (1,0) x 2,
  # (0,1) x 2, each 1/6) and the edge {a, b} (pair (1,1)): with p = E[X]
  # and q = E[XY], r = (q - p^2) / (p (1 - p)).
  p <- (n * 2 / 3 + 1) / (n + 1)
  q <- (n / 3 + 1) / (n + 1)
  expect_equal(got$assortativity, (q - p^2) / (p * (1 - p)), tolerance = 1e-10)
  # top_2 and top_bottom are genuinely constant here: still NA
  expect_true(all(is.na(hg_assortativity(hg, type = c("top_2", "top_bottom"),
                                         scale = "degree")$assortativity)))
})

test_that("A07: coefficients equal the enumerated weighted Pearson", {
  lapply(1:3, \(seed) {
    hg <- random_hypergraph("gnp", n = 10, m = 12, p = 0.35, seed = seed)
    lapply(c("uniform", "top_2", "top_bottom"), \(ty) lapply(
      c("degree", "rank"), \(sc) {
        expect_equal(hg_assortativity(hg, type = ty, scale = sc)$assortativity,
                     .as_enumerated(hg, ty, sc), tolerance = 1e-10)
      }))
  })
})

test_that("A07: scores are shift and scale invariant; regular stays NA", {
  hg <- random_hypergraph("gnp", n = 10, m = 12, p = 0.35, seed = 2)
  membership <- (as.matrix(hg$incidence) > 0) * 1
  members <- lapply(seq_len(ncol(membership)), \(j) which(membership[, j] > 0))
  members <- members[lengths(members) >= 2L]
  x <- rowSums(membership)
  lapply(c("uniform", "top_2", "top_bottom"), \(ty) {
    base <- .hg_assort_one(.hg_standardise_scores(x, members), members, ty)
    moved <- .hg_assort_one(.hg_standardise_scores(7 * x + 1e9, members),
                            members, ty)
    expect_equal(moved, base, tolerance = 1e-10)
  })
  regular <- random_hypergraph("regular", n = 12, m = 8, k = 3, seed = 1)
  expect_true(all(is.na(hg_assortativity(
    regular, type = c("uniform", "top_2", "top_bottom"),
    scale = "degree")$assortativity)))
})
