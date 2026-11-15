# Small exhaustive domains and independent formulations of numerical results.

test_that("every binary 3-node 2-edge incidence obeys set-based definitions", {
  for (mask in 0:63) {
    b <- matrix(as.integer(intToBits(mask))[1:6], 3, 2,
                dimnames = list(letters[1:3], c("e1", "e2")))
    h <- .thg_from_incidence(b, list())
    sp <- .thg_from_incidence(Matrix::Matrix(b, sparse = TRUE), list())
    sets <- lapply(1:2, function(j) letters[which(b[, j] == 1)])
    expected <- outer(letters[1:3], letters[1:3], Vectorize(function(a, z) {
      if (a == z) return(0)
      sum(vapply(sets, function(e) a %in% e && z %in% e, logical(1)))
    }))
    for (hg in list(h, sp)) {
      nodes <- hg_measures(hg)
      expect_equal(nodes$hyperdegree, vapply(letters[1:3], function(a) {
        sum(vapply(sets, function(e) a %in% e, logical(1)))
      }, integer(1)), ignore_attr = TRUE)
      expect_equal(nodes$strength, vapply(letters[1:3], function(a) {
        sum(lengths(sets)[vapply(sets, function(e) a %in% e, logical(1))])
      }, integer(1)), ignore_attr = TRUE)
      expect_equal(unname(pairwise_network(hg, weighted = FALSE)$weights),
                   unname(expected))
      overlap <- hg_measures(hg, what = "overlap")
      inter <- length(intersect(sets[[1]], sets[[2]]))
      union <- length(union(sets[[1]], sets[[2]]))
      expect_equal(overlap$overlap, inter)
      expect_equal(overlap$jaccard, if (union) inter / union else 0)
      for (s in 1:3) {
        line <- as.matrix(hg_line_graph(hg, s = s, what = "matrix"))
        expect_equal(unname(line), matrix(c(0, inter * (inter >= s),
                                            inter * (inter >= s), 0), 2))
      }
    }
    # Solve PageRank directly, including dangling rows, without using any
    # transition-matrix helper from the package.
    p <- matrix(0, 3, 3)
    for (i in 1:3) {
      incident <- which(b[i, ] == 1)
      if (!length(incident)) {
        p[i, ] <- 1 / 3
      } else for (e in incident) {
        p[i, ] <- p[i, ] + b[, e] / sum(b[, e]) / length(incident)
      }
    }
    expected_rank <- as.numeric(solve(diag(3) - 0.85 * t(p), rep(0.05, 3)))
    expect_equal(hg_pagerank(h)$pagerank, expected_rank, tolerance = 1e-10)
    expect_equal(hg_pagerank(sp)$pagerank, expected_rank, tolerance = 1e-10)
  }
})

test_that("Wasserstein agrees with exhaustive partial matchings", {
  # Enumerate every partial one-to-one matching. Points not matched to
  # another diagram are matched to their nearest point on the diagonal.
  exhaustive <- function(a, b, q, p) {
    diag_cost <- function(point) {
      if (is.infinite(p)) abs(diff(point)) / 2 else
        abs(diff(point)) * 2^(1 / p - 1)
    }
    ground <- function(x, y) {
      d <- abs(x - y)
      if (is.infinite(p)) max(d) else sum(d^p)^(1 / p)
    }
    visit <- function(i, remaining) {
      if (i > nrow(a)) {
        return(sum(vapply(remaining, function(j) diag_cost(b[j, ])^q, numeric(1))))
      }
      candidates <- c(diag_cost(a[i, ])^q + visit(i + 1L, remaining),
        vapply(remaining, function(j) {
          ground(a[i, ], b[j, ])^q + visit(i + 1L, setdiff(remaining, j))
        }, numeric(1)))
      min(candidates)
    }
    visit(1L, seq_len(nrow(b)))^(1 / q)
  }
  set.seed(1602)
  for (trial in 1:8) {
    a <- matrix(runif(6), 3, 2)
    b <- matrix(runif(4), 2, 2)
    a[, 2] <- a[, 1] + a[, 2]
    b[, 2] <- b[, 1] + b[, 2]
    da <- data.frame(dimension = 0L, birth = a[, 1], death = a[, 2])
    db <- data.frame(dimension = 0L, birth = b[, 1], death = b[, 2])
    for (q in c(1, 2, 5)) for (p in c(1, 2, Inf)) {
      expected <- exhaustive(a, b, q, p)
      expect_equal(unname(hg_wasserstein(da, db, order = q, internal_p = p)),
                   expected, tolerance = 1e-12)
      expect_equal(unname(hg_wasserstein(db, da, order = q, internal_p = p)),
                   expected, tolerance = 1e-12)
    }
  }
})

test_that("topic shares are invariant to reciprocal factor rescaling", {
  X <- Matrix::Matrix(matrix(1:12, 3, 4,
                             dimnames = list(letters[1:3], LETTERS[1:4])),
                       sparse = TRUE)
  fit <- list(W = matrix(c(1, 2, 4, 8, 2, 1), 3, 2),
              Ht = matrix(c(1, 2, 3, 4, 8, 4, 2, 1), 4, 2),
              iterations = 1L)
  result <- function(f) .tm_result(f, X, 2L, c(NA_real_, NA_real_), list(f),
    0, 1L, TRUE, list(list(mean = 1)), integer(), NULL, list())
  original <- result(fit)
  rescaled <- fit
  rescaled$W <- sweep(fit$W, 2, c(1000, 0.01), "*")
  rescaled$Ht <- sweep(fit$Ht, 2, c(1000, 0.01), "/")
  expect_equal(hg_get(original, what = "shares"), hg_get(result(rescaled), what = "shares"))
  expect_equal(hg_get(original, what = "words"), hg_get(result(rescaled), what = "words"))
  expected <- sweep(fit$W, 2, colSums(fit$Ht), "*")
  expected <- expected / rowSums(expected)
  expected <- expected[, order(-colMeans(expected)), drop = FALSE]
  actual <- hg_get(original, what = "shares")
  expect_equal(actual$share, as.vector(t(expected)), tolerance = 1e-12)
})
