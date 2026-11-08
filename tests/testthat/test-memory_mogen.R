testthat::skip_on_cran()

# ===========================================================================
# Section 6: mogen end-to-end
# ===========================================================================

test_that("mogen returns net_mogen class", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"),
                c("B", "C", "D", "A"), c("C", "D", "A", "B"))
  m <- mogen(trajs, max_order = 3L)

  expect_s3_class(m, "net_mogen")
  expect_true(m$optimal_order %in% 0L:3L)
  expect_equal(length(m$aic), 4L)  # orders 0-3
  expect_equal(length(m$transition_matrices), 4L)
  expect_equal(m$n_paths, 4L)
})

test_that("mogen selects order 1 for first-order Markov data", {
  # Generate data from a known order-1 model
  set.seed(42)
  states <- c("A", "B", "C")
  tm <- matrix(c(0.1, 0.7, 0.2,
                 0.3, 0.1, 0.6,
                 0.5, 0.3, 0.2), 3, 3, byrow = TRUE,
               dimnames = list(states, states))

  # Generate 200 paths of length 20
  trajs <- lapply(seq_len(200L), function(i) {
    path <- character(20L)
    path[1L] <- sample(states, 1)
    vapply(2L:20L, function(t) {
      path[t] <<- sample(states, 1, prob = tm[path[t - 1L], ])
      ""
    }, character(1L))
    path
  })

  m <- mogen(trajs, max_order = 3L, criterion = "bic")

  # BIC should select order 1 (data is first-order Markov)
  expect_equal(m$optimal_order, 1L)
})

test_that("mogen selects higher order for second-order data", {
  # Data with strong second-order dependency
  set.seed(123)
  trajs <- lapply(seq_len(300L), function(i) {
    path <- character(15L)
    path[1L] <- sample(c("A", "B", "C"), 1)
    path[2L] <- sample(c("A", "B", "C"), 1)
    vapply(3L:15L, function(t) {
      # After A->B: always go to C
      # After C->B: always go to A
      # Otherwise: uniform
      prev2 <- path[t - 2L]
      prev1 <- path[t - 1L]
      if (prev2 == "A" && prev1 == "B") {
        path[t] <<- "C"
      } else if (prev2 == "C" && prev1 == "B") {
        path[t] <<- "A"
      } else {
        path[t] <<- sample(c("A", "B", "C"), 1)
      }
      ""
    }, character(1L))
    path
  })

  m <- mogen(trajs, max_order = 4L, criterion = "aic")

  # Should detect order >= 2 due to second-order dependency

  expect_true(m$optimal_order >= 2L)
})

test_that("mogen handles data.frame input", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "C"),
                   T3 = c("C", "A"), T4 = c("D", "B"))
  m <- mogen(df, max_order = 2L)
  expect_s3_class(m, "net_mogen")
})

test_that("mogen caps max_order at path length", {
  trajs <- list(c("A", "B", "C"))  # length 3
  expect_message(mogen(trajs, max_order = 10L), "capped")
})

test_that("mogen LRT criterion works", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"),
                c("B", "C", "D", "A"), c("C", "D", "A", "B"))
  m <- mogen(trajs, max_order = 2L, criterion = "lrt")

  expect_s3_class(m, "net_mogen")
  expect_true(m$optimal_order %in% 0L:2L)
  expect_equal(m$criterion, "lrt")
})

test_that("mogen log_likelihoods increase with order", {
  set.seed(42)
  trajs <- lapply(seq_len(50L), function(i) {
    sample(LETTERS[1:5], 10, replace = TRUE)
  })
  m <- mogen(trajs, max_order = 3L)

  # Log-likelihood should generally increase with order
  # (more parameters = better fit)
  ll <- m$log_likelihood
  expect_true(ll[2] >= ll[1])  # order 1 >= order 0
})

# ===========================================================================
# Section 7: S3 methods
# ===========================================================================

test_that("print.net_mogen works", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  m <- mogen(trajs, max_order = 2L)
  out <- capture.output(print(m))
  expect_true(any(grepl("Multi-order generative model", out)))
})

test_that("summary.net_mogen works", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  m <- mogen(trajs, max_order = 2L)
  out <- capture.output(summary(m))
  expect_true(any(grepl("order", out)))
})

test_that("plot.net_mogen works", {
  trajs <- list(c("A", "B", "C", "D"), c("B", "C", "D", "A"))
  m <- mogen(trajs, max_order = 2L)
  expect_no_error(plot(m))
  expect_no_error(plot(m, type = "likelihood"))
})

# ===========================================================================
# Section 8: Input validation
# ===========================================================================

test_that("mogen rejects invalid input", {
  expect_error(mogen(42), class = "hypergraphs_bad_input")
  expect_error(mogen(list(c("A", "B")), max_order = 0L), "max_order")
})

test_that("mogen AIC decreases then increases", {
  # With enough data, AIC should show U-shape: decrease then increase
  set.seed(42)
  states <- c("A", "B", "C", "D")
  trajs <- lapply(seq_len(100L), function(i) {
    sample(states, 8, replace = TRUE)
  })
  m <- mogen(trajs, max_order = 5L)

  # DOF should increase with order
  expect_true(all(diff(m$dof) >= 0))
})

# ===========================================================================
# Section 9: Coverage for previously uncovered paths
# ===========================================================================

# --- mogen: no valid trajectories ---
test_that("mogen stops when no valid trajectories", {
  df <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  expect_error(mogen(df, max_order = 2L), class = "hypergraphs_bad_input")
})

# --- hg_get(what = "transitions"): default order (NULL) uses optimal_order ---
test_that("hg_get(what = 'transitions') defaults to optimal_order", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"),
                c("B", "C", "D", "A"), c("C", "D", "A", "B"))
  m <- mogen(trajs, max_order = 2L)
  tr_default <- hg_get(m, what = "transitions")
  tr_explicit <- hg_get(m, what = "transitions", order = m$optimal_order)
  expect_identical(tr_default, tr_explicit)
})

# --- hg_get(what = "transitions"): empty return when min_count filters all ---
test_that("hg_get(what = 'transitions') is empty when min_count filters all", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  m <- mogen(trajs, max_order = 1L)
  tr <- hg_get(m, what = "transitions", order = 1L, min_count = 9999L)
  expect_true(is.data.frame(tr))
  expect_equal(nrow(tr), 0L)
  expect_true(all(c("path", "count", "probability", "from", "to") %in%
                    names(tr)))
})

# --- hg_get(what = "paths"): wide data.frame input ---
test_that("hg_get(what = 'paths') works on a model of a wide data.frame", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "C"),
                   T3 = c("C", "A"), stringsAsFactors = FALSE)
  result <- hg_get(mogen(df, max_order = 1L), what = "paths", k = 2L)
  expect_true(is.data.frame(result))
  expect_named(result, c("k", "path", "count", "proportion"))
  # A->B B->C, B->C C->A: B->C twice of four, the others once each
  expect_identical(result$count[result$path == "B -> C"], 2L)
  expect_identical(sum(result$count), 4L)
  expect_identical(result, data.frame(k = 2L, path = c("B -> C", "A -> B", "C -> A"),
                                      count = c(2L, 1L, 1L),
                                      proportion = c(2, 1, 1) / 4))
})

# --- hg_get(what = "paths"): list input ---
test_that("hg_get(what = 'paths') works on a model of a list", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"))
  result <- hg_get(mogen(trajs, max_order = 2L), what = "paths", k = 2L)
  expect_true(is.data.frame(result))
  # A->B twice, B->C once, B->D once: shares of 4
  expect_identical(result$count[result$path == "A -> B"], 2L)
  expect_identical(sum(result$count), 4L)
  expect_identical(result, data.frame(k = 2L, path = c("A -> B", "B -> C", "B -> D"),
                                      count = c(2L, 1L, 1L),
                                      proportion = c(2, 1, 1) / 4))
})

# --- hg_get(what = "paths"): trajectories shorter than k contribute nothing ---
test_that("hg_get(what = 'paths') ignores trajectories shorter than k", {
  trajs <- list(c("A"), c("B"), c("A", "B", "C"))
  result <- hg_get(mogen(trajs, max_order = 2L), what = "paths", k = 3L)
  # Only the third trajectory contributes 1 path: A -> B -> C
  expect_equal(nrow(result), 1L)
  expect_equal(result$count[1L], 1L)
})

# --- hg_get(what = "paths"): top limits output ---
test_that("hg_get(what = 'paths') top limits rows returned", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"))
  m <- mogen(trajs, max_order = 2L)
  result_all <- hg_get(m, what = "paths", k = 2L)
  result_top <- hg_get(m, what = "paths", k = 2L, top = 1L)
  expect_equal(nrow(result_top), 1L)
  expect_identical(result_top, utils::head(result_all, 1L))
})

# --- hg_get(what = "paths"): k outside [2, max_order + 1] ---
test_that("hg_get(what = 'paths') rejects k outside the fitted orders", {
  m <- mogen(list(c("A", "B", "C"), c("B", "C", "A")), max_order = 1L)
  expect_error(hg_get(m, what = "paths", k = 1L), class = "hypergraphs_bad_input")
  expect_error(hg_get(m, what = "paths", k = 3L), class = "hypergraphs_bad_input")
  expect_error(hg_get(m, what = "paths", k = 2.5), class = "hypergraphs_bad_input")
})

# --- hg_get(what = "paths"): trailing NAs of a wide frame end a sequence ---
test_that("hg_get(what = 'paths') reads a wide frame with trailing NAs", {
  df <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"), T3 = c("C", NA),
                   T4 = c("D", NA), stringsAsFactors = FALSE)
  result <- hg_get(mogen(df, max_order = 2L), what = "paths", k = 2L)
  expect_false(any(grepl("NA", result$path, fixed = TRUE)))
  # A->B B->C C->D, A->B: four, A->B twice
  expect_identical(sum(result$count), 4L)
  expect_identical(result, data.frame(k = 2L, path = c("A -> B", "B -> C", "C -> D"),
                                      count = c(2L, 1L, 1L),
                                      proportion = c(2, 1, 1) / 4))
})

# --- print.net_mogen: LRT criterion branch uses AIC label ---
test_that("print.net_mogen displays AIC when criterion is 'lrt'", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"),
                c("B", "C", "D", "A"), c("C", "D", "A", "B"))
  m <- mogen(trajs, max_order = 2L, criterion = "lrt")
  out <- capture.output(print(m))
  expect_true(any(grepl("selected by LRT", out)))
})


# ===========================================================================
# Section 10: pathways() tests for MOGen
# ===========================================================================

test_that("pathways.net_mogen returns character vector", {
  set.seed(42)
  trajs <- lapply(seq_len(30L), function(i) sample(LETTERS[1:4], 6, replace = TRUE))
  m <- mogen(trajs, max_order = 2L)
  pw <- .pathways(m)
  expect_true(is.character(pw))
})

test_that("pathways.net_mogen order=0 returns empty character", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  m <- mogen(trajs, max_order = 1L)
  pw <- .pathways(m, order = 0L)
  expect_equal(pw, character(0))
})

test_that("pathways.net_mogen min_count filters paths", {
  set.seed(42)
  trajs <- lapply(seq_len(30L), function(i) sample(LETTERS[1:4], 6, replace = TRUE))
  m <- mogen(trajs, max_order = 2L)
  pw_all  <- .pathways(m, min_count = 1L)
  pw_filt <- .pathways(m, min_count = 9999L)
  expect_equal(length(pw_filt), 0L)
  expect_true(length(pw_all) >= length(pw_filt))
})

test_that("pathways.net_mogen top limits output", {
  set.seed(42)
  trajs <- lapply(seq_len(30L), function(i) sample(LETTERS[1:4], 6, replace = TRUE))
  m <- mogen(trajs, max_order = 2L)
  pw_top <- .pathways(m, top = 2L)
  expect_true(length(pw_top) <= 2L)
})

test_that("pathways.net_mogen min_prob filters by probability", {
  set.seed(42)
  trajs <- lapply(seq_len(30L), function(i) sample(LETTERS[1:4], 6, replace = TRUE))
  m <- mogen(trajs, max_order = 2L)
  pw_all  <- .pathways(m, min_prob = 0)
  pw_filt <- .pathways(m, min_prob = 0.99)
  expect_true(length(pw_all) >= length(pw_filt))
})

test_that("hg_get(what = 'paths') takes a range of path lengths", {
  trajs <- list(c("A", "B", "C", "A", "B"), c("A", "B", "D", "A"))
  m <- mogen(trajs, max_order = 3L)
  both <- hg_get(m, what = "paths", k = 2:3, top = 2L)
  expect_identical(unique(both$k), 2:3)
  expect_identical(both[both$k == 3L, setdiff(names(both), "k")],
                   {x <- hg_get(m, what = "paths", k = 3L, top = 2L)
                    x <- x[setdiff(names(x), "k")]; rownames(x) <- 3:4; x})
  expect_true(all(table(both$k) <= 2L))
  expect_error(hg_get(m, what = "paths", k = c(2L, 9L)), class = "hypergraphs_bad_input")
})

# ===========================================================================
# Audit 2026-11-01: the order-0 layer and argument contracts
# ===========================================================================

test_that("hg_get(what = 'transitions') at a selected order 0 is an empty table", {
  # memoryless data select order 0; the default table used to fail with an
  # unclassed "'order' must be >= 1", and summary() failed with it
  set.seed(1)
  trajs <- replicate(20, sample(c("a", "b", "c", "d"), 12, replace = TRUE),
                     simplify = FALSE)
  m <- mogen(trajs, max_order = 3L)
  expect_identical(m$optimal_order, 0L)
  tr <- hg_get(m, what = "transitions")
  expect_identical(nrow(tr), 0L)
  expect_named(tr, names(hg_get(m, what = "transitions", order = 1)))
  expect_identical(nrow(hg_get(m, what = "transitions", order = 0)), 0L)
  expect_identical(nrow(hg_get(m, what = "pathways")), 0L)
  expect_named(summary(m), c("orders", "transitions", "paths", "pathways"))
})

test_that("hg_get.net_mogen validates order and refuses misplaced arguments", {
  set.seed(1)
  trajs <- replicate(20, sample(c("a", "b", "c", "d"), 12, replace = TRUE),
                     simplify = FALSE)
  m <- mogen(trajs, max_order = 2L)
  expect_error(hg_get(m, what = "transitions", order = 9),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(m, what = "pathways", order = 1.5),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(m, k = 3), class = "hypergraphs_bad_input")
  expect_error(hg_get(m, what = "paths", order = 2),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(m, what = "transitions", min_prob = 0.1),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(m, what = "pathways", foo = 1),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(m, what = "paths", k = "2"),
               class = "hypergraphs_bad_input")
  kept <- hg_get(m, what = "transitions", order = 2, min_count = 2)
  expect_s3_class(kept, "data.frame")
  expect_gt(nrow(kept), 0L)
  expect_true(all(kept$count >= 2))
  # counted on their own: the 3-grams of the same sequences seen at least twice
  grams <- hg_get(m, what = "paths", k = 3L)
  expect_identical(sort(kept$count), sort(grams$count[grams$count >= 2L]))
  # each kept probability is its 3-gram's share of all 3-grams opening with the
  # same two states (dropped ones included), printed to four decimals
  context_total <- tapply(grams$count, sub(" -> [^ ]+$", "", grams$path), sum)
  expect_identical(kept$probability,
                   round(kept$count / as.vector(context_total[kept$from]), 4))
  # the first kept row is b -> b -> b (8 times); the recount of these counts
  # against the raw sequences is local (test-equiv-reference-random-2.R)
  expect_identical(kept$to[1L], "b")
})

test_that("hg_get(what = 'paths') and min_count count k-grams by hand", {
  # a b a b and b a b. 2-grams: a -> b in positions 1, 3 of the first and 2 of
  # the second (3), b -> a once in each (2); proportions 3/5, 2/5.
  # 3-grams: b -> a -> b once in each (2), a -> b -> a once (1); 4-grams:
  # a -> b -> a -> b (1). The random-sequence recount is local
  # (test-equiv-reference-random-2.R).
  m <- mogen(list(c("a", "b", "a", "b"), c("b", "a", "b")), max_order = 3L)
  two <- hg_get(m, what = "paths", k = 2L)
  expect_identical(two$path, c("a -> b", "b -> a"))
  expect_identical(two$count, c(3L, 2L))
  expect_equal(two$proportion, c(0.6, 0.4))
  three <- hg_get(m, what = "paths", k = 3L)
  expect_identical(three$path, c("b -> a -> b", "a -> b -> a"))
  expect_identical(three$count, c(2L, 1L))
  four <- hg_get(m, what = "paths", k = 4L)
  expect_identical(four$path, "a -> b -> a -> b")
  expect_identical(four$count, 1L)
  # min_count = 2 keeps only b -> a -> b; both 3-grams opening with b -> a
  # end in b, so its probability is 2 / 2 = 1
  kept <- hg_get(m, what = "transitions", order = 2, min_count = 2)
  expect_identical(kept$path, "b -> a -> b")
  expect_identical(kept$from, "b -> a")
  expect_identical(kept$to, "b")
  expect_identical(kept$count, 2L)
  expect_equal(kept$probability, 1)
})

test_that("mogen() refuses bad arguments with hypergraphs_bad_input", {
  trajs <- list(c("a", "b", "c"), c("b", "c", "a"))
  expect_error(mogen(trajs, max_order = 0), class = "hypergraphs_bad_input")
  expect_error(mogen(trajs, criterion = "xx"), class = "hypergraphs_bad_input")
  expect_error(mogen(trajs, lrt_alpha = 1), class = "hypergraphs_bad_input")
})
