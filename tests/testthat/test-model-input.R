# Networks reach the sequence verbs only as their sequences: a network's fields,
# a one-hot network's 0/1 codes, an association network's variables and a group of
# networks are never read as states.

.mi_onehot <- function() {
  data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
             C = rep(c(0, 0, 1, 1), 6L), actor = rep(paste0("p", 1:6), each = 4L))
}

test_that("window constructors read a sequence model as its sequences", {
  d <- data.frame(T1 = c("A", "B"), T2 = c("B", "C"), T3 = c("C", "A"))
  n <- Nestimate::build_network(d, method = "frequency")
  # windows of 2: {A,B} {B,C} from row one, {B,C} {C,A} from row two
  edges <- hg_get(window_hypergraph(n, window = 2L))
  expect_identical(edges$members, c("A, B", "A, C", "B, C"))
  expect_identical(edges$weight, c(1, 1, 2))
  expect_identical(edges, hg_get(window_hypergraph(d, window = 2L)))
  sc <- simplicial(n, type = "window", window = 2L)
  expect_identical(sort(unique(unlist(strsplit(hg_get(sc)$members, ", ")))),
                   c("A", "B", "C"))
  # three vertices and the three pairs AB, BC, CA; no window holds three states
  expect_equal(unname(sc$f_vector), c(3, 3))
})

test_that("window constructors refuse one-hot networks and groups", {
  hot <- Nestimate::build_network(.mi_onehot(), method = "frequency",
                                  format = "onehot", actor = "actor",
                                  window_size = 2L)
  expect_error(window_hypergraph(hot, window = 2L), class = "hypergraphs_bad_input")
  expect_error(simplicial(hot, type = "window", window = 2L),
               class = "hypergraphs_bad_input")
  d <- data.frame(T1 = c("A", "B"), T2 = c("B", "C"), T3 = c("C", "A"),
                  g = c("x", "y"))
  grp <- Nestimate::build_network(d, method = "frequency", group = "g",
                                  state_cols = c("T1", "T2", "T3"))
  expect_error(window_hypergraph(grp, window = 2L), class = "hypergraphs_bad_input")
})

test_that("group_hypergraph() refuses a group of one-hot networks", {
  q <- .mi_onehot()
  q$cohort <- rep(c("g1", "g2"), each = 12L)
  grp <- Nestimate::build_network(q, method = "frequency", format = "onehot",
                                  actor = "actor", group = "cohort",
                                  window_size = 2L)
  expect_error(group_hypergraph(grp, top = 3L), class = "hypergraphs_bad_input")
})

test_that("a bare cograph_network and an association network are not sequences", {
  hot <- Nestimate::build_network(.mi_onehot(), method = "frequency",
                                  format = "onehot", actor = "actor",
                                  window_size = 2L)
  cg <- hot
  class(cg) <- "cograph_network"
  expect_error(hon(cg, max_order = 1L), class = "hypergraphs_bad_input")
  scores <- data.frame(a = c(1, 3, 2, 5, 4), b = c(2, 1, 4, 3, 5),
                       c = c(5, 4, 3, 1, 2))
  expect_error(hon(Nestimate::build_network(scores, method = "cor"), max_order = 1L),
               class = "hypergraphs_bad_input")
})

test_that("duplicated column names and one 0/1 column are refused", {
  d <- data.frame(week_a = c("A", "B"), week_b = c("B", "A"), week_c = c("C", "B"))
  x <- Nestimate::extract_sequences(d, cols = names(d))
  x[[4]] <- c("NOVEL", "NOVEL")
  names(x)[4] <- "week_a"
  expect_null(.ho_marked_positions(x))
  expect_error(hon(x, max_order = 1L), class = "hypergraphs_bad_input")
  expect_error(hon(data.frame(T1 = c(1, 0, 1)), max_order = 1L),
               class = "hypergraphs_bad_input")
})

test_that("a model of declared 0/1 positions reads as the states 0 and 1", {
  u <- data.frame(T1 = c(1, 0), T2 = c(0, 0), T3 = c(1, 1))
  n <- Nestimate::build_network(u, method = "frequency", state_cols = names(u))
  # rows 1 0 1 and 0 0 1: transitions 1->0, 0->1, 0->0, 0->1
  tr <- hg_get(hon(n, max_order = 1L))
  expect_identical(tr$path, c("0 -> 0", "0 -> 1", "1 -> 0"))
  expect_identical(tr$count, c(1L, 2L, 1L))
  expect_identical(tr$from, c("0", "0", "1"))
  expect_identical(tr$to, c("0", "1", "0"))
  # from 0: one of three to 0, two to 1; from 1: always to 0
  expect_equal(tr$probability, c(1 / 3, 2 / 3, 1))
  mo <- mogen(n, max_order = 1L)
  expect_s3_class(mo, "net_mogen")
  # three 0s and three 1s; order 1 adds the transitions to the first states
  expect_equal(unname(mo$count_matrices[[1]]), c(3, 3))
  expect_equal(unname(mo$count_matrices[[2]]), matrix(c(1, 1, 2, 0), 2L))
  expect_equal(unname(mo$log_likelihood),
               c(6 * log(1 / 2), 2 * log(1 / 2) + log(1 / 3) + 2 * log(2 / 3)))
  hy <- hypa(n, order = 1L, min_count = 1L)
  expect_s3_class(hy, "net_hypa")
  sc <- hy$scores
  # the fitted expectations keep the four observed transitions, and the tails
  # are P(X <= f) and P(X > f) of one hypergeometric, so they add up to 1
  expect_equal(sum(sc$expected), 4)
  expect_identical(sc$path, c("0 -> 0", "1 -> 0", "0 -> 1"))
  expect_equal(sc$observed, c(1, 1, 2))
  # the fit preserves every strength, and here it reproduces the counts
  expect_equal(sc$expected, c(1, 1, 2), tolerance = 1e-8)
  expect_equal(sc$p_under + sc$p_over, rep(1, 3L))
  # xi = 4, 4, 8 of 16, four draws: P(X <= 1) for xi = 4 is
  # (C(12,4) + 4 C(12,3)) / C(16,4) = 1375 / 1820, P(X <= 2) for xi = 8 is
  # (C(8,4) + 8 C(8,3) + C(8,2)^2) / C(16,4) = 1302 / 1820
  expect_equal(sc$p_under, c(1375, 1375, 1302) / 1820)
})

test_that("a fitted result is not read as sequences", {
  expect_error(hon(stats::t.test(1:6, 2:7), max_order = 1L))
})
