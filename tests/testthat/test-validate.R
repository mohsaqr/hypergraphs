test_that(".ho_check_count accepts whole numbers and returns integers", {
  expect_identical(.ho_check_count(3, "n"), 3L)
  expect_identical(.ho_check_count(3L, "n"), 3L)
  expect_identical(.ho_check_count(0, "n", min = 0), 0L)
  expect_identical(.ho_check_count(Inf, "n", allow_inf = TRUE), Inf)
})

test_that(".ho_check_count refuses fractions, NA, infinities and overflow", {
  bad <- list(1.5, NA_real_, NaN, -Inf, Inf, c(1, 2), numeric(0), "3",
              factor(3), 0, 1e12)
  lapply(bad, \(x) expect_error(.ho_check_count(x, "n"),
                                class = "hypergraphs_bad_input"))
  expect_error(.ho_check_count(-Inf, "n", allow_inf = TRUE),
               class = "hypergraphs_bad_input")
  expect_error(.ho_check_count(19.5, "n", allow_inf = TRUE),
               class = "hypergraphs_bad_input")
  expect_error(.ho_check_count(5, "k", max = 4), class = "hypergraphs_bad_input")
})

test_that(".ho_check_number checks finiteness and range", {
  expect_identical(.ho_check_number(0.5, "tol", min = 0, max = 1), 0.5)
  lapply(list(NA_real_, Inf, -1, 2, c(0.1, 0.2), "a"),
         \(x) expect_error(.ho_check_number(x, "tol", min = 0, max = 1),
                           class = "hypergraphs_bad_input"))
})

test_that(".ho_check_ids refuses missing, empty and repeated identifiers", {
  expect_identical(.ho_check_ids(c("a", "b"), "nodes"), c("a", "b"))
  lapply(list(c("a", NA), c("a", ""), c("a", "a")),
         \(x) expect_error(.ho_check_ids(x, "nodes"),
                           class = "hypergraphs_bad_input"))
})

test_that(".ho_check_weights refuses negative, non-finite and factor weights", {
  expect_identical(.ho_check_weights(c(0, 1.5), "weight"), c(0, 1.5))
  lapply(list(c(1, -1), c(1, Inf), c(1, NA), factor(c(1, 2)), "1"),
         \(x) expect_error(.ho_check_weights(x, "weight"),
                           class = "hypergraphs_bad_input"))
})
