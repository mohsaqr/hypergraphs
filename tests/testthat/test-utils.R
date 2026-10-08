test_that("padding an empty sequence list is a classed error (R23)", {
  expect_error(.ho_wide_sequences(list()), class = "hypergraphs_bad_input")
  expect_error(.ho_wide_sequences(list(character(), character())),
               class = "hypergraphs_bad_input")
  expect_error(.ho_wide_sequences(list(c("a", "b"), list("c"))),
               class = "hypergraphs_bad_input")
  expect_error(.ho_wide_sequences(list(c("a", "b"), NULL)),
               class = "hypergraphs_bad_input")
  wide <- .ho_wide_sequences(list(c("a", "b"), "c", character()))
  expect_identical(dim(wide), c(3L, 2L))
  expect_identical(names(wide), c("T1", "T2"))
  expect_true(all(is.na(unlist(wide[3L, ]))))
  expect_error(.ho_sequence_input(list(), lists = "wide"),
               class = "hypergraphs_bad_input")
})

test_that(".ho_apply surfaces a failed parallel worker as a classed error", {
  skip_on_os("windows")
  expect_identical(.ho_apply(1:4, \(i) i^2, parallel = TRUE, n_cores = 2),
                   lapply(1:4, \(i) i^2))
  expect_error(.ho_apply(1:4, \(i) if (i == 3) stop("boom") else i,
                         parallel = TRUE, n_cores = 2),
               class = "hypergraphs_parallel_failed")
})

test_that(".ho_apply reruns a killed worker serially, with a warning", {
  skip_on_os("windows")
  parent <- Sys.getpid()
  # the worker for 3 kills its own forked process; the serial rerun in the
  # parent computes it normally
  fn <- \(i) {
    if (i == 3 && Sys.getpid() != parent) tools::pskill(Sys.getpid())
    i^2
  }
  expect_warning(out <- .ho_apply(1:4, fn, parallel = TRUE, n_cores = 2),
                 class = "hypergraphs_parallel_fallback")
  expect_identical(out, lapply(1:4, \(i) i^2))
})

test_that(".ho_kmeans equals stats::kmeans when every start finishes", {
  set.seed(21)
  for (r in 1:30) {
    n <- sample(c(20, 60, 150), 1)
    # rounded data has duplicate rows, which changes how kmeans draws
    x <- matrix(round(stats::rnorm(n * 3), sample(c(1, 6), 1)), n)
    k <- sample(2:5, 1)
    nstart <- sample(c(1L, 2L, 10L), 1)
    seed <- sample.int(1e6, 1)
    set.seed(seed)
    reference <- tryCatch(stats::kmeans(x, k, nstart = nstart, iter.max = 100L),
                          warning = function(w) NULL)
    if (is.null(reference)) next
    set.seed(seed)
    expect_identical(unclass(.ho_kmeans(x, k, nstart = nstart)),
                     unclass(reference))
  }
})

test_that(".ho_kmeans finishes a start that stops early", {
  set.seed(2)
  x <- matrix(stats::rnorm(400), 200)
  # iter.max = 1 stops every Hartigan-Wong start early
  expect_warning(stats::kmeans(x, 4, iter.max = 1L))
  set.seed(2)
  fit <- expect_no_warning(.ho_kmeans(x, 4, nstart = 3, iter.max = 1L))
  # a finished k-means partition is a fixed point: every centre is the mean
  # of its points, and every point is closest to its own centre
  means <- rowsum(x, fit$cluster) / as.vector(table(fit$cluster))
  expect_equal(unname(means), unname(fit$centers), tolerance = 1e-10)
  d2 <- outer(rowSums(x^2), rowSums(fit$centers^2), `+`) -
    2 * x %*% t(fit$centers)
  expect_identical(max.col(-d2, ties.method = "first"), unname(fit$cluster))
  # a start not finished even by Lloyd is reported
  set.seed(2)
  expect_warning(.ho_kmeans(x, 4, nstart = 2, iter.max = 1L,
                            max_continue = 0L, lloyd_iter = 1L),
                 class = "hypergraphs_no_converge")
})
