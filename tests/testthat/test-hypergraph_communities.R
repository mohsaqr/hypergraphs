testthat::skip_on_cran()

test_that("hg_get() filters, orders and truncates the runs of a community fit", {
  skip_if_not_installed("igraph")
  dat <- data.frame(
    member = c("a", "b", "c", "a", "b", "d", "x", "y", "z", "x", "y", "c"),
    edge = rep(paste0("e", 1:4), each = 3)
  )
  h <- group_hypergraph(dat, "member", "edge")
  expect_warning(
    fit <- hg_communities(h, type = "irmm", n_runs = 4, seeds = 1:4, max_iter = 1),
    class = "hypergraphs_no_converge")
  runs <- hg_get(fit, what = "runs")
  expect_identical(nrow(hg_get(fit, what = "runs", converged = TRUE)), 0L)
  expect_identical(hg_get(fit, what = "runs", converged = FALSE),
                   subset(runs, !converged, drop = FALSE) |> `rownames<-`(NULL))
  best <- hg_get(fit, what = "runs", sort_by = "modularity", top = 1)
  expect_identical(nrow(best), 1L)
  expect_equal(best$modularity, max(runs$modularity))
  # INVARIANT: sorting permutes the runs, it never drops or adds one
  sorted <- hg_get(fit, what = "runs", sort_by = "modularity")
  expect_setequal(sorted$run, runs$run)
  expect_false(is.unsorted(rev(sorted$modularity)))
  expect_identical(nrow(hg_get(fit, what = "sizes", top = 1)), 1L)
  expect_error(hg_get(fit, what = "sizes", converged = TRUE),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "runs", sort_by = "seed"),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "runs", converged = NA),
               class = "hypergraphs_bad_input")
  infomap <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
  expect_error(hg_get(infomap, what = "runs", converged = TRUE),
               class = "hypergraphs_bad_input")
})
