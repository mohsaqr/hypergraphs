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

# ---- Audit regressions (2026-10-06) ---------------------------------------

.cm_hg <- function(node, edge) {
  group_hypergraph(data.frame(node = node, hyperedge = edge),
                   node = "node", hyperedge = "hyperedge")
}

test_that("A13: fits sharing fewer than two nodes have undefined agreement", {
  skip_if_not_installed("igraph")
  path <- .cm_hg(c("a", "b", "b", "c"), c("X", "X", "Y", "Y"))
  xy <- .cm_hg(c("x", "y"), c("Z", "Z"))
  ay <- .cm_hg(c("a", "y"), c("Z", "Z"))
  f <- hg_communities(path, n_runs = 2, trials = 2)
  g <- hg_communities(xy, n_runs = 2, trials = 2)
  one <- hg_communities(ay, n_runs = 2, trials = 2)
  expect_warning(cmp <- hg_compare_communities(a = f, b = g),
                 class = "hypergraphs_undefined_statistic")
  sim <- hg_get(cmp, what = "similarity")
  expect_identical(sim$n_nodes, 0L)
  expect_true(all(is.na(unlist(sim[c("ami", "ari", "nmi")]))))
  # a single shared node: still no pair of nodes to agree on
  expect_warning(cmp1 <- hg_compare_communities(a = f, b = one),
                 class = "hypergraphs_undefined_statistic")
  expect_identical(hg_get(cmp1, what = "similarity")$n_nodes, 1L)
  expect_true(is.na(hg_get(cmp1, what = "similarity")$ami))
  # same universe: defined, and identical fits agree perfectly
  expect_no_warning(same <- hg_compare_communities(a = f, b = f))
  expect_equal(unlist(hg_get(same, what = "similarity")[c("ami", "ari", "nmi")]),
               c(ami = 1, ari = 1, nmi = 1))
})

test_that("A14: comparison quality scores each fit's saved projection", {
  skip_if_not_installed("igraph")
  h <- .cm_hg(c("a", "b", "c", "d", "e", "f", "c", "d"),
              c("E1", "E1", "E1", "E2", "E2", "E2", "E3", "E3"))
  irmm <- hg_communities(h, type = "irmm", n_runs = 1)
  info <- hg_communities(h, n_runs = 2, trials = 2)
  q <- hg_get(hg_compare_communities(irmm = irmm, info = info, hg = h),
              what = "quality")
  # independent: cograph's modularity on the IRMM reweighted graph
  direct <- cograph::cluster_quality(irmm$projection, irmm$medoid$community,
                                     weighted = TRUE, directed = FALSE)
  expect_equal(q$modularity[[1L]], direct$global$modularity, tolerance = 1e-12)
  # Infomap fit: its saved projection is the association graph, so the score
  # equals the public quality verb on the same hypergraph
  expect_equal(q[2L, -1L],
               hg_community_quality(h, info), tolerance = 1e-12,
               ignore_attr = TRUE)
  # a directed citation fit has no undirected quality score
  hs <- h
  hs$edge_data <- data.frame(edge = c("E1", "E2", "E3"),
                             source = c("a", "d", "c"))
  directed <- hg_communities(hs, n_runs = 2, trials = 2,
                             method = "citation", directed = TRUE)
  expect_warning(
    mixed <- hg_compare_communities(info = info, dir = directed, hg = hs),
    class = "hypergraphs_undefined_statistic"
  )
  mixed_quality <- hg_get(mixed, what = "quality")
  expect_identical(mixed_quality$model, c("info", "dir"))
  expect_true(all(is.na(unlist(subset(mixed_quality, model == "dir",
                                      select = -model)))))
  expect_false(anyNA(subset(mixed_quality, model == "info",
                            select = c(coverage, modularity))))
  # edge_source is deprecated: the saved projection carries the sources
  expect_warning(hg_compare_communities(irmm = irmm, info = info, hg = h,
                                        edge_source = "source"),
                 class = "hypergraphs_deprecated")
})

test_that("A15: conflicting repeated partition entries are refused", {
  path <- .cm_hg(c("a", "b", "b", "c"), c("X", "X", "Y", "Y"))
  conflicting <- list(
    c(a = "A", a = "B", b = "A", c = "B"),
    data.frame(node = c("a", "a", "b", "c"), community = c("A", "B", "A", "B")),
    c(a = "A", b = NA, c = "B"),
    stats::setNames(c("A", "B", "A"), c("a", "", "c"))
  )
  lapply(conflicting, \(p) {
    expect_error(hg_modularity(path, partition = p),
                 class = "hypergraphs_bad_input")
    expect_error(hg_community_quality(path, p),
                 class = "hypergraphs_bad_input")
  })
  # an identical repeat is one assignment; a shuffled mapping is the same
  base <- hg_modularity(path, partition = c(a = "A", b = "A", c = "B"))
  expect_identical(
    hg_modularity(path, partition = c(a = "A", a = "A", b = "A", c = "B")),
    base)
  expect_identical(
    hg_modularity(path, partition = c(c = "B", a = "A", b = "A")), base)
  expect_identical(
    hg_community_quality(path, c(c = "B", a = "A", a = "A", b = "A")),
    hg_community_quality(path, c(a = "A", b = "A", c = "B")))
})

test_that("A17: the agreement heatmap keeps negative, zero and missing apart", {
  skip_if_not_installed("igraph")
  h <- .cm_hg(c("a", "b", "c", "d", "e", "f", "c", "d"),
              c("E1", "E1", "E1", "E2", "E2", "E2", "E3", "E3"))
  f <- hg_communities(h, n_runs = 2, trials = 2)
  g <- hg_communities(h, n_runs = 2, trials = 2, duplicate_edges = "collapse")
  k <- hg_communities(h, n_runs = 2, trials = 2, seeds = 3:4)
  cmp <- hg_compare_communities(a = f, b = g, c = k)
  cmp$similarity$ami <- c(-0.5, 0, NA)
  cmp$similarity$ari <- c(0.5, -1, 1)
  tiles <- ggplot2::ggplot_build(plot(cmp, what = "similarity"))$data[[1L]]
  m <- .thg_similarity_matrix(cmp)
  # map each tile back to its matrix cell (y counts rows from the bottom)
  value <- m[cbind(nrow(m) + 1L - tiles$y, tiles$x)]
  fill_of <- function(v) unique(tiles$fill[!is.na(value) & value == v])
  expect_identical(fill_of(0), "#FFFFFF")
  expect_length(unique(c(fill_of(-0.5), fill_of(-1), fill_of(0),
                         fill_of(0.5), fill_of(1))), 5L)
  expect_identical(fill_of(-1), "#D33F6A")
  expect_identical(fill_of(1), "#4A6FE3")
  missing_fill <- unique(tiles$fill[is.na(value)])
  expect_identical(missing_fill, "#999999")
  expect_false(missing_fill %in% c(fill_of(-0.5), fill_of(0), fill_of(0.5)))
})

test_that("parallel community runs reproduce the serial fit exactly", {
  skip_if_not_installed("igraph")
  skip_on_os("windows")
  h <- group_hypergraph(utils::head(icsid_tribunals, 300),
                        node = "arbitrator", hyperedge = "case")
  serial <- hg_communities(h, n_runs = 6, trials = 5)
  forked <- hg_communities(h, n_runs = 6, trials = 5, parallel = TRUE,
                           n_cores = 2)
  expect_identical(forked, serial)
  serial <- suppressWarnings(hg_communities(h, type = "irmm", n_runs = 6))
  forked <- suppressWarnings(hg_communities(h, type = "irmm", n_runs = 6,
                                            parallel = TRUE, n_cores = 2))
  expect_identical(forked, serial)
})

test_that("parallel and n_cores are validated for both community types", {
  skip_if_not_installed("igraph")
  h <- group_hypergraph(utils::head(icsid_tribunals, 60),
                        node = "arbitrator", hyperedge = "case")
  expect_error(hg_communities(h, n_runs = 2, parallel = NA),
               class = "hypergraphs_bad_input")
  expect_error(hg_communities(h, n_runs = 2, parallel = TRUE, n_cores = 0),
               class = "hypergraphs_bad_input")
  expect_error(hg_communities(h, type = "irmm", n_runs = 2, parallel = "yes"),
               class = "hypergraphs_bad_input")
})

test_that("IRMM's sparse passes give the dense clique reduction exactly", {
  h <- group_hypergraph(utils::head(icsid_tribunals, 300),
                        node = "arbitrator", hyperedge = "case")
  dense <- .hg_mod_parts(h)
  sparse <- dense
  sparse$pattern <- methods::as(methods::as(dense$pattern, "CsparseMatrix"),
                                "generalMatrix")
  w <- seq(0.5, 2, length.out = length(dense$weights))
  expect_identical(unname(as.matrix(.hg_irmm_two_section(sparse, w))),
                   unname(.hg_irmm_two_section(dense, w)))
})
