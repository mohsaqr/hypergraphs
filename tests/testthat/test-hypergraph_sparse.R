# Sparse core: parity with the dense engines (the in-package oracle chain:
# dense == HyperNetX), plus scale behavior and guards.

sparse_corpus <- c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
)
sparse_stop <- c("the", "with", "and", "a", "this", "at", "on", "all")

both <- function(...) {
  list(
    dense = text_hypergraph(sparse_corpus, stop_words = sparse_stop, ...),
    sparse = text_hypergraph(sparse_corpus, stop_words = sparse_stop,
                             sparse = TRUE, ...)
  )
}

test_that("sparse construction reproduces the dense incidence exactly", {
  for (w in c("n", "tfidf")) {
    for (nd in c("doc", "word")) {
      hg <- both(weight = w, nodes = nd)
      expect_s4_class(hg$sparse$incidence, "sparseMatrix")
      expect_identical(as.matrix(hg$sparse$incidence), hg$dense$incidence)
      expect_identical(hg$sparse$nodes, hg$dense$nodes)
      expect_identical(hg$sparse$size_distribution,
                       hg$dense$size_distribution)
      sparse_table <- hg_get(hg$sparse)
      dense_table <- hg_get(hg$dense)
      expect_identical(sparse_table, dense_table)
    }
  }
})

test_that("sparse measures match the dense engine", {
  hg <- both(weight = "tfidf")
  for (what in c("nodes", "edges", "summary", "overlap")) {
    sparse_measures <- hg_measures(hg$sparse, what = what)
    dense_measures <- hg_measures(hg$dense, what = what)
    expect_equal(sparse_measures, dense_measures, tolerance = 1e-12)
  }
})

test_that("sparse pagerank matches dense to machine precision", {
  hg <- both(weight = "tfidf")
  for (d in c(0.85, 1)) {
    sparse_pr <- hg_pagerank(hg$sparse, damping = d, tol = 1e-14,
                             max_iter = 20000L)
    dense_pr <- hg_pagerank(hg$dense, damping = d, tol = 1e-14,
                            max_iter = 20000L)
    expect_equal(sparse_pr, dense_pr, tolerance = 1e-10)
  }
  sparse_personalized <- hg_pagerank(hg$sparse, personalized = c(cooking_1 = 1))
  dense_personalized <- hg_pagerank(hg$dense, personalized = c(cooking_1 = 1))
  expect_equal(sparse_personalized, dense_personalized, tolerance = 1e-10)
})

test_that("sparse transduction (CG) matches the dense closed form", {
  hg <- both(weight = "tfidf")
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  for (type in c("zhou", "random_walk")) {
    sparse_fit <- hg_classify(hg$sparse, labels = labels, type = type)
    dense_fit <- hg_classify(hg$dense, labels = labels, type = type)
    expect_identical(sparse_fit$predicted, dense_fit$predicted)
    expect_equal(sparse_fit$score, dense_fit$score, tolerance = 1e-8)
    expect_equal(sparse_fit$margin, dense_fit$margin, tolerance = 1e-8)
  }
})

test_that("sparse class_mass transduction matches the dense engine", {
  hg <- both(weight = "tfidf")
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  sparse_fit <- hg_classify(hg$sparse, labels = labels,
                            normalization = "class_mass")
  dense_fit <- hg_classify(hg$dense, labels = labels,
                           normalization = "class_mass")
  expect_identical(sparse_fit$predicted, dense_fit$predicted)
  expect_equal(sparse_fit$score, dense_fit$score, tolerance = 1e-8)
  expect_equal(sparse_fit$margin, dense_fit$margin, tolerance = 1e-8)
})

test_that("sparse clustering recovers the same planted partition", {
  events <- data.frame(
    person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
               "d", "e", "f", "c", "d"),
    meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
                "m4", "m4", "m4", "m5", "m5"),
    w = 1
  )
  dense_hg <- group_hypergraph(events, node = "person",
                                          hyperedge = "meeting", weight = "w")
  sparse_hg <- hypergraphs:::.thg_sparse_bipartite(
    events, node = "person", hyperedge = "meeting", weight = "w"
  )
  dense_cl <- hg_cluster(dense_hg, k = 2, seed = 1)
  sparse_cl <- hg_cluster(sparse_hg, k = 2, seed = 1)
  expect_identical(sparse_cl, dense_cl)

  fit <- hypergraphs:::.thg_sparse_cluster(
    sparse_hg, k = 2, type = "zhou", edge_weights = NULL, nstart = 25L,
    seed = 1
  )
  dense_fit <- .hg_cluster_fit(dense_hg, k = 2, seed = 1)
  expect_equal(fit$eigenvalues[seq_len(3)],
               dense_fit$eigenvalues[seq_len(3)], tolerance = 1e-8)
  expect_equal(fit$pi, dense_fit$pi, tolerance = 1e-10)
})

test_that("sparse duals round-trip and stay sparse", {
  hg <- both(weight = "tfidf")
  dual <- dual_hypergraph(hg$sparse)
  expect_s4_class(dual$incidence, "sparseMatrix")
  expect_identical(as.matrix(dual$incidence), t(hg$dense$incidence))
  back <- dual_hypergraph(dual)
  expect_identical(as.matrix(back$incidence), hg$dense$incidence)
})

test_that("unsupported sparse paths refuse with classed errors", {
  hg <- both()
  expect_error(hg_centrality(hg$sparse), class = "hypergraphs_sparse_unsupported")
  expect_error(hg_null_test(hg$sparse), class = "hypergraphs_sparse_unsupported")
  expect_error(
    text_hypergraph(sparse_corpus, construction = "window", sparse = TRUE),
    class = "hypergraphs_bad_input"
  )
  big_long <- data.frame(v = rep(c("x", "y"), each = 2100),
                         e = paste0("e", c(seq_len(2100), seq_len(2100))),
                         w = 1)
  big <- hypergraphs:::.thg_sparse_bipartite(big_long, node = "v",
                                                hyperedge = "e", weight = "w")
  expect_error(hg_measures(big, what = "overlap"),
               class = "hypergraphs_sparse_too_large")
})

test_that("sparse scale: thousands of documents classify in seconds", {
  skip_on_cran()
  set.seed(42)
  # letters-only tokens: the tokenizer treats digits as separators
  vocab <- head(apply(expand.grid(letters, letters, letters), 1, paste0,
                      collapse = ""), 1405)
  # two disjoint vocabulary blocks; only the first 100 documents carry two
  # tokens from a tiny shared bridge pool, so the corpus is connected but
  # the bridge hyperedges stay small (the zhou Laplacian is binary -- a
  # bridge word present in every document would wash all labels into ties)
  make_doc <- \(i) {
    pool <- if (i %% 2L == 0L) seq_len(700) else 701:1400
    body <- sample(vocab[pool], 30, replace = TRUE)
    if (i <= 100L) {
      body <- c(body, sample(vocab[1401:1405], 2))
    }
    paste(body, collapse = " ")
  }
  docs <- vapply(seq_len(3000), make_doc, character(1))
  names(docs) <- sprintf("d%04d", seq_len(3000))
  elapsed <- system.time({
    hg <- text_hypergraph(docs, weight = "tfidf", sparse = TRUE)
    labels <- c(d0002 = "even", d0004 = "even", d0006 = "even",
                d0001 = "odd", d0003 = "odd", d0005 = "odd")
    fit <- hg_classify(hg, labels = labels, type = "zhou")
    pr <- hg_pagerank(hg)
  })[["elapsed"]]
  expect_identical(nrow(fit), 3000L)
  expect_identical(nrow(pr), 3000L)
  expect_lt(elapsed, 60)
  # the two vocabulary populations are recovered from one label each
  truth <- rep(c("odd", "even"), length.out = 3000)
  accuracy <- mean(fit$predicted == truth)
  expect_gt(accuracy, 0.95)
})

# ---- R11: sparse measures use the dense engine's definitions ------------

.sparse_general <- function(m) {
  methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix")
}

.measure_fixtures <- function() {
  inc <- function(cols, nodes) {
    m <- vapply(cols, \(members) as.numeric(nodes %in% members),
                numeric(length(nodes)))
    m <- matrix(m, nrow = length(nodes),
                dimnames = list(nodes, names(cols) %||% character(0)))
    m
  }
  list(
    # the audit fixture: {a,b}, {b,c} plus the isolate z
    uniform_isolate = inc(list(e1 = c("a", "b"), e2 = c("b", "c")),
                          c("a", "b", "c", "z")),
    uniform_3 = inc(list(e1 = c("a", "b", "c"), e2 = c("b", "c", "d"),
                         e3 = c("a", "c", "d")), c("a", "b", "c", "d")),
    mixed = inc(list(e1 = c("a", "b", "c"), e2 = c("c", "d"),
                     e3 = c("d")), c("a", "b", "c", "d", "z")),
    with_empty_edge = inc(list(e1 = c("a", "b"), e2 = character(0),
                               e3 = c("b", "c", "d")), c("a", "b", "c", "d")),
    singleton = inc(list(e1 = "a"), "a"),
    no_edges = inc(list(), c("a", "b", "c"))
  )
}

test_that("REGRESSION R11: sparse measures equal dense measures on every shape", {
  for (nm in names(.measure_fixtures())) {
    m <- .measure_fixtures()[[nm]]
    dense <- hypergraphs:::.thg_from_incidence(m, list())
    sparse <- hypergraphs:::.thg_from_incidence(.sparse_general(m), list())
    for (what in c("nodes", "edges", "overlap", "summary")) {
      expect_equal(hg_measures(sparse, what = what),
                   hg_measures(dense, what = what),
                   tolerance = 1e-12, info = paste(nm, what))
    }
  }
})

test_that("REGRESSION R11: the documented values on the audit fixture", {
  m <- .measure_fixtures()$uniform_isolate
  sparse <- hypergraphs:::.thg_from_incidence(.sparse_general(m), list())
  summary_tab <- hg_measures(sparse, what = "summary")
  value <- stats::setNames(summary_tab$value, summary_tab$measure)
  # 2 of the choose(4, 2) = 6 vertex pairs co-occur: ab, bc
  expect_equal(value[["pairwise_participation"]], 2 / 6)
  # 2-uniform: m / choose(n, 2)
  expect_equal(value[["density"]], 2 / choose(4, 2))
  nodes <- hg_measures(sparse, what = "nodes")
  expect_identical(nodes$max_edge_size, c(2L, 2L, 2L, 0L))
  expect_identical(nodes$n_neighbors, c(1L, 2L, 1L, 0L))
  no_edges <- hypergraphs:::.thg_from_incidence(
    .sparse_general(.measure_fixtures()$no_edges), list()
  )
  empty_summary <- hg_measures(no_edges, what = "summary")
  expect_identical(empty_summary$value[empty_summary$measure == "density"], 0)
  expect_true(is.na(
    empty_summary$value[empty_summary$measure == "avg_edge_size"]
  ))
})

# ---- R25: sparse clustering honours the dense count contract ------------

test_that("REGRESSION R25: sparse hg_cluster refuses invalid k and nstart", {
  path <- data.frame(node = c("a", "b", "b", "c"),
                     edge = c("e1", "e1", "e2", "e2"))
  sparse <- group_hypergraph(path, node = "node", hyperedge = "edge",
                             sparse = TRUE)
  dense <- group_hypergraph(path, node = "node", hyperedge = "edge")
  bad_k <- list(1, 1.5, 3, Inf, NA_real_, -2, 3e9, c(2, 2))
  for (k in bad_k) {
    expect_error(hg_cluster(sparse, k = k), class = "hypergraphs_bad_input",
                 info = paste(k, collapse = ","))
  }
  for (nstart in list(0, 1.5, NA_real_, Inf, 3e9)) {
    expect_error(hg_cluster(sparse, k = 2, nstart = nstart),
                 class = "hypergraphs_bad_input", info = nstart)
  }
  # the dense engine refuses the same controls with the same class
  for (k in bad_k) {
    expect_error(hg_cluster(dense, k = k), class = "hypergraphs_bad_input",
                 info = paste(k, collapse = ","))
  }
  for (nstart in list(0, 1.5, NA_real_, Inf, 3e9)) {
    expect_error(hg_cluster(dense, k = 2, nstart = nstart),
                 class = "hypergraphs_bad_input", info = nstart)
  }
  # the smallest valid sparse case still clusters, as the dense one does
  expect_identical(hg_cluster(sparse, k = 2, seed = 1),
                   hg_cluster(dense, k = 2, seed = 1))
  expect_identical(nrow(hg_cluster(sparse, k = 2, nstart = 1, seed = 1)), 3L)
})

# ---- A04: empty hyperedges contribute nothing to the sparse walk ---------

test_that("REGRESSION A04: sparse PageRank ignores an empty hyperedge", {
  # seed 4 draws one empty column (h2) on a connected support
  h <- random_hypergraph("gnp", n = 3, m = 4, p = 0.5, seed = 4)
  expect_identical(unname(colSums(h$incidence))[2L], 0)
  sparse <- h
  sparse$incidence <- .sparse_general(h$incidence)
  without_empty <- hg_subset(h, edges = c("h1", "h3", "h4"),
                             drop_isolated = FALSE)
  # the default dispersion weight of an empty column is 1, not 0/0
  expect_identical(hypergraphs:::.thg_sparse_edge_weights(sparse$incidence),
                   c(1, 1, 1, 1))
  for (damping in c(0.85, 1)) {
    sparse_pr <- hg_pagerank(sparse, damping = damping)
    expect_true(all(is.finite(sparse_pr$pagerank)))
    expect_equal(sparse_pr, hg_pagerank(h, damping = damping),
                 tolerance = 1e-10)
    expect_equal(sparse_pr, hg_pagerank(without_empty, damping = damping),
                 tolerance = 1e-10)
  }
  # the Zhou similarity operator behind sparse transduction, likewise
  labels <- c(V1 = "x", V3 = "y")
  expect_equal(hg_classify(sparse, labels = labels),
               hg_classify(h, labels = labels), tolerance = 1e-8)
})

test_that("SymNMF on a sparse hypergraph matches the dense engine", {
  # repeated words and uneven document lengths, so the random-walk
  # (edge-dependent weights) and Zhou similarities genuinely differ
  docs <- c(
    d1 = "solar solar solar wind energy grid night",
    d2 = "solar panels energy energy night",
    d3 = "wind wind turbines grid energy night night",
    d4 = "coal plant emissions emissions emissions night",
    d5 = "coal coal mining emissions plant night",
    d6 = "gas plant emissions night",
    d7 = "solar wind turbines turbines night",
    d8 = "coal gas gas gas emissions night night")
  dense <- text_hypergraph(docs)
  sparse <- text_hypergraph(docs, sparse = TRUE)
  expect_true(.thg_is_sparse(sparse))
  similarity <- \(type) as.matrix(Matrix::Diagonal(length(docs)) -
                                   .hl_build(sparse, type, NULL)$L)
  expect_gt(max(abs(similarity("zhou") - similarity("random_walk"))), 0.01)
  for (type in c("zhou", "random_walk")) {
    a <- hg_cluster(dense, k = 2, type = type, algorithm = "symnmf",
                    seed = 3, nstart = 3, what = "membership")
    b <- hg_cluster(sparse, k = 2, type = type, algorithm = "symnmf",
                    seed = 3, nstart = 3, what = "membership")
    expect_equal(b, a, tolerance = 1e-6)
    expect_equal(as.numeric(tapply(b$membership, b$node, sum)),
                 rep(1, length(docs)), tolerance = 1e-12)
  }
  # leading eigenvalues agree with the dense spectrum
  ev_dense <- hg_cluster(dense, k = 2, algorithm = "symnmf", seed = 3,
                         what = "eigenvalues")
  ev_sparse <- hg_cluster(sparse, k = 2, algorithm = "symnmf", seed = 3,
                          what = "eigenvalues")
  expect_equal(ev_sparse$value, utils::head(ev_dense$value, 3),
               tolerance = 1e-8)
})

test_that("sparse SymNMF refuses a similarity too large to hold", {
  sparse <- text_hypergraph(c(a = "x y night", b = "y z night",
                              c = "z x night"), sparse = TRUE)
  expect_error(.thg_sparse_symnmf(sparse, k = 2, type = "zhou",
                                  edge_weights = NULL, nstart = 1, seed = 1,
                                  max_iter = 10, tol = 1e-6, max_nodes = 2),
               class = "hypergraphs_sparse_too_large")
})

test_that("hg_cluster runs SymNMF starts in parallel with the serial result", {
  skip_on_os("windows")
  docs <- c(d1 = "solar wind energy night", d2 = "solar panels energy night",
            d3 = "coal plant emissions night", d4 = "coal gas emissions night",
            d5 = "wind turbines grid night", d6 = "gas plant mining night")
  for (sparse in c(FALSE, TRUE)) {
    hg <- text_hypergraph(docs, sparse = sparse)
    serial <- hg_cluster(hg, k = 2, algorithm = "symnmf", seed = 4,
                         nstart = 3, what = "membership")
    forked <- suppressWarnings(  # a macOS fork may fall back to serial
      hg_cluster(hg, k = 2, algorithm = "symnmf", seed = 4, nstart = 3,
                 what = "membership", parallel = TRUE, n_cores = 2))
    expect_identical(forked, serial)
  }
  expect_error(hg_cluster(text_hypergraph(docs), k = 2, parallel = TRUE),
               class = "hypergraphs_bad_input")
})
