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

test_that("node order and the results built on it do not depend on the locale", {
  d <- data.frame(node = c("b", "B", "a", "A"), hyperedge = c(1, 1, 2, 2))
  e <- rbind(Q = c(1, 0), a = c(0, 1), B = c(0, 1), c = c(1, 1))
  results <- function() {
    g <- group_hypergraph(d, node = "node", hyperedge = "hyperedge")
    list(nodes = g$nodes,
         modularity = hg_modularity(g, partition = c(1, 1, 2, 2)),
         knn = hg_get(knn_hypergraph(e, k = 1, weight = "binary"),
                      what = "memberships"),
         sorted = .ho_sort(c("b", "B", "a", "A")))
  }
  in_c <- results()
  expect_identical(in_c$nodes, c("A", "B", "a", "b"))
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
  skip_if(identical(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"), ""),
          "no en_US.UTF-8 locale")
  expect_identical(results(), in_c)
})

test_that("an object that is not sequences is refused, not read as sequences", {
  hg <- group_hypergraph(data.frame(n = c("a", "b", "c"), e = c("x", "x", "y")),
                         node = "n", hyperedge = "e")
  # a hypergraph is a list: its fields used to be read as sequences, so
  # hon() built a network over states such as "group_hypergraph" and "1:2"
  expect_error(.ho_sequence_input(hg), class = "hypergraphs_bad_input")
  expect_error(hon(hg), class = "hypergraphs_bad_input")
  expect_error(mogen(hg), class = "hypergraphs_bad_input")
  expect_error(markov_order(hg), class = "hypergraphs_bad_input")
  expect_error(hg_bootstrap(hg), class = "hypergraphs_bad_input")
  expect_error(window_hypergraph(hg), class = "hypergraphs_bad_input")
  expect_error(memory(hg), class = "hypergraphs_bad_input")
  # plain lists and the model objects that carry sequences still pass
  seqs <- list(c("a", "b", "a"), c("b", "a", "b"))
  expect_identical(.ho_sequence_input(seqs), seqs)
  expect_identical(.ho_sequence_input(structure(seqs, class = "list")),
                   structure(seqs, class = "list"))
  net <- Nestimate::build_network(.ho_wide_sequences(seqs), method = "relative")
  expect_identical(.ho_sequence_input(net), net)
})

test_that("the long route checks time_threshold with a classed error", {
  long <- data.frame(action = c("a", "b", "a"), actor = "u", time = 1:3)
  expect_error(.ho_sequence_input(long, actor = "actor", time_threshold = -1),
               class = "hypergraphs_bad_input")
  expect_error(.ho_sequence_input(long, actor = "actor", time_threshold = "a"),
               class = "hypergraphs_bad_input")
  expect_error(.ho_sequence_input(long, actor = "actor",
                                  time_threshold = c(1, 2)),
               class = "hypergraphs_bad_input")
})

test_that(".ho_local_seed refuses a seed that set.seed() would bend", {
  draw <- function(seed) {
    .ho_local_seed(seed)
    stats::runif(1L)
  }
  # set.seed(1.7) is set.seed(1) and set.seed(c(1, 2)) is set.seed(1)
  lapply(list(1.7, c(1, 2), "1", NA, Inf, 1e12, factor(1)), \(s) {
    expect_error(draw(s), class = "hypergraphs_bad_input")
  })
  expect_identical(draw(1), draw(1L))
  expect_identical(draw(-5), draw(-5L))
  set.seed(9)
  before <- .Random.seed
  draw(3)
  expect_identical(.Random.seed, before)
})

test_that("a long table found by its column names must name its actor", {
  # build_network() never detects the actor, so this table used to become
  # one sequence of both actors' events interleaved, with only a message
  s <- data.frame(actor = rep(c("ana", "bo"), each = 4), time = rep(1:4, 2),
                  action = c("x", "x", "y", "y", "y", "x", "x", "y"))
  expect_error(.ho_sequence_input(s), class = "hypergraphs_long_format")
  expect_error(hon(s), class = "hypergraphs_long_format")
  expect_error(hon(s), class = "hypergraphs_bad_input")
  # naming the actor gives one sequence per actor, as build_network() does
  ours <- .ho_sequence_input(s, actor = "actor", time_threshold = Inf)
  expect_identical(ours, Nestimate::build_network(
    s, method = "relative", actor = "actor", time_threshold = Inf)$data)
  expect_identical(nrow(ours), 2L)
  # naming action is the deliberate one-sequence reading
  expect_message(one <- .ho_sequence_input(s, action = "action", time = "time",
                                           time_threshold = Inf),
                 class = "hypergraphs_single_sequence")
  expect_identical(nrow(one), 1L)
  # a session column detected as the session is not an undeclared actor
  sess <- data.frame(session_id = rep(c("s1", "s2"), each = 3),
                     action = c("a", "b", "a", "b", "a", "b"), time = 1:6)
  expect_no_error(suppressMessages(.ho_sequence_input(sess)))
})
