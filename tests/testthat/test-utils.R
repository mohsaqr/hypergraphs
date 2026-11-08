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

test_that(".ho_kmeans gives the hand partition of two separated pairs", {
  # (0, 0), (0, 1) and (10, 10), (10, 11): centres are the pair means, each
  # pair is 0.5 from its centre in squared distance, and the grand mean
  # (5, 5.5) is 201 away in total
  x <- rbind(c(0, 0), c(0, 1), c(10, 10), c(10, 11))
  set.seed(1)
  fit <- .ho_kmeans(x, 2, nstart = 1)
  expect_identical(fit$cluster[1], fit$cluster[2])
  expect_identical(fit$cluster[3], fit$cluster[4])
  expect_false(identical(fit$cluster[1], fit$cluster[3]))
  centre <- rbind(c(0, 0.5), c(10, 10.5))
  expect_equal(unname(fit$centers[fit$cluster[c(1, 3)], ]), centre)
  expect_equal(fit$withinss, c(0.5, 0.5))
  expect_equal(fit$totss, 201)
  expect_equal(fit$betweenss, 200)
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

test_that("a long table declares its roles", {
  s <- data.frame(actor = rep(c("ana", "bo"), each = 4), time = rep(1:4, 2),
                  action = c("x", "x", "y", "y", "y", "x", "x", "y"))
  expect_error(.ho_sequence_input(s), class = "hypergraphs_long_format")
  expect_error(hon(s), class = "hypergraphs_long_format")
  expect_error(hon(s), class = "hypergraphs_bad_input")
  # a column called action alone is enough to refuse
  expect_error(.ho_sequence_input(data.frame(action = c("a", "b"), T2 = "c")),
               class = "hypergraphs_long_format")
  # other roles mean nothing without action
  expect_error(.ho_sequence_input(s, actor = "actor"),
               class = "hypergraphs_bad_input")
  # declaring the roles gives one sequence per actor
  ours <- .ho_sequence_input(s, action = "action", actor = "actor",
                             time = "time", time_threshold = Inf)
  expect_identical(nrow(ours), 2L)
  expect_identical(unlist(ours[1L, paste0("T", 1:4)], use.names = FALSE),
                   s$action[1:4])
  expect_identical(unlist(ours[2L, paste0("T", 1:4)], use.names = FALSE),
                   s$action[5:8])
  # the positions are the columns, in time order: ana x x y y, bo y x x y
  expect_identical(names(ours), paste0("T", 1:4))
  expect_identical(unname(as.matrix(ours)),
                   matrix(c("x", "x", "y", "y", "y", "x", "x", "y"), 2L, byrow = TRUE))
  # naming action alone is the deliberate one-sequence reading
  expect_message(one <- .ho_sequence_input(s, action = "action", time = "time",
                                           time_threshold = Inf),
                 class = "hypergraphs_single_sequence")
  expect_identical(nrow(one), 1L)
})

test_that("an undeclared time or session column beside action is named, not used", {
  sess <- data.frame(session_id = rep(c("s1", "s2"), each = 3),
                     action = c("a", "b", "c", "c", "a", "b"), time = 6:1,
                     who = "p")
  expect_message(
    undeclared <- .ho_sequence_input(sess, action = "action", actor = "who"),
    class = "hypergraphs_column_role")
  # a column declared as the actor is in use, whatever it is called
  expect_no_message(
    .ho_sequence_input(sess, action = "action", actor = "session_id",
                       time = "time", time_threshold = Inf))
  # row order, one sequence
  expect_identical(nrow(undeclared), 1L)
  expect_identical(unlist(undeclared[paste0("T", 1:6)], use.names = FALSE),
                   sess$action)
  # declared, nothing is announced and the sessions split
  expect_no_message(
    declared <- .ho_sequence_input(sess, action = "action", actor = "who",
                                   session = "session_id", time = "time",
                                   time_threshold = Inf))
  expect_identical(nrow(declared), 2L)
  # within each session, by time: s1 at 4, 5, 6 is c b a; s2 at 1, 2, 3 is b a c
  expect_identical(unlist(declared[1L, paste0("T", 1:3)], use.names = FALSE),
                   c("c", "b", "a"))
  expect_identical(unlist(declared[2L, paste0("T", 1:3)], use.names = FALSE),
                   c("b", "a", "c"))
  # the same by position: the columns are T1, T2, T3 in that order
  expect_identical(names(declared), paste0("T", 1:3))
  expect_identical(unname(as.matrix(declared)),
                   matrix(c("c", "b", "a", "b", "a", "c"), 2L, byrow = TRUE))
  expect_identical(names(undeclared), paste0("T", 1:6))
})

test_that("a column declared under any role is never called unused", {
  d <- data.frame(session_id = c(1, 2, 3, 4), Time = 1:4,
                  who = "p", code = c("a", "b", "a", "b"))
  says_unused <- function(...) {
    said <- character()
    withCallingHandlers(
      .ho_sequence_input(d, action = "code", time_threshold = Inf, ...),
      hypergraphs_column_role = function(m) {
        said <<- c(said, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    said
  }
  # session_id declared as the time: only the other clock is unused
  said <- says_unused(actor = "who", time = "session_id")
  expect_length(said, 0L)
  # nothing declared but the action: both are named
  said <- says_unused()
  expect_match(said, "`Time`", fixed = TRUE)
  expect_match(said, "`session_id`", fixed = TRUE)
  expect_no_message(.ho_sequence_input(d, action = "code", actor = "session_id",
                                       time = "Time", time_threshold = Inf))
})

test_that("the one-sequence notice follows what was built, not the arguments", {
  d <- data.frame(session_id = c("s1", "s1", "s2", "s2"), code = c("a", "b", "b", "a"),
                  order = 1:4)
  # two sessions with no actor: two sequences, no notice
  built <- expect_no_message(
    .ho_sequence_input(d, action = "code", session = "session_id", time = "order",
                       time_threshold = Inf),
    class = "hypergraphs_single_sequence")
  expect_identical(nrow(built), 2L)
  # one actor is one sequence: the notice is true
  expect_message(.ho_sequence_input(transform(d, who = "p"), action = "code",
                                    actor = "who", time = "order",
                                    time_threshold = Inf),
                 class = "hypergraphs_single_sequence")
})

test_that("a frame that is not a run of positions is refused, not read", {
  guard <- hypergraphs:::.ho_stop_undeclared_frame
  expect_error(guard(data.frame(first = c("a", "b"), second = c("b", "a"))),
               class = "hypergraphs_long_format")
  # a run with something else in it would read the something else as states
  expect_error(guard(data.frame(T1 = "a", T2 = "b", outcome = "pass")),
               class = "hypergraphs_long_format")
  expect_null(guard(data.frame(T1 = "a", T2 = "b")))
  expect_null(guard(data.frame(`1` = "a", `2` = "b", check.names = FALSE)))
  # 0/1 could be one-hot codes or positions; no reading is guessed
  expect_error(guard(data.frame(reading = c(0, 1), coding = c(1, 0))),
               class = "hypergraphs_long_format")
  expect_error(guard(data.frame(T1 = c(0, 1), T2 = c(1, 0))),
               class = "hypergraphs_long_format")
  expect_null(guard(data.frame(code = c("a", "b"))))
  # positions a Nestimate verb already resolved read back in
  months <- data.frame(Sep.93 = c("a", "b"), Oct.93 = c("b", "a"), Nov.93 = c("a", "a"))
  expect_error(guard(months), class = "hypergraphs_long_format")
  expect_null(guard(Nestimate::extract_sequences(months, cols = names(months))))
  expect_null(guard(list(c("a", "b"))))
  # bundled log, nothing declared, on a verb
  expect_error(hon(utils::head(Nestimate::human_long, 20L)),
               class = "hypergraphs_long_format")
})

test_that("the unused-role message says what holds for the call", {
  d <- data.frame(actor = c("p", "p", "q", "q"), code = c("B", "A", "B", "A"),
                  time = c(2, 1, 2, 1), session_id = c("s1", "s1", "s2", "s2"))
  says <- function(...) {
    said <- character()
    withCallingHandlers(
      .ho_sequence_input(d, action = "code", actor = "actor", time_threshold = Inf, ...),
      hypergraphs_column_role = function(m) {
        said <<- c(said, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    said
  }
  # time declared: events follow it, so the message must not say row order
  with_time <- says(time = "time")
  expect_match(with_time, "`session_id` not used", fixed = TRUE)
  expect_false(grepl("row order", with_time, fixed = TRUE))
  # session declared: the message must not say that no session column splits them
  with_session <- says(session = "session_id")
  expect_match(with_session, "`time` not used", fixed = TRUE)
  expect_false(grepl("no session column", with_session, fixed = TRUE))
  expect_length(says(time = "time", session = "session_id"), 0L)
})

test_that("a declared wide table reaches every estimator, as the refusal advises", {
  weeks <- data.frame(week_a = c("a", "b", "a", "b"), week_b = c("b", "c", "c", "a"),
                      week_c = c("c", "a", "b", "c"), week_d = c("a", "b", "c", "a"))
  expect_error(hon(weeks), class = "hypergraphs_long_format")
  # the refusal's advice: the positions alone, as a matrix or a list
  numbered0 <- weeks
  names(numbered0) <- paste0("T", 1:4)
  expect_identical(hon(as.matrix(weeks), max_order = 2L)$n_edges,
                   hon(numbered0, max_order = 2L)$n_edges)
  expect_s3_class(hon(unname(as.list(as.data.frame(t(as.matrix(weeks))))), max_order = 2L),
                  "net_hon")
  net <- Nestimate::build_network(weeks, method = "relative",
                                  state_cols = names(weeks))
  numbered <- weeks
  names(numbered) <- paste0("T", 1:4)
  expect_identical(hon(net, max_order = 2L)$n_edges,
                   hon(numbered, max_order = 2L)$n_edges)
  expect_s3_class(mogen(net, max_order = 1L), "net_mogen")
  expect_s3_class(hypa(net, min_count = 1L), "net_hypa")
  expect_s3_class(markov_order(net, max_order = 1L, n_perm = 2L, seed = 1L),
                  "net_markov_order")
  expect_s3_class(memory(net, min_count = 1L), "net_path_dependence")
  # positions a Nestimate verb resolved read back in, a column added to them does not
  action_named <- weeks
  names(action_named)[1L] <- "action"
  resolved <- Nestimate::extract_sequences(action_named, cols = names(action_named))
  expect_s3_class(hon(resolved, max_order = 2L), "net_hon")
  resolved$grade <- "p"
  expect_error(hon(resolved, max_order = 2L), class = "hypergraphs_long_format")
})

test_that("a character matrix of positions can be split by group", {
  weeks <- matrix(c("a", "b", "c", "a", "b", "c", "a", "b",
                    "b", "c", "a", "b", "c", "a", "b", "c"), 4L,
                  dimnames = list(NULL, c("wk_a", "wk_b", "wk_c", "wk_d")))
  grouped <- hon(weeks, group = c("g1", "g1", "g2", "g2"), max_order = 2L)
  as_rows <- lapply(seq_len(nrow(weeks)), function(i) unname(weeks[i, ]))
  expect_identical(names(grouped), c("g1", "g2"))
  expect_identical(grouped[["g1"]]$n_edges,
                   hon(as_rows[1:2], max_order = 2L)$n_edges)
})

test_that("window columns are one run of positions, and one-hot networks are not read", {
  windows <- data.frame(W0_T1 = c("a", "b"), W0_T2 = c("b", "a"),
                        W1_T1 = c("a", "a"), W1_T2 = c("b", "b"))
  expect_null(hypergraphs:::.ho_stop_undeclared_frame(windows))
  expect_s3_class(hon(windows, max_order = 1L), "net_hon")
  hot <- Nestimate::build_network(
    data.frame(A = rep(c(1, 0, 0, 1), 6L), B = rep(c(0, 1, 0, 0), 6L),
               C = rep(c(0, 0, 1, 1), 6L), who = rep(paste0("p", 1:6), each = 4L)),
    method = "frequency", format = "onehot", actor = "who", window_size = 2)
  invisible(lapply(list(hon, mogen, hypa, markov_order, memory), function(read) {
    expect_error(read(hot), class = "hypergraphs_bad_input")
  }))
})

