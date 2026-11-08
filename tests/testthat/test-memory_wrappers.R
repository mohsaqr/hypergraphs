testthat::skip_on_cran()

# ---- Memory-network wrappers: one answer on every input form ----
#
# hon(), mogen(), hypa() on sequences, markov_order(), memory() and
# hg_markov_stability() coerce the input and call the estimator. The tests below
# assert that every input form the wrappers accept gives the same model: a list,
# a wide frame, a long event table, a netobject and a tna model. Agreement with
# the direct Nestimate call is local only:
# local_testing_and_equivalence/test-equiv-wrappers-nestimate.R.

.mw_fixture <- function() {
  set.seed(3)
  seqs <- lapply(seq_len(12L), function(i) {
    sample(c("a", "b", "c", "d"), sample(6:10, 1L), replace = TRUE)
  })
  names(seqs) <- sprintf("u%02d", seq_len(12L))
  long <- do.call(rbind, lapply(names(seqs), function(n) {
    data.frame(who = n, when = seq_along(seqs[[n]]) * 10, what = seqs[[n]],
               stringsAsFactors = FALSE)
  }))
  long <- long[sample(nrow(long)), ]          # row order must not matter
  rownames(long) <- NULL
  wide <- .ho_wide_sequences(seqs)
  net <- Nestimate::build_network(long, method = "relative", actor = "who",
                                  action = "what", time = "when",
                                  time_threshold = Inf)
  list(seqs = seqs, long = long, wide = wide, net = net)
}

test_that("the long route builds one sequence per declared unit", {
  same <- function(data, ...) nrow(.ho_sequence_input(data, ...))
  sessions <- human_long
  n <- c(
    position = same(sessions, action = "code", actor = "session_id",
                    time = "order_in_session"),
    timestamp = same(sessions, action = "code", actor = "session_id",
                     time = "timestamp"),
    no_gaps = same(sessions, action = "code", actor = "session_id",
                   time = "timestamp", time_threshold = FALSE),
    project = same(sessions, action = "code", actor = "project",
                   time = "order_in_session"),
    project_session = same(sessions, action = "code", actor = "project",
                           session = "session_id",
                           time = "order_in_session"),
    two_actors = same(sessions, action = "code",
                      actor = c("project", "session_id"),
                      time = "order_in_session")
  )
  # a gap of more than 900 seconds starts a new sequence; a session column
  # counts only when it is declared, so `project` alone joins its sessions
  expect_identical(unname(n), c(429L, 526L, 429L, 34L, 429L, 429L))
  # positions are columns T1, T2, ... in time order: the first session in the
  # data, 0086cabebd15, opens Specify, Command, Specify, Interrupt, Verify,
  # Specify (order_in_session 1, 2, 3, 4, 7, 10)
  by_session <- .ho_sequence_input(sessions, action = "code", actor = "session_id",
                                   time = "order_in_session")
  expect_identical(names(by_session), paste0("T", seq_len(ncol(by_session))))
  expect_identical(unlist(by_session[1L, 1:6], use.names = FALSE),
                   c("Specify", "Command", "Specify", "Interrupt", "Verify",
                     "Specify"))
  # every sequence starts at T1: no row is empty there
  expect_false(anyNA(by_session[[1L]]))
  # names give no role: an `action` column has to be declared
  renamed <- human_long
  names(renamed)[names(renamed) == "code"] <- "action"
  expect_error(.ho_sequence_input(renamed), class = "hypergraphs_long_format")
  expect_identical(.ho_sequence_input(renamed, action = "action",
                                      actor = "session_id",
                                      time = "order_in_session"),
                   .ho_sequence_input(human_long, action = "code",
                                      actor = "session_id",
                                      time = "order_in_session"))
  fx <- .mw_fixture()
  expect_identical(same(fx$long, action = "what", actor = "who",
                        time = "when"), 12L)
})

test_that("session = FALSE means no session, the same as NULL", {
  by_project <- .ho_sequence_input(human_long, action = "code",
                                   actor = "project",
                                   time = "order_in_session",
                                   session = FALSE)
  # without the session column the call is the same
  without <- human_long[setdiff(names(human_long), "session_id")]
  expect_identical(by_project,
                   .ho_sequence_input(without, action = "code", actor = "project",
                                      time = "order_in_session"))
  expect_identical(nrow(by_project), length(unique(human_long$project)))
  # positions are columns T1, T2, ... and every project's sequence starts at T1
  expect_identical(names(by_project), paste0("T", seq_len(ncol(by_project))))
  expect_false(anyNA(by_project[[1L]]))
})

test_that("missing actors and sessions are refused, by class", {
  fx <- .mw_fixture()
  long <- fx$long
  long$who[3L] <- NA
  expect_error(hon(long, action = "what", actor = "who", time = "when"),
               class = "hypergraphs_bad_input")
  expect_error(hon(fx$long, action = "what", actor = "who", time = "when",
                   time_threshold = -1))
})

test_that("hon() gives one model on every input form", {
  fx <- .mw_fixture()
  ref <- hon(fx$seqs, max_order = 2L)
  expect_s3_class(ref, "net_hon")
  expect_identical(hon(fx$wide, max_order = 2L), ref)
  expect_identical(hon(fx$long, action = "what", actor = "who", time = "when",
                       max_order = 2L), ref)
  expect_identical(hon(fx$net, max_order = 2L), ref)
  # no state repeats itself once repeats are collapsed
  collapsed <- hon(fx$seqs, max_order = 3L, min_freq = 2L,
                   collapse_repeats = TRUE, method = "hon")
  expect_s3_class(collapsed, "net_hon")
  expect_false(any(grepl("\\b([a-d]) -> \\1\\b",
                         hg_get(collapsed, what = "pathways")$pathway, perl = TRUE)))
  # x's next state depends on what came before it: second-order paths, counted by hand
  memory_seqs <- c(rep(list(c("a", "x", "c")), 20L), rep(list(c("b", "x", "d")), 20L))
  edges <- hg_get(hon(memory_seqs, max_order = 2L))
  expect_setequal(edges$path, c("a -> x", "a -> x -> c", "b -> x", "b -> x -> d",
                                "x -> c", "x -> d"))
  expect_equal(edges$count, rep(20, 6L))
  expect_setequal(hg_get(hon(memory_seqs, max_order = 1L))$path,
                  c("a -> x", "b -> x", "x -> c", "x -> d"))
})

# A tna model as the tna package stores it: integer-coded sequences with
# their labels (built here so the test needs no extra package).
.mw_tna <- function(wide) {
  labels <- sort(unique(stats::na.omit(unlist(wide, use.names = FALSE))))
  codes <- vapply(wide, function(col) match(col, labels), integer(nrow(wide)))
  codes <- matrix(codes, nrow(wide), dimnames = list(NULL, names(wide)))
  attr(codes, "labels") <- labels
  structure(list(data = codes, labels = labels), class = "tna")
}

test_that("hon() reads a tna model", {
  fx <- .mw_fixture()
  tm <- .mw_tna(fx$wide)
  from_tna <- hon(tm, max_order = 2L)
  expect_identical(from_tna, hon(fx$seqs, max_order = 2L))
  # the weights are transition probabilities: every row that leaves sums to 1
  expect_equal(unique(round(rowSums(from_tna$weights), 12)), 1)
  # a tna model of a model worked by hand: a x c and b x d, 20 times each
  memory_wide <- .ho_wide_sequences(c(rep(list(c("a", "x", "c")), 20L),
                                      rep(list(c("b", "x", "d")), 20L)))
  memory_hon <- hon(.mw_tna(memory_wide), max_order = 2L)
  expect_equal(memory_hon$weights["x", c("c", "d")], c(c = 0.5, d = 0.5))
  expect_equal(memory_hon$weights["a -> x", "c"], 1)
  expect_equal(memory_hon$weights["b -> x", "d"], 1)
  expect_equal(unname(rowSums(memory_hon$weights)), c(1, 1, 1, 1, 0, 0, 1))
  expect_equal(sort(memory_hon$edges$weight), c(0.5, 0.5, 1, 1, 1, 1))
})

test_that("honem() embeds a HON", {
  fx <- .mw_fixture()
  h <- hon(fx$seqs, max_order = 2L)
  embedding <- honem(h, dim = 3L, max_power = 5L)
  expect_s3_class(embedding, "net_honem")
  expect_identical(ncol(embedding$embeddings), 3L)
  # the columns are the singular vectors scaled by the root of the singular values
  expect_equal(unname(colSums(embedding$embeddings^2)),
               unname(embedding$singular_values[1:3]), tolerance = 1e-6)
  # each node keeps its own row: squared row lengths, which no sign flip of a
  # singular vector changes, pinned for this fixture
  expect_identical(rownames(embedding$embeddings),
                   c("a", "a -> d", "b", "b -> a", "b -> c", "b -> d", "c",
                     "c -> a", "d"))
  expect_equal(unname(rowSums(embedding$embeddings^2)),
               c(0.227463542349195, 0.503770192637955, 0.136863767515392,
                 0.243353693417727, 0.343849760076545, 0.238502081498088,
                 0.151407284367499, 0.465326842812449, 0.187074812542932),
               tolerance = 1e-10)
})

test_that("mogen() gives one model on every input form", {
  fx <- .mw_fixture()
  ref <- mogen(fx$seqs, max_order = 2L)
  expect_s3_class(ref, "net_mogen")
  expect_identical(mogen(fx$wide, max_order = 2L), ref)
  expect_identical(mogen(fx$long, action = "what", actor = "who", time = "when",
                         max_order = 2L), ref)
  expect_identical(mogen(fx$net, max_order = 2L), ref)
  by_bic <- mogen(fx$wide, max_order = 2L, criterion = "bic")
  expect_s3_class(by_bic, "net_mogen")
  expect_identical(by_bic$criterion, "bic")
  # pinned for this fixture: BIC keeps order 0
  expect_equal(unname(by_bic$log_likelihood),
               c(-134.75097465654, -131.627451439646, -103.266842399967),
               tolerance = 1e-10)
  expect_equal(unname(by_bic$bic),
               c(283.256851749092, 332.029415059351, 495.386635956179),
               tolerance = 1e-10)
  expect_identical(by_bic$optimal_order, 0L)
  # worked by hand: a x c and b x d, 20 times each, 120 observations.
  # order 0: a, b, c, d 20 each and x 40 of 120;
  # order 1: the first state at its order-0 rate 1/6, x then c or d at 1/2;
  # order 2: the first state at 1/6 and everything after it certain
  memory_seqs <- c(rep(list(c("a", "x", "c")), 20L), rep(list(c("b", "x", "d")), 20L))
  orders <- hg_get(mogen(memory_seqs, max_order = 2L))
  expect_equal(orders$log_likelihood,
               c(80 * log(1 / 6) + 40 * log(1 / 3),
                 40 * log(1 / 6) + 40 * log(1 / 2),
                 40 * log(1 / 6)))
  expect_equal(orders$aic, 2 * orders$df - 2 * orders$log_likelihood)
  expect_equal(orders$bic, orders$df * log(120) - 2 * orders$log_likelihood)
})

test_that("hypa() on sequences gives one model on every input form", {
  fx <- .mw_fixture()
  ref <- hypa(fx$seqs, min_count = 1L)
  expect_s3_class(ref, "net_hypa")
  expect_identical(hypa(fx$wide, min_count = 1L), ref)
  expect_identical(hypa(fx$long, action = "what", actor = "who", time = "when",
                        min_count = 1L), ref)
  expect_identical(hypa(fx$net, min_count = 1L), ref)
  both <- hypa(fx$wide, order = c(1L, 2L), min_count = 2L, p_adjust = "holm")
  expect_s3_class(both, "net_hypa")
  expect_setequal(unique(hg_get(both, what = "scores")$order), c(1L, 2L))
  # worked by hand at order 2: a x and b x each send 20 into x c and x d, which
  # each take 20, so all four routes expect 20 * 20 / 40 = 10
  memory_seqs <- c(rep(list(c("a", "x", "c")), 20L), rep(list(c("b", "x", "d")), 20L))
  scored <- hg_get(hypa(memory_seqs, order = c(1L, 2L), min_count = 2L,
                        p_adjust = "holm"), what = "scores")
  # order 1 scores only the four transitions seen, and expects each 20 times,
  # as often as it was seen (pinned)
  first <- scored[scored$order == 1L, , drop = FALSE]
  expect_setequal(first$path, c("a -> x", "b -> x", "x -> c", "x -> d"))
  expect_equal(first$expected, rep(20, 4L))
  second <- scored[scored$order == 2L, , drop = FALSE]
  expect_setequal(second$path, c("a -> x -> c", "a -> x -> d", "b -> x -> c",
                                 "b -> x -> d"))
  expect_equal(second$expected, rep(10, 4L))
  expect_equal(second$count[second$path %in% c("a -> x -> c", "b -> x -> d")],
               c(20, 20))
})

test_that("markov_order() gives one test on every input form", {
  fx <- .mw_fixture()
  ref <- markov_order(fx$seqs, max_order = 2L, n_perm = 20L, seed = 1L)
  expect_s3_class(ref, "net_markov_order")
  expect_identical(markov_order(fx$wide, max_order = 2L, n_perm = 20L, seed = 1L), ref)
  expect_identical(markov_order(fx$long, action = "what", actor = "who", time = "when",
                                max_order = 2L, n_perm = 20L, seed = 1L), ref)
  expect_identical(markov_order(fx$net, max_order = 2L, n_perm = 20L, seed = 1L), ref)
  # order 1: G2 of the table of consecutive pairs, worked from its counts
  memory_seqs <- c(rep(list(c("a", "x", "c")), 20L), rep(list(c("b", "x", "d")), 20L))
  observed <- matrix(c(20, 0, 0, 20, 0, 0, 0, 20, 20), 3L, byrow = TRUE,
                     dimnames = list(c("a", "b", "x"), c("x", "c", "d")))
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  seen <- observed > 0
  by_hand <- 2 * sum(observed[seen] * log(observed[seen] / expected[seen]))
  tested <- hg_get(markov_order(memory_seqs, max_order = 1L, n_perm = 5L, seed = 1L))
  expect_equal(tested$g2[tested$order == 1L], by_hand, tolerance = 1e-8)
  # four cells at 20 against 10 expected: G2 = 160 log 2 on (3 - 1)(3 - 1) = 4 df
  expect_equal(by_hand, 160 * log(2))
  expect_identical(tested$df[tested$order == 1L], 4L)
  expect_equal(tested$p_asymptotic[tested$order == 1L],
               stats::pchisq(160 * log(2), df = 4, lower.tail = FALSE))
  # no shuffle of the next states reaches the perfect dependence: (0 + 1) / (5 + 1)
  expect_equal(tested$p_value[tested$order == 1L], 1 / 6)
})

test_that("markov_order() decodes a tna model to its sequences", {
  fx <- .mw_fixture()
  expect_identical(markov_order(.mw_tna(fx$wide), max_order = 2L,
                                n_perm = 20L, seed = 1L),
                   markov_order(fx$seqs, max_order = 2L, n_perm = 20L,
                                seed = 1L))
})

test_that("memory() gives one estimate on every input form", {
  fx <- .mw_fixture()
  ref <- memory(fx$wide, order = 2L, min_count = 1L)
  expect_s3_class(ref, "net_path_dependence")
  # the estimator reads only frames: a list is padded into the same frame
  expect_identical(memory(fx$seqs, order = 2L, min_count = 1L), ref)
  from_long <- memory(fx$long, action = "what", actor = "who", time = "when",
                      min_count = 1L)
  expect_identical(from_long, memory(fx$net, min_count = 1L))
  expect_s3_class(memory(fx$net, min_count = 1L, base = exp(1)), "net_path_dependence")
  # x goes to c or d evenly, one bit, and knowing what came before removes it all
  memory_seqs <- c(rep(list(c("a", "x", "c")), 20L), rep(list(c("b", "x", "d")), 20L))
  in_bits <- hg_get(memory(memory_seqs, order = 2L, min_count = 1L))
  expect_equal(in_bits$entropy_first_order, c(1, 1))
  expect_equal(in_bits$entropy_order_k, c(0, 0))
  in_nats <- hg_get(memory(memory_seqs, order = 2L, min_count = 1L, base = exp(1)))
  expect_equal(in_nats$entropy_first_order, c(log(2), log(2)), tolerance = 1e-6)
  # KL of the context's next state (c or d for certain) from x's (1/2, 1/2):
  # 1 * log2(1 / 0.5) = 1 bit, log 2 nats; 20 observations each
  expect_equal(in_bits$kl, c(1, 1))
  expect_equal(in_nats$kl, c(log(2), log(2)), tolerance = 1e-6)
  expect_equal(in_bits$entropy_drop, c(1, 1))
  expect_equal(in_bits$count, c(20, 20))
  expect_identical(in_bits$context, c("a -> x", "b -> x"))
})

test_that("hg_markov_stability() gives one estimate on every input form", {
  fx <- .mw_fixture()
  P <- matrix(c(0.5, 0.3, 0.2, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), 3, 3,
              byrow = TRUE, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_s3_class(hg_markov_stability(P), "net_markov_stability")
  # the balance equations of P give 13/55, 23/55, 19/55; staying has probability 0.5
  by_state <- hg_get(hg_markov_stability(P))
  expect_equal(by_state$stationary, c(13, 23, 19) / 55, tolerance = 1e-3)
  expect_equal(by_state$persistence, rep(0.5, 3L))
  expect_equal(by_state$sojourn_time, rep(2, 3L))
  # mean passage times solve m_ij = 1 + sum_{k != j} P_ik m_kj; into A:
  # 0.5 m_B - 0.3 m_C = 1 and -0.4 m_B + 0.5 m_C = 1 give 80/13 and 90/13
  # (into B: 70/23, 60/23; into C: 80/19, 70/19); the diagonal is 1 / pi
  steps <- hg_get(hg_markov_stability(P), what = "passage_time")
  expect_identical(steps$from, rep(c("A", "B", "C"), 3L))
  expect_identical(steps$to, rep(c("A", "B", "C"), each = 3L))
  expect_equal(steps$steps, c(55 / 13, 80 / 13, 90 / 13, 70 / 23, 55 / 23,
                              60 / 23, 80 / 19, 70 / 19, 55 / 19))
  ref <- hg_markov_stability(fx$wide)
  expect_identical(hg_markov_stability(fx$seqs), ref)
  from_long <- hg_markov_stability(fx$long, action = "what", actor = "who",
                                   time = "when")
  expect_identical(from_long, hg_markov_stability(fx$net))
})

test_that("session splits an actor's events into separate sequences", {
  fx <- .mw_fixture()
  long <- fx$long
  long$sess <- ifelse(long$when > 50, "s2", "s1")
  seqs <- .ho_sequence_input(long, action = "what", actor = "who",
                             session = "sess", time = "when")
  # every actor with events on both sides of the split gives two sequences
  expect_identical(nrow(seqs), 24L)
  expect_identical(sum(!is.na(as.matrix(seqs))), nrow(long))
  # the model is the model of those sequences
  expect_identical(hon(long, action = "what", actor = "who", session = "sess",
                       time = "when", max_order = 2L), hon(seqs, max_order = 2L))
  # worked by hand: 20 people each do a x c in one session and b x d in the
  # next. Split, nothing follows c, and x goes to c or d half the time; joined,
  # c would lead into b
  split_long <- data.frame(who = rep(sprintf("p%02d", 1:20), each = 6L),
                           sess = rep(rep(c("s1", "s2"), each = 3L), 20L),
                           when = rep(1:6, 20L),
                           what = rep(c("a", "x", "c", "b", "x", "d"), 20L))
  by_session <- hon(split_long, action = "what", actor = "who", session = "sess",
                    time = "when", max_order = 2L)
  expect_equal(unname(rowSums(by_session$weights[c("c", "d"), , drop = FALSE])), c(0, 0))
  expect_equal(by_session$weights["x", c("c", "d")], c(c = 0.5, d = 0.5))
  expect_equal(by_session$weights["a -> x", "c"], 1)
  expect_equal(sort(by_session$edges$weight), c(0.5, 0.5, 1, 1, 1, 1))
})

test_that("a long table passed without action is refused, by class", {
  fx <- .mw_fixture()
  canonical <- data.frame(who = c("p", "p", "q"), state = c("a", "b", "a"),
                          timestamp = 1:3)
  lapply(list(
    function(d) hon(d),
    function(d) mogen(d),
    function(d) hypa(d),
    function(d) markov_order(d),
    function(d) memory(d),
    function(d) hg_markov_stability(d),
    function(d) hg_bootstrap(d, n_boot = 2L)
  ), function(f) {
    expect_error(f(canonical), class = "hypergraphs_long_format")
    expect_error(f(human_long), class = "hypergraphs_bad_input")
  })
  # a wide frame of states passes the guard
  expect_s3_class(hon(fx$wide, max_order = 2L), "net_hon")
})

test_that("long-format column arguments are validated", {
  fx <- .mw_fixture()
  expect_error(hon(fx$seqs, actor = "who"), class = "hypergraphs_bad_input")
  expect_error(hon(fx$wide, time = "when"), class = "hypergraphs_bad_input")
  expect_error(hon(fx$long, action = "nope"), class = "hypergraphs_bad_input")
  expect_error(hon(fx$long, action = "what", actor = "nope"),
               class = "hypergraphs_bad_input")
  expect_error(hon(fx$seqs, action = "what"), class = "hypergraphs_bad_input")
})

test_that("the wrappers are invariant to the row order of a long table", {
  fx <- .mw_fixture()
  shuffled <- fx$long[rev(seq_len(nrow(fx$long))), ]
  expect_identical(
    hon(shuffled, action = "what", actor = "who", time = "when",
        max_order = 2L),
    hon(fx$long, action = "what", actor = "who", time = "when",
        max_order = 2L))
})

test_that("hg_markov_stability() refuses a chain that is not irreducible", {
  absorbing <- matrix(c(0.5, 0.4, 0.1,
                        0.3, 0.6, 0.1,
                        0.0, 0.0, 1.0), 3, byrow = TRUE,
                      dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_error(hg_markov_stability(absorbing),
               class = "hypergraphs_not_ergodic")
  err <- tryCatch(hg_markov_stability(absorbing), error = identity)
  expect_match(conditionMessage(err), "transient: A, B")
  expect_match(conditionMessage(err), "absorbing: C")
  # an irreducible chain is estimated
  regular <- matrix(c(0.5, 0.3, 0.2, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), 3,
                    byrow = TRUE,
                    dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  estimated <- hg_markov_stability(regular)
  expect_s3_class(estimated, "net_markov_stability")
  # staying put has probability P_ii = 0.5 in every state, so 2 steps per visit
  by_state <- hg_get(estimated)
  expect_equal(by_state$persistence, c(0.5, 0.5, 0.5))
  expect_equal(by_state$sojourn_time, c(2, 2, 2))
})

test_that("a state with no outgoing transition is refused by class", {
  dead_end <- matrix(c(0.5, 0.5, 0,
                       0.3, 0.3, 0.4,
                       0, 0, 0), 3, byrow = TRUE,
                     dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  err <- tryCatch(hg_markov_stability(dead_end), error = identity)
  expect_s3_class(err, "hypergraphs_not_ergodic")
  expect_match(conditionMessage(err), "no outgoing transition from: C")
})

# ---- Missing states are gaps (M01) ----------------------------------------
# Each estimator sees a sequence with a gap as its runs of observed states.

test_that("mogen(), hypa() and markov_order() read a gap as a break", {
  gap <- list(c(NA, "a", "b", NA, "c", "d", "a", "b"),
              c("b", "c", NA, NA, "d", "a", "b", NA),
              c("a", "b", "c", "d", "a"))
  runs <- list(c("a", "b"), c("c", "d", "a", "b"), c("b", "c"),
               c("d", "a", "b"), c("a", "b", "c", "d", "a"))
  expect_identical(mogen(gap, max_order = 2L), mogen(runs, max_order = 2L))
  expect_identical(hg_get(hypa(gap, order = 2L, min_count = 1L),
                          what = "scores"),
                   hg_get(hypa(runs, order = 2L, min_count = 1L),
                          what = "scores"))
  expect_identical(markov_order(gap, max_order = 2L, n_perm = 10L, seed = 1L),
                   markov_order(runs, max_order = 2L, n_perm = 10L,
                                seed = 1L))
  scores <- hg_get(hypa(gap, order = 2L, min_count = 1L), what = "scores")
  expect_false(any(grepl("NA", scores$path, fixed = TRUE)))
})

test_that("memory() already ends a context at a gap and counts none across", {
  gap <- list(c("a", "b", NA, "c", "a", "b", "c"), c("a", "b", "c", NA))
  runs <- list(c("a", "b"), c("c", "a", "b", "c"), c("a", "b", "c"))
  expect_identical(hg_get(memory(gap, order = 2L, min_count = 1L)),
                   hg_get(memory(runs, order = 2L, min_count = 1L)))
})

# ---- Labels that collide with the memory-node notation (M07) ---------------

test_that("path-naming estimators refuse state labels holding the arrow", {
  seqs <- rep(list(c("a -> b", "c", "a -> b", "c")), 2L)
  expect_error(mogen(seqs, max_order = 1L), class = "hypergraphs_bad_input")
  expect_error(hypa(seqs, order = 1L, min_count = 1L),
               class = "hypergraphs_bad_input")
  expect_error(memory(seqs, order = 2L, min_count = 1L),
               class = "hypergraphs_bad_input")
  expect_error(hg_bootstrap(seqs, n_boot = 3L, max_order = 1L),
               class = "hypergraphs_bad_input")
  # the order test names no path, so such labels are harmless there
  expect_s3_class(markov_order(seqs, max_order = 1L, n_perm = 5L, seed = 1L),
                  "net_markov_order")
})

# ---- memory(): order is the number of conditioning states (M16) ------------

test_that("memory(order = k) conditions on k states", {
  seqs <- list(c("a", "b", "c", "a", "b", "c", "a", "b", "c"))
  states_in <- \(context) lengths(strsplit(context, " -> ", fixed = TRUE))
  two <- hg_get(memory(seqs, order = 2L, min_count = 1L))
  three <- hg_get(memory(seqs, order = 3L, min_count = 1L))
  expect_true(all(states_in(two$context) == 2L))
  expect_true(all(states_in(three$context) == 3L))
})
