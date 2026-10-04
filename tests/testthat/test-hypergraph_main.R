# hypergraph() hands every input to its constructor unchanged.

test_that("a data frame of groups goes to group_hypergraph()", {
  meetings <- data.frame(person = c("a", "b", "c", "a", "b", "d"),
                         meeting = c("m1", "m1", "m1", "m2", "m2", "m2"))
  expect_identical(hypergraph(meetings, node = "person", hyperedge = "meeting"),
                   group_hypergraph(meetings, node = "person",
                                    hyperedge = "meeting"))
  expect_identical(hypergraph(tutoring_events, node = "event", hyperedge = "step",
                              group = "outcome", top = 4),
                   group_hypergraph(tutoring_events, node = "event",
                                    hyperedge = "step", group = "outcome", top = 4))
})

test_that("sequences go to window_hypergraph()", {
  sessions <- list(c("a", "b", "c", "a"), c("b", "c", "d"))
  expect_identical(hypergraph(sessions, window = 2),
                   window_hypergraph(sessions, window = 2))
  long <- data.frame(user = c(1, 1, 1, 2, 2), act = c("a", "b", "c", "b", "c"),
                     at = c(1, 2, 3, 1, 2))
  expect_identical(hypergraph(long, action = "act", actor = "user",
                              time = "at", window = 2),
                   window_hypergraph(long, action = "act", actor = "user",
                                     time = "at", window = 2))
})

test_that("data with a clock go to temporal_hypergraph()", {
  contacts <- data.frame(actor = c("a", "b", "a", "c"),
                         group = c("g1", "g1", "g2", "g2"),
                         time = c(1, 1, 2, 2))
  expect_identical(hypergraph(contacts, node = "actor", hyperedge = "group",
                              time = "time"),
                   temporal_hypergraph(contacts, node = "actor",
                                       hyperedge = "group", time = "time"))
})

test_that("a network goes to network_hypergraph()", {
  weights <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3,
                    dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
  expect_identical(hypergraph(weights), network_hypergraph(weights))
  sparse <- Matrix::Matrix(weights, sparse = TRUE)
  expect_identical(hypergraph(sparse), network_hypergraph(weights))
  expect_identical(hypergraph(weights, include_pairwise = FALSE),
                   network_hypergraph(weights, include_pairwise = FALSE))
})

test_that("a topic model goes to group_hypergraph()", {
  hg <- text_hypergraph(c(a = "soup salt onion soup broth",
                          b = "salt soup broth onion",
                          c = "stars sky moon night",
                          d = "sky stars night moon moon",
                          e = "soup stars salt sky night broth"))
  fit <- hg_topics(hg, k = 2, nstart = 2)
  expect_identical(hypergraph(fit, threshold = 0.2),
                   group_hypergraph(fit, threshold = 0.2))
})

test_that("an unreadable input is refused by class", {
  expect_error(hypergraph(42), class = "hypergraphs_bad_input")
  expect_error(hypergraph("a"), class = "hypergraphs_bad_input")
})

test_that("event data without a window give one hyperedge per session", {
  visits <- data.frame(user = c("u1", "u1", "u1", "u1", "u2", "u2"),
                       visit = c(1, 1, 2, 2, 1, 1),
                       page = c("home", "cart", "home", "help", "home", "cart"))
  # sessions within their actor: u1's visit 1 and u2's visit 1 differ
  nested <- hypergraph(visits, action = "page", actor = "user",
                       session = "visit")
  expect_identical(nested$n_hyperedges, 3L)
  keyed <- transform(visits, visit = paste(user, visit, sep = "."))
  expect_identical(nested, group_hypergraph(keyed, node = "page",
                                            hyperedge = "visit"))
  # a session alone, and an actor alone
  expect_identical(hypergraph(visits, action = "page", session = "visit")$n_hyperedges,
                   2L)
  per_user <- hypergraph(visits, action = "page", actor = "user")
  expect_setequal(hg_get(per_user)$members, c("cart, help, home", "cart, home"))
  # counted sets with a comparison group
  visits$device <- c("phone", "phone", "laptop", "laptop", "phone", "phone")
  counted <- hypergraph(visits, action = "page", actor = "user",
                        session = "visit", group = "device", top = 2)
  expect_setequal(unique(hg_get(counted, what = "sets")$group),
                  c("laptop", "phone"))
  # a clock without a window, or no unit at all, is refused
  expect_error(hypergraph(visits, action = "page", session = "visit",
                          time = "visit"), class = "hypergraphs_bad_input")
  expect_error(hypergraph(visits, action = "page"),
               class = "hypergraphs_bad_input")
  expect_error(hypergraph(visits, action = "nope", session = "visit"),
               class = "hypergraphs_bad_input")
})

test_that("event data with a window go to window_hypergraph()", {
  visits <- data.frame(user = c("u1", "u1", "u1", "u2", "u2", "u2"),
                       at = c(1, 2, 3, 1, 2, 3),
                       page = c("home", "cart", "pay", "home", "help", "cart"))
  expect_identical(hypergraph(visits, action = "page", actor = "user",
                              time = "at", window = 2),
                   window_hypergraph(visits, action = "page", actor = "user",
                                     time = "at", window = 2))
})
