# hypergraph() hands every input to its constructor unchanged.

test_that("a data frame of groups goes to group_hypergraph()", {
  meetings <- data.frame(person = c("a", "b", "c", "a", "b", "d"),
                         meeting = c("m1", "m1", "m1", "m2", "m2", "m2"))
  expect_identical(hypergraph(meetings, actor = "person", group = "meeting"),
                   group_hypergraph(meetings, actor = "person",
                                    group = "meeting"))
  expect_identical(hypergraph(debug_events, actor = "event", group = "session",
                              by = "group", top = 4),
                   group_hypergraph(debug_events, actor = "event",
                                    group = "session", by = "group", top = 4))
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
  expect_identical(hypergraph(contacts, actor = "actor", group = "group",
                              time = "time"),
                   temporal_hypergraph(contacts, actor = "actor",
                                       group = "group", time = "time"))
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
