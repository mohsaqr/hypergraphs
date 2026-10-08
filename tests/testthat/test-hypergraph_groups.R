# ---- group_hypergraph() tests --------------------------------------------

# Helpers ------------------------------------------------------------------

.bg_sample_data <- function() {
  data.frame(
    member = c("Alice", "Bob", "Carol", "Alice", "Bob",
               "Dave", "Carol", "Dave", "Eve"),
    session = c("S1", "S1", "S1", "S2", "S2",
                "S3", "S3", "S3", "S3"),
    stringsAsFactors = FALSE
  )
}

# Structure ----------------------------------------------------------------

test_that("returns a net_hg with required fields", {
  hg <- group_hypergraph(.bg_sample_data(), node = "member", hyperedge = "session")
  expect_s3_class(hg, "net_hg")
  expect_named(hg, c("hyperedges", "incidence", "nodes", "n_nodes",
                     "n_hyperedges", "size_distribution", "params"))
  expect_equal(hg$n_nodes, 5L)
  expect_equal(hg$n_hyperedges, 3L)
  expect_setequal(hg$nodes, c("Alice", "Bob", "Carol", "Dave", "Eve"))
})

# Hyperedge content --------------------------------------------------------

test_that("each group becomes a hyperedge spanning its members", {
  hg <- group_hypergraph(.bg_sample_data(), node = "member", hyperedge = "session")
  # Map back: hyperedges are integer indices into hg$nodes
  members_by_session <- lapply(hg$hyperedges, function(idx) sort(hg$nodes[idx]))
  names(members_by_session) <- colnames(hg$incidence)
  expect_equal(members_by_session[["S1"]], c("Alice", "Bob", "Carol"))
  expect_equal(members_by_session[["S2"]], c("Alice", "Bob"))
  expect_equal(members_by_session[["S3"]], c("Carol", "Dave", "Eve"))
})

# Incidence matrix ---------------------------------------------------------

test_that("incidence is binary by default with correct dimensions and sums", {
  hg <- group_hypergraph(.bg_sample_data(), node = "member", hyperedge = "session")
  expect_equal(dim(hg$incidence), c(5L, 3L))
  expect_true(all(hg$incidence %in% c(0L, 1L)))
  # column sums = group sizes
  expect_equal(unname(colSums(hg$incidence)), c(3L, 2L, 3L))
  # row sums = member hyperdegree (num groups they appeared in)
  expect_equal(hg$incidence["Alice", , drop = TRUE], c(S1 = 1L, S2 = 1L, S3 = 0L))
  expect_equal(rowSums(hg$incidence)[["Carol"]], 2L)
})

# Size distribution --------------------------------------------------------

test_that("size_distribution counts hyperedges by size", {
  hg <- group_hypergraph(.bg_sample_data(), node = "member", hyperedge = "session")
  expect_equal(hg$size_distribution[["size_2"]], 1L)
  expect_equal(hg$size_distribution[["size_3"]], 2L)
})

# Weighted incidence -------------------------------------------------------

test_that("weight column produces weighted incidence (sum per cell)", {
  d <- data.frame(
    member = c("A", "B", "A", "A", "B"),
    grp    = c("g1", "g1", "g1", "g2", "g2"),
    n      = c(2, 5, 3, 1, 4),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, node = "member", hyperedge = "grp", weight = "n")
  # A in g1: 2 + 3 = 5; B in g1: 5; A in g2: 1; B in g2: 4
  expect_equal(hg$incidence["A", "g1"], 5)
  expect_equal(hg$incidence["B", "g1"], 5)
  expect_equal(hg$incidence["A", "g2"], 1)
  expect_equal(hg$incidence["B", "g2"], 4)
  # hyperedges still based on membership (any nonzero), not weight
  expect_setequal(hg$nodes[hg$hyperedges[[1]]], c("A", "B"))
})

# NA handling --------------------------------------------------------------

test_that("rows with NA in member or group are dropped silently", {
  d <- data.frame(
    member = c("A", "B", NA,  "C", "D"),
    grp    = c("g1", "g1", "g1", NA,  "g2"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, node = "member", hyperedge = "grp")
  expect_equal(hg$n_nodes, 3L)         # A, B, D (C dropped because group NA)
  expect_setequal(hg$nodes, c("A", "B", "D"))
  expect_equal(hg$n_hyperedges, 2L)    # g1, g2
  # n_observations records post-drop count
  expect_equal(hg$params$n_observations, 3L)
})

test_that("all-NA data raises an error", {
  d <- data.frame(member = NA, grp = NA, stringsAsFactors = FALSE)
  expect_error(group_hypergraph(d, "member", "grp"),
               "No complete observations")
})

# Repeated rows are deduplicated in binary mode ----------------------------

test_that("duplicate rows do not duplicate membership in binary mode", {
  d <- data.frame(
    member = c("A", "A", "A", "B"),
    grp    = c("g1", "g1", "g1", "g1"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "grp")
  expect_equal(hg$incidence["A", "g1"], 1L)
  expect_equal(hg$incidence["B", "g1"], 1L)
  expect_equal(hg$n_hyperedges, 1L)
  expect_equal(length(hg$hyperedges[[1]]), 2L)
})

# Single-member groups -----------------------------------------------------

test_that("singleton groups become size-1 hyperedges", {
  d <- data.frame(
    member = c("A", "B", "C"),
    grp    = c("g1", "g2", "g3"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "grp")
  expect_equal(hg$n_hyperedges, 3L)
  expect_true(all(vapply(hg$hyperedges, length, integer(1)) == 1L))
})

test_that("an explicit node universe preserves isolates", {
  d <- data.frame(member = c("A", "B"), grp = c("g1", "g1"))
  hg <- group_hypergraph(d, "member", "grp", nodes = c("A", "B", "C"))
  expect_equal(hg$nodes, c("A", "B", "C"))
  expect_equal(unname(rowSums(hg$incidence)), c(1, 1, 0))
  expect_equal(hg$n_nodes, 3L)
  expect_error(group_hypergraph(d, "member", "grp", nodes = c("A", "C")),
               "Every observed member")
  expect_error(group_hypergraph(d, "member", "grp", nodes = c("A", "B", "B")),
               "unique")
})

test_that("sparse group incidence is identical to dense incidence", {
  d <- rbind(.bg_sample_data(), .bg_sample_data()[1, , drop = FALSE])
  dense <- group_hypergraph(d, "member", "session")
  sparse <- group_hypergraph(d, "member", "session", sparse = TRUE)
  expect_s4_class(sparse$incidence, "sparseMatrix")
  expect_equal(as.matrix(sparse$incidence), dense$incidence)
  expect_equal(sparse$hyperedges, dense$hyperedges)

  weighted <- transform(d, n = seq_len(nrow(d)))
  dense_w <- group_hypergraph(weighted, "member", "session", "n")
  sparse_w <- group_hypergraph(weighted, "member", "session", "n",
                               sparse = TRUE)
  expect_equal(as.matrix(sparse_w$incidence), dense_w$incidence)
})

# Validation ---------------------------------------------------------------

test_that("missing member column raises error", {
  d <- data.frame(member = c("A"), grp = c("g1"), stringsAsFactors = FALSE)
  expect_error(group_hypergraph(d, node = "playor", hyperedge = "grp"))
})

test_that("missing group column raises error", {
  d <- data.frame(member = c("A"), grp = c("g1"), stringsAsFactors = FALSE)
  expect_error(group_hypergraph(d, node = "member", hyperedge = "guppy"))
})

test_that("non-data.frame input rejected", {
  expect_error(group_hypergraph(list(member = "A", grp = "g"), "member", "grp"))
})

# print and summary work --------------------------------------------------

test_that("print and summary work via shared net_hg methods", {
  hg <- group_hypergraph(.bg_sample_data(), node = "member", hyperedge = "session")
  expect_invisible(print(hg))
  # summary now returns a tidy node-degree data.frame (visible)
  s <- summary(hg)
  expect_s3_class(s, "data.frame")
  expect_setequal(names(s), c("node", "degree"))
})

# Integration with bundled dataset ----------------------------------------

test_that("works on bundled human_long dataset (long-format event data)", {
  data("human_long", package = "hypergraphs")
  hg <- group_hypergraph(human_long, node = "code", hyperedge = "session_id")
  expect_s3_class(hg, "net_hg")
  expect_gt(hg$n_nodes, 0L)
  expect_gt(hg$n_hyperedges, 0L)
  # Each session is a hyperedge; members = codes appearing in that session
  expect_equal(ncol(hg$incidence), hg$n_hyperedges)
  expect_equal(nrow(hg$incidence), hg$n_nodes)
})

# Hyperedge attributes -----------------------------------------------------

test_that("columns constant within a group become hyperedge attributes", {
  dat <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "d", "e"),
    event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3"),
    kind = c("x", "x", "x", "x", "x", "x", "y", "y"),
    seat = c("p", "q", "r", "p", "q", "r", "p", "q"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(dat, node = "member", hyperedge = "event")
  expect_identical(names(hg$edge_data), c("edge", "kind"))
  expect_identical(hg$edge_data$edge, c("e1", "e2", "e3"))
  expect_identical(hg$edge_data$kind, c("x", "x", "y"))
  kept <- hg_subset(hg, where = c(kind = "y"))
  expect_identical(kept$nodes, c("d", "e"))
  # an NA inside a group does not make the attribute vary
  dat$kind[2L] <- NA
  with_na <- group_hypergraph(dat, node = "member", hyperedge = "event")
  expect_identical(with_na$edge_data$kind, c("x", "x", "y"))
  # sparse and dense agree
  sparse <- group_hypergraph(dat, node = "member", hyperedge = "event",
                             sparse = TRUE)
  expect_identical(sparse$edge_data, with_na$edge_data)
})

test_that("INVARIANT: a plain member/group table keeps the original layout", {
  hg <- group_hypergraph(.bg_sample_data(), node = "member", hyperedge = "session")
  expect_false("edge_data" %in% names(hg))
  weighted <- data.frame(member = c("a", "b"), session = c("s", "s"),
                         w = c(1, 2))
  hg_w <- group_hypergraph(weighted, node = "member", hyperedge = "session",
                           weight = "w")
  expect_false("edge_data" %in% names(hg_w))
})

test_that("an edge list carries its attributes onto the size-2 hyperedges", {
  contacts <- data.frame(from = c("a", "b"), to = c("b", "c"),
                         kind = c("call", "mail"), stringsAsFactors = FALSE)
  hg <- group_hypergraph(contacts, from = "from", to = "to")
  expect_identical(hg$edge_data,
                   data.frame(edge = c("e1", "e2"), kind = c("call", "mail"),
                              stringsAsFactors = FALSE))
})

test_that("an unaddressable dense incidence is refused, not attempted", {
  skip_on_cran()
  # 50000 members x 50000 groups = 2.5e9 cells: past the integer cell index,
  # and ~20 Gb if it were allocated.
  n <- 50000L
  d <- data.frame(member = sprintf("m%d", seq_len(n)),
                  group = sprintf("g%d", seq_len(n)))
  expect_error(group_hypergraph(d, node = "member", hyperedge = "group"),
               class = "hypergraphs_dense_too_large")
  expect_s3_class(group_hypergraph(d, node = "member", hyperedge = "group",
                                   sparse = TRUE),
                  "net_hg")
})

test_that("group_hypergraph(group =) counts each group's set within by", {
  events <- data.frame(
    item = c("a", "b", "a", "b", "a", "c", "a", "b", "b", "c"),
    basket = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    shop = c("x", "x", "x", "x", "x", "x", "y", "y", "y", "y"))
  sets <- group_hypergraph(events, node = "item", hyperedge = "basket",
                           group = "shop", top = Inf)
  table_sets <- hg_get(sets, what = "sets")
  x_sets <- subset(table_sets, group == "x")
  expect_identical(x_sets$set, c("a + b", "a + c"))
  expect_identical(x_sets$count, c(2L, 1L))
  expect_equal(x_sets$share, c(2, 1) / 3)
  y_sets <- subset(table_sets, group == "y")
  expect_setequal(y_sets$set, c("a + b", "b + c"))
  expect_output(print(sets), "sets of item per basket, counted within each shop")
  # top keeps the most frequent; min_size drops the small sets
  top_one <- hg_get(group_hypergraph(events, node = "item", hyperedge = "basket",
                                     group = "shop", top = 1), what = "sets")
  expect_identical(subset(top_one, group == "x")$set, "a + b")
  expect_s3_class(plot(sets, group = "x"), "ggplot")
})

test_that("group_hypergraph(group =) refuses a by that varies within a group", {
  events <- data.frame(item = c("a", "b"), basket = c(1, 1),
                       shop = c("x", "y"))
  expect_error(group_hypergraph(events, node = "item", hyperedge = "basket",
                                group = "shop"),
               class = "hypergraphs_bad_input")
  expect_error(group_hypergraph(events, node = "item", hyperedge = "basket",
                                group = "nope"),
               class = "hypergraphs_bad_input")
  expect_error(group_hypergraph(events, node = "item", hyperedge = "basket",
                                min_size = 2),
               class = "hypergraphs_bad_input")
})

test_that("tutoring_events trials follow the after-Incorrect rule", {
  second <- subset(tutoring_events, step == 2L)
  after_incorrect <- c(FALSE, head(second$event, -1L) == "Incorrect")
  trial_number <- as.integer(sub("^[0-9]+[.]", "", second$trial))
  expect_identical(trial_number, cumsum(after_incorrect) + 1L)
  groups <- tapply(tutoring_events$outcome, tutoring_events$step,
                   \(v) length(unique(v)))
  expect_true(all(groups == 1L))
  expect_setequal(unique(tutoring_events$outcome), c("completed", "stopped"))
  expect_setequal(unique(tutoring_events$event),
                  c("Begin", "Task", "Attempt", "Reattempt", "Correct",
                    "Incorrect", "Completed", "Stopped", "Guidance",
                    "Question", "Order", "Refute", "Reflect", "Comfort",
                    "GiveUp"))
  expect_identical(anyNA(tutoring_events), FALSE)
  expect_identical(nrow(tutoring_events), 84356L)
  expect_length(unique(tutoring_events$step), 13309L)
})

test_that("group_hypergraph(top =) without by counts over all groups", {
  events <- data.frame(
    item = c("a", "b", "a", "b", "a", "c", "a", "b"),
    basket = c(1, 1, 2, 2, 3, 3, 4, 4))
  sets <- hg_get(group_hypergraph(events, node = "item", hyperedge = "basket",
                                  top = Inf), what = "sets")
  expect_identical(unique(sets$group), "All baskets")
  expect_identical(sets$set, c("a + b", "a + c"))
  expect_identical(sets$count, c(3L, 1L))
  # neither by nor top: one hyperedge per basket, as before
  plain <- group_hypergraph(events, node = "item", hyperedge = "basket")
  expect_identical(plain$n_hyperedges, 4L)
})

test_that("min_share keeps the sets that reach a minimum support", {
  events <- data.frame(
    item = c("a", "b", "a", "b", "a", "b", "a", "c", "a", "b", "d"),
    basket = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6)
  )
  all_sets <- hg_get(group_hypergraph(events, node = "item", hyperedge = "basket",
                                      top = Inf), what = "sets")
  kept <- group_hypergraph(events, node = "item", hyperedge = "basket",
                           min_share = 0.2)
  kept_sets <- hg_get(kept, what = "sets")
  # the minimum support is the share of the baskets that hold the set
  expect_identical(kept_sets$set, subset(all_sets, share >= 0.2)$set)
  expect_true(all(kept_sets$share >= 0.2))
  expect_identical(kept$params$min_share, 0.2)
  expect_output(print(kept), "sets of at least 20% per group")
  # top still caps the number kept
  capped <- hg_get(group_hypergraph(events, node = "item", hyperedge = "basket",
                                    min_share = 0.1, top = 1), what = "sets")
  expect_identical(nrow(capped), 1L)
  expect_error(group_hypergraph(events, node = "item", hyperedge = "basket",
                                min_share = 0), class = "hypergraphs_bad_input")
  expect_error(group_hypergraph(events, node = "item", hyperedge = "basket",
                                min_share = 2), class = "hypergraphs_bad_input")
})

test_that("a set's members, not its display label, identify it (R02)", {
  d <- data.frame(node = c("a + b", "a", "b"), edge = c("x", "y", "y"))
  hg <- group_hypergraph(d, node = "node", hyperedge = "edge", top = Inf)
  sets <- hg_get(hg, what = "sets")
  expect_identical(nrow(sets), 2L)
  expect_identical(sort(sets$size), c(1L, 2L))
  expect_identical(sets$count, c(1L, 1L))
  expect_false(anyDuplicated(sets$hyperedge) > 0L)
  expect_setequal(hg_get(hg)$members, c("a + b", "a, b"))
  # grouped prefixes: a group name holding ": " cannot merge two hyperedges
  g <- data.frame(node = c("z", "y: z"), edge = c("t1", "t2"),
                  grp = c("x: y", "x"))
  grouped <- group_hypergraph(g, node = "node", hyperedge = "edge",
                              group = "grp", top = Inf)
  expect_identical(grouped$n_hyperedges, 2L)
  expect_identical(nrow(hg_get(grouped, what = "sets")), 2L)
  # the plot's distinct-set view keeps them apart too
  repeated <- data.frame(node = c("a + b", "a + b", "a", "b"),
                         edge = c("x1", "x2", "y", "y"))
  distinct <- .thg_distinct_sets(group_hypergraph(repeated, node = "node",
                                                  hyperedge = "edge"))
  expect_identical(distinct$n_hyperedges, 2L)
  expect_setequal(hg_get(distinct, what = "sets")$count, c(2L, 1L))
})

test_that("metadata columns never replace structural columns (R03)", {
  d <- data.frame(from = c("a", "b"), to = c("b", "c"), edge = c("same", "same"),
                  actor = c("p", "q"), weight = c(5, 6))
  hg <- group_hypergraph(d, from = "from", to = "to")
  expect_identical(hg$n_hyperedges, 2L)
  expect_setequal(hg_get(hg)$members, c("a, b", "b, c"))
  # the metadata stay as attributes, under a name that does not clash
  expect_identical(hg_get(hg, what = "edge_data")$edge, c("e1", "e2"))
  expect_true("edge_1" %in% names(hg_get(hg, what = "edge_data")))
  expect_identical(hg$params$node, "actor")
  # weights still come from the named weight column, not from `weight`
  w <- data.frame(from = c("a", "b"), to = c("b", "c"), w = c(2, 3),
                  weight = c(100, 100))
  weighted <- group_hypergraph(w, from = "from", to = "to", weight = "w")
  expect_equal(unname(colSums(weighted$incidence)), c(4, 6))
  # membership data with a metadata column called `edge`
  m <- data.frame(person = c("a", "b", "c"), meeting = c("m1", "m1", "m2"),
                  edge = c("k", "k", "l"))
  mh <- group_hypergraph(m, node = "person", hyperedge = "meeting")
  ed <- hg_get(mh, what = "edge_data")
  expect_identical(ed$edge, c("m1", "m2"))
  expect_identical(ed$edge_1, c("k", "l"))
})

test_that("membership weights are validated, zero membership is absent (R04)", {
  base <- data.frame(node = c("a", "b"), edge = "x")
  bad <- list(c(1, -1), c(1, Inf), factor(c("2", "3")), c("1", "2"))
  lapply(bad, function(w) {
    d <- base
    d$w <- w
    expect_error(group_hypergraph(d, node = "node", hyperedge = "edge",
                                  weight = "w"),
                 class = "hypergraphs_bad_input")
    expect_error(group_hypergraph(d, node = "node", hyperedge = "edge",
                                  weight = "w", sparse = TRUE),
                 class = "hypergraphs_bad_input")
  })
  zero <- data.frame(node = c("a", "b", "c"), edge = c("x", "x", "y"),
                     w = c(1, 0, 2))
  dense <- group_hypergraph(zero, node = "node", hyperedge = "edge", weight = "w")
  sparse <- group_hypergraph(zero, node = "node", hyperedge = "edge", weight = "w",
                             sparse = TRUE)
  expect_identical(hg_get(dense), hg_get(sparse))
  expect_identical(hg_get(dense)$members, c("a", "c"))
  expect_identical(hg_get(dense, what = "nodes")$degree, c(1L, 0L, 1L))
})

test_that("edge attributes do not depend on row order (R05)", {
  d <- data.frame(node = c("a", "b", "c", "d"), edge = c("x", "x", "y", "y"),
                  color = c(NA, "red", "blue", NA),
                  when = as.Date(c(NA, "2020-01-02", NA, NA)))
  orders <- list(1:4, 4:1, c(2, 1, 4, 3))
  tables <- lapply(orders, function(o) {
    hg_get(group_hypergraph(d[o, ], node = "node", hyperedge = "edge"),
           what = "edge_data")
  })
  lapply(tables, function(t) expect_identical(t, tables[[1L]]))
  expect_identical(tables[[1L]]$color, c("red", "blue"))
  expect_s3_class(tables[[1L]]$when, "Date")
  expect_true(is.na(tables[[1L]]$when[2L]))
  sparse <- hg_get(group_hypergraph(d[4:1, ], node = "node", hyperedge = "edge",
                                    sparse = TRUE), what = "edge_data")
  expect_identical(sparse, tables[[1L]])
  picked <- hg_subset(group_hypergraph(d[4:1, ], node = "node", hyperedge = "edge"),
                      where = list(color = "red"))
  expect_identical(picked$n_hyperedges, 1L)
})

test_that("the canonical node and hyperedge columns are detected (R20)", {
  d <- data.frame(node = c("a", "b", "c"), hyperedge = c("e", "e", "f"))
  detected <- group_hypergraph(d)
  expect_identical(detected, group_hypergraph(d, node = "node",
                                              hyperedge = "hyperedge"))
  mixed <- data.frame(Node = c("a", "b"), HyperEdge = c("e", "e"))
  expect_identical(group_hypergraph(mixed)$n_hyperedges, 1L)
  # an explicit name still wins over detection
  d$team <- c("t1", "t1", "t1")
  expect_identical(group_hypergraph(d, node = "node", hyperedge = "team")$n_hyperedges,
                   1L)
  td <- data.frame(node = c("a", "b"), hyperedge = c("e", "e"), time = c(1, 1))
  expect_identical(temporal_hypergraph(td)$edges, "e")
})
