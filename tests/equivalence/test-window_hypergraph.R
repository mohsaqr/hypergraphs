# window_hypergraph(): hyperedges from windows slid over sequences. Each window of `window` consecutive events,
# started every `step` events, gives the set of distinct states in it; equal sets are one hyperedge whose weight
# is the number of windows that gave it, and a member's weight is its number of occurrences over those windows.
# Reference: the windows enumerated in base R. Nestimate's windowed co-occurrence (internal: hypergraphs imports
# Nestimate) for the pairs of two-event tumbling windows.

#' Windows of a list of sequences, as one row per (hyperedge set, member) with the window count and occurrences
by_definition <- function(sequences, window, step) {
  windows <- unlist(lapply(sequences, \(events) {
    starts <- seq(1L, length(events) - window + 1L, by = step)
    lapply(starts, \(start) events[start:(start + window - 1L)])
  }), recursive = FALSE)
  sets <- vapply(windows, \(events) paste(sort(unique(events)), collapse = ", "), character(1))
  edges <- data.frame(members = names(table(sets)), weight = as.vector(table(sets)), stringsAsFactors = FALSE)
  edges$size <- lengths(strsplit(edges$members, ", ", fixed = TRUE))
  occurrences <- do.call(rbind, Map(\(events, set) data.frame(members = set, node = events), windows, sets))
  members <- stats::aggregate(list(weight = rep(1, nrow(occurrences))), occurrences[c("members", "node")], sum)
  list(edges = edges, members = members)
}

sequences <- list(s1 = c("a", "b", "a", "c", "c", "b", "d", "a", "b", "b"),
                  s2 = c("b", "c", "d", "a", "a", "c", "b"),
                  s3 = c("d", "d", "c", "a", "b", "c", "a", "d", "b"))

test_that("hyperedges, their window counts and member occurrences equal the enumerated windows", {
  settings <- list(c(window = 3L, step = 1L), c(window = 2L, step = 2L), c(window = 4L, step = 2L))
  invisible(lapply(settings, \(setting) {
    expected <- by_definition(sequences, setting[["window"]], setting[["step"]])
    ours <- window_hypergraph(sequences, window = setting[["window"]], step = setting[["step"]])
    label <- sprintf("window %d, step %d", setting[["window"]], setting[["step"]])
    # Counts: integers on both sides; a hyperedge is keyed by its member set.
    expect_reference(hg_get(ours, what = "edges"), expected$edges,
                     reference = "independent: windows enumerated over the sequences, base R",
                     covers = c("window_hypergraph()", "hg_get.net_hg()"), quantity = paste(label, "hyperedge weight and size"),
                     tolerance = 0, key = "members")
    edges <- hg_get(ours, what = "edges")
    memberships <- merge(hg_get(ours, what = "memberships"), data.frame(hyperedge = edges$hyperedge, members = edges$members))
    expect_reference(memberships, expected$members,
                     reference = "independent: windows enumerated over the sequences, base R",
                     covers = c("window_hypergraph()", "hg_get.net_hg()"), quantity = paste(label, "member occurrences"),
                     tolerance = 0, key = c("members", "node"))
  }))
})

test_that("pairs of two-event tumbling windows equal Nestimate's windowed co-occurrence (internal)", {
  skip_without_reference("Nestimate")
  set.seed(515)
  events <- sample(letters[1:5], 40L, replace = TRUE)
  states <- sort(unique(events))
  onehot <- as.data.frame(outer(events, states, `==`) * 1L)
  names(onehot) <- states
  cooccurrence <- Nestimate::wtna(onehot, method = "cooccurrence", type = "frequency", window_size = 2L,
                                  mode = "non-overlapping")$weights
  pairs <- ref_pairs_upper(cooccurrence[states, states], name = "weight")
  edges <- hg_get(window_hypergraph(list(events), window = 2L, step = 2L), what = "edges")
  ours <- subset(edges, size == 2L)
  ours <- cbind(do.call(rbind, lapply(strsplit(ours$members, ", ", fixed = TRUE), \(pair) data.frame(from = pair[[1L]], to = pair[[2L]]))),
                weight = ours$weight)
  expect_reference(ours, pairs, reference = "Nestimate", covers = "window_hypergraph()",
                   quantity = "two-event tumbling windows, pair counts", tolerance = 0, key = c("from", "to"),
                   note = "Nestimate::wtna(method = 'cooccurrence', mode = 'non-overlapping'); off-diagonal cells only.")
})
