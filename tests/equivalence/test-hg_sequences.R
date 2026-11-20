# hg_sequences(): a clustered corpus as per-actor state sequences, one row per document with its actor, its time
# and its cluster as the state. Reference: the documents table joined to the clusters in base R, ordered by actor
# and time.

test_that("the sequence of each student equals the documents ordered by turn with their cluster", {
  hg <- text_hypergraph(forum_posts, column = "text", id = "post", stop_words = stop_words_en(), min_count = 3L)
  topics <- hg_cluster(hg, k = 4L, seed = 1L)
  posts <- forum_posts[order(forum_posts$student, forum_posts$turn), ]
  expected <- data.frame(actor = posts$student, time = posts$turn, action = topics$cluster[match(posts$post, topics$node)])
  expected$position <- stats::ave(expected$time, expected$actor, FUN = seq_along)
  ours <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
  ours$position <- stats::ave(ours$time, ours$actor, FUN = seq_along)
  expect_reference(ours, expected, reference = "independent: the documents joined to their clusters, ordered by student and turn",
                   covers = "hg_sequences()", quantity = "forum posts, four clusters", tolerance = 0, key = c("actor", "position"))
})
