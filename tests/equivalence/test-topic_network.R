# topic_network(): the network of the topics of a partition of the documents, the strength of two topics being the
# sum over the words of the product of the numbers of documents in each topic that contain the word (full counting,
# the co-occurrence aggregation of bibliometrics). Reference: N N' in base R, N the topic-by-word matrix of
# document counts.

test_that("the raw strength between two topics equals the off-diagonal of N N'", {
  abstracts <- utils::head(covid_abstracts, 80L)
  hg <- text_hypergraph(abstracts, column = "abstract", id = "doc", stop_words = stop_words_en(), min_count = 3L)
  topics <- hg_cluster(hg, k = 5L, seed = 1L)
  cells <- hg_get(hg, what = "weights")
  documents <- unclass(table(topics$cluster[match(cells$doc, topics$node)], cells$word))
  strength <- tcrossprod(documents)
  expected <- ref_pairs_upper(strength, name = "weight")
  names(expected) <- c("source", "target", "weight")
  # Products and sums of document counts: integers.
  expect_reference(topic_network(hg, topics, similarity = "none"), expected,
                   reference = "independent: N N' of the topic-by-word document counts, base R", covers = "topic_network()",
                   quantity = "five clusters, raw strength", tolerance = 0, key = c("source", "target"))
})
