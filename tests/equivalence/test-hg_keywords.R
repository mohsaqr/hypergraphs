# hg_keywords(): the words of each cluster of a clustered corpus, scored by "mass" (the summed incidence weight the
# cluster's documents place on the word), "frequency" (the word's token count in the cluster's documents) or
# "ctfidf" (class-based tf-idf, Grootendorst 2022: the cluster's L1-normalised term frequency times
# log(A / f_w + 1), f_w the word's corpus count and A the average token count per cluster, truncated), with the
# share of the word's summed score that falls in the cluster and the cluster's documents containing it.
# Reference: the definitions on the document-word table in base R.

abstracts <- utils::head(covid_abstracts, 80L)
hg <- text_hypergraph(abstracts, column = "abstract", id = "doc", stop_words = stop_words_en(), min_count = 3L, weight = "tfidf")
topics <- hg_cluster(hg, k = 4L, seed = 1L)

test_that("mass, frequency and c-TF-IDF scores, shares and document counts equal their definitions", {
  cells <- hg_get(hg, what = "weights")
  cells$cluster <- topics$cluster[match(cells$doc, topics$node)]
  per_cluster <- stats::aggregate(list(mass = cells$weight, frequency = cells$count, n_docs = rep(1, nrow(cells))),
                                  list(cluster = cells$cluster, word = cells$word), sum)
  tokens <- stats::aggregate(list(total = per_cluster$frequency), list(cluster = per_cluster$cluster), sum)
  corpus <- stats::aggregate(list(corpus = per_cluster$frequency), list(word = per_cluster$word), sum)
  per_cluster <- merge(merge(per_cluster, tokens), corpus)
  average <- floor(mean(tokens$total))
  per_cluster$ctfidf <- per_cluster$frequency / per_cluster$total * log(average / per_cluster$corpus + 1)
  expected <- do.call(rbind, lapply(c("mass", "frequency", "ctfidf"), \(type) {
    score <- per_cluster[[type]]
    summed <- stats::ave(score, per_cluster$word, FUN = sum)
    data.frame(type = type, cluster = per_cluster$cluster, word = per_cluster$word, score = score, share = score / summed,
               n_docs = per_cluster$n_docs)
  }))
  ours <- hg_keywords(hg, topics, n = Inf, type = c("mass", "frequency", "ctfidf"))
  # Sums of tf-idf weights and one logarithm per cell.
  expect_reference(as.data.frame(subset(ours, select = c(type, cluster, word, score, share, n_docs))), subset(expected, score > 0),
                   reference = "independent: per-cluster sums of the document-word table, base R", covers = "hg_keywords()",
                   quantity = "four clusters, every word", tolerance = 1e-10, key = c("type", "cluster", "word"))
})
