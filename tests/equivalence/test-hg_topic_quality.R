# hg_topic_quality(), hg_topic_sizes(): the quality and size of the topics of a clustered corpus. NPMI coherence
# (Bouma 2009; Lau et al. 2014) of a topic's top words on document co-occurrence, over the whole corpus ("npmi")
# or the topic's own documents ("npmi_cluster"); topic sizes as counts and shares of the clustered documents.
# References: text2vec::coherence(metrics = "mean_npmi") and the counts in base R. text2vec smooths every joint
# probability by 1e-12, so it is compared on the topics whose word pairs all co-occur in some but not all of the
# reference documents (where the smoothing moves a value by at most about 1e-12 / p).

abstracts <- utils::head(covid_abstracts, 80L)
hg <- text_hypergraph(abstracts, column = "abstract", id = "doc", stop_words = stop_words_en(), min_count = 3L)
topics <- hg_cluster(hg, k = 5L, seed = 1L)
documents <- hg_get(hg, what = "weights")

test_that("topic sizes equal the counts of the clustered documents", {
  counted <- table(topics$cluster)
  expect_reference(hg_topic_sizes(hg, topics),
                   data.frame(topic = names(counted), n = as.vector(counted), share = as.vector(counted) / nrow(topics)),
                   reference = "independent: counts of the clustered documents, base R", covers = "hg_topic_sizes()",
                   quantity = "five clusters", tolerance = 1e-12, key = "topic")
})

test_that("corpus and in-cluster NPMI of given top words equal text2vec", {
  skip_without_reference("text2vec")
  presence <- unclass(table(documents$doc, documents$word)) > 0
  presence <- presence * 1
  # the five words of each topic with the largest count in its documents
  counts <- stats::aggregate(list(count = documents$count), list(topic = topics$cluster[match(documents$doc, topics$node)], word = documents$word), sum)
  counts <- counts[order(counts$topic, -counts$count, counts$word), ]
  words <- do.call(rbind, lapply(split(counts, counts$topic), \(t) data.frame(topic = t$topic[1:5], word = t$word[1:5], rank = 1:5)))
  # the topics whose word pairs co-occur in some but not all of the reference documents
  interior <- \(present, top) {
    joint <- crossprod(present[, top, drop = FALSE])
    all(joint[upper.tri(joint)] > 0 & joint[upper.tri(joint)] < nrow(present))
  }
  npmi <- \(present, top) unname(text2vec::coherence(matrix(top, ncol = 1L), crossprod(present), metrics = "mean_npmi",
                                                     n_doc_tcm = nrow(present))[1L, "mean_npmi"])
  by_topic <- split(words$word, words$topic)
  corpus_topics <- Filter(\(topic) interior(presence, by_topic[[topic]]), names(by_topic))
  corpus <- hg_topic_quality(hg, words = words, n = 5L, coherence = "npmi", exclusivity = "none")
  # Logarithms of ratios of document counts; text2vec's smoothing moves them by about 1e-12 / p.
  expect_reference(subset(corpus, topic %in% corpus_topics, select = c(topic, coherence)),
                   data.frame(topic = corpus_topics, coherence = vapply(corpus_topics, \(topic) npmi(presence, by_topic[[topic]]), numeric(1))),
                   reference = "text2vec", covers = "hg_topic_quality()", quantity = "corpus NPMI", tolerance = 1e-8, key = "topic")
  own <- \(topic) presence[topics$node[topics$cluster == topic], , drop = FALSE]
  cluster_topics <- Filter(\(topic) interior(own(topic), by_topic[[topic]]), names(by_topic))
  in_cluster <- hg_topic_quality(hg, topics, words = words, n = 5L, coherence = "npmi_cluster", exclusivity = "none")
  expect_reference(subset(in_cluster, topic %in% cluster_topics, select = c(topic, coherence)),
                   data.frame(topic = cluster_topics, coherence = vapply(cluster_topics, \(topic) npmi(own(topic), by_topic[[topic]]), numeric(1))),
                   reference = "text2vec", covers = "hg_topic_quality()", quantity = "in-cluster NPMI", tolerance = 1e-8, key = "topic")
})
