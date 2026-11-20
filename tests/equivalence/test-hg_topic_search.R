# hg_topic_search(): fits hg_topics() for every number of topics in `k` and reports, per k, the mean UMass
# coherence (Mimno et al. 2011) and the mean FREX exclusivity (Bischof and Airoldi 2012) of the topics' top `n`
# words, the divergence of the kept start, and whether the model is on the frontier of the two criteria (no other
# k is higher on both; Roberts et al. 2014). The fits draw their starts inside from `seed`; the same call to
# hg_topics() with the same seed reproduces each fit, whose topics are then scored here: coherence from its
# definition, exclusivity by stm::exclusivity(), the divergence from its definition, the frontier from its
# definition. hg_topic_quality(topics = ) is scored the same way.

corpus <- ref_topic_corpus(5L, len = 60L)
# FREX ranks the words of a topic: the anchor words of a topic must not tie (equal counts give mathematically equal
# probabilities, whose order then rests on rounding). Seed 5 is a corpus whose anchor counts are distinct per topic.
stopifnot("the anchor words of a topic have distinct counts" =
            all(tapply(colSums(corpus$X)[seq_len(18L)], rep(1:3, each = 6L), \(counts) anyDuplicated(counts) == 0L)))
hg <- group_hypergraph(corpus$memberships, node = "doc", hyperedge = "word", weight = "count")
ks <- 2:5
search <- hg_topic_search(hg, k = ks, nstart = 2L, n = 8L, tol = 0, max_iter = 1500L, seed = 3L)
fits <- lapply(ks, \(k) hg_topics(hg, k = k, nstart = 2L, tol = 0, max_iter = 1500L, depth = 8L, seed = 3L))
present <- (corpus$X > 0) * 1

#' The word probabilities of a fit as a topics x words matrix, topics in label order
topic_matrix <- function(fit) {
  words <- hg_get(fit, what = "words", n = Inf)
  ref_wide(words, "topic", "word", "probability", sort(unique(words$topic)), colnames(corpus$X))
}

#' A fitted topic model in the form stm's diagnostics read: one log-beta matrix
stm_model <- function(probabilities) {
  structure(list(beta = list(logbeta = list(log(probabilities))),
                 settings = list(dim = list(K = nrow(probabilities), V = ncol(probabilities))),
                 vocab = colnames(probabilities)), class = "STM")
}

#' The top n words of every topic of a fit, in rank order
top_words <- function(fit, n = 8L) {
  words <- hg_get(fit, what = "words", n = n)
  split(words$word[order(words$topic, words$rank)], words$topic[order(words$topic, words$rank)])
}

test_that("the frontier is the set of k that no other k beats on both coherence and exclusivity", {
  expect_reference(subset(search, select = c(k, frontier)), data.frame(k = search$k, frontier = ref_frontier(search$coherence, search$exclusivity)),
                   reference = "independent: Pareto frontier of mean coherence and mean exclusivity (Roberts et al. 2014)",
                   covers = "hg_topic_search()", quantity = "on the frontier, k = 2 to 5", tolerance = 0, key = "k")
})

test_that("the coherence of every k is the mean UMass coherence of its topics' top words", {
  umass <- lapply(fits, \(fit) vapply(top_words(fit), \(top) ref_umass(present, top), numeric(1)))
  expect_reference(subset(search, select = c(k, coherence)), data.frame(k = ks, coherence = vapply(umass, mean, numeric(1))),
                   reference = "independent: UMass coherence of the top 8 words, Mimno et al. (2011)", covers = "hg_topic_search()",
                   quantity = "mean coherence, k = 2 to 5", tolerance = 1e-10, key = "k")
  quality <- hg_topic_quality(hg, topics = fits[[3L]], n = 8L)
  expect_reference(subset(quality, select = c(topic, coherence)), data.frame(topic = names(umass[[3L]]), coherence = unname(umass[[3L]])),
                   reference = "independent: UMass coherence of the top 8 words, Mimno et al. (2011)", covers = "hg_topic_quality()",
                   quantity = "coherence of each topic of a 4-topic model", tolerance = 1e-10, key = "topic")
})

test_that("the exclusivity of every k is the mean FREX of its topics, as stm computes it", {
  skip_without_reference("stm")
  frex <- lapply(fits, \(fit) stm::exclusivity(stm_model(topic_matrix(fit)), M = 8L, frexw = 0.7))
  expect_reference(subset(search, select = c(k, exclusivity)), data.frame(k = ks, exclusivity = vapply(frex, mean, numeric(1))),
                   reference = "stm", covers = "hg_topic_search()", quantity = "mean exclusivity (FREX, frexw = 0.7), k = 2 to 5",
                   tolerance = 1e-10, key = "k")
  quality <- hg_topic_quality(hg, topics = fits[[3L]], n = 8L)
  expect_reference(subset(quality, select = c(topic, exclusivity)), data.frame(topic = rownames(topic_matrix(fits[[3L]])), exclusivity = frex[[3L]]),
                   reference = "stm", covers = "hg_topic_quality()", quantity = "FREX of each topic of a 4-topic model",
                   tolerance = 1e-10, key = "topic")
})

test_that("stm's semantic coherence orients a word pair by the order of the words across all topics", {
  skip_without_reference("stm")
  model <- stm_model(topic_matrix(fits[[3L]]))
  documents <- lapply(seq_len(nrow(corpus$X)), \(d) {
    used <- which(corpus$X[d, ] > 0)
    rbind(as.integer(used), as.integer(corpus$X[d, used]))
  })
  quality <- hg_topic_quality(hg, topics = fits[[3L]], n = 8L)
  expect_known_difference(subset(quality, select = c(topic, coherence)),
                          data.frame(topic = rownames(topic_matrix(fits[[3L]])), coherence = stm::semanticCoherence(model, documents, M = 8L)),
                          reference = "stm", covers = "hg_topic_quality()", quantity = "UMass coherence of each topic, stm",
                          tolerance = 1e-10, key = "topic",
                          reason = paste("stm's semCoh1beta() conditions each pair on the word that comes later in the list of all",
                                         "topics' top words (first appearance across topics), not on the higher-ranked word of the",
                                         "topic as Mimno et al. (2011) define it; the two agree on a topic whose top words first",
                                         "appear in that topic, and differ on a topic that shares top words with an earlier one."),
                          agrees_with = "Mimno et al. (2011): each pair conditioned on the higher-ranked word (the independent UMass above)")
})

test_that("the divergence of every k is the generalised KL divergence of its kept model", {
  divergence <- vapply(fits, \(fit) ref_kl(corpus$X, ref_plsa_fitted(corpus$X, hg_get(fit, what = "shares"), hg_get(fit, what = "words", n = Inf))),
                       numeric(1))
  expect_reference(subset(search, select = c(k, divergence)), data.frame(k = ks, divergence = divergence),
                   reference = "independent: generalised KL divergence of the PLSA fit from its shares and word probabilities (Lee and Seung 2001)",
                   covers = "hg_topic_search()", quantity = "divergence of the kept start, k = 2 to 5", tolerance = 1e-8,
                   measure = "relative", key = "k")
})
