# hg_topics(): the KL non-negative factorization topic model (Lee & Seung
# 2001; equivalent to PLSA, Gaussier & Goutte 2005), its restarts and the
# agreement of topics across starts (Greene et al. 2014).

.tm_toy <- function() {
  c(a = "soup salt onion soup broth", b = "salt soup broth onion",
    c = "stars sky moon night", d = "sky stars night moon moon",
    e = "soup stars salt sky night broth")
}

test_that("one update is the Lee-Seung KL update, computed densely", {
  set.seed(3)
  X <- matrix(rpois(6 * 5, 2), 6, 5)
  X[X == 0] <- 0
  Xs <- methods::as(X, "CsparseMatrix")
  k <- 2L
  W <- matrix(runif(6 * k, 0.5, 1.5), 6, k)
  H <- matrix(runif(k * 5, 0.5, 1.5), k, 5)
  fit <- .tm_fit(Xs, t(W), H, max_iter = 1L, tol = 0)
  # the update written with dense matrices, straight from the paper: W
  # first, then H with the new W; the model floored where X is non-zero
  floored <- function(W, H) {
    wh <- W %*% H
    wh[X > 0] <- pmax(wh[X > 0], .TM_EPSILON)
    wh
  }
  q <- ifelse(X > 0, X / floored(W, H), 0)
  W1 <- W * (q %*% t(H)) / matrix(rowSums(H), 6, k, byrow = TRUE)
  q1 <- ifelse(X > 0, X / floored(W1, H), 0)
  H1 <- H * (t(W1) %*% q1) / matrix(colSums(W1), k, 5)
  H1[H1 < .Machine$double.eps] <- 0
  expect_equal(fit$W, W1, tolerance = 1e-12)
  expect_equal(t(fit$Ht), H1, tolerance = 1e-12)
})

test_that("the divergence never increases and reaches zero on exact data", {
  set.seed(5)
  W0 <- matrix(rexp(20 * 3), 20, 3)
  H0 <- matrix(rexp(3 * 15), 3, 15)
  X <- methods::as(W0 %*% H0, "CsparseMatrix")
  start <- local({
    set.seed(9)
    .tm_start(X, 3L)
  })
  divergences <- vapply(c(1L, 5L, 20L, 80L, 400L), \(steps) {
    .tm_fit(X, start$Wt, start$Htt, max_iter = steps, tol = 0)$divergence
  }, numeric(1L))
  expect_false(is.unsorted(rev(divergences)))
  expect_lt(utils::tail(divergences, 1L), 1e-3 * divergences[1L])
})

test_that("two planted vocabularies give pure documents and a mixed one", {
  fit <- hg_topics(text_hypergraph(.tm_toy()), k = 2, nstart = 3)
  shares <- hg_get(fit, what = "shares")
  mixed <- subset(shares, node == "e")
  pure <- subset(shares, node %in% c("a", "b", "c", "d"))
  expect_true(all(pure$share < 1e-6 | pure$share > 1 - 1e-6))
  expect_equal(sort(mixed$share), c(0.5, 0.5), tolerance = 1e-3)
  # the documents of one vocabulary share their topic
  dominant <- hg_get(fit, what = "documents")
  expect_identical(subset(dominant, node == "a")$topic,
                   subset(dominant, node == "b")$topic)
  expect_false(identical(subset(dominant, node == "a")$topic,
                         subset(dominant, node == "c")$topic))
})

test_that("shares and word probabilities are distributions", {
  fit <- hg_topics(text_hypergraph(.tm_toy()), k = 2, nstart = 2)
  shares <- hg_get(fit, what = "shares")
  per_document <- tapply(shares$share, shares$node, sum)
  expect_equal(as.vector(per_document), rep(1, 5), tolerance = 1e-12)
  words <- hg_get(fit, what = "words", n = Inf)
  per_topic <- tapply(words$probability, words$topic, sum)
  expect_equal(as.vector(per_topic), rep(1, 2), tolerance = 1e-12)
  topics <- hg_get(fit)
  expect_equal(sum(topics$documents), 5, tolerance = 1e-10)
  expect_false(is.unsorted(rev(topics$prevalence)))
  expect_named(topics, c("topic", "prevalence", "documents", "agreement",
                         "top_words"))
})

test_that("the average Jaccard and the Hungarian assignment are exact", {
  # Greene et al. (2014): depth 1 {a}/{a} = 1, depth 2 {a,b}/{a,c} = 1/3,
  # depth 3 {a,b,c}/{a,c,b} = 1
  expect_equal(.tm_average_jaccard(c("a", "b", "c"), c("a", "c", "b")),
               mean(c(1, 1 / 3, 1)))
  set.seed(11)
  permutations <- function(n) {
    if (n == 1L) return(matrix(1L))
    smaller <- permutations(n - 1L)
    do.call(rbind, lapply(seq_len(n), \(first) {
      cbind(first, matrix(setdiff(seq_len(n), first)[smaller],
                          nrow(smaller)))
    }))
  }
  lapply(2:6, \(n) {
    cost <- matrix(runif(n * n), n, n)
    all <- permutations(n)
    brute <- min(apply(all, 1L, \(p) sum(cost[cbind(seq_len(n), p)])))
    assignment <- .tm_hungarian(cost)
    expect_setequal(assignment, seq_len(n))
    expect_equal(sum(cost[cbind(seq_len(n), assignment)]), brute,
                 tolerance = 1e-12)
  })
})

test_that("a start's agreement with itself is one", {
  fit <- hg_topics(text_hypergraph(.tm_toy()), k = 2, nstart = 3)
  restarts <- hg_get(fit, what = "restarts")
  expect_identical(subset(restarts, best)$agreement, 1)
  expect_true(all(restarts$agreement >= 0 & restarts$agreement <= 1))
  expect_identical(sum(restarts$best), 1L)
  expect_identical(subset(restarts, best)$divergence, min(restarts$divergence))
})

test_that("the seed reproduces the fit and leaves the caller's stream alone", {
  hg <- text_hypergraph(.tm_toy())
  set.seed(42)
  before <- .Random.seed
  a <- hg_topics(hg, k = 2, nstart = 2, seed = 7)
  expect_identical(.Random.seed, before)
  b <- hg_topics(hg, k = 2, nstart = 2, seed = 7)
  expect_identical(a, b)
  skip_on_os("windows")
  p <- hg_topics(hg, k = 2, nstart = 2, seed = 7, parallel = TRUE)
  expect_identical(p, a)
})

test_that("any hypergraph is factorized on its incidence", {
  hg <- group_hypergraph(
    data.frame(actor = c("a", "b", "c", "a", "b", "d", "e", "f", "e", "f"),
               group = c("g1", "g1", "g1", "g2", "g2", "g2", "g3", "g3",
                         "g4", "g4")),
    node = "actor", hyperedge = "group")
  fit <- hg_topics(hg, k = 2, nstart = 2)
  expect_identical(fit$n_documents, 6L)
  expect_identical(fit$n_words, 4L)
})

test_that("bad input is refused by class and non-convergence is reported", {
  hg <- text_hypergraph(.tm_toy())
  expect_error(hg_topics(hg, k = 1), class = "hypergraphs_bad_input")
  expect_error(hg_topics(hg, k = 5), class = "hypergraphs_bad_input")
  expect_error(hg_topics(hg, k = 2, nstart = 0), class = "hypergraphs_bad_input")
  expect_error(hg_topics(hg, k = 2, tol = -1), class = "hypergraphs_bad_input")
  expect_error(hg_topics(list(), k = 2))
  expect_warning(hg_topics(hg, k = 2, nstart = 1, max_iter = 2, tol = 1e-12),
                 class = "hypergraphs_no_converge")
  fit <- hg_topics(hg, k = 2, nstart = 2)
  expect_error(hg_get(fit, what = "words", topic = "Topic 9"),
               class = "hypergraphs_bad_input")
})

test_that("the readers, print, summary and plot", {
  fit <- hg_topics(text_hypergraph(.tm_toy()), k = 2, nstart = 2)
  expect_identical(nrow(hg_get(fit, what = "words", n = 3)), 6L)
  expect_error(hg_get(fit, what = "words", n = 2.5),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "words", n = -Inf),
               class = "hypergraphs_bad_input")
  expect_identical(unique(hg_get(fit, what = "shares", topic = "Topic 1")$topic),
                   "Topic 1")
  expect_output(print(fit), "Topic model \\(KL factorization\\): 2 topics")
  s <- summary(fit)
  expect_s3_class(s, "hypergraphs_summary")
  expect_identical(s$shares, hg_get(fit, what = "shares"))
  expect_s3_class(plot(fit, n = 3), "ggplot")
})

test_that("SymNMF memberships are distributions whose maximum is the cluster", {
  hg <- text_hypergraph(.tm_toy())
  membership <- hg_cluster(hg, k = 2, algorithm = "symnmf", seed = 1,
                           what = "membership")
  expect_named(membership, c("node", "cluster", "membership"))
  per_node <- tapply(membership$membership, membership$node, sum)
  expect_equal(as.vector(per_node), rep(1, 5), tolerance = 1e-12)
  clusters <- hg_cluster(hg, k = 2, algorithm = "symnmf", seed = 1)
  strongest <- do.call(rbind, lapply(split(membership, membership$node),
                                     \(d) d[which.max(d$membership), ]))
  expect_identical(strongest$cluster[match(clusters$node, strongest$node)],
                   clusters$cluster)
  expect_error(hg_cluster(hg, k = 2, what = "membership"),
               class = "hypergraphs_bad_input")
})

test_that("hg_topic_quality() scores a topic model from its distributions", {
  hg <- text_hypergraph(.tm_toy())
  fit <- hg_topics(hg, k = 2, nstart = 2)
  quality <- hg_topic_quality(hg, topics = fit, n = 3)
  expect_identical(quality$topic, c("Topic 1", "Topic 2"))
  expect_equal(quality$size, hg_get(fit)$documents)
  # the same top words scored through `words` give the same coherence
  words <- hg_get(fit, what = "words", n = 3)
  by_words <- hg_topic_quality(hg, words = words, n = 3, exclusivity = "none")
  expect_equal(quality$coherence, by_words$coherence)
  # FREX by hand from P(w | z) (Bischof & Airoldi 2012, as stm computes it)
  all_words <- hg_get(fit, what = "words", n = Inf)
  beta <- rbind(
    subset(all_words, topic == "Topic 1")$probability[
      order(subset(all_words, topic == "Topic 1")$word)],
    subset(all_words, topic == "Topic 2")$probability[
      order(subset(all_words, topic == "Topic 2")$word)])
  vocabulary <- sort(unique(all_words$word))
  keep <- colSums(beta) > 0
  beta <- beta[, keep]
  vocabulary <- vocabulary[keep]
  share <- sweep(beta, 2L, colSums(beta), `/`)
  frex <- vapply(1:2, \(t) {
    ex <- rank(share[t, ]) / ncol(beta)
    fr <- rank(beta[t, ]) / ncol(beta)
    value <- 1 / (0.7 / ex + 0.3 / fr)
    sum(value[match(subset(words, topic == sprintf("Topic %d", t))$word,
                    vocabulary)])
  }, numeric(1L))
  expect_equal(quality$exclusivity, frex)
  expect_error(hg_topic_quality(hg, topics = list()),
               class = "hypergraphs_bad_input")
  expect_error(hg_topic_quality(hg, topics = fit, words = words),
               class = "hypergraphs_bad_input")
  expect_error(hg_topic_quality(hg, topics = fit, coherence = "npmi_cluster"),
               class = "hypergraphs_bad_input")
  other <- text_hypergraph(c(x = "one two three", y = "two three four",
                             z = "four five one"))
  expect_error(hg_topic_quality(other, topics = fit),
               class = "hypergraphs_bad_input")
})

test_that("hg_topic_search() scores every k and marks the frontier", {
  corpus <- c(
    a = "soup salt onion soup broth", b = "salt soup broth onion",
    c = "stars sky moon night", d = "sky stars night moon moon",
    e = "soup stars salt sky night broth", f = "onion salt stars broth")
  hg <- text_hypergraph(corpus)
  search <- hg_topic_search(hg, k = 3:2, nstart = 2, n = 4)
  expect_s3_class(search, "hypergraphs_topic_search")
  expect_identical(search$k, 2:3)
  expect_named(search, c("k", "coherence", "exclusivity", "divergence",
                         "agreement", "converged", "frontier"))
  # each row is the quality of the single fit with the same settings
  fit <- hg_topics(hg, k = 3, nstart = 2, depth = 4)
  quality <- hg_topic_quality(hg, topics = fit, n = 4)
  expect_equal(subset(search, k == 3)$coherence, mean(quality$coherence))
  expect_equal(subset(search, k == 3)$exclusivity, mean(quality$exclusivity))
  # the frontier: no row on it is beaten on both measures
  on <- subset(search, frontier)
  expect_gte(nrow(on), 1L)
  expect_false(any(vapply(seq_len(nrow(on)), \(r) {
    any(search$coherence > on$coherence[r] &
          search$exclusivity > on$exclusivity[r])
  }, logical(1L))))
  skip_on_os("windows")
  expect_identical(hg_topic_search(hg, k = 2:3, nstart = 2, n = 4,
                                   parallel = TRUE), search)
  expect_s3_class(plot(search), "ggplot")
  expect_error(hg_topic_search(hg, k = 2), class = "hypergraphs_bad_input")
})

.tm_network_corpus <- function() {
  c(a = "soup salt onion soup broth", b = "salt soup broth onion",
    c = "stars sky moon night", d = "sky stars night moon moon",
    e = "soup stars salt sky night broth", f = "onion salt stars broth",
    g = "train rail track station", h = "rail station train ticket",
    i = "soup train salt rail")
}

test_that("topic_network() on a topic model: thresholded co-occurrence counts", {
  hg <- text_hypergraph(.tm_network_corpus())
  fit <- hg_topics(hg, k = 3, nstart = 3)
  shares <- hg_get(fit, what = "shares")
  edges <- topic_network(hg, topics = fit, threshold = 0.2)
  # by hand: documents in which both topics reach the threshold
  present <- subset(shares, share >= 0.2)
  pair_count <- function(t1, t2) {
    length(intersect(subset(present, topic == t1)$node,
                     subset(present, topic == t2)$node))
  }
  labels <- hg_get(fit)$topic
  pairs <- t(combn(labels, 2))
  expected <- data.frame(source = pairs[, 1], target = pairs[, 2],
                         weight = mapply(pair_count, pairs[, 1], pairs[, 2],
                                         USE.NAMES = FALSE))
  expected <- subset(expected, weight > 0)
  rownames(expected) <- NULL
  expect_equal(edges, expected)
  # cosine normalisation divides by the topics' document counts
  cosine <- topic_network(hg, topics = fit, threshold = 0.2,
                       similarity = "cosine")
  n_present <- table(factor(present$topic, levels = labels))
  expect_equal(cosine$weight,
               expected$weight / sqrt(as.numeric(n_present[expected$source]) *
                                        as.numeric(n_present[expected$target])))
})

test_that("topic_network() on a topic model: stm's simple topic correlation", {
  hg <- text_hypergraph(.tm_network_corpus())
  fit <- hg_topics(hg, k = 3, nstart = 3)
  shares <- hg_get(fit, what = "shares")
  labels <- hg_get(fit)$topic
  theta <- vapply(labels, \(t) subset(shares, topic == t)$share,
                  numeric(length(unique(shares$node))))
  correlation <- cor(theta)
  edges <- topic_network(hg, topics = fit)
  expected <- subset(
    data.frame(source = labels[row(correlation)[upper.tri(correlation)]],
               target = labels[col(correlation)[upper.tri(correlation)]],
               weight = correlation[upper.tri(correlation)]),
    weight > 0.01)
  rownames(expected) <- NULL
  expect_equal(edges, expected)
})

test_that("topic_network() refuses inconsistent input by class", {
  hg <- text_hypergraph(.tm_network_corpus())
  fit <- hg_topics(hg, k = 3, nstart = 2)
  clusters <- hg_cluster(hg, k = 3, seed = 1)
  expect_error(topic_network(hg), class = "hypergraphs_bad_input")
  expect_error(topic_network(hg, clusters = clusters, topics = fit),
               class = "hypergraphs_bad_input")
  expect_error(topic_network(hg, clusters = clusters, threshold = 0.2),
               class = "hypergraphs_bad_input")
  expect_error(topic_network(hg, topics = fit, similarity = "cosine"),
               class = "hypergraphs_bad_input")
  expect_error(topic_network(hg, topics = fit, threshold = 2),
               class = "hypergraphs_bad_input")
  expect_error(topic_network(text_hypergraph(.tm_toy()), topics = fit),
               class = "hypergraphs_bad_input")
  net <- topic_network(hg, topics = fit, threshold = 0.2, what = "network")
  expect_s3_class(net, "cograph_network")
})

test_that("group_hypergraph() counts the topic combinations of documents", {
  hg <- text_hypergraph(.tm_network_corpus())
  fit <- hg_topics(hg, k = 3, nstart = 3)
  shares <- hg_get(fit, what = "shares")
  combos <- group_hypergraph(fit, threshold = 0.2, top = Inf)
  sets <- hg_get(combos, what = "sets")
  # by hand: each document's set of topics at the threshold
  present <- subset(shares, share >= 0.2)
  by_doc <- lapply(split(present$topic, present$node), \(v) {
    paste(v[order(as.integer(sub("Topic ", "", v)))], collapse = " + ")
  })
  counted <- table(unlist(by_doc))
  expect_setequal(sets$set, names(counted))
  expect_equal(sets$count, as.integer(counted[sets$set]))
  expect_equal(sets$share, sets$count / 9)
  expect_false(is.unsorted(rev(sets$count)))
  # min_size keeps the multi-topic documents, shares still over all nine
  multi <- hg_get(group_hypergraph(fit, threshold = 0.2, top = Inf,
                                   min_size = 2), what = "sets")
  expect_true(all(multi$size >= 2))
  expect_equal(multi$share, multi$count / 9)
  # grouped by a column of the documents table
  meta_hg <- text_hypergraph(
    data.frame(id = names(.tm_network_corpus()),
               text = unname(.tm_network_corpus()),
               half = rep(c("first", "second"), c(5, 4))),
    column = "text", id = "id")
  meta_fit <- hg_topics(meta_hg, k = 3, nstart = 2)
  halves <- group_hypergraph(meta_fit, threshold = 0.2, group = "half",
                             top = Inf)
  expect_setequal(unique(hg_get(halves, what = "sets")$group),
                  c("first", "second"))
  expect_output(print(halves), "within each half")
})

test_that("group_hypergraph() refuses bad topic-combination arguments", {
  hg <- text_hypergraph(.tm_network_corpus())
  fit <- hg_topics(hg, k = 3, nstart = 2)
  expect_error(group_hypergraph(fit), class = "hypergraphs_bad_input")
  expect_error(group_hypergraph(fit, threshold = 0), class = "hypergraphs_bad_input")
  expect_error(group_hypergraph(fit, threshold = 0.2, group = "nope"),
               class = "hypergraphs_bad_input")
  expect_error(group_hypergraph(fit, threshold = 0.2, min_size = 0),
               class = "hypergraphs_bad_input")
  expect_error(group_hypergraph(data.frame(actor = "a", group = "g"),
                                node = "actor", hyperedge = "group",
                                threshold = 0.2),
               class = "hypergraphs_bad_input")
})

test_that("hg_sequences(topics =) uses each document's main topic", {
  corpus <- data.frame(
    id = names(.tm_network_corpus()), text = unname(.tm_network_corpus()),
    person = rep(c("p1", "p2", "p3"), each = 3), turn = rep(1:3, 3))
  hg <- text_hypergraph(corpus, column = "text", id = "id")
  fit <- hg_topics(hg, k = 3, nstart = 2)
  seqs <- hg_sequences(hg, topics = fit, actor = "person", order_by = "turn")
  main <- hg_get(fit, what = "documents")
  expect_identical(seqs, hg_sequences(hg, main, actor = "person",
                                      order_by = "turn", state = "topic"))
  expect_identical(nrow(seqs), 9L)
  expect_error(hg_sequences(hg, main, topics = fit, actor = "person",
                            order_by = "turn"),
               class = "hypergraphs_bad_input")
})

test_that("hg_get(what = \"prevalence\") gives the mean topic share by group", {
  corpus <- c(
    a = "soup salt onion soup broth", b = "salt soup broth onion",
    c = "stars sky moon night", d = "sky stars night moon moon",
    e = "soup stars broth sky", f = "onion moon salt night"
  )
  hg <- text_hypergraph(corpus)
  fit <- hg_topics(hg, k = 2, nstart = 2L, seed = 1L)
  # groups of three, so a mean and a median of the shares differ
  groups <- c(a = "food", b = "food", e = "food", c = "sky", d = "sky",
              f = "sky")
  prev <- hg_get(fit, what = "prevalence", group = groups)
  expect_named(prev, c("topic", "group", "prevalence", "documents", "n",
                       "top_words"))
  expect_identical(nrow(prev), 2L * 2L)
  # the formula: the mean of the documents' shares, against a hand count
  shares <- hg_get(fit, what = "shares")
  food_t1 <- mean(subset(shares, topic == "Topic 1" & node %in% c("a", "b", "e"))$share)
  expect_equal(subset(prev, topic == "Topic 1" & group == "food")$prevalence,
               food_t1, tolerance = 1e-12)
  # invariant: each group's prevalences sum to one, sums = mean * size
  expect_equal(as.numeric(tapply(prev$prevalence, prev$group, sum)), rep(1, 2),
               tolerance = 1e-12)
  expect_equal(prev$documents, prev$prevalence * prev$n, tolerance = 1e-12)
  # a tidy data.frame works as group, and missing documents are left out
  tidy <- data.frame(node = c("a", "b", "c"), label = c("x", "x", "y"))
  expect_setequal(hg_get(fit, what = "prevalence", group = tidy)$group,
                  c("x", "y"))
  expect_error(hg_get(fit, what = "prevalence"),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "shares", group = groups),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "prevalence", group = c(z = "q")),
               class = "hypergraphs_bad_input")
})

# ---- audit regressions (2026-10-06) -----------------------------------------

test_that("explicit sparse zeros do not change the KL-NMF fit (TXT-10)", {
  dense <- matrix(c(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1), 3, 3, byrow = TRUE)
  implicit <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1, dims = c(3, 3))
  explicit <- Matrix::sparseMatrix(i = c(1, 2, 3, 1), j = c(1, 2, 3, 2),
                                   x = c(1, 1, 1, 0), dims = c(3, 3))
  expect_identical(length(explicit@x), 4L)   # the zero is really stored
  set.seed(1)
  start <- .tm_start(implicit, 2L)
  fits <- lapply(list(methods::as(dense, "CsparseMatrix"), implicit, explicit),
                 \(X) .tm_fit(X, start$Wt, start$Htt, 50L, 1e-5))
  divergences <- vapply(fits, \(f) f$divergence, numeric(1L))
  expect_true(all(is.finite(divergences)))
  expect_equal(fits[[3L]]$divergence, fits[[2L]]$divergence, tolerance = 1e-12)
  expect_equal(fits[[3L]]$W, fits[[2L]]$W, tolerance = 1e-12)
  expect_equal(fits[[1L]]$Ht, fits[[3L]]$Ht, tolerance = 1e-12)
  # the divergence skips zero cells: 0 log 0 = 0
  expect_equal(.tm_divergence(c(1, 0), c(0.5, 0.25), matrix(1), matrix(1)),
               log(2) - 1 + 1)
})

test_that("a topic model is checked against the corpus it was fitted on (TXT-11)", {
  docs <- .tm_network_corpus()
  hg <- text_hypergraph(docs)
  fit <- hg_topics(hg, k = 2, nstart = 1, max_iter = 50L, tol = 0)
  expect_s3_class(hg_topic_quality(hg, topics = fit), "hypergraphs_topic_quality")
  # same texts and words under other document ids
  renamed <- text_hypergraph(stats::setNames(unname(docs),
                                             paste0("w", seq_along(docs))))
  expect_error(hg_topic_quality(renamed, topics = fit),
               class = "hypergraphs_bad_input")
  expect_error(topic_network(renamed, topics = fit),
               class = "hypergraphs_bad_input")
  # same ids and words, different counts
  recounted <- docs
  recounted[["a"]] <- paste(recounted[["a"]], "soup soup")
  expect_error(hg_topic_quality(text_hypergraph(recounted), topics = fit),
               class = "hypergraphs_bad_input")
  # the same corpus with its documents in another order is the same corpus
  reordered <- text_hypergraph(rev(docs))
  expect_identical(hg_topic_quality(reordered, topics = fit),
                   hg_topic_quality(hg, topics = fit))
  # a model saved before the signature is still checked on ids and words
  legacy <- fit
  legacy$corpus <- NULL
  expect_s3_class(hg_topic_quality(hg, topics = legacy),
                  "hypergraphs_topic_quality")
  expect_error(hg_topic_quality(renamed, topics = legacy),
               class = "hypergraphs_bad_input")
  # counts that are negative are refused before fitting
  plain <- group_hypergraph(data.frame(node = c("a", "a", "b", "b", "c", "c"),
                                       hyperedge = c("x", "y", "x", "z", "y", "z")),
                            node = "node", hyperedge = "hyperedge")
  plain$incidence["a", "y"] <- -1
  expect_error(hg_topics(plain, k = 2, nstart = 1),
               class = "hypergraphs_bad_input")
})

test_that("constant topic shares give isolated topics, never NA edges (TXT-13)", {
  hg <- text_hypergraph(.tm_network_corpus())
  fit <- hg_topics(hg, k = 3, nstart = 1, max_iter = 50L, tol = 0)
  constant <- fit
  constant$shares$share <- 1 / 3
  expect_warning(edges <- topic_network(hg, topics = constant),
                 class = "hypergraphs_constant_topics")
  expect_identical(nrow(edges), 0L)
  expect_warning(net <- topic_network(hg, topics = constant, what = "network"),
                 class = "hypergraphs_constant_topics")
  expect_identical(net$nodes$name, fit$topics$topic)
  expect_identical(nrow(net$edges), 0L)
  # one constant topic among varying ones: the varying pair keeps its value
  mixed <- fit
  t3 <- mixed$shares$topic == "Topic 3"
  mixed$shares$share[t3] <- 0.2
  expect_warning(mixed_edges <- topic_network(hg, topics = mixed),
                 class = "hypergraphs_constant_topics")
  expect_false(anyNA(mixed_edges))
  expect_false(any(c(mixed_edges$source, mixed_edges$target) == "Topic 3"))
  shares <- hg_get(mixed, what = "shares")
  wide <- stats::reshape(shares, idvar = "node", timevar = "topic",
                         direction = "wide")
  expected <- stats::cor(wide$`share.Topic 1`, wide$`share.Topic 2`)
  pair <- subset(mixed_edges, source == "Topic 1" & target == "Topic 2")
  if (expected > 0.01) {
    expect_equal(pair$weight, expected)
  } else {
    expect_identical(nrow(pair), 0L)
  }
})

test_that("prevalence refuses a document given two groups (TXT-07)", {
  hg <- text_hypergraph(.tm_network_corpus())
  fit <- hg_topics(hg, k = 2, nstart = 1, max_iter = 50L, tol = 0)
  groups <- c(a = "x", a = "y", b = "x", c = "y")
  expect_error(hg_get(fit, what = "prevalence", group = groups),
               class = "hypergraphs_bad_input")
})

test_that("the expected counts are the floored per-cell dot products", {
  set.seed(4)
  X <- Matrix::rsparsematrix(40, 70, density = 0.2,
                             rand.x = \(n) stats::rpois(n, 3) + 1)
  i <- X@i + 1L
  j <- rep.int(seq_len(ncol(X)), diff(X@p))
  st <- .tm_start(X, 6L)
  # spread the scale so some cells fall under the floor
  st$Wt[, 1:5] <- st$Wt[, 1:5] * 1e-9
  by_cell <- vapply(seq_along(i), \(cell)
    max(sum(st$Wt[, i[cell]] * st$Htt[, j[cell]]), .TM_EPSILON), numeric(1))
  expect_equal(.tm_fitted(st$Wt, st$Htt, i, j), by_cell)
  expect_true(any(by_cell == .TM_EPSILON))
})

test_that("all-pairs topic similarity equals the pairwise average Jaccard", {
  set.seed(8)
  vocabulary <- sprintf("w%02d", 1:40)
  for (k in c(1L, 3L, 9L)) for (depth in c(1L, 4L, 10L)) {
    a <- replicate(k, sample(vocabulary, depth), simplify = FALSE)
    b <- replicate(k, sample(vocabulary, depth), simplify = FALSE)
    pairwise <- outer(seq_len(k), seq_len(k), Vectorize(
      \(x, y) .tm_average_jaccard(a[[x]], b[[y]])))
    expect_identical(.tm_similarity(a, b), pairwise)
  }
  # identical lists agree perfectly; disjoint lists not at all
  same <- list(c("a", "b", "c"), c("d", "e", "f"))
  expect_identical(diag(.tm_similarity(same, same)), c(1, 1))
  expect_identical(.tm_similarity(same, rev(same))[1, 1], 0)
})
