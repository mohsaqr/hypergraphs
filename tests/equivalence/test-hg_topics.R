# hg_topics() and the methods of its result: the mixed-membership topic model fitted by the multiplicative updates
# that minimise the generalised Kullback-Leibler divergence (Lee and Seung 1999, 2001), the maximum-likelihood fit
# of probabilistic latent semantic analysis (Hofmann 1999; Gaussier and Goutte 2005). Its starting factors are drawn
# inside the function from `seed`, with no argument to supply them, so it is checked in three ways:
#   - properties that hold for any start: the shares and the word probabilities are distributions, the reported
#     divergence is the divergence of the reported model recomputed from its definition, and the divergence never
#     increases from one update to the next (Lee and Seung 2001, Theorem 2);
#   - recovery of planted topics on a separable corpus generated here (helper-topics.R), across 12 seeds;
#   - the distribution of the final divergence over 20 random starts against scikit-learn's NMF with the same
#     divergence and solver (beta_loss = "kullback-leibler", solver = "mu") on the same corpus.
# The tables of the fit are read with hg_get() and summary(); hypergraph() of a fit counts the documents' topic sets.

recovery_seeds <- 1:12
planted <- lapply(recovery_seeds, \(seed) {
  corpus <- ref_topic_corpus(seed)
  hg <- group_hypergraph(corpus$memberships, node = "doc", hyperedge = "word", weight = "count")
  list(seed = seed, corpus = corpus, fit = hg_topics(hg, k = 3L, nstart = 5L, seed = seed))
})

test_that("the topic shares and the word probabilities are distributions, for every start", {
  shares <- do.call(rbind, lapply(planted, \(p) transform(hg_get(p$fit, what = "shares"), seed = p$seed)))
  words <- do.call(rbind, lapply(planted, \(p) transform(hg_get(p$fit, what = "words", n = Inf), seed = p$seed)))
  expect_reference(subset(shares, select = c(seed, node, topic, share)), transform(subset(shares, select = c(seed, node, topic, share)), share = pmax(share, 0)),
                   reference = "independent: P(z | d) is non-negative", covers = c("hg_topics()", "hg_get.net_hg_topics()"),
                   quantity = "12 seeds, topic share of every document", tolerance = 0, key = c("seed", "node", "topic"))
  totals <- stats::aggregate(share ~ seed + node, data = shares, FUN = sum)
  expect_reference(totals, transform(totals, share = 1), reference = "independent: P(z | d) sums to one over the topics",
                   covers = c("hg_topics()", "hg_get.net_hg_topics()"), quantity = "12 seeds, sum of a document's topic shares",
                   tolerance = 1e-12, key = c("seed", "node"))
  expect_reference(subset(words, select = c(seed, topic, word, probability)), transform(subset(words, select = c(seed, topic, word, probability)), probability = pmax(probability, 0)),
                   reference = "independent: P(w | z) is non-negative", covers = c("hg_topics()", "hg_get.net_hg_topics()"),
                   quantity = "12 seeds, word probability of every topic", tolerance = 0, key = c("seed", "topic", "word"))
  mass <- stats::aggregate(probability ~ seed + topic, data = words, FUN = sum)
  expect_reference(mass, transform(mass, probability = 1), reference = "independent: P(w | z) sums to one over the words",
                   covers = c("hg_topics()", "hg_get.net_hg_topics()"), quantity = "12 seeds, sum of a topic's word probabilities",
                   tolerance = 1e-12, key = c("seed", "topic"))
})

test_that("the reported divergence is the generalised KL divergence of the reported model", {
  rows <- lapply(1:5, \(seed) {
    corpus <- ref_topic_corpus(seed)
    hg <- group_hypergraph(corpus$memberships, node = "doc", hyperedge = "word", weight = "count")
    # tol = 0 runs all 2000 updates: the row totals of W H equal the documents' token counts at the fixed point,
    # which the PLSA form of the fit below assumes
    fit <- hg_topics(hg, k = 3L, nstart = 1L, seed = seed, tol = 0, max_iter = 2000L)
    fitted <- ref_plsa_fitted(corpus$X, hg_get(fit, what = "shares"), hg_get(fit, what = "words", n = Inf))
    data.frame(seed = seed, ours = subset(hg_get(fit, what = "restarts"), best)$divergence, theirs = ref_kl(corpus$X, fitted))
  })
  rows <- do.call(rbind, rows)
  expect_reference(data.frame(seed = rows$seed, divergence = rows$ours), data.frame(seed = rows$seed, divergence = rows$theirs),
                   reference = "independent: generalised KL divergence of the PLSA fit from its shares and word probabilities (Lee and Seung 2001)",
                   covers = c("hg_topics()", "hg_get.net_hg_topics()"), quantity = "5 seeds, divergence of the kept start",
                   tolerance = 1e-8, measure = "relative", key = "seed")
})

test_that("the divergence never increases from one update to the next", {
  corpus <- ref_topic_corpus(1L)
  hg <- group_hypergraph(corpus$memberships, node = "doc", hyperedge = "word", weight = "count")
  # one start (seed 1), stopped after a growing number of updates: every fit is a truncation of the same path
  updates <- c(1L, 2L, 5L, 10L, 20L, 50L, 100L, 200L, 500L)
  path <- vapply(updates, \(m) hg_get(hg_topics(hg, k = 3L, nstart = 1L, seed = 1L, tol = 0, max_iter = m), what = "restarts")$divergence,
                 numeric(1))
  expect_reference(data.frame(max_iter = updates, divergence = path), data.frame(max_iter = updates, divergence = cummin(path)),
                   reference = "independent: the KL divergence is non-increasing under the multiplicative updates (Lee and Seung 2001, Theorem 2)",
                   covers = "hg_topics()", quantity = "divergence after a given number of updates, running minimum",
                   tolerance = 1e-12, measure = "relative", key = "max_iter")
})

test_that("the topics of a separable corpus recover the planted topics, on every one of 12 seeds", {
  skip_without_reference("clue")
  matched <- lapply(planted, \(p) {
    words <- hg_get(p$fit, what = "words", n = Inf)
    recovered <- ref_wide(words, "topic", "word", "probability", sort(unique(words$topic)), p$corpus$words)
    # each recovered topic matched to one planted topic by the Hungarian method on the L1 distance of the word
    # distributions (Kuhn 1955)
    cost <- outer(seq_len(nrow(recovered)), seq_len(nrow(p$corpus$H)),
                  Vectorize(\(a, b) sum(abs(recovered[a, ] - p$corpus$H[b, ]))))
    assignment <- as.integer(clue::solve_LSAP(cost))
    list(seed = p$seed, recovered = ref_long(`rownames<-`(recovered, assignment), name = "probability", row = "topic", column = "word"),
         planted = ref_long(p$corpus$H, name = "probability", row = "topic", column = "word"),
         tokens = min(colSums(p$corpus$theta)) * 1000)
  })
  invisible(lapply(matched, \(m) {
    # Tolerance fixed from the sampling error before the run: a planted probability p estimated from N tokens has
    # standard error sqrt(p (1 - p) / N) <= 0.5 / sqrt(N); N is the expected token count of the smallest topic
    # (about 19,000). Six such bounds cover the largest of the 90 errors of a seed and the extra error of the mixed
    # documents (study of 40 seeds before this test: largest error 0.0087 against a bound of 0.022; 40 of 40 within).
    expect_reference(m$recovered, m$planted, reference = "independent: the planted topics of the simulated corpus",
                     covers = c("hg_topics()", "hg_get.net_hg_topics()"),
                     quantity = sprintf("Monte Carlo summary: seed %d, word probabilities of the topics matched to the planted ones", m$seed),
                     tolerance = 6 * 0.5 / sqrt(m$tokens), key = c("topic", "word"))
  }))
  # a negative control: the same topics against the planted ones in another order must differ, on each topic's own
  # anchor words (a word shared by all topics, or an anchor of a third topic, is near zero in both orders)
  first <- matched[[1L]]
  shifted <- transform(first$planted, topic = as.character((as.integer(topic) %% 3L) + 1L))
  own_anchors <- \(table) subset(table, startsWith(word, paste0("t", topic, "w")))
  expect_differs(own_anchors(first$recovered), own_anchors(shifted), reference = "independent: the planted topics of the simulated corpus", covers = "hg_topics()",
                 quantity = "seed 1, recovered topics against the planted topics in another order", tolerance = 6 * 0.5 / sqrt(first$tokens),
                 key = c("topic", "word"))
})

test_that("the final divergence over 20 starts is distributed as scikit-learn's KL NMF", {
  skip_without_python(c("scikit-learn" = "sklearn", "numpy"))
  corpus <- ref_topic_corpus(1L)
  hg <- group_hypergraph(corpus$memberships, node = "doc", hyperedge = "word", weight = "count")
  # 20 starts (seeds 1 to 20), each run for all 1000 updates
  ours <- hg_get(hg_topics(hg, k = 3L, nstart = 20L, seed = 1L, tol = 0, max_iter = 1000L), what = "restarts")$divergence
  script <- c(
    "import sys, json, warnings, numpy as np",
    "from sklearn.decomposition import NMF",
    "spec = json.load(open(sys.argv[1]))",
    "X = np.array(spec['X'], dtype=float)",
    "out = []",
    "for seed in spec['seeds']:",
    "    with warnings.catch_warnings():",
    "        warnings.simplefilter('ignore')  # ConvergenceWarning: tol = 0 runs every update by design",
    "        model = NMF(n_components=3, beta_loss='kullback-leibler', solver='mu', init='random', random_state=int(seed),",
    "                    max_iter=int(spec['max_iter']), tol=0.0)",
    "        W = model.fit_transform(X)",
    "    WH = W @ model.components_",
    "    positive = X > 0",
    "    out.append(float(np.sum(X[positive] * np.log(X[positive] / WH[positive])) - X.sum() + WH.sum()))",
    "json.dump({'divergence': out}, open(sys.argv[2], 'w'))")
  theirs <- ref_python_run(script, list(X = unname(corpus$X), seeds = 1:20, max_iter = 1000L))$divergence
  probs <- c(minimum = 0, lower_quartile = 0.25, median = 0.5)
  # Monte Carlo summary of two samples of 20 starts. The starts that reach the maximum-likelihood fit end on the same
  # divergence, so the Monte Carlo error of these quantiles is below 1e-12 unless a quarter of the starts miss the
  # optimum. In a study before this test (3 corpora x 5 blocks of 20 seeds per side) 0 to 2 starts per block stopped
  # at a local optimum on either side, and the three quantiles agreed in 15 of 15 blocks, the largest relative
  # difference 1.1e-14. 1e-9 relative allows for what is left of the convergence after 1000 updates and for the order
  # of summation in R and numpy.
  # The upper quantiles are not compared: they hold the local optima, whose frequency depends on each solver's start.
  expect_reference(stats::quantile(ours, probs, names = FALSE), stats::quantile(theirs, probs, names = FALSE),
                   reference = "scikit-learn", covers = c("hg_topics()", "hg_get.net_hg_topics()"),
                   quantity = "Monte Carlo summary: quantile of the final KL divergence over 20 random starts", tolerance = 1e-9,
                   measure = "relative", key = names(probs),
                   note = sprintf("%s; NMF(beta_loss = 'kullback-leibler', solver = 'mu', init = 'random', tol = 0, max_iter = 1000)",
                                  ref_python_version("scikit-learn")))
})

test_that("summary() of a fit holds the topic prevalence and the dominant topics of its shares", {
  fit <- planted[[1L]]$fit
  tables <- summary(fit)
  shares <- hg_get(tables, what = "shares")
  prevalence <- stats::aggregate(share ~ topic, data = shares, FUN = mean)
  expected <- data.frame(topic = prevalence$topic, prevalence = prevalence$share,
                         documents = stats::aggregate(share ~ topic, data = shares, FUN = sum)$share)
  expect_reference(subset(hg_get(tables, what = "topics"), select = c(topic, prevalence, documents)), expected,
                   reference = "independent: mean and sum of the topic shares over the documents (Roberts et al. 2014)",
                   covers = c("summary.net_hg_topics()", "hg_get.net_hg_topics()"), quantity = "seed 1, topic prevalence",
                   tolerance = 1e-12, key = "topic")
  dominant <- stats::aggregate(share ~ node, data = shares, FUN = max)
  expect_reference(subset(hg_get(tables, what = "documents"), select = c(node, share)), dominant,
                   reference = "independent: the largest topic share of each document", covers = c("summary.net_hg_topics()", "hg_get.net_hg_topics()"),
                   quantity = "seed 1, share of the dominant topic", tolerance = 0, key = "node")
})

test_that("hypergraph() of a fit counts the sets of topics each document holds at the threshold", {
  fit <- planted[[1L]]$fit
  shares <- hg_get(fit, what = "shares")
  # threshold 0.3 < 1/3: every document holds at least one topic, so every document is counted
  held <- subset(shares, share >= 0.3)
  sets <- vapply(split(held$topic, held$node), \(topics) paste(sort(topics), collapse = " + "), character(1))
  counted <- as.data.frame(table(set = sets), stringsAsFactors = FALSE)
  expected <- data.frame(group = "All documents", set = counted$set, count = counted$Freq, share = counted$Freq / length(unique(shares$node)),
                         size = lengths(strsplit(counted$set, " + ", fixed = TRUE)))
  expect_reference(hg_get(hypergraph(fit, threshold = 0.3, top = Inf), what = "sets"), expected,
                   reference = "independent: sets of topics with share at least the threshold, counted over the documents (Cassi et al. 2017)",
                   covers = c("hypergraph()", "hypergraph.net_hg_topics()", "group_hypergraph()", "hg_get.net_hg()"),
                   quantity = "seed 1, topic sets at threshold 0.3", tolerance = 1e-12, key = c("group", "set"))
})
