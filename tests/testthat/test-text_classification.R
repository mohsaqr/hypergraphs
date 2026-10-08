# hg_classify(labels = <column>, holdout =) and the hg_classification result

corpus <- data.frame(
  text = c("soup salt onion broth", "salt soup onion", "broth soup salt",
           "onion salt broth", "salt broth night", "soup onion warm",
           "stars sky moon night", "sky stars moon", "night sky stars",
           "moon night sky", "stars moon bright", "sky night dark"),
  theme = rep(c("food", "sky"), each = 6)
)

test_that("labels can name a column of the document table", {
  hg <- text_hypergraph(corpus, column = "text")
  by_column <- hg_classify(hg, labels = "theme")
  # the documents are doc_1, doc_2, ... in the order of the corpus rows
  by_vector <- hg_classify(hg, labels = stats::setNames(
    corpus$theme, paste0("doc_", seq_len(nrow(corpus)))))
  expect_identical(by_column, by_vector)
  expect_error(hg_classify(hg, labels = "nope"), class = "hypergraphs_bad_input")
})

test_that("holdout hides a stratified share and scores it", {
  hg <- text_hypergraph(corpus, column = "text")
  fit <- hg_classify(hg, labels = "theme", holdout = 1 / 3, seed = 2L)
  expect_s3_class(fit, "hg_classification")
  expect_s3_class(fit, "data.frame")
  test <- subset(hg_get(fit), split == "test")
  # stratified: two of the six documents of each theme are held out
  expect_identical(as.integer(table(test$label)), c(2L, 2L))
  # the held-out labels were hidden from the fit, so a fit on the known
  # labels alone gives the same predictions
  known <- subset(hg_get(fit), split == "train")
  direct <- hg_classify(hg, labels = stats::setNames(known$label, known$node))
  expect_identical(subset(direct, node %in% test$node)$predicted,
                   test$predicted)
  # the scores, which shift whenever a hidden label leaks into the fit
  expect_equal(subset(direct, node %in% test$node)$score, test$score,
               tolerance = 1e-12)
  leaked <- hg_classify(hg, labels = "theme")
  expect_false(isTRUE(all.equal(subset(leaked, node %in% test$node)$score,
                                test$score)))
  # the evaluation tables agree with the per-document table
  accuracy <- hg_get(fit, what = "accuracy")
  expect_identical(accuracy$n_test, 4L)
  expect_equal(accuracy$accuracy, mean(test$predicted == test$label))
  classes <- hg_get(fit, what = "classes")
  expect_equal(accuracy$balanced_accuracy, mean(classes$recall))
  confusion <- hg_get(fit, what = "confusion")
  expect_identical(sum(confusion$n), 4L)
  expect_identical(accuracy$chance, 0.5)
  # the caller's random stream is untouched
  set.seed(9); before <- runif(1)
  set.seed(9); invisible(hg_classify(hg, labels = "theme", holdout = 0.3))
  expect_identical(runif(1), before)
})

test_that("holdout output prints, summarises and plots", {
  hg <- text_hypergraph(corpus, column = "text")
  fit <- hg_classify(hg, labels = "theme", holdout = 1 / 3)
  expect_output(print(fit), "balanced accuracy")
  expect_s3_class(summary(fit), "data.frame")
  expect_s3_class(plot(fit), "ggplot")
  expect_error(hg_classify(hg, labels = "theme", holdout = 1.5))
  expect_error(hg_classify(hg, labels = "theme", holdout = 0.3,
                           seed = NA))
})

test_that("hg_hypergat() takes a label column and a holdout", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  docs <- transform(corpus, text = paste0(text, ". ", text, "."))
  fit <- hg_hypergat(docs, labels = "theme", column = "text", epochs = 2L,
                     embed_dim = 8L, hidden = 4L, validation = 0,
                     holdout = 1 / 3, seed = 3L)
  expect_s3_class(fit, "hg_classification")
  test <- subset(hg_get(fit), split == "test")
  expect_identical(as.integer(table(test$label)), c(2L, 2L))
  expect_identical(hg_get(fit, what = "accuracy")$n_test, 4L)
  expect_error(hg_hypergat(docs, labels = "nope", column = "text",
                           epochs = 1L), class = "hypergraphs_bad_input")
})

test_that("keywords, topic sizes and prevalence take a document column", {
  hg <- text_hypergraph(corpus, column = "text")
  tidy <- data.frame(node = paste0("doc_", seq_len(nrow(corpus))),
                     label = corpus$theme)
  expect_identical(hg_keywords(hg, "theme", n = 3),
                   hg_keywords(hg, tidy, n = 3))
  expect_identical(hg_topic_sizes(hg, "theme"), hg_topic_sizes(hg, tidy))
  topics <- hg_topics(hg, k = 2, nstart = 2L, seed = 1L)
  expect_identical(hg_get(topics, what = "prevalence", group = "theme"),
                   hg_get(topics, what = "prevalence", group = tidy))
  expect_error(hg_get(topics, what = "prevalence", group = "nope"),
               class = "hypergraphs_bad_input")
})

test_that("hg_get() filters, sorts and truncates the predictions", {
  hg <- text_hypergraph(corpus, column = "text")
  fit <- hg_classify(hg, labels = "theme", holdout = 1 / 2, seed = 4L)
  all_rows <- hg_get(fit)
  test <- hg_get(fit, split = "test")
  expect_identical(test$node, all_rows$node[all_rows$split %in% "test"])
  right <- hg_get(fit, split = "test", correct = TRUE)
  expect_true(all(right$correct))
  # top applies within each true label, after the ordering
  best <- hg_get(fit, split = "test", sort_by = "margin", top = 1)
  by_hand <- do.call(rbind, lapply(split(test, test$label), \(d) {
    utils::head(d[order(-d$margin, d$node), , drop = FALSE], 1L)
  }))
  rownames(by_hand) <- NULL
  expect_identical(best[order(best$label), ], by_hand)
  expect_identical(nrow(best), 2L)
  # a table with a node column selects those documents
  expect_identical(hg_get(fit, node = best)$node, sort(best$node))
  expect_error(hg_get(fit, what = "accuracy", split = "test"),
               class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, correct = NA), class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "hyperedges"),
               class = "hypergraphs_bad_input")
})

test_that("a hg_hypergat() holdout keeps the attention of its own fit", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  docs <- transform(corpus, text = paste0(text, ". ", text, " night."))
  settings <- list(column = "text", epochs = 2L, embed_dim = 8L, hidden = 4L,
                   validation = 0, seed = 3L)
  fit <- do.call(hg_hypergat, c(list(docs, labels = "theme", holdout = 1 / 3),
                                settings))
  known <- hg_get(fit, split = "train")
  # the same network refitted on the known labels alone gives the same
  # attention tables, so the stored ones belong to this fit
  direct_fit <- do.call(hg_hypergat, c(list(docs, labels = stats::setNames(
    known$label, known$node)), settings))
  direct <- hg_get(direct_fit, what = "hyperedges")
  stored <- hg_get(fit, what = "hyperedges")
  expect_equal(stored[setdiff(names(stored), c("label", "predicted"))],
               direct[setdiff(names(direct), c("label", "predicted"))],
               tolerance = 1e-6)
  words <- hg_get(fit, what = "hyperedge_words")
  expect_equal(stats::aggregate(word_weight ~ node + hyperedge, words,
                                sum)$word_weight,
               rep(1, nrow(stored)), tolerance = 1e-5)
  example <- hg_get(fit, split = "test", sort_by = "margin", top = 1)
  heaviest <- hg_get(fit, what = "hyperedges", node = example, top = 1)
  expect_identical(nrow(heaviest), nrow(example))
  one_edge <- hg_get(fit, what = "hyperedge_words", node = heaviest)
  expect_identical(unique(paste(one_edge$node, one_edge$hyperedge)),
                   paste(heaviest$node, heaviest$hyperedge))
  expect_s3_class(plot(fit, type = "attention", top = 2), "ggplot")
})

test_that("plot() draws topic prevalence by a document column", {
  hg <- text_hypergraph(corpus, column = "text")
  topics <- hg_topics(hg, k = 2, nstart = 2L, seed = 1L)
  expect_s3_class(plot(topics, type = "prevalence", group = "theme"),
                  "ggplot")
  expect_error(plot(topics, group = "theme"), class = "hypergraphs_bad_input")
})

# ---- audit regressions (2026-10-06) -----------------------------------------

test_that("the confusion table's no-prediction column never takes a class name (TXT-19)", {
  predictions <- data.frame(
    node = c("a", "b", "c", "d", "e", "f"),
    predicted = c("(unscored)", "x", "x", NA, "(unscored)", "x"),
    score = 1, margin = 1, stringsAsFactors = FALSE)
  labels <- c(a = "(unscored)", b = "x", c = "(unscored)", d = "x",
              e = "(unscored)", f = "x")
  hidden <- labels[c("a", "c", "d", "f")]
  result <- .thg_classification(predictions, labels, hidden,
                                method = "hg_classify", holdout = 0.5)
  confusion <- hg_get(result, what = "confusion")
  expect_identical(sum(confusion$n), 4L)
  expect_setequal(unique(confusion$label), c("(unscored)", "x"))
  expect_setequal(unique(confusion$predicted),
                  c("(unscored)", "x", "(unscored).1"))
  missing_row <- subset(confusion, label == "x" & predicted == "(unscored).1")
  expect_identical(missing_row$n, 1L)
  real_row <- subset(confusion, label == "(unscored)" & predicted == "(unscored)")
  expect_identical(real_row$n, 1L)
  # without a colliding class the column keeps its plain name
  plain_labels <- c(a = "y", b = "x", c = "y", d = "x", e = "y", f = "x")
  plain <- .thg_classification(transform(predictions, predicted = c("y", "x", "x", NA, "y", "x")),
                               plain_labels, plain_labels[c("a", "c", "d", "f")],
                               method = "hg_classify", holdout = 0.5)
  expect_true("(unscored)" %in% hg_get(plain, what = "confusion")$predicted)
})

test_that("hg_classify refuses conflicting labels before the held-out draw (TXT-07)", {
  hg <- text_hypergraph(c(a = "apple pear night", b = "apple peach",
                          c = "star moon night", d = "star sky"))
  conflict <- c(a = "food", a = "space", b = "food", c = "space", d = "space")
  expect_error(hg_classify(hg, conflict, holdout = 0.5),
               class = "hypergraphs_bad_input")
  expect_error(hg_classify(hg, data.frame(node = c("a", "a", "c"),
                                          label = c("food", "space", "space"))),
               class = "hypergraphs_bad_input")
})
