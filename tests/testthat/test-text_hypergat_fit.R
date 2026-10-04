test_that("a fitted network predicts new text without training or changing vocabulary", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  docs <- c(a = "Soup salt onion. Broth warm soup.",
            b = "Salt soup broth. Onion warm soup.",
            c = "Stars sky moon. Night bright stars.",
            d = "Moon stars sky. Bright night moon.")
  fit <- hg_hypergat(docs, labels = c(a = "food", b = "food", c = "sky", d = "sky"),
    embed_dim = 8L, hidden = 4L, epochs = 2L, validation = 0, seed = 9L)
  expect_s3_class(fit, "hg_hypergat")
  expect_s3_class(fit, "hg_classification")
  expect_output(print(fit), "No held-out evaluation")
  expect_true(is.na(hg_get(fit, what = "accuracy")$accuracy))
  expect_error(plot(fit), "held-out labels")
  expect_s3_class(plot(fit, type = "history"), "ggplot")
  # The raw reader contains neither model pointers nor attached fit tables.
  expect_named(attributes(hg_get(fit)), c("names", "row.names", "class"))
  expect_identical(hg_get(fit, what = "history"), attr(fit, "history"))
  expect_null(attr(fit, "hypergat")$cache$words)
  weights <- lapply(attr(fit, "hypergat")$model$state_dict(), as.array)
  # Fitting is unavailable during every reader and prediction below.
  testthat::local_mocked_bindings(.thg_hypergat_net = function(...) stop("refit"),
    .thg_hypergat_lda = function(...) stop("refit"))
  scored <- predict(fit, newdata = docs)
  expect_equal(scored$score, fit$score, tolerance = 1e-6)
  expect_identical(scored$predicted, fit$predicted)
  p <- predict(fit, newdata = docs, type = "probabilities")
  expect_equal(rowSums(p), rep(1, length(docs)), ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(p, predict(fit, type = "probabilities"), tolerance = 1e-6)
  # OOV text cannot create a new embedding or change another prediction.
  mixed <- c(docs[1], novel = "qwertyxyz zzzunknown", missing = NA_character_)
  expect_warning(result <- predict(fit, newdata = mixed),
    class = "hypergraphs_dropped_documents")
  expect_identical(result$node, names(mixed))
  expect_true(all(is.na(result$predicted[2:3])))
  expect_equal(result$score[1], scored$score[1], tolerance = 1e-6)
  expect_warning(probabilities <- predict(fit, newdata = mixed,
    type = "probabilities"), class = "hypergraphs_dropped_documents")
  expect_true(all(is.na(probabilities[2:3, ])))
  expect_identical(nrow(predict(fit, newdata = character())), 0L)
  expect_error(predict(fit, newdata = setNames(docs, rep("same", 4))),
    class = "hypergraphs_bad_input")
  words <- hg_get(fit, what = "hyperedge_words", node = "a")
  expect_identical(unique(words$node), "a")
  expect_null(attr(fit, "hypergat")$cache$words)
  expect_equal(lapply(attr(fit, "hypergat")$model$state_dict(), as.array), weights)
})

test_that("prediction keeps fitted columns, semantic edges and observed test labels", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  train <- data.frame(doc = c("a", "b", "c", "d"),
    text = c("Soup salt onion. Broth soup.", "Broth salt soup. Onion soup.",
             "Stars sky moon. Night sky.", "Moon stars night. Bright sky."),
    subject = c("food", "food", "sky", "sky"))
  fit <- hg_hypergat(train, labels = "subject", column = "text", id = "doc",
    semantic = "lda", lda_keywords = list(food = c("soup", "salt"), sky = "moon"),
    embed_dim = 8L, hidden = 4L, epochs = 2L, validation = 0)
  testthat::local_mocked_bindings(.thg_hypergat_net = function(...) stop("refit"),
    .thg_hypergat_lda = function(...) stop("refit"))
  predicted <- predict(fit, newdata = train, labels = "subject")
  expect_s3_class(predicted, "hg_classification")
  expect_equal(predicted$score, fit$score, tolerance = 1e-6)
  expect_true(all(predicted$split == "test"))
  expect_equal(hg_get(predicted, what = "accuracy")$accuracy,
               mean(predicted$predicted == train$subject))
  expect_identical(hg_get(predicted, what = "documents")$text, train$text)
  expect_s3_class(plot(predicted), "ggplot")
  new <- data.frame(doc = "new", text = "zzzzunseen", subject = "food")
  expect_warning(unscored <- predict(fit, newdata = new, labels = "subject"),
    class = "hypergraphs_dropped_documents")
  expect_identical(hg_get(unscored, what = "accuracy")$n_test, 1L)
  expect_equal(hg_get(unscored, what = "accuracy")$accuracy, 0)
  expect_identical(sum(hg_get(unscored, what = "confusion")$n), 1L)
  expect_true("(unscored)" %in% hg_get(unscored, what = "confusion")$predicted)
  expect_error(predict(fit, newdata = transform(train, subject = "newclass"),
    labels = "subject"), class = "hypergraphs_bad_input")
})

test_that("held-out vocabulary is frozen and unscorable held-out rows count as errors", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  docs <- c(a = "alpha uniquealpha", b = "beta uniquebeta",
            c = "gamma uniquegamma", d = "delta uniquedelta")
  labels <- c(a = "one", b = "one", c = "two", d = "two")
  expect_warning(fit <- hg_hypergat(docs, labels = labels, holdout = 0.5,
    embed_dim = 8L, hidden = 4L, epochs = 1L, validation = 0),
    class = "hypergraphs_dropped_documents")
  heldout <- hg_get(fit, split = "test")
  expect_identical(nrow(heldout), 2L)
  expect_true(all(is.na(heldout$predicted)))
  expect_true(all(!heldout$correct))
  expect_equal(hg_get(fit, what = "accuracy")$accuracy, 0)
  expect_identical(sum(hg_get(fit, what = "confusion")$n), 2L)
  expect_true(all(is.na(predict(fit, type = "probabilities")[heldout$node, ])))
})

test_that("exclusive words have edge weight one and the diagnostic makes it explicit", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  fit <- hg_hypergat(c(a = "Soup onion. Broth salt.", b = "Moon sky. Stars night."),
    labels = c(a = "food", b = "sky"), embed_dim = 8L, hidden = 4L,
    epochs = 1L, validation = 0)
  edges <- hg_get(fit, what = "hyperedges")
  words <- hg_get(fit, what = "hyperedge_words")
  expect_equal(words$edge_weight, rep(1, nrow(words)), tolerance = 1e-6)
  expect_equal(words$edge_uniform, rep(1, nrow(words)))
  expect_equal(edges$weight, edges$uniform_weight, tolerance = 1e-6)
  expect_equal(edges$exclusive_fraction, rep(1, nrow(edges)))
  expect_s3_class(plot(fit, type = "attention", node = "a"), "ggplot")
  hypergraph <- plot(fit, type = "hyperedges", node = "a")
  polygons <- Filter(function(l) inherits(l$geom, "GeomPolygon"), hypergraph$layers)
  expect_length(polygons, 1L)
  expect_setequal(as.character(unique(polygons[[1L]]$data$hyperedge)),
                  c("sentence 1", "sentence 2"))
  expect_error(plot(fit, type = "hyperedges", node = c("a", "b")),
    class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "hyperedges", split = "test"),
    class = "hypergraphs_bad_input")
  expect_error(hg_get(fit, what = "attention", node = "absent"),
    class = "hypergraphs_bad_input")
})


test_that("vocabulary filtering cannot silently remove a labelled class", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  docs <- c(a = "Soup salt. Soup salt.", b = "Soup salt. Soup salt.",
            c = "Moon sky. Moon sky.", d = "Moon sky. Moon sky.",
            e = "rareterm", f = "rareterm")
  expect_error(suppressWarnings(hg_hypergat(docs,
    labels = c(a = "food", b = "food", c = "sky", d = "sky", e = "rare", f = "rare"),
    holdout = 0.5, min_count = 2, epochs = 1L)),
    "Every labelled class needs usable training text",
    class = "hypergraphs_bad_input")
})
