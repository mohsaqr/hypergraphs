.hgat_fixture <- function() {
  nodes <- c("d1", "d2", "t1", "e1")
  A <- matrix(0, 4, 4, dimnames = list(nodes, nodes))
  A["d1", c("t1", "e1")] <- 1
  A["d2", "t1"] <- 1
  A <- A + t(A)
  types <- c(d1 = "document", d2 = "document", t1 = "topic", e1 = "entity")
  features <- list(
    document = matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE,
                      dimnames = list(c("d1", "d2"), NULL)),
    topic = matrix(c(0.8, 0.2, 0), 1, 3,
                   dimnames = list("t1", NULL)),
    entity = matrix(c(0.1, 0.9, 0.2, 0.8), 1, 4,
                    dimnames = list("e1", NULL))
  )
  list(A = A, types = types, features = features)
}

test_that("HGAT aligns heterogeneous feature spaces and normalizes A+I", {
  x <- .hgat_fixture()
  input <- .thg_hgat_inputs(x$A, x$types, x$features)
  expect_identical(input$type_names, c("document", "entity", "topic"))
  expect_equal(input$features$document,
               x$features$document[c("d1", "d2"), , drop = FALSE])
  expect_true(isSymmetric(input$adjacency))
  expect_true(all(diag(input$adjacency) > 0))
})

test_that("heterogeneous HGAT trains across unequal feature widths", {
  skip_if_not_installed("torch")
  x <- .hgat_fixture()
  fit <- heterogeneous_hgat(
    x$A, x$types, x$features, labels = c(d1 = "a", d2 = "b"),
    hidden = 6, layers = 2, epochs = 5, validation = 0, seed = 1
  )
  expect_s3_class(fit, "data.frame")
  expect_equal(nrow(fit), 4L)
  expect_identical(attr(fit, "node_types"), x$types)
  expect_true(all(is.finite(attr(fit, "history")$loss)))
})

test_that("HGAT validates type-specific feature rows", {
  x <- .hgat_fixture()
  rownames(x$features$topic) <- "wrong"
  expect_error(.thg_hgat_inputs(x$A, x$types, x$features),
               "feature rows for type")
})

# ---- audit regressions (2026-10-06) -----------------------------------------

test_that("HGAT refuses duplicated feature rows and node types (TXT-15)", {
  x <- .hgat_fixture()
  dup_rows <- x$features
  dup_rows$document <- rbind(dup_rows$document,
                             d1 = c(100, 100))
  expect_error(.thg_hgat_inputs(x$A, x$types, dup_rows),
               class = "hypergraphs_bad_input")
  dup_types <- c(x$types, d1 = "topic")
  expect_error(.thg_hgat_inputs(x$A, dup_types, x$features),
               class = "hypergraphs_bad_input")
  # a permutation of unique rows aligns to the adjacency order
  shuffled <- x$features
  shuffled$document <- shuffled$document[c("d2", "d1"), , drop = FALSE]
  expect_identical(.thg_hgat_inputs(x$A, x$types, shuffled)$features,
                   .thg_hgat_inputs(x$A, x$types, x$features)$features)
})

test_that("heterogeneous HGAT refuses conflicting labels (TXT-07)", {
  skip_if_not_installed("torch")
  x <- .hgat_fixture()
  expect_error(heterogeneous_hgat(x$A, x$types, x$features,
                                  labels = c(d1 = "a", d1 = "b", d2 = "b"),
                                  epochs = 1L, hidden = 4L),
               class = "hypergraphs_bad_input")
})
