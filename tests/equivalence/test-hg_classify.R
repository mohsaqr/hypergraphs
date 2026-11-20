# hg_classify(): transductive label spreading on a hypergraph (Zhou, Huang and Scholkopf 2007). With S = I - L
# the similarity operator of the Zhou Laplacian and Y the indicators of the given labels, the scores are
# F = (1 - xi) (I - xi S)^-1 Y; a node takes the class of its largest score, `score` is that value and `margin`
# its distance to the second. With `holdout`, a share of the labels is hidden and predicted from the others.
# Reference: the closed form written from the paper in base R.

set.seed(61)
blocks <- rep(1:3, each = 6L)
nodes <- sprintf("n%02d", seq_len(18L))
edges <- c(lapply(seq_len(30L), \(i) sample(nodes[blocks == ((i - 1L) %% 3L + 1L)], 3L)), list(c("n01", "n07"), c("n08", "n13")))
memberships <- data.frame(member = unlist(edges), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges))))
hg <- group_hypergraph(memberships, node = "member", hyperedge = "edge")

#' Label spreading by the closed form: one row per node with the predicted class, its score and margin
by_definition <- function(labels, xi = 0.99) {
  incidence <- ref_incidence(memberships, node = "member", edge = "edge") > 0
  degree_root <- diag(1 / sqrt(rowSums(incidence)))
  similarity <- degree_root %*% incidence %*% diag(1 / colSums(incidence)) %*% t(incidence) %*% degree_root
  given <- labels[rownames(incidence)]
  classes <- sort(unique(stats::na.omit(given)))
  indicators <- vapply(classes, \(class) as.numeric(!is.na(given) & given == class), numeric(nrow(incidence)))
  scores <- (1 - xi) * solve(diag(nrow(incidence)) - xi * similarity, indicators)
  ordered <- t(apply(scores, 1L, sort, decreasing = TRUE))
  data.frame(node = rownames(incidence), predicted = classes[max.col(scores, ties.method = "first")],
             score = ordered[, 1L], margin = ordered[, 1L] - ordered[, 2L])
}

test_that("predicted classes, scores and margins equal the closed form of label spreading", {
  labels <- stats::setNames(rep(NA_character_, 18L), nodes)
  labels[c("n01", "n02", "n07", "n08", "n13", "n14")] <- c("x", "x", "y", "y", "z", "z")
  # A linear solve of order 18 with xi = 0.99 (condition number near 100): 1e-10 on scores of order 0.01 to 1.
  expect_reference(hg_classify(hg, labels = stats::na.omit(labels), xi = 0.99), by_definition(labels),
                   reference = "independent: closed form F = (1 - xi) (I - xi S)^-1 Y of Zhou et al. (2007), base R",
                   covers = "hg_classify()", quantity = "six seed labels", tolerance = 1e-10, key = "node")
})

test_that("held-out predictions and accuracy equal the closed form from the training labels", {
  truth <- stats::setNames(c("x", "y", "z")[blocks], nodes)
  found <- hg_classify(hg, labels = truth, xi = 0.9, holdout = 0.5, seed = 1L)
  table <- hg_get(found)
  # the draw of the held-out labels is the package's; the reference spreads the labels it kept
  training <- stats::setNames(ifelse(table$split == "train", table$label, NA_character_), table$node)
  expected <- by_definition(training, xi = 0.9)
  tested <- subset(merge(table, expected, by = "node", suffixes = c("", "_reference")), split == "test")
  expect_reference(subset(table, split == "test", select = c(node, predicted, score, margin)),
                   subset(expected, node %in% tested$node), reference = "independent: closed form of Zhou et al. (2007), base R",
                   covers = c("hg_classify()", "hg_get.hg_classification()"), quantity = "held-out half, predictions",
                   tolerance = 1e-10, key = "node")
  accuracy <- mean(tested$predicted_reference == tested$label)
  expect_reference(hg_get(found, what = "accuracy")$accuracy, accuracy,
                   reference = "independent: closed form of Zhou et al. (2007), base R",
                   covers = c("hg_classify()", "hg_get.hg_classification()"), quantity = "held-out accuracy", tolerance = 1e-12)
  # summary() is the table of the classes: test documents, correct predictions, recall and precision
  classes <- sort(unique(tested$label))
  by_class <- data.frame(class = classes,
                         n_test = vapply(classes, \(class) sum(tested$label == class), numeric(1)),
                         correct = vapply(classes, \(class) sum(tested$label == class & tested$predicted_reference == class), numeric(1)),
                         recall = vapply(classes, \(class) mean(tested$predicted_reference[tested$label == class] == class), numeric(1)),
                         precision = vapply(classes, \(class) mean(tested$label[tested$predicted_reference == class] == class), numeric(1)))
  expect_reference(summary(found), by_class, reference = "independent: closed form of Zhou et al. (2007), base R",
                   covers = "summary.hg_classification()", quantity = "held-out half, per class", tolerance = 1e-12, key = "class")
})
