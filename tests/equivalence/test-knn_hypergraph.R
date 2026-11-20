# knn_hypergraph(): one hyperedge per item, holding the item and its k nearest neighbours by cosine similarity,
# with the similarity as the member's weight (the item itself weighs 1). References: the cosine similarity
# matrix in base R, and HyperG::knn_hypergraph() for the binary incidence (on unit-length rows Euclidean and
# cosine order the neighbours alike).

set.seed(42)
ids <- sprintf("it%02d", seq_len(30L))
embeddings <- matrix(abs(stats::rnorm(30L * 5L)), nrow = 30L, dimnames = list(ids, NULL))
unit <- embeddings / sqrt(rowSums(embeddings^2))

test_that("members and cosine weights equal the k nearest neighbours by cosine similarity", {
  k <- 3L
  similarity <- tcrossprod(unit)
  expected <- do.call(rbind, lapply(ids, \(center) {
    others <- setdiff(ids, center)
    nearest <- utils::head(others[order(-similarity[center, others])], k)
    data.frame(node = c(center, nearest), hyperedge = center, weight = c(1, similarity[center, nearest]))
  }))
  # Cosine similarities of five-dimensional rows: rounding of a few sums.
  expect_reference(hg_get(knn_hypergraph(embeddings, k = k, weight = "cosine"), what = "memberships"), expected,
                   reference = "independent: cosine similarity matrix, base R", covers = c("knn_hypergraph()", "hg_get.net_hg()"),
                   quantity = "k = 3 memberships and cosine weights", tolerance = 1e-12, key = c("node", "hyperedge"))
})

test_that("the binary incidence equals HyperG::knn_hypergraph()", {
  skip_without_reference("HyperG")
  k <- 3L
  reference <- HyperG::knn_hypergraph(unit, k = k)
  incidence <- t(as.matrix(HyperG::incidence_matrix(reference)))
  dimnames(incidence) <- list(ids, ids)
  expected <- subset(ref_long(incidence > 0, name = "member", row = "node", column = "hyperedge"), member)
  ours <- transform(hg_get(knn_hypergraph(unit, k = k, weight = "binary"), what = "memberships"), member = TRUE)
  expect_reference(ours, expected, reference = "HyperG", covers = "knn_hypergraph()",
                   quantity = "k = 3 binary memberships", tolerance = 0, key = c("node", "hyperedge"),
                   note = "HyperG builds hyperedge i from item i and its k nearest neighbours (Euclidean), in input order.")
})
