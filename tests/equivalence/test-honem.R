# honem(): the higher-order network embedding of Saebi et al. (2020). With D the row-normalised matrix of the
# higher-order network, the neighbourhood matrix is S = sum_{k = 0..L} exp(-k) D^(k + 1) / sum_k exp(-k), and the
# embedding is U sqrt(Sigma) of its truncated singular value decomposition. Reference: the formula with base
# svd(), the network read from the rules of hon(). The embedding is compared through its Gram matrix E E' =
# U Sigma U', which does not depend on the signs of the singular vectors.

sequences <- ref_second_order(n_sequences = 40L, length = 20L, seed = 9L)

#' The node-by-node matrix of the higher-order network, from its rules: a rule leads from the node `from` to the
#' node made of the last `to_order` states of its path
network_matrix <- function(rules) {
  states <- strsplit(rules$path, " -> ", fixed = TRUE)
  target <- mapply(\(path, order) paste(utils::tail(path, order), collapse = " -> "), states, rules$to_order)
  nodes <- sort(unique(c(rules$from, target)))
  counts <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  counts[cbind(rules$from, target)] <- rules$count
  counts
}

test_that("singular values and the embedding's Gram matrix equal the SVD of the neighbourhood matrix", {
  network <- hon(sequences, max_order = 2L, min_freq = 1L)
  counts <- network_matrix(hg_get(network, what = "rules"))
  transition <- counts / pmax(rowSums(counts), 1)
  max_power <- 5L
  weights <- exp(-(0:max_power))
  powers <- Reduce(\(power, k) power %*% transition, seq_len(max_power), transition, accumulate = TRUE)
  neighbourhood <- Reduce(`+`, Map(`*`, weights, powers)) / sum(weights)
  decomposition <- svd(neighbourhood)
  dim <- 3L
  embedding <- decomposition$u[, seq_len(dim)] %*% diag(sqrt(decomposition$d[seq_len(dim)]))
  gram <- tcrossprod(embedding)
  dimnames(gram) <- list(rownames(counts), rownames(counts))
  fit <- honem(network, dim = dim, max_power = max_power)
  # Dense SVD of a matrix of order about ten on both sides.
  expect_reference(subset(hg_get(fit, what = "variance"), select = c(dimension, singular_value)),
                   data.frame(dimension = seq_len(dim), singular_value = decomposition$d[seq_len(dim)]),
                   reference = "independent: SVD of the neighbourhood matrix of Saebi et al. (2020), base R",
                   covers = c("honem()", "hg_get.net_honem()", "hg_get.hypergraphs_result()"), quantity = "singular values",
                   tolerance = 1e-10, key = "dimension")
  coordinates <- hg_get(fit, what = "embeddings")
  ours <- tcrossprod(as.matrix(subset(coordinates, select = -node)))
  dimnames(ours) <- list(coordinates$node, coordinates$node)
  expect_reference(ours, gram, reference = "independent: SVD of the neighbourhood matrix of Saebi et al. (2020), base R",
                   covers = c("honem()", "hg_get.net_honem()"), quantity = "Gram matrix of the embedding", tolerance = 1e-10)
})
