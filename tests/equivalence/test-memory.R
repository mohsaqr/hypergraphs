# memory(): the path dependence of a chain. For each context of `order` states seen at least `min_count` times,
# the next-state distribution after the context is set against the first-order distribution after its last state:
# their entropies (base 2), the drop between them, the Kullback-Leibler divergence of the context's distribution
# from the first-order one, the most likely next state under each and whether it flips. Reference: the
# definitions in base R on the k-gram counts.

sequences <- ref_second_order(n_sequences = 40L, length = 20L, seed = 9L)

#' One row per context from the k-gram counts
by_definition <- function(order) {
  continuation <- \(k) {
    grams <- ref_kgrams(sequences, k + 1L)
    parts <- strsplit(names(grams), " -> ", fixed = TRUE)
    data.frame(context = vapply(parts, \(p) paste(utils::head(p, -1L), collapse = " -> "), character(1)),
               following = vapply(parts, \(p) utils::tail(p, 1L), character(1)), count = unname(grams))
  }
  first <- continuation(1L)
  higher <- continuation(order)
  states <- sort(unique(unlist(sequences)))
  distribution <- \(table, context) {
    counts <- stats::setNames(numeric(length(states)), states)
    rows <- table[table$context == context, ]
    counts[rows$following] <- rows$count
    counts / sum(counts)
  }
  entropy <- \(p) -sum(p[p > 0] * log2(p[p > 0]))
  contexts <- unique(higher$context)
  do.call(rbind, lapply(contexts, \(context) {
    p_k <- distribution(higher, context)
    p_1 <- distribution(first, utils::tail(strsplit(context, " -> ", fixed = TRUE)[[1L]], 1L))
    kept <- p_k > 0
    data.frame(context = context, count = sum(higher$count[higher$context == context]),
               entropy_first_order = entropy(p_1), entropy_order_k = entropy(p_k),
               entropy_drop = entropy(p_1) - entropy(p_k), kl = sum(p_k[kept] * log2(p_k[kept] / p_1[kept])),
               top_first_order = names(which.max(p_1)), top_order_k = names(which.max(p_k)),
               flips = names(which.max(p_1)) != names(which.max(p_k)))
  }))
}

test_that("entropies, divergence and flips of every context equal their definitions", {
  fit <- memory(sequences, order = 2L, min_count = 1L, base = 2)
  # Sums of a few logarithms.
  expect_reference(hg_get(fit, what = "contexts"), by_definition(2L), reference = "independent: entropies and KL divergence of the k-gram distributions, base R",
                   covers = c("memory()", "hg_get.net_path_dependence()", "hg_get.hypergraphs_result()"),
                   quantity = "order 2 contexts", tolerance = 1e-12, key = "context")
  expect_reference(hg_get(summary(fit), what = "contexts"), by_definition(2L),
                   reference = "independent: entropies and KL divergence of the k-gram distributions, base R",
                   covers = c("summary.hypergraphs_result()", "hg_get.hypergraphs_summary()"),
                   quantity = "order 2 contexts through summary()", tolerance = 1e-12, key = "context")
})
