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

# The same contexts from pathpy's higher-order networks (the k-gram counts of order 1 and of `order`) and the
# entropy package (entropy.plugin() and KL.plugin() in bits): no definition of the test's own beyond reading
# the distributions out of pathpy's edge weights.
from_pathpy <- function(sequences, order) {
  grams <- ref_pathpy_kgrams(sequences, unique(c(1L, order)))
  states <- sort(unique(unlist(sequences)))
  distribution <- \(k, context) {
    rows <- grams[grams$k == k & grams$context == context, ]
    counts <- stats::setNames(numeric(length(states)), states)
    counts[rows$following] <- rows$count
    counts
  }
  contexts <- unique(grams$context[grams$k == order])
  do.call(rbind, lapply(contexts, \(context) {
    n_k <- distribution(order, context)
    n_1 <- distribution(1L, utils::tail(strsplit(context, " -> ", fixed = TRUE)[[1L]], 1L))
    support <- n_1 > 0
    data.frame(context = context, count = sum(n_k),
               entropy_first_order = entropy::entropy.plugin(n_1, unit = "log2"),
               entropy_order_k = entropy::entropy.plugin(n_k, unit = "log2"),
               entropy_drop = entropy::entropy.plugin(n_1, unit = "log2") - entropy::entropy.plugin(n_k, unit = "log2"),
               kl = entropy::KL.plugin(n_k[support], n_1[support], unit = "log2"),
               top_first_order = names(which.max(n_1)), top_order_k = names(which.max(n_k)))
  }))
}

test_that("contexts, counts, entropies and divergences equal pathpy's k-gram counts and the entropy package", {
  skip_without_python("pathpy")
  skip_without_reference("entropy")
  version <- ref_python_version("pathpy")
  inputs <- list(second_order = sequences,
                 random = ref_random_sequences(30L, 20L, c("a", "b", "c", "d"), seed = 4L))
  invisible(Map(\(input, name) {
    theirs <- from_pathpy(input, 2L)
    ours <- hg_get(memory(input, order = 2L, min_count = 1L, base = 2), what = "contexts")
    expect_reference(ours, subset(theirs, select = c(context, count, top_first_order, top_order_k)), reference = "pathpy",
                     covers = c("memory()", "hg_get.net_path_dependence()"), quantity = paste(name, "contexts and counts"),
                     tolerance = 0, key = "context", note = paste(version, "HigherOrderNetwork(k) edge weights"))
    expect_reference(ours, subset(theirs, select = -c(count, top_first_order, top_order_k)), reference = "entropy",
                     covers = c("memory()", "hg_get.net_path_dependence()"),
                     quantity = paste(name, "entropies and divergence (bits)"), tolerance = 1e-12, key = "context",
                     note = "entropy::entropy.plugin() and KL.plugin() of pathpy's k-gram counts")
    expect_reference(hg_get(summary(memory(input, order = 2L, min_count = 1L, base = 2)), what = "contexts"),
                     subset(theirs, select = -c(top_first_order, top_order_k)), reference = "entropy",
                     covers = c("summary.hypergraphs_result()", "hg_get.hypergraphs_summary()"),
                     quantity = paste(name, "counts, entropies and divergence through summary()"), tolerance = 1e-12,
                     key = "context", note = "pathpy's k-gram counts, entropy::entropy.plugin() and KL.plugin()")
  }, inputs, names(inputs)))
})
