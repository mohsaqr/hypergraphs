# hon(): the higher-order network of Xu, Wickramarathne and Chawla (2016). Each rule leads from a context (a path
# of one or more states) to a next state; its count is the number of times the context is followed by that state
# in the sequences, and its probability that count over all continuations of the context. References: the
# k-gram counts in base R for every rule the network keeps, and markovchain's maximum-likelihood transition matrix
# for the rules out of single states. Which contexts become higher-order (the KL rule extraction) is not compared:
# no established implementation of BuildHON is installed.

sequences <- ref_second_order()

#' The count and probability of each rule, from the k-grams of the sequences
by_definition <- function(rules) {
  states <- strsplit(rules$path, " -> ", fixed = TRUE)
  contexts <- strsplit(rules$from, " -> ", fixed = TRUE)
  counts <- unlist(lapply(sort(unique(lengths(states))), \(k) ref_kgrams(sequences, k)))
  count <- unname(counts[rules$path])
  # all continuations of a context: the k-grams one longer that begin with it
  continued <- vapply(contexts, \(context) {
    longer <- ref_kgrams(sequences, length(context) + 1L)
    starts <- vapply(strsplit(names(longer), " -> ", fixed = TRUE), \(gram) paste(utils::head(gram, -1L), collapse = " -> "), character(1))
    sum(longer[starts == paste(context, collapse = " -> ")])
  }, numeric(1))
  data.frame(path = rules$path, count = count, probability = count / continued)
}

test_that("rule counts and probabilities equal the k-gram counts of the sequences", {
  network <- hon(sequences, max_order = 2L, min_freq = 1L)
  rules <- hg_get(network, what = "rules")
  # Counts are integers; a probability is one division.
  expect_reference(rules, by_definition(rules), reference = "independent: k-gram counts of the sequences, base R",
                   covers = c("hon()", "hg_get.net_hon()", "hg_get.hypergraphs_result()"),
                   quantity = "second-order chain, rules", tolerance = 1e-12, key = "path")
  summary_rules <- hg_get(summary(network), what = "rules")
  expect_reference(summary_rules, by_definition(summary_rules), reference = "independent: k-gram counts of the sequences, base R",
                   covers = c("summary.hypergraphs_result()", "hg_get.hypergraphs_summary()"),
                   quantity = "second-order chain, rules through summary()", tolerance = 1e-12, key = "path")
})

test_that("on first-order sequences the rules out of single states equal markovchain's transition matrix", {
  skip_without_reference("markovchain")
  random <- ref_random_sequences(25L, 20L, c("a", "b", "c", "d"), seed = 11L)
  fit <- markovchain::markovchainFit(data = random, method = "mle")$estimate
  matrix <- methods::slot(fit, "transitionMatrix")
  rules <- subset(hg_get(hon(random, max_order = 1L, min_freq = 1L), what = "rules"), select = c(from, to, probability))
  expect_reference(rules, ref_long(matrix, name = "probability"), reference = "markovchain", covers = c("hon()", "hg_get.net_hon()"),
                   quantity = "first-order transition probabilities", tolerance = 1e-12, key = c("from", "to"),
                   note = "markovchain::markovchainFit(method = 'mle') on the list of sequences.")
})
