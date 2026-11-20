# hg_bootstrap(): resampling of sequences for the higher-order rules. The estimate of each rule (its count and
# probability in the sequences themselves) is compared with the k-gram counts in base R. The intervals come from
# resamples drawn inside the function and are not compared.

test_that("the estimates of the first- and second-order rules equal the k-gram counts", {
  sequences <- ref_second_order(n_sequences = 40L, length = 20L, seed = 9L)
  found <- hg_get(hg_bootstrap(sequences, n_boot = 20L, max_order = 2L, seed = 1L), what = "edges")
  # every rule as a path: the states of `from` followed by `to`
  path <- paste(found$from, found$to, sep = " -> ")
  counts <- c(ref_kgrams(sequences, 2L), ref_kgrams(sequences, 3L))
  contexts <- vapply(strsplit(path, " -> ", fixed = TRUE), \(p) paste(utils::head(p, -1L), collapse = " -> "), character(1))
  continued <- vapply(contexts, \(context) {
    k <- length(strsplit(context, " -> ", fixed = TRUE)[[1L]]) + 1L
    grams <- ref_kgrams(sequences, k)
    starts <- vapply(strsplit(names(grams), " -> ", fixed = TRUE), \(p) paste(utils::head(p, -1L), collapse = " -> "), character(1))
    sum(grams[starts == context])
  }, numeric(1))
  expect_reference(data.frame(path = path, count = found$count, probability = found$probability),
                   data.frame(path = path, count = unname(counts[path]), probability = unname(counts[path]) / continued),
                   reference = "independent: k-gram counts of the sequences, base R",
                   covers = c("hg_bootstrap()", "hg_get.hypergraphs_bootstrap()"), quantity = "rule estimates",
                   tolerance = 1e-12, key = "path")
  summarised <- hg_get(summary(hg_bootstrap(sequences, n_boot = 20L, max_order = 2L, seed = 1L)), what = "edges")
  summarised_path <- paste(summarised$from, summarised$to, sep = " -> ")
  expect_reference(data.frame(path = summarised_path, count = summarised$count),
                   data.frame(path = summarised_path, count = unname(counts[summarised_path])),
                   reference = "independent: k-gram counts of the sequences, base R",
                   covers = c("summary.hypergraphs_bootstrap()", "hg_get.hypergraphs_summary()"), quantity = "rule counts through summary()",
                   tolerance = 0, key = "path")
})
