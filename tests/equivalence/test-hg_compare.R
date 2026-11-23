# hg_compare(): the permutation comparison of the higher-order rules of two groups. The rule set comes from the
# pooled sequences; for each rule the two groups' counts and probabilities, their difference and the global statistic
# (the pooled-count-weighted mean absolute difference over the rules whose context both groups observe) are no draws
# and are written here from the documented definition with the k-gram counts of each group's sequences. The
# Benjamini-Hochberg adjustment is compared with stats::p.adjust(). The permutation p-values are drawn inside the
# function with only a seed: they are checked by the bounds of any permutation p-value (count + 1) / (B + 1)
# (Phipson & Smyth 2010).

first <- ref_second_order(n_sequences = 30L, length = 15L, seed = 9L)
second <- ref_random_sequences(30L, 15L, c("A", "B", "C"), seed = 21L)
labels <- rep(c("g1", "g2"), each = 30L)
n_perm <- 99L

# Count of a rule (context, then next state) and of its context continued by anything, in some sequences.
hc_counts <- function(sequences, path) {
  vapply(strsplit(path, " -> ", fixed = TRUE), \(p) {
    grams <- ref_kgrams(sequences, length(p))
    starts <- vapply(strsplit(names(grams), " -> ", fixed = TRUE), \(g) paste(utils::head(g, -1L), collapse = " -> "),
                     character(1))
    full <- paste(p, collapse = " -> ")
    c(rule = if (full %in% names(grams)) grams[[full]] else 0,
      context = sum(grams[starts == paste(utils::head(p, -1L), collapse = " -> ")]))
  }, numeric(2))
}

comparison <- hg_compare(hon(c(first, second), group = labels, max_order = 2L), n_perm = n_perm, seed = 1L)

test_that("each group's rule counts, probabilities and their difference equal the k-gram counts", {
  rules <- hg_get(comparison)
  path <- paste(rules$from, rules$to, sep = " -> ")
  one <- hc_counts(first, path)
  two <- hc_counts(second, path)
  probability <- \(counts) ifelse(counts["context", ] > 0, counts["rule", ] / counts["context", ], NA_real_)
  expected <- data.frame(path = path, count_g1 = one["rule", ], count_g2 = two["rule", ],
                         probability_g1 = probability(one), probability_g2 = probability(two),
                         diff = probability(one) - probability(two))
  expect_reference(cbind(path = path, subset(rules, select = c(count_g1, count_g2, probability_g1, probability_g2, diff))),
                   expected, reference = "independent: k-gram counts of each group's sequences, base R",
                   covers = c("hg_compare()", "hg_compare.hypergraphs_memory_group()", "hg_get.hypergraphs_comparison()"), quantity = "rule estimates by group",
                   tolerance = 1e-12, key = "path")
  summarised <- hg_get(summary(comparison))
  expect_reference(data.frame(path = paste(summarised$from, summarised$to, sep = " -> "), diff = summarised$diff),
                   subset(expected, select = c(path, diff)),
                   reference = "independent: k-gram counts of each group's sequences, base R",
                   covers = c("summary.hypergraphs_comparison()", "hg_get.hypergraphs_summary()"),
                   quantity = "rule differences through summary()", tolerance = 1e-12, key = "path")
  # The global statistic: pooled-count-weighted mean |diff| over the rules both groups observe.
  pooled <- hc_counts(c(first, second), path)["rule", ]
  comparable <- !is.na(expected$diff)
  expect_reference(cbind(test = "global", subset(hg_get(comparison, what = "global"), select = c(statistic, n_comparable, n_rules))),
                   data.frame(test = "global", statistic = sum(pooled[comparable] * abs(expected$diff[comparable])) / sum(pooled[comparable]),
                              n_comparable = sum(comparable), n_rules = length(path)),
                   reference = "independent: k-gram counts of each group's sequences, base R", covers = c("hg_compare()", "hg_compare.hypergraphs_memory_group()"),
                   quantity = "global statistic", tolerance = 1e-12, key = "test")
})

# Both checks below take the package's own p-values as input: they are range and adjustment checks, declared as
# independent computations, not outside numbers (the permutation p-values themselves have no outside reference).
test_that("adjusted p-values equal stats::p.adjust(BH); p-values lie where any permutation p-value lies", {
  rules <- hg_get(comparison)
  path <- paste(rules$from, rules$to, sep = " -> ")
  expect_reference(data.frame(path = path, p_adj = rules$p_adj),
                   data.frame(path = path, p_adj = stats::p.adjust(rules$p_value, method = "BH")),
                   reference = "independent: stats::p.adjust() of the package's own p-values",
                   covers = c("hg_compare()", "hg_compare.hypergraphs_memory_group()", "hg_get.hypergraphs_comparison()"), quantity = "Benjamini-Hochberg adjustment",
                   tolerance = 1e-12, key = "path")
  p_values <- data.frame(path = path, p_value = rules$p_value)
  expect_reference(p_values, within(p_values, p_value <- pmin(pmax(p_value, 1 / (rules$n_perm_used + 1)), 1)),
                   reference = "independent: range of a permutation p-value, [1/(B+1), 1] (Phipson & Smyth 2010)",
                   covers = c("hg_compare()", "hg_compare.hypergraphs_memory_group()", "hg_get.hypergraphs_comparison()"),
                   quantity = "Monte Carlo summary: rule p-values within their bounds", tolerance = 0, key = "path")
  global <- subset(hg_get(comparison, what = "global"), select = c(p_value, n_perm_used))
  expect_reference(global$p_value, min(max(global$p_value, 1 / (global$n_perm_used + 1)), 1),
                   reference = "independent: range of a permutation p-value, [1/(B+1), 1] (Phipson & Smyth 2010)",
                   covers = c("hg_compare()", "hg_compare.hypergraphs_memory_group()"), quantity = "Monte Carlo summary: global p-value within its bounds", tolerance = 0)
})
