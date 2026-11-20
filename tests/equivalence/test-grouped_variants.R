# The grouped variants of the memory-network verbs: hon(group =), hg_bootstrap() of a group model, markov_order()
# and hg_markov_stability() of a grouped network (a Nestimate netobject_group). Each grouped table stacks the
# ungrouped table of every group under a `group` column, so the reference is the ungrouped comparison of the
# matching test file (test-hg_bootstrap.R, test-markov_order.R, test-hg_markov_stability.R) repeated on each
# group's own sequences: the k-gram counts in base R, the G2 of the conditional independence tables, and
# markovchain's stationary distribution.

group_one <- ref_second_order(n_sequences = 30L, length = 15L, seed = 9L)
group_two <- ref_random_sequences(30L, 15L, c("A", "B", "C"), seed = 21L)
sequences <- c(group_one, group_two)
labels <- rep(c("g1", "g2"), each = 30L)
by_group <- list(g1 = group_one, g2 = group_two)

# The count and probability of each rule `path` (the context followed by the next state) in some sequences.
gv_rule_estimates <- function(sequences, path) {
  states <- strsplit(path, " -> ", fixed = TRUE)
  counts <- c(ref_kgrams(sequences, 2L), ref_kgrams(sequences, 3L))
  continued <- vapply(states, \(p) {
    context <- paste(utils::head(p, -1L), collapse = " -> ")
    grams <- ref_kgrams(sequences, length(p))
    starts <- vapply(strsplit(names(grams), " -> ", fixed = TRUE), \(g) paste(utils::head(g, -1L), collapse = " -> "),
                     character(1))
    sum(grams[starts == context])
  }, numeric(1))
  data.frame(path = path, count = unname(counts[path]), probability = unname(counts[path]) / continued)
}

# A grouped table's rules as group, path, count, probability, and the same from the k-grams of each group.
gv_compare_rules <- function(found, covers, quantity) {
  path <- paste(found$from, found$to, sep = " -> ")
  expected <- do.call(rbind, lapply(names(by_group), \(g) {
    cbind(group = g, gv_rule_estimates(by_group[[g]], path[found$group == g]))
  }))
  expect_reference(data.frame(group = found$group, path = path, count = found$count, probability = found$probability),
                   expected, reference = "independent: k-gram counts of each group's sequences, base R", covers = covers,
                   quantity = quantity, tolerance = 1e-12, key = c("group", "path"))
}

test_that("a group model's rules equal the k-gram counts of each group", {
  model <- hon(sequences, group = labels, max_order = 2L)
  gv_compare_rules(hg_get(model), c("hon()", "hg_get.hypergraphs_memory_group()"), "rules of a group model")
  gv_compare_rules(hg_get(summary(model)), c("summary.hypergraphs_memory_group()", "hg_get.hypergraphs_summary()"),
                   "rules of a group model through summary()")
})

test_that("the bootstrap of a group model estimates each group's rules by its k-gram counts", {
  boot <- hg_bootstrap(hon(sequences, group = labels, max_order = 2L), n_boot = 20L, seed = 1L)
  gv_compare_rules(hg_get(boot), c("hg_bootstrap()", "hg_get.hypergraphs_bootstrap_group()"), "bootstrap estimates by group")
  gv_compare_rules(hg_get(summary(boot)), c("summary.hypergraphs_bootstrap_group()", "hg_get.hypergraphs_summary()"),
                   "bootstrap estimates by group through summary()")
})

# G2 and df of order k against order k - 1 (Anderson and Goodman 1957), as in test-markov_order.R.
gv_order_test <- function(sequences, k) {
  tuples <- do.call(rbind, lapply(Filter(\(s) length(s) >= k + 1L, sequences), \(s) {
    at <- seq_len(length(s) - k)
    data.frame(back = s[at],
               between = vapply(at, \(p) if (k == 1L) "" else paste(s[(p + 1L):(p + k - 1L)], collapse = "\r"), character(1)),
               following = s[at + k])
  }))
  parts <- lapply(split(tuples, tuples$between), \(part) {
    observed <- table(part$back, part$following)
    if (nrow(observed) < 2L || ncol(observed) < 2L) return(c(g2 = 0, df = 0))
    expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
    kept <- observed > 0
    c(g2 = 2 * sum(observed[kept] * log(observed[kept] / expected[kept])), df = (nrow(observed) - 1) * (ncol(observed) - 1))
  })
  total <- Reduce(`+`, parts)
  data.frame(order = k, g2 = total[["g2"]], df = total[["df"]])
}

grouped_network <- function() {
  skip_without_reference("Nestimate")
  wide <- as.data.frame(do.call(rbind, sequences))
  wide$grp <- labels
  Nestimate::build_network(wide, method = "relative", group = "grp")
}

test_that("the Markov order test of a grouped network gives each group's G2 and df", {
  fit <- markov_order(grouped_network(), max_order = 2L, n_perm = 19L, seed = 1L)
  expected <- do.call(rbind, lapply(names(by_group), \(g) {
    cbind(group = g, do.call(rbind, lapply(1:2, \(k) gv_order_test(by_group[[g]], k))))
  }))
  expect_reference(subset(hg_get(fit, what = "orders"), order >= 1, select = c(group, order, g2, df)), expected,
                   reference = "independent: G2 of the conditional independence tables of each group, base R",
                   covers = c("markov_order()", "hg_get.net_markov_order_group()"), quantity = "G2 and df by group",
                   tolerance = 1e-9, key = c("group", "order"), note = "the network is built with Nestimate::build_network(group =)")
})

test_that("the stability of a grouped network gives markovchain's stationary distribution of each group", {
  skip_without_reference("markovchain")
  stability <- hg_markov_stability(grouped_network())
  expected <- do.call(rbind, lapply(names(by_group), \(g) {
    chain <- markovchain::markovchainFit(data = by_group[[g]], method = "mle")$estimate
    stationary <- markovchain::steadyStates(chain)[1L, ]
    data.frame(group = g, state = names(stationary), stationary = unname(stationary),
               return_time = unname(markovchain::meanRecurrenceTime(chain)[names(stationary)]))
  }))
  expect_reference(hg_get(stability, what = "stationary"), expected, reference = "markovchain",
                   covers = c("hg_markov_stability()", "hg_get.net_markov_stability_group()"),
                   quantity = "stationary distribution and return time by group", tolerance = 1e-10, key = c("group", "state"),
                   note = "markovchainFit(method = 'mle') on each group's sequences, steadyStates(), meanRecurrenceTime().")
})
