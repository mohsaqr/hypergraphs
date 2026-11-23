# markov_order(): the test of order k against order k - 1 (Anderson and Goodman 1957). For each context of the
# k - 1 states between, the table of the state k steps back against the next state is tested for independence;
# G2 = 2 sum O log(O / E) over the contexts and df = sum (rows - 1)(columns - 1) over the contexts whose table has
# at least two rows and two columns; the asymptotic p-value is the chi-square tail. Reference: the formula in base
# R with stats::pchisq(). markovchain::verifyMarkovProperty() is not a reference here: it computes a Pearson
# statistic of a first-order fit on triples, another test.

#' G2 and df of order k against order k - 1, from the definition
by_definition <- function(sequences, k) {
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
  data.frame(order = k, g2 = total[["g2"]], df = total[["df"]], p_asymptotic = stats::pchisq(total[["g2"]], total[["df"]], lower.tail = FALSE))
}

test_that("G2, df and asymptotic p-values of each order equal the conditional independence test", {
  inputs <- list(second_order = ref_second_order(n_sequences = 40L, length = 20L, seed = 9L),
                 first_order = ref_random_sequences(30L, 40L, c("a", "b", "c", "d"), seed = 21L))
  invisible(Map(\(sequences, name) {
    fit <- markov_order(sequences, max_order = 3L, n_perm = 19L, seed = 1L)
    expected <- do.call(rbind, lapply(1:3, \(k) by_definition(sequences, k)))
    # Sums of a few hundred logarithms; p-values are compared relative, as they run down to 1e-30.
    expect_reference(subset(hg_get(fit, what = "orders"), order >= 1, select = c(order, g2, df)), subset(expected, select = c(order, g2, df)),
                     reference = "independent: G2 of the conditional independence tables, base R",
                     covers = c("markov_order()", "hg_get.net_markov_order()", "hg_get.hypergraphs_result()"),
                     quantity = paste(name, "G2 and df"), tolerance = 1e-9, key = "order")
    expect_reference(subset(hg_get(fit, what = "orders"), order >= 1, select = c(order, p_asymptotic)),
                     subset(expected, select = c(order, p_asymptotic)),
                     reference = "independent: chi-squared tail, stats::pchisq() of the G2 written out above",
                     covers = c("markov_order()", "hg_get.net_markov_order()"), quantity = paste(name, "asymptotic p-value"),
                     tolerance = 1e-8, measure = "relative", key = "order",
                     note = "stats::pchisq(G2, df, lower.tail = FALSE) of the G2 written out above.")
  }, inputs, names(inputs)))
})
