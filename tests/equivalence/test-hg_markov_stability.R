# hg_markov_stability(): the first-order chain fitted to the sequences (pooled maximum likelihood) and its
# stability: the stationary distribution, the mean return time 1 / pi, the mean first passage times between
# states, the persistence P[i, i] and the sojourn time 1 / (1 - P[i, i]). References: markovchain
# (markovchainFit(), steadyStates(), meanFirstPassageTime(), meanRecurrenceTime()) and the definitions in base R.

sequences <- ref_random_sequences(30L, 25L, c("A", "B", "C", "D"), seed = 13L)
stability <- hg_markov_stability(sequences)

test_that("the stationary distribution, return times and first passage times equal markovchain", {
  skip_without_reference("markovchain")
  chain <- markovchain::markovchainFit(data = sequences, method = "mle")$estimate
  stationary <- markovchain::steadyStates(chain)[1L, ]
  expect_reference(hg_get(stability, what = "stationary"),
                   data.frame(state = names(stationary), stationary = unname(stationary),
                              return_time = unname(markovchain::meanRecurrenceTime(chain)[names(stationary)])),
                   reference = "markovchain", covers = c("hg_markov_stability()", "hg_get.net_markov_stability()", "hg_get.hypergraphs_result()"),
                   quantity = "stationary distribution and return time", tolerance = 1e-10, key = "state",
                   note = "markovchainFit(method = 'mle') on the list of sequences, steadyStates(), meanRecurrenceTime().")
  passage <- markovchain::meanFirstPassageTime(chain)
  # meanFirstPassageTime() leaves the diagonal at 0; the package puts the return time there
  diag(passage) <- markovchain::meanRecurrenceTime(chain)[rownames(passage)]
  expect_reference(hg_get(stability, what = "passage_time"), ref_long(passage, name = "steps"), reference = "markovchain",
                   covers = c("hg_markov_stability()", "hg_get.net_markov_stability()"), quantity = "mean first passage time",
                   tolerance = 1e-10, key = c("from", "to"), note = "meanFirstPassageTime(), its diagonal set to meanRecurrenceTime().")
})

test_that("persistence and sojourn times equal P[i, i] and 1 / (1 - P[i, i]) of the fitted chain", {
  transitions <- table(unlist(lapply(sequences, utils::head, -1L)), unlist(lapply(sequences, utils::tail, -1L)))
  probability <- unclass(prop.table(transitions, 1L))
  persistence <- diag(probability)
  # The states table rounds persistence and the stationary share to 4 decimals and the times to 2 (the package's
  # display precision): the tolerance is half a unit of the last decimal shown.
  expect_reference(subset(hg_get(stability, what = "states"), select = c(state, persistence)),
                   data.frame(state = names(persistence), persistence = unname(persistence)),
                   reference = "independent: diagonal of the maximum-likelihood transition matrix, base R",
                   covers = c("hg_markov_stability()", "hg_get.net_markov_stability()"), quantity = "persistence (4 decimals)",
                   tolerance = 5e-5, key = "state")
  expect_reference(subset(hg_get(stability, what = "states"), select = c(state, sojourn_time)),
                   data.frame(state = names(persistence), sojourn_time = unname(1 / (1 - persistence))),
                   reference = "independent: 1 / (1 - P[i, i]) of the maximum-likelihood transition matrix, base R",
                   covers = c("hg_markov_stability()", "hg_get.net_markov_stability()"), quantity = "sojourn time (2 decimals)",
                   tolerance = 5e-3, key = "state")
  expect_reference(subset(hg_get(summary(stability), what = "states"), select = c(state, persistence)),
                   data.frame(state = names(persistence), persistence = unname(persistence)),
                   reference = "independent: diagonal of the maximum-likelihood transition matrix, base R",
                   covers = c("summary.hypergraphs_result()", "hg_get.hypergraphs_summary()"),
                   quantity = "persistence through summary() (4 decimals)", tolerance = 5e-5, key = "state")
})
