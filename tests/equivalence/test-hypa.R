# hypa(): the hypergeometric path anomalies of LaRock et al. (2020). For each path of k + 1 states the count of
# its last transition is compared with a null in which the m observed transitions of the k-th order network are
# drawn without replacement from an urn of sum(Xi) balls, Xi[v, w] of them for the transition v -> w, where Xi is
# fitted so that Xi / m reproduces the observed strengths (expected count = Xi / m). The p-value of a path is
# P(X <= count) with X ~ Wallenius(m1 = Xi, m2 = sum(Xi) - Xi, n = m, odds 1), the upper p-value P(X > count). References: BiasedUrn's
# pWNCHypergeo() and stats::phyper() (the two coincide at odds 1), the counts by k-gram in base R.
#
# No outside software computes HYPA here: pathpy 2.2.0 (in ~/equiv-py) has no HYPA module, and the authors' code
# (tlarock/hypa) is not installed. So Xi is fitted again in this file from the definition, on the de Bruijn support
# built from the sequences (Sinkhorn & Knopp scaling of s_out(v) s_in(w)), and the hypergeometric tails are taken
# of that Xi. Every row is an independent computation (Indirect in Prova), none is an outside number. Tails taken of
# the package's own Xi are labelled as such.

# The HYPA null of order 2 from the definition (LaRock et al. 2020): the contexts v = (a, b) with an outgoing path
# and the targets w = (b, c) with an incoming one, joined on b, are the possible paths a -> b -> c; Xi[v, w] =
# x_v y_w on them, scaled (Sinkhorn & Knopp 1967) until the expected counts m Xi / sum(Xi) reproduce the observed
# out- and in-strengths, then Xi is put on sum(Xi) = m^2. With `round_sweeps = TRUE` Xi is kept on that scale and
# rounded to integers after every row and every column step: the rounding that Nestimate's documentation
# (R/hypa.R, .hypa_fit_xi) says the original Python fitXi()/xifix_row() does. That variant is written from that
# description, not from the authors' code.
hypa_reference <- function(sequences, round_sweeps = FALSE, sweeps = 5000L) {
  grams <- ref_kgrams(sequences, 3L)
  parts <- strsplit(names(grams), " -> ", fixed = TRUE)
  observed <- data.frame(v = vapply(parts, \(p) paste(p[1:2], collapse = " -> "), character(1)),
                         w = vapply(parts, \(p) paste(p[2:3], collapse = " -> "), character(1)), count = unname(grams))
  s_out <- tapply(observed$count, observed$v, sum)
  s_in <- tapply(observed$count, observed$w, sum)
  possible <- merge(data.frame(v = names(s_out), b = sub("^.* -> ", "", names(s_out))),
                    data.frame(w = names(s_in), b = sub(" -> .*$", "", names(s_in))), by = "b")
  possible <- merge(possible, observed, by = c("v", "w"), all.x = TRUE)
  possible$count[is.na(possible$count)] <- 0
  m <- sum(observed$count)
  rows <- match(possible$v, names(s_out))
  cols <- match(possible$w, names(s_in))
  target_out <- as.numeric(s_out)[rows]
  target_in <- as.numeric(s_in)[cols]
  scale_to <- \(xi, group, target, total) {
    current <- stats::ave(xi, group, FUN = sum)
    fitted <- xi * target * total / current
    if (round_sweeps) round(fitted) else fitted
  }
  start <- target_out * target_in * if (round_sweeps) 1 else 1 / m^2
  total <- if (round_sweeps) m else 1 / m
  xi <- Reduce(\(xi, i) scale_to(scale_to(xi, rows, target_out, total), cols, target_in, total), seq_len(sweeps), start)
  xi <- xi * m^2 / sum(xi)
  big_k <- round(xi)
  data.frame(path = paste(possible$v, sub("^.* -> ", "", possible$w), sep = " -> "), count = possible$count,
             expected = m * xi / sum(xi),
             p_under = stats::phyper(possible$count, big_k, round(m^2) - big_k, m),
             p_over = stats::phyper(possible$count, big_k, round(m^2) - big_k, m, lower.tail = FALSE))
}

sequences <- ref_second_order(n_sequences = 40L, length = 20L, seed = 9L)

test_that("path counts equal the k-grams of the sequences, and expected counts reproduce the observed strengths", {
  scores <- hg_get(hypa(sequences, order = 2L, p_adjust = "none"), what = "scores")
  grams <- ref_kgrams(sequences, 3L)
  expect_reference(subset(scores, count > 0, select = c(path, count)), data.frame(path = names(grams), count = unname(grams)),
                   reference = "independent: k-gram counts of the sequences, base R",
                   covers = c("hypa()", "hypa.default()", "hg_get.net_hypa()", "hg_get.hypergraphs_hypa()"),
                   quantity = "observed path counts", tolerance = 0, key = "path")
  # the expected counts of the transitions out of a context sum to the observed transitions out of it
  observed_out <- stats::aggregate(list(count = scores$count), list(from = scores$from), sum)
  expected_out <- stats::aggregate(list(count = scores$expected), list(from = scores$from), sum)
  expect_reference(expected_out, observed_out, reference = "independent: observed out-strengths of the second-order network",
                   covers = c("hypa()", "hg_get.net_hypa()"), quantity = "expected out-strengths", tolerance = 1e-8,
                   key = "from", note = "Xi is fitted iteratively (LaRock et al. 2020); its strength constraints hold to the fit's tolerance.")
})

test_that("p-values equal the Wallenius distribution at odds 1 (BiasedUrn) and the hypergeometric (stats)", {
  skip_without_reference("BiasedUrn")
  scores <- hg_get(hypa(sequences, order = 2L, p_adjust = "none"), what = "scores")
  m <- sum(scores$count)
  xi <- scores$expected * m
  m1 <- round(xi)
  m2 <- round(m^2) - m1
  wallenius <- mapply(\(count, a, b) BiasedUrn::pWNCHypergeo(count, a, b, min(m, a + b), odds = 1, precision = 1e-12),
                      scores$count, m1, m2)
  # Wallenius probabilities to BiasedUrn's precision of 1e-12, summed over the lower tail.
  expect_known_difference(subset(scores, select = c(path, p_under)), data.frame(path = scores$path, p_under = wallenius),
                          reference = "BiasedUrn", covers = c("hypa()", "hypa.default()", "hg_get.net_hypa()"),
                          quantity = "lower-tail p-values", tolerance = 1e-9, key = "path",
                          # BiasedUrn computes the tails of the package's own fitted Xi: a check of the tail, not of the
                          # null, so its agreeing rows are no outside evidence for hypa().
                          authority = "independent",
                          note = "BiasedUrn::pWNCHypergeo(count, round(Xi), round(sum Xi) - round(Xi), min(m, sum Xi), odds = 1, precision = 1e-12).",
                          reason = "At odds 1 the Wallenius distribution is the hypergeometric; BiasedUrn's pWNCHypergeo() (precision = 1e-12 requested) differs from the exact stats::phyper() by up to 1.6e-9 on these inputs, and hypergraphs equals stats::phyper() (compared below).")
  expect_reference(subset(scores, select = c(path, p_under)),
                   data.frame(path = scores$path, p_under = stats::phyper(scores$count, m1, m2, pmin(m, m1 + m2))),
                   reference = "independent: stats::phyper() of the package's own fitted Xi",
                   covers = c("hypa()", "hypa.default()", "hg_get.net_hypa()", "hg_get.hypergraphs_hypa()"),
                   quantity = "lower-tail p-values", tolerance = 1e-9, key = "path",
                   note = "The urn is the package's fitted Xi: this checks the tail, not the null. The null is checked against an independent fit below.")
  expect_reference(subset(scores, select = c(path, p_over)),
                   data.frame(path = scores$path, p_over = stats::phyper(scores$count, m1, m2, pmin(m, m1 + m2), lower.tail = FALSE)),
                   reference = "independent: stats::phyper() of the package's own fitted Xi", covers = c("hypa()", "hypa.default()", "hg_get.net_hypa()"),
                   quantity = "upper-tail p-values", tolerance = 1e-9, key = "path",
                   note = "The strict upper tail P(X > count), as hypergraphs documents it (a path is over-represented when P(X > f) <= alpha).")
})

test_that("adjusted p-values equal stats::p.adjust()", {
  scores <- hg_get(hypa(sequences, order = 2L, p_adjust = "BH"), what = "scores")
  expect_reference(subset(scores, select = c(path, p_adj_under, p_adj_over)),
                   data.frame(path = scores$path, p_adj_under = stats::p.adjust(scores$p_under, method = "BH"),
                              p_adj_over = stats::p.adjust(scores$p_over, method = "BH")),
                   reference = "independent: stats::p.adjust() of the package's own p-values", covers = c("hypa()", "hg_get.net_hypa()"),
                   quantity = "Benjamini-Hochberg adjustment", tolerance = 1e-12, key = "path")
  expect_reference(subset(hg_get(summary(hypa(sequences, order = 2L, p_adjust = "BH")), what = "scores"), select = c(path, p_adj_under)),
                   data.frame(path = scores$path, p_adj_under = stats::p.adjust(scores$p_under, method = "BH")),
                   reference = "independent: stats::p.adjust() of the package's own p-values",
                   covers = c("summary.hypergraphs_hypa()", "hg_get.hypergraphs_summary()"),
                   quantity = "Benjamini-Hochberg adjustment through summary()", tolerance = 1e-12, key = "path")
})

test_that("expected counts and tails equal an independent fit of the HYPA null on the de Bruijn support", {
  scores <- hg_get(hypa(sequences, order = 2L, p_adjust = "none"), what = "scores")
  continuous <- hypa_reference(sequences)
  covers <- c("hypa()", "hypa.default()", "hg_get.net_hypa()", "hg_get.hypergraphs_hypa()")
  # Possible paths and their counts: a key in one table and not the other fails.
  expect_reference(subset(scores, select = c(path, count)), subset(continuous, select = c(path, count)),
                   reference = "independent: possible paths of the de Bruijn graph and their counts, base R", covers = covers,
                   quantity = "possible paths and counts", tolerance = 0, key = "path")
  # The package stops its scaling when the largest strength error is below 1e-10 m (m = 720 paths here, so
  # 7.2e-8 on a strength); the reference scales 5000 times. A cell may carry the strength error: 1e-6.
  expect_reference(subset(scores, select = c(path, expected)), subset(continuous, select = c(path, expected)),
                   reference = "independent: Sinkhorn-Knopp fit of Xi on the de Bruijn support, base R", covers = covers,
                   quantity = "expected counts", tolerance = 1e-6, key = "path")
  expect_reference(subset(scores, select = c(path, p_under, p_over)), subset(continuous, select = c(path, p_under, p_over)),
                   reference = "independent: Sinkhorn-Knopp fit of Xi on the de Bruijn support, base R", covers = covers,
                   quantity = "hypergeometric tails of the independent Xi", tolerance = 1e-9, key = "path",
                   note = "stats::phyper(count, round(Xi), m^2 - round(Xi), m) of the Xi fitted in this file.")
  rounded <- hypa_reference(sequences, round_sweeps = TRUE, sweeps = 200L)
  expect_known_difference(subset(scores, select = c(path, expected)), subset(rounded, select = c(path, expected)),
                          reference = "independent: Xi rounded to integers in every sweep, as described for tlarock/hypa", covers = covers,
                          quantity = "expected counts with Xi rounded in every sweep", tolerance = 1e-6, key = "path",
                          agrees_with = "the continuous Sinkhorn-Knopp fit (compared above)",
                          reason = "Nestimate (which hypa() calls) fits Xi without rounding inside the loop and documents this as a departure from the original Python fitXi()/xifix_row() of LaRock et al., which round Xi to integers in every sweep. The rounded variant here is written from that description, not from the authors' code, which is not installed.")
})

# hypa.net_hg(): the hypergraph method. Written from its documented definition (?hypa.net_hg): c_ij the hyperedges
# both nodes share (from the incidence, base R), m = sum c_ij, strengths s_i; the Chung-Lu fitness x with
# x_i sum_{j != i} x_j = s_i solved here by Newton's method (the package uses a damped fixed point to 1e-10 m);
# E_ij = x_i x_j; p_under = phyper(c_ij, round(m E_ij), m^2 - round(m E_ij), m) and p_over its strict upper tail;
# the scored pairs those with c_ij >= min_count or E_ij >= min_count.
hn_reference <- function(memberships, min_count) {
  incidence <- ref_incidence(memberships, node = "member", edge = "edge")
  shared <- tcrossprod(incidence)
  diag(shared) <- 0
  strength <- rowSums(shared)
  active <- strength > 0
  shared <- shared[active, active]
  strength <- strength[active]
  m <- sum(shared) / 2
  newton <- function(x, step) {
    if (step == 0L) return(x)
    residual <- x * (sum(x) - x) - strength
    # d f_i / d x_j = x_i (j != i), d f_i / d x_i = sum_{j != i} x_j
    jacobian <- outer(x, rep(1, length(x)))
    diag(jacobian) <- sum(x) - x
    newton(x - solve(jacobian, residual), step - 1L)
  }
  fitness <- newton(strength / sqrt(2 * m), 60L)
  expected <- outer(fitness, fitness)
  keep <- \(cells) upper.tri(cells) & (shared >= min_count | expected >= min_count)
  pairs <- ref_pairs_upper(shared, keep, name = "observed")
  e <- expected[cbind(pairs$from, pairs$to)]
  big_k <- round(m * e)
  data.frame(pair = paste(pmin(pairs$from, pairs$to), pmax(pairs$from, pairs$to)), observed = pairs$observed,
             expected = e, p_under = stats::phyper(pairs$observed, big_k, m^2 - big_k, m),
             p_over = stats::phyper(pairs$observed, big_k, m^2 - big_k, m, lower.tail = FALSE))
}

test_that("hypa() of a hypergraph equals the documented hypergeometric null of its co-occurrences", {
  set.seed(31)
  nodes <- sprintf("v%02d", seq_len(14L))
  sets <- c(lapply(seq_len(40L), \(i) sample(nodes, sample(2:4, 1L), prob = rep(c(3, 1), c(4L, 10L)))),
            rep(list(c("v05", "v06")), 6L))
  memberships <- data.frame(member = unlist(sets), edge = sprintf("e%02d", rep(seq_along(sets), lengths(sets))))
  found <- hypa(group_hypergraph(memberships, node = "member", hyperedge = "edge"), min_count = 3L)
  expected <- hn_reference(memberships, min_count = 3L)
  found <- data.frame(pair = paste(pmin(found$from, found$to), pmax(found$from, found$to)),
                      subset(found, select = c(observed, expected, p_under, p_over)))
  expect_reference(subset(found, select = c(pair, observed)), subset(expected, select = c(pair, observed)),
                   reference = "independent: shared hyperedges of the incidence, base R", covers = "hypa.net_hg()",
                   quantity = "scored pairs and their co-occurrence counts", tolerance = 0, key = "pair")
  # The package's fit stops at a strength error of 1e-10 m; Newton's method here converges to rounding.
  expect_reference(subset(found, select = c(pair, expected)), subset(expected, select = c(pair, expected)),
                   reference = "independent: Chung-Lu fitness solved by Newton's method, base R", covers = "hypa.net_hg()",
                   quantity = "expected co-occurrences", tolerance = 1e-8, key = "pair")
  expect_reference(subset(found, select = c(pair, p_under, p_over)), subset(expected, select = c(pair, p_under, p_over)),
                   reference = "independent: stats::phyper() of the Chung-Lu null solved above", covers = "hypa.net_hg()",
                   quantity = "hypergeometric tails", tolerance = 1e-9,
                   key = "pair", note = "stats::phyper() with K = round(m E_ij), N = m^2, n = m")
})
