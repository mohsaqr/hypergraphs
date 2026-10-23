# Hypergeometric anomaly detection for hypergraph co-occurrence.
#
# The memory family's sequence HYPA (hypa() on sequences) scores the paths of
# a De Bruijn graph against a hypergeometric null whose propensity matrix is
# fitted so that the expected counts preserve every node's observed strength
# (LaRock et al. 2020). This file is the undirected analogue for a
# hypergraph: a pair of nodes co-occurs in some number of hyperedges, and the
# question is whether that count departs from what the two nodes' strengths
# imply.
#
# The null is a fitness model, E_ij = x_i x_j over every unordered pair of
# nodes that take part in a co-occurrence, with x fitted so that the expected
# strengths equal the observed ones. The fit needs only the strength vector,
# so each sweep is O(n) and no n x n matrix is ever built. The scored pairs
# are those that could be classified at all: co-occurring at least
# `min_count` times, or expected to co-occur at least `min_count` times. The
# second set is found from the sorted fitness values and holds at most
# m / min_count pairs, because the expected counts sum to m.

#' Hypergeometric anomaly detection for co-occurring node pairs
#'
#' The hypergraph method of [hypa()]. Scores pairs of nodes against a
#' hypergeometric null that keeps every node's co-occurrence strength and
#' otherwise pairs the nodes at random. A pair that shares far more
#' hyperedges than the null predicts is over-represented; a pair that shares
#' far fewer, possibly none, is under-represented. This is the undirected
#' counterpart of [hypa.default()], which applies the same null to the paths
#' of a higher-order network, and it reports the same columns so the two
#' read alike.
#'
#' The test is analytic, so it needs no resampling. That is what makes it
#' usable where [hg_null_test()] is not: [hg_null_test()] answers whether the
#' hypergraph as a whole carries more repetition than chance, by Monte Carlo,
#' and returns one row; this returns one row per pair.
#'
#' @details
#' **Observed counts.** For nodes \eqn{i < j}, \eqn{c_{ij}} is the number of
#' hyperedges that contain both (hyperedge weights are ignored; a repeated
#' hyperedge counts each time). The total is \eqn{m = \sum_{i<j} c_{ij}} and
#' the strength of node \eqn{i} is \eqn{s_i = \sum_{j \ne i} c_{ij}}, so
#' \eqn{\sum_i s_i = 2m}.
#'
#' **Null model.** Every unordered pair of distinct nodes with positive
#' strength is a possible cell, whether it co-occurs or not. The expected
#' count is \eqn{E_{ij} = x_i x_j}, the product form of Chung and Lu (2002),
#' with the fitness \eqn{x} fitted so that \eqn{\sum_{j \ne i} E_{ij} = s_i}
#' for every node: the null preserves each node's observed strength, as the
#' propensity fit of LaRock et al. (2020) preserves the in- and out-strengths
#' of a De Bruijn graph. The fit is the damped fixed point
#' \eqn{x_i \leftarrow \sqrt{x_i s_i / (X - x_i)}}, \eqn{X = \sum_k x_k},
#' iterated until the largest strength error falls below \eqn{10^{-10} m} or
#' 10,000 sweeps have run. Each sweep costs \eqn{O(n)}; no
#' \eqn{n \times n} matrix is formed. With three or more nodes the fit has a
#' solution with positive fitness exactly when no node takes part in every
#' co-occurring pair. In a star-like hypergraph, where one node is in every
#' co-occurring pair, the fit cannot converge, and the result is returned
#' with a warning of class `hypergraphs_no_converge`.
#'
#' **Test.** The propensity of a pair is \eqn{\Xi_{ij} = m E_{ij}}, which
#' sums to \eqn{m^2}. The observed count \eqn{f = c_{ij}} is referred to a
#' \eqn{\mathrm{Hypergeometric}(N, K, n)} law with \eqn{N = m^2},
#' \eqn{K = \mathrm{round}(\Xi_{ij})} and \eqn{n = m} draws, whose mean is
#' \eqn{E_{ij}}. `p_under` is \eqn{P(X \le f)} and `p_over` is the strict
#' upper tail \eqn{P(X > f) = 1 - } `p_under`, so a pair is over-represented
#' when its HYPA score \eqn{P(X \le f)} exceeds \eqn{1 - \alpha}, the rule of
#' LaRock et al. (2020). With small counts the strict upper tail rejects more
#' often than \eqn{\alpha} under the null; `min_count` limits the problem.
#'
#' **Scored pairs and classification.** A pair is scored when it could be
#' classified: when it co-occurs at least `min_count` times (eligible for
#' over-representation) or is expected to co-occur at least `min_count`
#' times (eligible for under-representation, including pairs that never
#' co-occur). There are at most \eqn{m /} `min_count` pairs of the second
#' kind, because the expected counts sum to \eqn{m}. A scored pair is
#' `"under"` when `expected >= min_count` and `p_adjusted_under < alpha`,
#' otherwise `"over"` when `observed >= min_count` and
#' `p_adjusted_over < alpha`, otherwise `"normal"`. `min_count` selects the
#' hypotheses by count, never by p-value, and the null is fitted on every
#' pair whatever the floor, so raising the floor changes which pairs are
#' tested but never the expected counts.
#'
#' **Multiple testing.** LaRock et al. (2020) classify against `alpha`
#' without correction, and so does the default `p_adjust = "none"`. Any
#' [stats::p.adjust] method can be applied instead. The under- and
#' over-representation p-values are then adjusted separately, and each
#' family is every possible pair (\eqn{n(n-1)/2} for the \eqn{n} nodes with
#' positive strength), as in the sequence method, which adjusts over every
#' possible path. Only the scored pairs are computed; the others enter
#' [stats::p.adjust()] through its `n` argument, which treats them as larger
#' p-values. Adjusting within the scored pairs alone would under-correct,
#' because a pair enters the scored set partly through its observed count.
#'
#' Time is handled by scoring a cumulative snapshot rather than by slicing
#' into periods. Disjoint periods split a pair's evidence across cells while
#' multiplying the number of tests, so both terms move the wrong way; a
#' sequence of `hg_snapshot(hg, at = ..., mode = "cumulative")`
#' accumulates evidence instead and shows when a pair becomes anomalous.
#'
#' @param x A `net_hg`, or a snapshot of a temporal one.
#' @param min_count Single whole number \eqn{\ge 1} (default `5L`). A pair
#'   is over-represented only if it co-occurs at least `min_count` times and
#'   under-represented only if it is expected to co-occur at least
#'   `min_count` times. Pairs eligible for neither are not scored.
#' @param alpha Significance threshold for the `anomaly` column, a single
#'   number in `(0, 0.5)` (default `0.05`), applied to the adjusted
#'   p-values.
#' @param p_adjust Multiplicity correction, any [stats::p.adjust] method or
#'   `"none"` (default `"none"`, as in LaRock et al. 2020). `"BH"` controls
#'   the false discovery rate within each direction.
#' @param top Optional number of rows to keep, the most significantly
#'   over-represented first. `NULL` (default) keeps every scored pair.
#' @return A base `data.frame`, one row per scored pair, most significantly
#'   over-represented first (by `p_adjusted_over`, ties broken by `ratio`,
#'   then by name), with columns `from`, `to` (the pair, in the node order of
#'   the hypergraph), `observed` (co-occurrence count, `0` for a pair that
#'   never co-occurs), `expected` (the fitted \eqn{E_{ij}}), `ratio`
#'   (observed / expected), `p_under`, `p_over`, `p_adjusted_under`,
#'   `p_adjusted_over` and `anomaly` (`"over"`, `"under"` or `"normal"`).
#'   The attribute `"hypa_fit"` records the fit of the null: `converged`,
#'   `iterations`, `max_error` (largest absolute strength error) and
#'   `n_cooccurrences` (\eqn{m}).
#'
#'   Errors of class `hypergraphs_bad_input` are raised for an invalid
#'   argument, and of class `hypergraphs_empty_result` when no pair
#'   co-occurs or no pair reaches `min_count` in either sense. A warning of
#'   class `hypergraphs_no_converge` is raised when the null cannot be
#'   fitted to tolerance.
#' @references
#' LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
#' and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies in
#' time series data on networks. *Proceedings of the 2020 SIAM International
#' Conference on Data Mining*, 460-468. \doi{10.1137/1.9781611976236.52}
#'
#' Chung, F., and Lu, L. (2002). Connected components in random graphs with
#' given expected degree sequences. *Annals of Combinatorics*, 6(2),
#' 125-145. \doi{10.1007/PL00012580}
#' @seealso [hypa.default()] for the same null on the paths of sequences,
#'   [hg_null_test()] for the Monte Carlo whole-hypergraph tests.
#' @examples
#' memberships <- data.frame(
#'   actor = c("a","b","c", "a","b","d", "a","b","e", "c","d","e", "a","b","f"),
#'   group = rep(c("g1","g2","g3","g4","g5"), each = 3)
#' )
#' hg <- group_hypergraph(memberships, node = "actor", hyperedge = "group")
#' hypa(hg, min_count = 2L)
#' @param ... Must be empty: an argument that only the sequence method
#'   takes (`order`, `action`, ...) raises `hypergraphs_bad_input`.
#' @export
hypa.net_hg <- function(x, min_count = 5L, alpha = 0.05, p_adjust = "none",
                        top = NULL, ...) {
  .ho_no_dots(..., .for = "a hypergraph")
  .thg_check_hg(x)
  min_count <- .ho_check_count(min_count, "min_count")
  if (!(is.numeric(alpha) && !is.factor(alpha) && length(alpha) == 1L &&
        !is.na(alpha) && alpha > 0 && alpha < 0.5)) {
    .ho_input_error("`alpha` must be one number in (0, 0.5)")
  }
  if (!is.null(top)) top <- .ho_check_count(top, "top")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !p_adjust %in% c(stats::p.adjust.methods, "none")) {
    .ho_input_error(
      "`p_adjust` must be one of stats::p.adjust.methods or \"none\"")
  }

  co <- .hg_hypa_cooccurrence(x$incidence)
  if (nrow(co$pairs) == 0L) {
    stop(errorCondition("no pair of nodes shares a hyperedge",
                        class = "hypergraphs_empty_result", call = NULL))
  }
  active <- which(co$strength > 0)
  fit <- .hg_hypa_fit(co$strength[active])
  scored <- .hg_hypa_scored_pairs(co$pairs, active, fit$x, min_count)
  if (nrow(scored) == 0L) {
    stop(errorCondition(
      sprintf(paste0("no pair co-occurs, or is expected to co-occur, at ",
                     "least %d times"), min_count),
      class = "hypergraphs_empty_result", call = NULL
    ))
  }
  tails <- .hg_hypa_tails(scored$observed, scored$expected, fit$m)

  out <- data.frame(
    from = x$nodes[active[scored$i]],
    to = x$nodes[active[scored$j]],
    observed = as.integer(scored$observed),
    expected = scored$expected,
    ratio = ifelse(scored$expected > 0, scored$observed / scored$expected,
                   NA_real_),
    p_under = tails$p_under,
    p_over = tails$p_over,
    stringsAsFactors = FALSE
  )
  # The family is every possible pair, as in the sequence HYPA, not only the
  # scored ones: the scored set is chosen partly by observed count, so
  # adjusting within it alone under-corrects (see @details). p.adjust(n =)
  # treats the unscored pairs as larger p-values, which is conservative
  # relative to adjusting the full family.
  n_possible <- choose(length(active), 2)
  adjust <- function(p) {
    if (identical(p_adjust, "none")) p else
      stats::p.adjust(p, p_adjust, n = n_possible)
  }
  out$p_adjusted_under <- adjust(out$p_under)
  out$p_adjusted_over <- adjust(out$p_over)
  is_under <- out$expected >= min_count & out$p_adjusted_under < alpha
  is_over <- out$observed >= min_count & out$p_adjusted_over < alpha
  out$anomaly <- ifelse(is_under, "under", ifelse(is_over, "over", "normal"))

  # Ordered by significance, not by `ratio`. Ratio is maximised by the rarest
  # pairs sitting just above `min_count` (a pair seen 5 times against an
  # expectation of 0.005 has a ratio near 1000 and says little), so ranking
  # by it buries the heavily co-occurring pairs that depart from the null.
  out <- out[.ho_order(out$p_adjusted_over, -out$ratio, out$from, out$to), ,
             drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(top)) out <- utils::head(out, top)
  attr(out, "hypa_fit") <- list(converged = fit$converged,
                                iterations = fit$iterations,
                                max_error = fit$max_error,
                                n_cooccurrences = fit$m)
  out
}

#' Observed co-occurrence counts and node strengths of a hypergraph
#'
#' @param incidence Node x hyperedge incidence matrix, dense or sparse; any
#'   non-zero entry is a membership.
#' @return A list: `pairs`, a data.frame with one row per co-occurring pair
#'   (`i < j`, row indices of `incidence`) and its count `co`; `strength`,
#'   the numeric vector \eqn{s_i = \sum_e B_{ie}(|e| - 1)}, one entry per
#'   node, which equals the row sums of the co-occurrence matrix without
#'   its diagonal.
#' @noRd
.hg_hypa_cooccurrence <- function(incidence) {
  b <- Matrix::Matrix(.thg_binary(incidence), sparse = TRUE)
  co <- Matrix::triu(Matrix::tcrossprod(b), k = 1L)
  trip <- Matrix::summary(methods::as(co, "TsparseMatrix"))
  pairs <- data.frame(i = as.integer(trip$i), j = as.integer(trip$j),
                      co = as.numeric(trip$x))
  pairs <- pairs[pairs$co > 0, , drop = FALSE]
  edge_size <- Matrix::colSums(b)
  strength <- as.numeric(b %*% pmax(edge_size - 1, 0))
  stopifnot(
    "co-occurrence pairs must be unique with i < j" =
      all(pairs$i < pairs$j) && !anyDuplicated(pairs[c("i", "j")]),
    "node strengths must sum to twice the co-occurrence total" =
      isTRUE(all.equal(sum(strength), 2 * sum(pairs$co)))
  )
  list(pairs = pairs, strength = strength)
}

#' Fit the fitness null of hypergraph HYPA
#'
#' Finds \eqn{x > 0} with \eqn{x_i (X - x_i) = s_i}, \eqn{X = \sum x}, by the
#' damped fixed point \eqn{x_i \leftarrow \sqrt{x_i s_i / (X - x_i)}} (the
#' geometric mean of the current value and the undamped update), started from
#' the Chung-Lu value \eqn{s_i / \sqrt{2m}}. Then \eqn{E_{ij} = x_i x_j}
#' preserves every strength: \eqn{\sum_{j \ne i} E_{ij} = s_i}.
#'
#' @param strength Positive node strengths (nodes without co-occurrence
#'   removed).
#' @param tol Relative tolerance: converged when the largest absolute
#'   strength error is below `tol * m`.
#' @param max_iter Maximum number of sweeps.
#' @return A list with `x` (fitness, one per node of `strength`), `m`
#'   (\eqn{\sum s / 2}), `max_error`, `iterations` and `converged`. A
#'   warning of class `hypergraphs_no_converge` is raised when the tolerance
#'   is not reached.
#' @noRd
.hg_hypa_fit <- function(strength, tol = 1e-10, max_iter = 10000L) {
  stopifnot(
    "`strength` must hold at least two positive, finite values" =
      is.numeric(strength) && length(strength) >= 2L &&
      all(is.finite(strength)) && all(strength > 0)
  )
  m <- sum(strength) / 2
  strength_error <- function(x) max(abs(x * (sum(x) - x) - strength))
  start <- strength / sqrt(2 * m)
  # Reduce over a bounded sweep count; once converged the state is carried
  # through unchanged, so the cost after convergence is one closure call.
  state <- Reduce(function(st, i) {
    if (st$done) return(st)
    x <- sqrt(st$x * strength / (sum(st$x) - st$x))
    err <- strength_error(x)
    list(x = x, err = err, iter = st$iter + 1L, done = err < tol * m)
  }, seq_len(max_iter),
  list(x = start, err = strength_error(start), iter = 0L,
       done = strength_error(start) < tol * m))
  if (!state$done) {
    warning(warningCondition(
      sprintf(paste0(
        "the HYPA null did not converge in %d sweeps (largest strength ",
        "error %.3g of %s co-occurrences); one node takes part in ",
        "(almost) every co-occurring pair, which leaves the null without ",
        "a fitted solution"),
        as.integer(max_iter), state$err, format(m)),
      class = "hypergraphs_no_converge", call = NULL))
  }
  list(x = state$x, m = m, max_error = state$err, iterations = state$iter,
       converged = state$done)
}

#' The pairs hypergraph HYPA scores
#'
#' The union of the pairs that co-occur at least `min_count` times and the
#' pairs whose expected count \eqn{x_i x_j} is at least `min_count`. The
#' second set is enumerated from the fitness sorted in decreasing order: for
#' the node at sorted position \eqn{p}, its partners at later positions with
#' \eqn{x_q \ge \mathrm{min\_count} / x_p} form a contiguous run, located by
#' [findInterval()]. The candidates are generated with a relative margin and
#' then filtered on the exact product, so rounding at the boundary cannot
#' drop a pair.
#'
#' @param pairs Co-occurring pairs (`i`, `j`, `co`) in original node indices.
#' @param active Original indices of the nodes with positive strength.
#' @param x Fitness of the active nodes.
#' @param min_count Count floor.
#' @return A data.frame (`i`, `j` in active-node indices with `i < j`,
#'   `observed`, `expected`), one row per scored pair, sorted by `(i, j)`.
#' @noRd
.hg_hypa_scored_pairs <- function(pairs, active, x, min_count) {
  n <- length(x)
  pos <- match(seq_len(max(active)), active)
  obs_i <- pos[pairs$i]
  obs_j <- pos[pairs$j]
  stopifnot("co-occurring nodes must have positive strength" =
              !anyNA(obs_i) && !anyNA(obs_j))
  pair_key <- function(i, j) (i - 1) * n + j

  over_set <- pairs$co >= min_count
  ord <- .ho_order(-x, seq_len(n))
  xs <- x[ord]
  threshold <- min_count * (1 - 1e-9) / xs
  n_at_least <- n - findInterval(threshold, rev(xs), left.open = TRUE)
  n_partner <- pmax(n_at_least - seq_len(n), 0L)
  p_idx <- rep(seq_len(n), n_partner)
  q_idx <- sequence(n_partner, from = seq_len(n) + 1L)
  under_i <- pmin(ord[p_idx], ord[q_idx])
  under_j <- pmax(ord[p_idx], ord[q_idx])
  under_set <- x[under_i] * x[under_j] >= min_count

  key <- .ho_sort(unique(c(pair_key(obs_i[over_set], obs_j[over_set]),
                       pair_key(under_i[under_set], under_j[under_set]))))
  i <- as.integer((key - 1) %/% n) + 1L
  j <- as.integer(key - (i - 1) * n)
  observed <- pairs$co[match(key, pair_key(obs_i, obs_j))]
  observed[is.na(observed)] <- 0
  out <- data.frame(i = i, j = j, observed = observed,
                    expected = x[i] * x[j])
  stopifnot(
    "scored pairs must have i < j" = all(out$i < out$j),
    "every scored pair must reach min_count in one sense" =
      all(out$observed >= min_count | out$expected >= min_count)
  )
  out
}

#' Hypergeometric tails of hypergraph HYPA
#'
#' @param observed Observed co-occurrence counts.
#' @param expected Fitted expected counts \eqn{E_{ij}}.
#' @param m Total co-occurrence count (the number of draws).
#' @return A list with `p_under = P(X <= f)` and the strict
#'   `p_over = P(X > f)`, \eqn{X \sim \mathrm{Hypergeometric}(m^2,
#'   \mathrm{round}(m E), m)}.
#' @noRd
.hg_hypa_tails <- function(observed, expected, m) {
  n_total <- round(m^2)
  k_white <- pmin(pmax(round(m * expected), 0), n_total)
  list(
    p_under = stats::phyper(observed, k_white, n_total - k_white, m),
    p_over = stats::phyper(observed, k_white, n_total - k_white, m,
                           lower.tail = FALSE)
  )
}
