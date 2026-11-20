# hg_mmsbm(): the hypergraph mixed-membership stochastic block model of Ruggeri, Contisciani, Battiston and De Bacco
# (2023, Science Advances 9, eadg9159), fitted by their EM from random starts drawn inside the function. The
# authors' Hy-MMSBM code is not an installed package, so the fit is checked by what holds for any start and by
# recovery of a planted structure:
#   - the memberships of a node sum to one, and the parameters u and w are non-negative;
#   - the reported log-likelihood equals the model's likelihood (Eq. 5) recomputed here from the returned u and w;
#   - on a planted two-block hypergraph the hard communities are the planted blocks, on each of 12 seeds;
#   - summary() reads its communities table off the fitted memberships and affinity.

#' A planted two-block hypergraph: 300 hyperedges of 2 or 3 nodes, 95% of them inside one block
planted_blocks <- function(seed, n = 30L, m = 300L, cross = 0.05) {
  set.seed(seed)
  block <- rep(1:2, each = n / 2L)
  edges <- lapply(seq_len(m), \(i) {
    size <- sample(2:3, 1L)
    pool <- if (stats::runif(1L) < cross) seq_len(n) else which(block == sample(2L, 1L))
    sample(pool, size)
  })
  list(block = stats::setNames(block, sprintf("n%02d", seq_len(n))),
       memberships = data.frame(node = sprintf("n%02d", unlist(edges)), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges)))))
}

#' The Hy-MMSBM log-likelihood, Ruggeri et al. (2023) Eq. 5, up to the constants free of u and w:
#' L = -C sum_{i<j} u_i' w u_j + sum_e A_e log(lambda_e), lambda_e = sum_{i<j in e} u_i' w u_j,
#' C = sum_{n=2}^{D} 2 / (n (n - 1)), D the largest hyperedge size. Every A_e is 1 here.
mmsbm_loglik <- function(u, w, edges) {
  affinity <- u %*% w %*% t(u)
  lambda <- vapply(edges, \(members) {
    block <- affinity[members, members, drop = FALSE]
    sum(block[upper.tri(block)])
  }, numeric(1))
  sizes <- lengths(edges)
  C <- sum(vapply(seq(2L, max(sizes)), \(s) 2 / (s * (s - 1)), numeric(1)))
  -C * sum(affinity[upper.tri(affinity)]) + sum(log(lambda))
}

seeds <- 1:12
fits <- lapply(seeds, \(seed) {
  planted <- planted_blocks(seed)
  hg <- group_hypergraph(planted$memberships, node = "node", hyperedge = "edge")
  list(seed = seed, planted = planted, hg = hg, fit = hg_mmsbm(hg, k = 2L, seed = seed))
})

test_that("memberships sum to one over the communities and u, w are non-negative, for every start", {
  membership <- do.call(rbind, lapply(fits, \(f) transform(hg_get(f$fit, what = "membership"), seed = f$seed)))
  totals <- stats::aggregate(membership ~ seed + node, data = membership, FUN = sum)
  expect_reference(totals, transform(totals, membership = 1), reference = "independent: memberships u_i / sum_q u_iq sum to one",
                   covers = c("hg_mmsbm()", "hg_get.net_hg_mmsbm()"), quantity = "12 seeds, sum of a node's memberships",
                   tolerance = 1e-12, key = c("seed", "node"))
  weights <- subset(membership, select = c(seed, node, community, membership_weight))
  expect_reference(weights, transform(weights, membership_weight = pmax(membership_weight, 0)),
                   reference = "independent: the membership parameters u are non-negative", covers = c("hg_mmsbm()", "hg_get.net_hg_mmsbm()"),
                   quantity = "12 seeds, membership parameter u", tolerance = 0, key = c("seed", "node", "community"))
  affinity <- do.call(rbind, lapply(fits, \(f) transform(hg_get(f$fit, what = "affinity"), seed = f$seed)))
  expect_reference(affinity, transform(affinity, affinity = pmax(affinity, 0)), reference = "independent: the affinity w is non-negative",
                   covers = c("hg_mmsbm()", "hg_get.net_hg_mmsbm()"), quantity = "12 seeds, affinity w", tolerance = 0,
                   key = c("seed", "from", "to"))
})

test_that("the reported log-likelihood is the Hy-MMSBM likelihood of the returned parameters", {
  rows <- lapply(fits, \(f) {
    membership <- hg_get(f$fit, what = "membership")
    nodes <- unique(membership$node)
    communities <- unique(membership$community)
    u <- ref_wide(membership, "node", "community", "membership_weight", nodes, communities)
    affinity <- hg_get(f$fit, what = "affinity")
    w <- ref_wide(affinity, "from", "to", "affinity", communities, communities)
    w <- w + t(w) - diag(diag(w))
    incidence <- hg_get(f$hg, what = "memberships")
    edges <- lapply(split(incidence$node, incidence$hyperedge), \(members) match(members, nodes))
    data.frame(seed = f$seed, ours = subset(hg_get(f$fit, what = "restarts"), best)$log_likelihood, theirs = mmsbm_loglik(u, w, edges))
  })
  rows <- do.call(rbind, rows)
  expect_reference(data.frame(seed = rows$seed, log_likelihood = rows$ours), data.frame(seed = rows$seed, log_likelihood = rows$theirs),
                   reference = "independent: Hy-MMSBM log-likelihood, Ruggeri et al. (2023) Eq. 5",
                   covers = c("hg_mmsbm()", "hg_get.net_hg_mmsbm()"), quantity = "12 seeds, log-likelihood of the kept start",
                   tolerance = 1e-10, measure = "relative", key = "seed")
})

test_that("the communities of a planted two-block hypergraph are the planted blocks, on every one of 12 seeds", {
  skip_without_reference("clue")
  matched <- lapply(fits, \(f) {
    nodes <- hg_get(f$fit, what = "nodes")
    blocks <- f$planted$block[nodes$node]
    agreement <- unclass(table(factor(nodes$community, levels = c("Community 1", "Community 2")), factor(blocks, levels = 1:2)))
    # each community matched to one block by the Hungarian method, maximising the nodes they share
    assignment <- as.integer(clue::solve_LSAP(agreement, maximum = TRUE))
    labels <- c("Community 1", "Community 2")[match(1:2, assignment)]
    list(ours = data.frame(seed = f$seed, node = nodes$node, community = nodes$community),
         planted = data.frame(seed = f$seed, node = nodes$node, community = labels[blocks]))
  })
  # Recovery is exact or not, node by node. Before this test, a study of 30 seeds of the same design recovered every
  # node on 30 of 30 (default nstart = 10); 12 seeds are declared here.
  expect_reference(do.call(rbind, lapply(matched, `[[`, "ours")), do.call(rbind, lapply(matched, `[[`, "planted")),
                   reference = "independent: the planted blocks of the simulated hypergraph", covers = c("hg_mmsbm()", "hg_get.net_hg_mmsbm()"),
                   quantity = "Monte Carlo summary: 12 seeds, community of every node matched to the planted block", tolerance = 0,
                   key = c("seed", "node"))
})

test_that("summary() of a fit tabulates the communities from the memberships and the affinity", {
  fit <- fits[[1L]]$fit
  nodes <- hg_get(fit, what = "nodes")
  membership <- hg_get(fit, what = "membership")
  affinity <- hg_get(fit, what = "affinity")
  communities <- sort(unique(membership$community))
  expected <- data.frame(
    community = communities,
    size = vapply(communities, \(c) sum(nodes$community == c), numeric(1)),
    mass = vapply(communities, \(c) sum(membership$membership[membership$community == c]), numeric(1)),
    mean_membership = vapply(communities, \(c) mean(nodes$membership[nodes$community == c]), numeric(1)),
    w_within = vapply(communities, \(c) affinity$affinity[affinity$from == c & affinity$to == c], numeric(1)))
  expect_reference(hg_get(summary(fit), what = "communities"), expected,
                   reference = "independent: community sizes, membership mass, mean membership and within-affinity from the fitted tables",
                   covers = c("summary.net_hg_mmsbm()", "hg_get.net_hg_mmsbm()"), quantity = "seed 1, communities", tolerance = 1e-12,
                   key = "community")
})
