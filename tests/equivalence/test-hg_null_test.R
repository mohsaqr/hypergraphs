# hg_null_test(): the observed statistics that are tested against a null model: the share of node pairs that
# share a hyperedge (pairwise participation), the mean Jaccard similarity over all pairs of hyperedges, the
# hyperedges whose member set already occurred, and the member pairs counted with multiplicity minus the distinct
# pairs. Reference: the definitions in base R on the incidence matrix. The null draws (swaps, configuration model)
# are made inside the function and are not compared.

test_that("the observed statistics equal their definitions", {
  set.seed(5)
  edges <- lapply(seq_len(25L), \(j) sample(12L, sample(2:5, 1L)))
  edges <- c(edges, edges[c(1L, 2L)])
  memberships <- data.frame(member = sprintf("n%02d", unlist(edges)), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges))))
  incidence <- ref_incidence(memberships, node = "member", edge = "edge")
  sizes <- colSums(incidence)
  shared <- tcrossprod(incidence)
  pairs <- sum(shared[upper.tri(shared)] > 0)
  overlap <- crossprod(incidence)
  union <- outer(sizes, sizes, `+`) - overlap
  sets <- vapply(split(memberships$member, memberships$edge), \(set) paste(sort(set), collapse = ","), character(1))
  expected <- data.frame(
    statistic = c("pairwise_participation", "avg_jaccard", "repeated_edges", "repeated_pairs"),
    observed = c(pairs / choose(nrow(incidence), 2), mean((overlap / union)[upper.tri(overlap)]),
                 length(sets) - length(unique(sets)), sum(choose(sizes, 2)) - pairs))
  found <- hg_null_test(group_hypergraph(memberships, node = "member", hyperedge = "edge"), n = 19L, seed = 1L)
  # Ratios and sums of small integers.
  expect_reference(found, expected, reference = "independent: the definitions on the incidence matrix, base R",
                   covers = "hg_null_test()", quantity = "observed statistics", tolerance = 1e-12, key = "statistic")
})
