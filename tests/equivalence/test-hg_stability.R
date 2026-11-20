# hg_stability(): cluster-wise subsampling stability of hg_cluster() (Hennig 2007; the subsetting scheme of Ben-Hur
# et al. 2002). The subsamples are drawn inside the function from `seed`. Its documentation states the procedure:
#   - the nodes are clustered once in full (hg_cluster() with the same k, type, seed and nstart);
#   - n_boot times, floor(fraction * n) nodes are drawn without replacement, the drawn subsamples being the ones
#     fpc::clusterboot(bootmethod = "subset", subtuning = floor(fraction * n), seed = seed) draws;
#   - the sub-hypergraph restricts every hyperedge to the drawn nodes, drops the hyperedges left empty and keeps the
#     incidence weights, and is clustered again; a subsample whose hypergraph is disconnected is not scored;
#   - each full cluster, restricted to the drawn nodes, is matched to its most similar subsample cluster by the
#     Jaccard coefficient (0 when it has no drawn member); its stability is the mean over the subsamples
#     (clusterboot's subsetmean), with the counts of subsamples at or below 0.5 (subsetbrd) and above 0.75
#     (subsetrecover).
# That rule is written below as the clustering method fpc::clusterboot() calls, and compared on three seeds. A Monte
# Carlo comparison over ten further seeds checks the stability summary against the same procedure written in base R
# with its own draws.

#' A planted hypergraph: three blocks of 8 nodes, 70 hyperedges of 3 to 6 nodes, every fifth one drawn from all nodes
planted_hypergraph <- function(seed, n_nodes = 24L, n_edges = 70L, k = 3L) {
  set.seed(seed)
  block <- rep(seq_len(k), length.out = n_nodes)
  edges <- lapply(seq_len(n_edges), \(i) {
    pool <- if (i %% 5L == 0L) seq_len(n_nodes) else which(block == sample(k, 1L))
    sample(pool, min(sample(3:6, 1L), length(pool)))
  })
  data.frame(member = sprintf("n%02d", unlist(edges)), edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges))))
}

memberships <- planted_hypergraph(7L)
hg <- group_hypergraph(memberships, node = "member", hyperedge = "edge")
nodes <- hg_get(hg, what = "nodes")$node
ks <- 2:4
n_boot <- 30L
drawn_size <- floor(0.5 * length(nodes))

#' The documented restriction: every hyperedge restricted to the drawn nodes, hyperedges left empty dropped
restrict <- function(drawn) group_hypergraph(memberships[memberships$member %in% drawn, ], node = "member", hyperedge = "edge")

#' The partition of the drawn nodes by hg_cluster(), clusters numbered by first appearance in the drawn order
partition <- function(drawn, k, seed) {
  clusters <- hg_cluster(restrict(drawn), k = k, seed = seed, nstart = 25L)
  label <- clusters$cluster[match(drawn, clusters$node)]
  match(label, unique(label))
}

#' The clustering method fpc::clusterboot() calls on the rows it draws: the rows are node indices. clusterboot()
#' keeps its own `seed` argument (the draws), so the k-means seed of hg_cluster() is passed as `kmeans_seed`.
clusterboot_method <- function(rows, k, kmeans_seed) {
  groups <- partition(nodes[rows[, 1L]], k, kmeans_seed)
  list(result = groups, nc = max(groups), clusterlist = lapply(seq_len(max(groups)), \(j) groups == j), partition = groups,
       clustermethod = "hg_cluster() on the restricted hypergraph")
}

#' fpc::clusterboot()'s cluster-wise subsampling statistics for one k and seed
clusterboot_table <- function(k, seed) {
  boot <- fpc::clusterboot(cbind(seq_along(nodes), seq_along(nodes)), B = n_boot, distances = FALSE, bootmethod = "subset",
                           subtuning = drawn_size, clustermethod = clusterboot_method, k = k, kmeans_seed = seed, seed = seed, count = FALSE)
  data.frame(k = k, cluster = paste("Cluster", seq_len(boot$nc)), jaccard = boot$subsetmean, n_dissolved = boot$subsetbrd,
             n_recovered = boot$subsetrecover, n_runs = n_boot)
}

test_that("the cluster-wise stability equals fpc::clusterboot(bootmethod = 'subset') on the documented subsamples", {
  skip_without_reference("fpc")
  invisible(lapply(1:3, \(seed) {
    expected <- do.call(rbind, lapply(ks, clusterboot_table, seed = seed))
    expect_reference(hg_stability(hg, k = ks, n_boot = n_boot, seed = seed, what = "clusters"), expected, reference = "fpc",
                     covers = "hg_stability()", quantity = sprintf("seed %d, stability of each cluster, k = 2 to 4", seed),
                     tolerance = 1e-12, key = c("k", "cluster"))
    summary_k <- do.call(rbind, lapply(split(expected, expected$k), \(t) {
      data.frame(k = t$k[[1L]], n_runs = n_boot, n_failed = 0, mean_jaccard = mean(t$jaccard), min_jaccard = min(t$jaccard),
                 n_stable = sum(t$jaccard >= 0.75), n_dissolved = sum(t$jaccard <= 0.5))
    }))
    expect_reference(hg_stability(hg, k = ks, n_boot = n_boot, seed = seed), summary_k, reference = "fpc",
                     covers = "hg_stability()", quantity = sprintf("seed %d, stability of each k", seed), tolerance = 1e-12, key = "k")
  }))
})

test_that("the eigengap of every k is the gap of the Zhou Laplacian's spectrum", {
  incidence <- ref_incidence(memberships, node = "member", edge = "edge") > 0
  degree_root <- diag(1 / sqrt(rowSums(incidence)))
  laplacian <- diag(nrow(incidence)) - degree_root %*% incidence %*% diag(1 / colSums(incidence)) %*% t(incidence) %*% degree_root
  values <- sort(eigen(laplacian, symmetric = TRUE, only.values = TRUE)$values)
  # eigenvalues of a symmetric matrix of order 24 by LAPACK on both sides
  expect_reference(subset(hg_stability(hg, k = ks, n_boot = 5L, seed = 1L), select = c(k, eigengap)),
                   data.frame(k = ks, eigengap = values[ks + 1L] - values[ks]),
                   reference = "independent: eigen() of the Zhou Laplacian, base R", covers = "hg_stability()",
                   quantity = "eigengap after k", tolerance = 1e-10, key = "k")
})

test_that("the mean stability over ten seeds agrees with the documented procedure run on its own draws", {
  #' The procedure of the documentation, in base R with its own subsamples: mean over the clusters of each cluster's
  #' mean best Jaccard over the scored subsamples
  procedure <- function(k, seed, draw_seed) {
    full <- partition(nodes, k, seed)
    set.seed(draw_seed)
    draws <- lapply(seq_len(n_boot), \(b) sample(length(nodes), drawn_size))
    scores <- lapply(draws, \(drawn) {
      # a disconnected subsample is not scored (the documented rule): hg_cluster() raises this condition for it
      sub <- tryCatch(partition(nodes[drawn], k, seed), hypergraphs_hypergraph_disconnected = \(e) NULL)
      if (is.null(sub)) return(NULL)
      vapply(seq_len(max(full)), \(j) {
        members <- full[drawn] == j
        if (!any(members)) return(0)
        max(vapply(unique(sub), \(l) sum(members & sub == l) / sum(members | sub == l), numeric(1)))
      }, numeric(1))
    })
    scored <- do.call(cbind, Filter(Negate(is.null), scores))
    mean(rowMeans(scored))
  }
  ours <- vapply(1:10, \(seed) {
    # a subsample that disconnects the hypergraph is reported by this documented warning; it is the rule compared
    stability <- withCallingHandlers(hg_stability(hg, k = ks, n_boot = n_boot, seed = seed),
                                     hypergraphs_hypergraph_disconnected = \(w) invokeRestart("muffleWarning"))
    stability$mean_jaccard
  }, numeric(length(ks)))
  theirs <- vapply(1:10, \(seed) vapply(ks, procedure, numeric(1), seed = seed, draw_seed = 1000L + seed), numeric(length(ks)))
  invisible(lapply(seq_along(ks), \(i) {
    # Tolerance fixed before the run: four standard errors of the difference of two means of 10 independent seeds,
    # each standard error from that side's own spread over its seeds. In a study before this test (6 blocks of 10
    # seeds, seeds 1 to 60) the three k passed in 6 of 6 blocks; the largest difference was 0.37 of its tolerance.
    tolerance <- 4 * sqrt(stats::var(ours[i, ]) / 10 + stats::var(theirs[i, ]) / 10)
    expect_reference(data.frame(k = ks[[i]], mean_jaccard = mean(ours[i, ])), data.frame(k = ks[[i]], mean_jaccard = mean(theirs[i, ])),
                     reference = "independent: cluster-wise subsampling stability, Hennig (2007), as documented, base R",
                     covers = "hg_stability()",
                     quantity = sprintf("Monte Carlo summary: mean over 10 seeds of the mean cluster stability, k = %d", ks[[i]]),
                     tolerance = tolerance, key = "k", note = sprintf("tolerance %.4g = 4 standard errors of the difference", tolerance))
  }))
})
