# hg_cocluster(): spectral co-clustering of nodes and hyperedges (Dhillon 2001, Algorithm Multipartition), against
# scikit-learn's SpectralCoclustering, which implements the same algorithm. The embedding (what = "embedding") is
# deterministic: the scaled incidence D1^-1/2 A D2^-1/2, its singular vectors 2 to l + 1 (l = ceiling(log2 k)), and
# the rows [D1^-1/2 U; D2^-1/2 V]. scikit-learn's are read from the estimator's own steps (_scale_normalize() and
# _svd(), the private methods its fit() calls, with randomized SVD given n_svd_vecs = min(dim): the sketch then spans
# the whole row space, so it is exact to rounding; with the default 10 oversamples it is not, and differed by 3e-4) and compared up to the sign of each vector. The k-means
# starts of both sides are drawn inside, so the partitions are compared as "same co-cluster or not" on a planted
# incidence whose co-clusters every start finds.

cc_script <- c(
  "import sys, json, numpy as np",
  "from sklearn.cluster import SpectralCoclustering",
  "from sklearn.cluster._bicluster import _scale_normalize",
  "spec = json.load(open(sys.argv[1]))",
  "x = np.array(spec['incidence'], dtype=float)",
  "k = int(spec['k'])",
  "model = SpectralCoclustering(n_clusters=k, svd_method='randomized', n_svd_vecs=min(x.shape), random_state=0)",
  "normalized, row_diag, col_diag = _scale_normalize(x)",
  "n_sv = 1 + int(np.ceil(np.log2(k)))",
  "u, v = model._svd(normalized, n_sv, n_discard=1)",
  "z = np.vstack((row_diag[:, np.newaxis] * u, col_diag[:, np.newaxis] * v))",
  "model.fit(x)",
  "json.dump({'z': z.tolist(), 'rows': model.row_labels_.tolist(), 'columns': model.column_labels_.tolist()},",
  "          open(sys.argv[2], 'w'))")

# Four planted co-clusters: 6 nodes each, 8 hyperedges of 3 drawn inside a block, and 3 hyperedges across blocks.
cc_memberships <- function(seed) {
  set.seed(seed)
  nodes <- sprintf("n%02d", seq_len(24L))
  block <- rep(1:4, each = 6L)
  sets <- c(lapply(seq_len(32L), \(i) sample(nodes[block == (i - 1L) %% 4L + 1L], 3L)),
            list(c("n01", "n07"), c("n08", "n13"), c("n14", "n19")))
  data.frame(member = unlist(sets), edge = sprintf("e%02d", rep(seq_along(sets), lengths(sets))))
}

# Each row labelled by the first row of its co-cluster, so two partitions compare label by label.
cc_canonical <- function(row, cluster) {
  first <- tapply(row, cluster, min)
  data.frame(row = row, cocluster = unname(first[as.character(cluster)]))
}

test_that("the co-clustering embedding equals scikit-learn's SpectralCoclustering, up to the sign of each vector", {
  skip_without_python(c("scikit-learn" = "sklearn"))
  invisible(lapply(c(41L, 42L), \(seed) {
    memberships <- cc_memberships(seed)
    incidence <- ref_incidence(memberships, node = "member", edge = "edge")
    theirs <- ref_python_run(cc_script, list(incidence = unname(incidence), k = 4L))
    found <- hg_cocluster(group_hypergraph(memberships, node = "member", hyperedge = "edge"), k = 4L, seed = 1L,
                          what = "embedding")
    row <- paste(found$role, found$node)
    reference_rows <- c(paste("node", rownames(incidence)), paste("hyperedge", colnames(incidence)))
    ours <- as.matrix(subset(found, select = c(dim1, dim2)))
    z <- theirs$z[match(row, reference_rows), , drop = FALSE]
    # A singular vector is defined up to its sign: each of scikit-learn's is turned to the sign of ours.
    z <- sweep(z, 2L, sign(colSums(z * ours)), "*")
    expect_reference(data.frame(row = row, ours), data.frame(row = row, dim1 = z[, 1L], dim2 = z[, 2L]),
                     reference = "scikit-learn", covers = "hg_cocluster()",
                     quantity = sprintf("embedding, planted seed %d", seed), tolerance = 1e-8, key = "row",
                     note = paste(ref_python_version("scikit-learn"), "- randomized SVD with n_svd_vecs = min(dim), exact to rounding"))
    labels <- c(theirs$rows, theirs$columns)[match(row, reference_rows)]
    expect_reference(cc_canonical(row, found$cluster), cc_canonical(row, labels), reference = "scikit-learn",
                     covers = "hg_cocluster()", quantity = sprintf("co-cluster of each node and hyperedge, planted seed %d", seed),
                     tolerance = 0, key = "row", note = ref_python_version("scikit-learn"))
  }))
})
