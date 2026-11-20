# hg_agreement(): the agreement of two labellings of the same nodes: the adjusted Rand index (Hubert and Arabie
# 1985), adjusted and normalised mutual information, and the cross-tabulation. References: scikit-learn's
# adjusted_rand_score, adjusted_mutual_info_score and normalized_mutual_info_score (arithmetic mean), and the
# adjusted Rand index and contingency table written out in base R.

labellings <- lapply(1:3, \(seed) {
  set.seed(seed)
  list(x = data.frame(node = sprintf("n%02d", seq_len(40L)), cluster = sample(1:3, 40L, replace = TRUE)),
       y = data.frame(node = sprintf("n%02d", seq_len(40L)), cluster = sample(1:4, 40L, replace = TRUE)))
})

test_that("the adjusted Rand index and the cross-tabulation equal their definitions", {
  invisible(Map(\(pair, index) {
    table_xy <- table(pair$x$cluster, pair$y$cluster)
    pairs <- \(n) n * (n - 1) / 2
    index_sum <- sum(pairs(table_xy))
    expected_index <- sum(pairs(rowSums(table_xy))) * sum(pairs(colSums(table_xy))) / pairs(sum(table_xy))
    maximum <- (sum(pairs(rowSums(table_xy))) + sum(pairs(colSums(table_xy)))) / 2
    ari <- (index_sum - expected_index) / (maximum - expected_index)
    expect_reference(hg_agreement(pair$x, pair$y, method = "ari")$ari, ari,
                     reference = "independent: adjusted Rand index of Hubert and Arabie (1985), base R", covers = "hg_agreement()",
                     quantity = sprintf("labellings %d, ARI", index), tolerance = 1e-12)
    crossed <- subset(as.data.frame(table_xy, stringsAsFactors = FALSE), Freq > 0)
    expect_reference(hg_agreement(pair$x, pair$y, what = "table"),
                     data.frame(label_x = crossed$Var1, label_y = crossed$Var2, n = crossed$Freq),
                     reference = "independent: contingency table, base R", covers = "hg_agreement()",
                     quantity = sprintf("labellings %d, cross-tabulation", index), tolerance = 0, key = c("label_x", "label_y"))
  }, labellings, seq_along(labellings)))
})

test_that("ARI, AMI and NMI equal scikit-learn", {
  skip_without_python(c("scikit-learn" = "sklearn"))
  script <- c(
    "import sys, json",
    "from sklearn.metrics import adjusted_mutual_info_score as ami, adjusted_rand_score as ari, normalized_mutual_info_score as nmi",
    "spec = json.load(open(sys.argv[1]))",
    "json.dump({'ari': ari(spec['a'], spec['b']), 'ami': ami(spec['a'], spec['b']), 'nmi': nmi(spec['a'], spec['b'])}, open(sys.argv[2], 'w'))")
  version <- ref_python_version("scikit-learn")
  invisible(Map(\(pair, index) {
    theirs <- ref_python_run(script, list(a = pair$x$cluster, b = pair$y$cluster))
    # Entropies and the expected mutual information of a small contingency table.
    expect_reference(unlist(subset(hg_agreement(pair$x, pair$y, method = c("ari", "ami", "nmi")), select = c(ari, ami, nmi))),
                     unlist(theirs), reference = "scikit-learn", covers = "hg_agreement()",
                     quantity = sprintf("labellings %d, ARI, AMI, NMI", index), tolerance = 1e-12, note = version)
  }, labellings, seq_along(labellings)))
})
