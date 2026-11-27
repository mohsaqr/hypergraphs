# mogen(): the multi-order generative model of Scholtes (2017). For each highest order k it reports the
# log-likelihood of the sequences under the model of orders 0 to k, its degrees of freedom sum_{l <= k} s^l (s - 1)
# (s states), AIC = 2 df - 2 LL and BIC = log(n) df - 2 LL, n the number of events; hg_get(what = "paths") counts
# the k-grams. References: pathpy's MultiOrderModel (the author's implementation) for the log-likelihoods and its
# HigherOrderNetwork for the k-gram counts, and the definitions in base R.

pathpy_preamble <- c(
  "import sys, json, scipy, numpy",
  "# pathpy 2.2.0 imports names that SciPy has since removed; they are taken from NumPy, pathpy is not changed",
  "for name in dir(numpy):",
  "    if not hasattr(scipy, name):",
  "        try: setattr(scipy, name, getattr(numpy, name))",
  "        except Exception: pass",
  "import pathpy as pp",
  "spec = json.load(open(sys.argv[1]))",
  "paths = pp.Paths()",
  "for s in spec['sequences']:",
  "    if len(s) >= 2: paths.add_path(','.join(s), frequency=1)")

inputs <- list(second_order = ref_second_order(), random = ref_random_sequences(20L, 25L, c("a", "b", "c", "d"), seed = 3L),
               longer = ref_second_order(n_sequences = 60L, length = 30L, seed = 8L))

test_that("the log-likelihood of each highest order equals pathpy's MultiOrderModel", {
  skip_without_python("pathpy")
  script <- c(pathpy_preamble,
              "out = {'order': [], 'log_likelihood': []}",
              "for k in range(0, spec['max_order'] + 1):",
              "    model = pp.MultiOrderModel(paths, max_order=k)",
              "    out['order'].append(k); out['log_likelihood'].append(float(model.likelihood(paths, log=True)))",
              "json.dump(out, open(sys.argv[2], 'w'))")
  version <- ref_python_version("pathpy")
  invisible(Map(\(sequences, name) {
    theirs <- as.data.frame(ref_python_run(script, list(sequences = sequences, max_order = 3L)))
    # Sums of a few hundred logarithms on each side.
    expect_reference(hg_get(mogen(sequences, max_order = 3L), what = "orders"), theirs, reference = "pathpy",
                     covers = c("mogen()", "hg_get.net_mogen()", "hg_get.hypergraphs_result()"),
                     quantity = paste(name, "log-likelihood"), tolerance = 1e-8, key = "order", note = version)
  }, inputs, names(inputs)))
})

test_that("k-gram counts equal pathpy's higher-order networks and the sequences' k-grams", {
  skip_without_python("pathpy")
  script <- c(pathpy_preamble,
              "out = {'k': [], 'path': [], 'count': []}",
              "for k in (1, 2):",
              "    network = pp.HigherOrderNetwork(paths, k=k)",
              "    for (u, v), attributes in network.edges.items():",
              "        weight = attributes['weight']",
              "        states = u.split(',') + [v.split(',')[-1]]",
              "        out['k'].append(k + 1); out['path'].append(' -> '.join(states)); out['count'].append(float(weight[0] + weight[1]))",
              "json.dump(out, open(sys.argv[2], 'w'))")
  version <- ref_python_version("pathpy")
  sequences <- inputs$second_order
  theirs <- as.data.frame(ref_python_run(script, list(sequences = sequences)))
  fit <- mogen(sequences, max_order = 2L)
  ours <- rbind(hg_get(fit, what = "paths", k = 2L), hg_get(fit, what = "paths", k = 3L))
  expect_reference(subset(ours, select = c(k, path, count)), theirs, reference = "pathpy", covers = c("mogen()", "hg_get.net_mogen()"),
                   quantity = "second-order chain, k-gram counts", tolerance = 0, key = c("k", "path"),
                   note = paste(version, "HigherOrderNetwork(k) edge weights, subpath plus longest-path counts"))
  grams <- c(ref_kgrams(sequences, 2L), ref_kgrams(sequences, 3L))
  expect_reference(subset(ours, select = c(path, count)), data.frame(path = names(grams), count = unname(grams)),
                   reference = "independent: k-gram counts of the sequences, base R", covers = c("mogen()", "hg_get.net_mogen()"),
                   quantity = "second-order chain, k-gram counts", tolerance = 0, key = "path")
})

test_that("degrees of freedom, AIC and BIC equal their definitions", {
  invisible(Map(\(sequences, name) {
    orders <- hg_get(mogen(sequences, max_order = 3L), what = "orders")
    s <- length(unique(unlist(sequences)))
    df <- cumsum(s^(0:3) * (s - 1))
    events <- length(unlist(sequences))
    expected <- data.frame(order = 0:3, df = df, aic = 2 * df - 2 * orders$log_likelihood,
                           bic = log(events) * df - 2 * orders$log_likelihood)
    expect_reference(orders, expected, reference = "independent: df = sum s^l (s - 1), AIC and BIC from the log-likelihood, base R",
                     covers = c("mogen()", "hg_get.net_mogen()"), quantity = paste(name, "df, AIC, BIC"), tolerance = 1e-9,
                     key = "order", note = "The log-likelihood is the package's (compared with pathpy above); df and n are counted here.")
  }, inputs, names(inputs)))
})

# Order selection (Scholtes 2017, Sec. 4.1). Layer k of the multi-order model has one parameter per path of k
# transitions the first-order network allows, less one per row of its transition matrix (a path of k - 1 transitions
# that can be extended): Eq. 7 as its footnote states it, which on a complete network gives |V|^k (|V| - 1). pathpy
# 2.2.0 codes Eq. 7 as printed, subtracting the non-zero rows of A^k, at most |V|, so its 'paths' count exceeds that
# from order 2 on, and on a complete network exceeds |V|^k (|V| - 1), the count its own docstring says it coincides
# with. The p-value of order k against k - 1 is the chi-square tail of -2 (LL_{k-1} - LL_k) on layer k's count (Eq. 6);
# the selected order is the highest one rejected at the level, as in pathpy's estimate_order().

#' A sequence model on a sparse first-order network: from each state only two (or one) states follow
#'
#' @return A list of character vectors.
ref_sparse_sequences <- function(n_sequences = 40L, length = 15L, seed = 5L) {
  following <- list(a = c("b", "c"), b = c("c", "d"), c = c("a", "d"), d = "a")
  set.seed(seed)
  lapply(seq_len(n_sequences), \(i) {
    Reduce(\(so_far, t) {
      options <- following[[utils::tail(so_far, 1L)]]
      c(so_far, options[sample.int(length(options), 1L)])
    }, seq_len(length - 1L), "a")
  })
}

#' Layer degrees of freedom by enumerating the paths of the observed first-order network
#'
#' @param sequences A list of character vectors.
#' @param max_order The highest order.
#' @return A data frame `order`, `layer_df`.
ref_layer_df <- function(sequences, max_order) {
  states <- sort(unique(unlist(sequences)))
  edges <- unique(do.call(rbind, lapply(sequences, \(s) {
    data.frame(from = utils::head(s, -1L), to = utils::tail(s, -1L))
  })))
  # every path of k transitions, as a data frame with one column per state
  walks <- Reduce(\(paths, k) {
    joined <- merge(paths, edges, by.x = paste0("s", k - 1L), by.y = "from")
    names(joined)[names(joined) == "to"] <- paste0("s", k)
    joined[, paste0("s", 0:k), drop = FALSE]
  }, seq_len(max_order), init = data.frame(s0 = states), accumulate = TRUE)
  layer <- vapply(seq_len(max_order), \(k) {
    paths <- walks[[k + 1L]]
    rows <- unique(paths[, paste0("s", 0:(k - 1L)), drop = FALSE])
    nrow(paths) - nrow(rows)
  }, numeric(1))
  data.frame(order = 0:max_order, layer_df = c(length(states) - 1, layer))
}

#' One side of a paired table, its value column renamed
#'
#' @param table A data frame with `input`, optionally `order`, and the columns `ours` and `pathpy`.
#' @param side `"ours"` or `"pathpy"`.
#' @param name The name of the value column.
#' @return A data frame of the key columns and `name`.
ref_side <- function(table, side, name) {
  keys <- intersect(c("input", "order"), names(table))
  out <- table[keys]
  out[[name]] <- table[[side]]
  out
}

selection_inputs <- c(inputs, list(sparse = ref_sparse_sequences()))

pathpy_selection <- c(pathpy_preamble,
  "K = spec['max_order']",
  "model = pp.MultiOrderModel(paths, max_order=K)",
  "out = {'order': list(range(K + 1)),",
  "       'df': [model.degrees_of_freedom(max_order=k, assumption='paths') for k in range(K + 1)],",
  "       'p_value': [None, None] + [float(model.likelihood_ratio_test(paths, max_order_null=k - 1, max_order=k)[1])",
  "                                  for k in range(2, K + 1)],",
  "       'selected': model.estimate_order(significance_threshold=spec['alpha'])}",
  "json.dump(out, open(sys.argv[2], 'w'))")

test_that("each layer's degrees of freedom count the allowed paths less the extendable rows", {
  invisible(Map(\(sequences, name) {
    orders <- hg_get(mogen(sequences, max_order = 3L), what = "orders")
    expect_reference(subset(orders, select = c(order, layer_df)), ref_layer_df(sequences, 3L),
                     reference = "independent: paths of the observed first-order network enumerated in base R, less the distinct prefixes (Scholtes 2017, Eq. 7 and footnote 1)",
                     covers = c("mogen()", "hg_get.net_mogen()"), quantity = paste(name, "layer degrees of freedom"),
                     tolerance = 0, key = "order")
  }, selection_inputs, names(selection_inputs)))
})

test_that("the likelihood-ratio p-value is the chi-square tail of the log-likelihood ratio", {
  invisible(Map(\(sequences, name) {
    orders <- hg_get(mogen(sequences, max_order = 3L), what = "orders")
    tested <- subset(orders, order >= 2L)
    statistic <- -2 * (utils::head(orders$log_likelihood, -1L) - utils::tail(orders$log_likelihood, -1L))
    statistic <- utils::tail(statistic, -1L)
    expected <- data.frame(order = tested$order,
                           p_value = stats::pchisq(statistic, df = tested$layer_df, lower.tail = FALSE))
    expect_reference(subset(tested, select = c(order, p_value)), expected,
                     reference = "independent: chi-square tail of -2 (LL_{k-1} - LL_k) on the layer's degrees of freedom, base R (Scholtes 2017, Eq. 6)",
                     covers = c("mogen()", "hg_get.net_mogen()"), quantity = paste(name, "likelihood-ratio p-value"),
                     tolerance = 1e-12, key = "order",
                     note = "The log-likelihoods and layer counts are the package's, each compared above.")
  }, selection_inputs, names(selection_inputs)))
})

test_that("the selected order and the counts against pathpy's MultiOrderModel", {
  skip_without_python("pathpy")
  version <- ref_python_version("pathpy")
  dof_reason <- paste("pathpy 2.2.0 codes Scholtes (2017) Eq. 7 as printed, subtracting the non-zero rows of A^k (at most |V|)",
                      "rather than one per row of the order-k transition matrix as the paper's footnote 1 states; from order 2",
                      "its count is larger, and on a complete network it exceeds |V|^k (|V| - 1).")
  runs <- Map(\(sequences, name) {
    theirs <- ref_python_run(pathpy_selection, list(sequences = sequences, max_order = 3L, alpha = 0.01))
    orders <- hg_get(mogen(sequences, max_order = 3L, criterion = "lrt", lrt_alpha = 0.01), what = "orders")
    tested <- subset(orders, order >= 2L)
    list(df = data.frame(input = name, order = orders$order, ours = orders$df, pathpy = theirs$df),
         p_value = data.frame(input = name, order = tested$order, ours = tested$p_value,
                              pathpy = utils::tail(as.numeric(theirs$p_value), 2L)),
         selected = data.frame(input = name, ours = subset(orders, optimal)$order, pathpy = as.integer(theirs$selected)))
  }, selection_inputs, names(selection_inputs))
  pooled <- lapply(c(df = "df", p_value = "p_value", selected = "selected"), \(part) {
    do.call(rbind, lapply(runs, \(run) run[[part]]))
  })
  expect_known_difference(ref_side(pooled$df, "ours", "df"),
                          ref_side(pooled$df, "pathpy", "df"),
                          reference = "pathpy", covers = c("mogen()", "hg_get.net_mogen()"),
                          quantity = "degrees of freedom of each highest order", tolerance = 0, key = c("input", "order"),
                          reason = dof_reason, agrees_with = "orders 0 and 1", note = version)
  # Where one test rejects overwhelmingly (p near 0) or not at all (p near 1), the two p-values agree although
  # their degrees of freedom differ; the difference is judged over all inputs together.
  expect_known_difference(ref_side(pooled$p_value, "ours", "p_value"),
                          ref_side(pooled$p_value, "pathpy", "p_value"),
                          reference = "pathpy", covers = c("mogen()", "hg_get.net_mogen()"),
                          quantity = "likelihood-ratio p-value of each order against the one below", tolerance = 1e-8,
                          key = c("input", "order"),
                          reason = paste("The p-values are taken on different degrees of freedom:", dof_reason),
                          note = version)
  selected <- pooled$selected
  expect_reference(ref_side(selected, "ours", "selected"),
                   ref_side(selected, "pathpy", "selected"),
                   reference = "pathpy", covers = c("mogen()", "hg_get.net_mogen()"),
                   quantity = "order selected by the likelihood-ratio test at 0.01", tolerance = 0, key = "input",
                   note = paste(version, "MultiOrderModel.estimate_order(significance_threshold = 0.01)"))
})
