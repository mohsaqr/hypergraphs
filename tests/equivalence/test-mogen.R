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
