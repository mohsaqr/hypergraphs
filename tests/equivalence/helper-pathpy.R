# pathpy's higher-order networks as k-gram counts. Nothing here calls hypergraphs.

#' The transition counts of pathpy's HigherOrderNetwork of each order k: one row per edge
#'
#' The weight of an edge is pathpy's subpath count plus its longest-path count; for a k-th order network the
#' edge (u, v) is the context u (k states) followed by the last state of v.
#'
#' @param sequences A list of character vectors.
#' @param orders The orders k.
#' @return A data frame `k`, `context` (states joined by " -> "), `following`, `count`.
ref_pathpy_kgrams <- function(sequences, orders) {
  script <- c(
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
    "    if len(s) >= 2: paths.add_path(','.join(s), frequency=1)",
    "out = {'k': [], 'context': [], 'following': [], 'count': []}",
    "for k in spec['orders']:",
    "    network = pp.HigherOrderNetwork(paths, k=k)",
    "    for (u, v), attributes in network.edges.items():",
    "        weight = attributes['weight']",
    "        out['k'].append(k); out['context'].append(' -> '.join(u.split(',')))",
    "        out['following'].append(v.split(',')[-1]); out['count'].append(float(weight[0] + weight[1]))",
    "json.dump(out, open(sys.argv[2], 'w'))")
  as.data.frame(ref_python_run(script, list(sequences = sequences, orders = as.list(orders))))
}
