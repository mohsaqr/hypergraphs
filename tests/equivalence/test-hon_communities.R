# hg_communities() on a higher-order network: Infomap's map equation for a memory (state-node) network (Rosvall et
# al. 2014), each higher-order node a state node of its last first-order state, with teleportation 0.15. Given a
# partition, the codelength and its index and module parts, the one-level codelength and the flow of every state
# node. Reference: Python infomap evaluating the same partition of the same state network (no search).

script <- c(
  "import sys, json, infomap",
  "spec = json.load(open(sys.argv[1]))",
  "nodes = spec['node']",
  "index = {n: i + 1 for i, n in enumerate(nodes)}",
  "physical = {s: i + 1 for i, s in enumerate(sorted(set(spec['state'])))}",
  "im = infomap.Infomap(directed=True, silent=True, two_level=True, teleportation_probability=spec['teleportation'],",
  "                     no_infomap=True, num_trials=1)",
  "for n, s in zip(nodes, spec['state']): im.add_state_node(index[n], physical[s])",
  "for a, b, w in zip(spec['from'], spec['to'], spec['weight']): im.add_link(index[a], index[b], float(w))",
  "im.run(initial_partition={index[n]: int(m) for n, m in zip(nodes, spec['module'])})",
  "flow = {n.state_id: n.flow for n in im.tree if n.is_leaf}",
  "json.dump({'codelength': im.codelength, 'index_codelength': im.index_codelength,",
  "           'module_codelength': im.module_codelength, 'one_level_codelength': im.one_level_codelength,",
  "           'node': nodes, 'flow': [float(flow[index[n]]) for n in nodes]}, open(sys.argv[2], 'w'))")

test_that("the codelength and state flows of a given partition equal Infomap's", {
  skip_without_python("infomap")
  sequences <- ref_second_order(n_sequences = 40L, length = 20L, seed = 9L)
  network <- hon(sequences, max_order = 2L, min_freq = 1L)
  rules <- hg_get(network, what = "rules")
  paths <- strsplit(rules$path, " -> ", fixed = TRUE)
  target <- mapply(\(path, order) paste(utils::tail(path, order), collapse = " -> "), paths, rules$to_order)
  states <- hg_get(hg_communities(network, trials = 2L, seed = 1L), what = "states")
  partition <- stats::setNames(c(1, 1, 2, 2, 3, 3)[seq_len(nrow(states))], states$node)
  found <- hg_communities(network, partition = partition)
  theirs <- ref_python_run(script, list(node = states$node, state = states$state, from = rules$from, to = target,
                                        weight = rules$probability, module = unname(partition), teleportation = 0.15))
  version <- ref_python_version("infomap")
  # Sums of flow-weighted logarithms; the flows are a power iteration on both sides.
  expect_reference(subset(hg_get(found, what = "codelength"), model == "memory",
                          select = c(model, codelength, index_codelength, module_codelength, one_level_codelength)),
                   data.frame(model = "memory", codelength = theirs$codelength, index_codelength = theirs$index_codelength,
                              module_codelength = theirs$module_codelength, one_level_codelength = theirs$one_level_codelength),
                   reference = "infomap", covers = c("hg_communities()", "hg_communities.net_hon()", "hg_get.hypergraphs_memory_communities()"),
                   quantity = "given partition, codelengths", tolerance = 1e-8, key = "model", note = version)
  expect_reference(subset(hg_get(found, what = "states"), select = c(node, flow)), data.frame(node = theirs$node, flow = theirs$flow),
                   reference = "infomap", covers = c("hg_communities.net_hon()", "hg_get.hypergraphs_memory_communities()"),
                   quantity = "given partition, state flows", tolerance = 1e-8, key = "node", note = version)
  expect_reference(subset(hg_get(summary(found), what = "states"), select = c(node, flow)), data.frame(node = theirs$node, flow = theirs$flow),
                   reference = "infomap", covers = c("summary.hypergraphs_memory_communities()", "hg_get.hypergraphs_summary()"),
                   quantity = "given partition, state flows through summary()", tolerance = 1e-8, key = "node", note = version)
})
