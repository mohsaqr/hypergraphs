# hg_representations(): the four representations of a hypergraph (binary graph, multigraph, binary hypergraph,
# multi-hypergraph) with their edge counts and mean and median degrees, against Table 2 of the paper that defines
# them on the bundled data: Coupette, Hartung and Katz (2024), Legal hypergraphs, Phil. Trans. R. Soc. A 382:
# 20230141. The table prints the means to one decimal: they are compared within 0.05, the half-unit of the last
# decimal shown; counts and medians exactly.

table_2 <- "published: Coupette, Hartung and Katz (2024), Legal hypergraphs, Phil. Trans. R. Soc. A 382: 20230141, Table 2"

test_that("the ICSID tribunals reproduce Table 2b", {
  tribunals <- temporal_hypergraph(icsid_tribunals, node = "arbitrator", hyperedge = "case", start = "constituted",
                                   end = "concluded", observation_end = as.Date("2023-06-15"))
  ours <- hg_representations(hg_snapshot(tribunals, mode = "cumulative"))
  paper <- data.frame(representation = c("bg", "mg", "bh", "mh"), n_edges = c(1869, 2226, 722, 742), median_degree = c(4, 4, 2, 2))
  expect_reference(ours, paper, reference = table_2, covers = c("hg_representations()", "temporal_hypergraph()", "hg_snapshot()"),
                   quantity = "ICSID, edges and median degree", tolerance = 0, key = "representation")
  # printed to one decimal in the paper: half a unit of the last decimal
  expect_reference(ours, data.frame(representation = c("bg", "mg", "bh", "mh"), mean_degree = c(8.5, 10.1, 4.9, 5.0)),
                   reference = table_2, covers = c("hg_representations()", "temporal_hypergraph()", "hg_snapshot()"),
                   quantity = "ICSID, mean degree (one decimal)", tolerance = 0.05, key = "representation")
  expect_reference(unlist(subset(summary(tribunals), select = n_hyperedges)), c(n_hyperedges = 742), reference = table_2,
                   covers = "summary.net_temporal_hypergraph()", quantity = "ICSID, hyperedges (cases)", tolerance = 0)
})

test_that("the GFCC citation blocks reproduce Table 2a, the multigraph and binary hypergraph short by the repeated citations", {
  blocks <- temporal_hypergraph(gfcc_citations, node = "cited", hyperedge = "block", time = "date_citing", nodes = gfcc_decisions,
                                sparse = TRUE)
  ours <- hg_representations(hg_snapshot(blocks, mode = "cumulative"), graph = "citation", edge_source = "citing")
  covers <- c("hg_representations()", "temporal_hypergraph()", "hg_snapshot()")
  expect_reference(subset(ours, representation %in% c("bg", "mh")), data.frame(representation = c("bg", "mh"), n_edges = c(39428, 46257)),
                   reference = table_2, covers = covers, quantity = "GFCC, edges", tolerance = 0, key = "representation")
  expect_known_difference(subset(ours, representation %in% c("mg", "bh")), data.frame(representation = c("mg", "bh"), n_edges = c(77284, 13541)),
                          reference = table_2, covers = covers, quantity = "GFCC, edges", tolerance = 0, key = "representation",
                          reason = "Five citations cite the same decision twice inside one block. The paper's multigraph keeps each repeat as an edge, and its hypergraph keeps a block that differs from another only by such a repeat as a distinct member multiset; an incidence matrix holds a membership once, so hypergraphs reports 5 fewer multigraph edges and 2 fewer binary hyperedges (counted from the raw rows in the standing battery).")
  expect_reference(ours, data.frame(representation = c("bg", "mg", "bh", "mh"), mean_degree = c(10.9, 21.4, 10.2, 21.4)),
                   reference = table_2, covers = covers, quantity = "GFCC, mean degree (one decimal)", tolerance = 0.05,
                   key = "representation")
  expect_reference(subset(ours, type == "graph"), data.frame(representation = c("bg", "mg"), median_out_degree = c(6, 9), median_in_degree = c(7, 10)),
                   reference = table_2, covers = covers, quantity = "GFCC, median out- and in-degree", tolerance = 0, key = "representation")
  expect_reference(subset(ours, type == "hypergraph"), data.frame(representation = c("bh", "mh"), median_degree = c(6, 10)),
                   reference = table_2, covers = covers, quantity = "GFCC, median degree", tolerance = 0, key = "representation")
})
