# simplicial(), hg_betti(), hg_euler(), hg_degree(), hg_qanalysis(): the clique complex of a network and its
# measures. Simplices are the cliques; Betti numbers the ranks of the homology, b_k = dim C_k - rank d_k -
# rank d_(k+1) from the boundary matrices; the Euler characteristic the alternating sum of the f-vector; the
# simplicial degree of a node the number of k-simplices that contain it; q-analysis (Atkin 1974) the number of
# q-connected components of the maximal simplices of dimension at least q and the largest q of each node.
# References: igraph's cliques, and the definitions in base R.

#' A network of 9 nodes: a 4-clique, two triangles around a hole of four nodes, and a pendant node
network <- function() {
  nodes <- letters[1:9]
  adjacency <- matrix(0, 9L, 9L, dimnames = list(nodes, nodes))
  links <- rbind(c("a", "b"), c("a", "c"), c("a", "d"), c("b", "c"), c("b", "d"), c("c", "d"), c("d", "e"), c("e", "f"),
                 c("f", "g"), c("g", "d"), c("e", "g"), c("g", "h"), c("h", "i"), c("i", "f"))
  adjacency[links] <- 1
  adjacency[links[, 2:1]] <- 1
  adjacency
}

#' Every clique of the network, by igraph, as sorted character vectors
igraph_simplices <- function(adjacency) {
  graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected")
  lapply(igraph::cliques(graph), \(clique) sort(igraph::V(graph)$name[clique]))
}

#' The rank of the boundary map from k-simplices to (k - 1)-simplices, over the reals
boundary_rank <- function(simplices, k) {
  upper <- Filter(\(s) length(s) == k + 1L, simplices)
  lower <- Filter(\(s) length(s) == k, simplices)
  if (!length(upper) || !length(lower) || k == 0L) return(0)
  keys <- vapply(lower, paste, character(1), collapse = ",")
  boundary <- vapply(upper, \(s) {
    column <- numeric(length(lower))
    faces <- vapply(seq_along(s), \(i) paste(s[-i], collapse = ","), character(1))
    column[match(faces, keys)] <- (-1)^(seq_along(s) - 1L)
    column
  }, numeric(length(lower)))
  qr(boundary)$rank
}

test_that("simplices, f-vector, Betti numbers, Euler characteristic and degrees equal igraph's cliques and the definitions", {
  skip_without_reference("igraph")
  adjacency <- network()
  simplices <- igraph_simplices(adjacency)
  complex <- simplicial(adjacency)
  expect_reference(transform(hg_get(complex, what = "simplices"), present = TRUE),
                   data.frame(members = vapply(simplices, paste, character(1), collapse = ", "), present = TRUE),
                   reference = "igraph", covers = c("simplicial()", "hg_get.simplicial_complex()", "hg_get.hypergraphs_simplicial()"),
                   quantity = "simplices", tolerance = 0, key = "members")
  sizes <- lengths(simplices)
  f_vector <- as.vector(table(sizes))
  expect_reference(hg_get(complex, what = "f_vector"), data.frame(dimension = seq_along(f_vector) - 1L, count = f_vector),
                   reference = "igraph", covers = c("simplicial()", "hg_get.simplicial_complex()"), quantity = "f-vector",
                   tolerance = 0, key = "dimension")
  dims <- seq_along(f_vector) - 1L
  ranks <- vapply(c(dims, max(dims) + 1L), \(k) boundary_rank(simplices, k), numeric(1))
  betti <- f_vector - utils::head(ranks, -1L) - utils::tail(ranks, -1L)
  expect_reference(hg_betti(complex), stats::setNames(betti, paste0("b", dims)),
                   reference = "independent: ranks of the boundary matrices of igraph's cliques, base R", covers = "hg_betti()",
                   quantity = "Betti numbers", tolerance = 0)
  expect_reference(hg_euler(complex), sum((-1)^dims * f_vector),
                   reference = "independent: alternating sum of the f-vector of igraph's cliques", covers = "hg_euler()",
                   quantity = "Euler characteristic", tolerance = 0)
  nodes <- rownames(adjacency)
  degree <- vapply(dims, \(k) vapply(nodes, \(node) sum(vapply(simplices, \(s) length(s) == k + 1L && node %in% s, logical(1))), numeric(1)),
                   numeric(length(nodes)))
  colnames(degree) <- paste0("d", dims)
  # `total` counts the simplices of dimension one and above (documented in ?hg_degree)
  expect_reference(hg_degree(complex), data.frame(node = nodes, degree, total = rowSums(degree[, -1L, drop = FALSE])),
                   reference = "independent: simplices of igraph's cliques containing each node", covers = "hg_degree()",
                   quantity = "simplicial degree", tolerance = 0, key = "node")
})

test_that("q-connectivity components and the largest q of each node equal Atkin's definition on igraph's maximal cliques", {
  skip_without_reference("igraph")
  adjacency <- network()
  graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected")
  maximal <- lapply(igraph::max_cliques(graph), \(clique) igraph::V(graph)$name[clique])
  top <- max(lengths(maximal)) - 1L
  components <- vapply(0:top, \(q) {
    kept <- Filter(\(s) length(s) - 1L >= q, maximal)
    shared <- outer(seq_along(kept), seq_along(kept), Vectorize(\(i, j) length(intersect(kept[[i]], kept[[j]])) >= q + 1L))
    igraph::components(igraph::graph_from_adjacency_matrix(shared * 1, mode = "undirected"))$no
  }, numeric(1))
  analysis <- hg_qanalysis(simplicial(adjacency))
  expect_reference(hg_get(analysis, what = "q_levels"), data.frame(q = 0:top, components = components),
                   reference = "independent: q-connectivity of igraph's maximal cliques (Atkin 1974)",
                   covers = c("hg_qanalysis()", "hg_get.q_analysis()"), quantity = "components at each q", tolerance = 0, key = "q")
  nodes <- rownames(adjacency)
  expect_reference(hg_get(analysis, what = "nodes"),
                   data.frame(node = nodes, max_q = vapply(nodes, \(node) max(lengths(Filter(\(s) node %in% s, maximal))) - 1, numeric(1))),
                   reference = "independent: q-connectivity of igraph's maximal cliques (Atkin 1974)",
                   covers = c("hg_qanalysis()", "hg_get.q_analysis()"), quantity = "largest q of each node", tolerance = 0, key = "node")
  expect_reference(hg_get(summary(simplicial(adjacency)), what = "f_vector"),
                   data.frame(dimension = 0:top, count = as.vector(table(lengths(igraph_simplices(adjacency))))),
                   reference = "igraph", covers = c("summary.hypergraphs_simplicial()", "hg_get.hypergraphs_summary()"),
                   quantity = "f-vector through summary()", tolerance = 0, key = "dimension")
})
