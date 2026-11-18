library(igraph)
memberships <- read.csv("../data/memberships.csv")
bipartite <- graph_from_data_frame(memberships, directed = FALSE)
V(bipartite)$type <- V(bipartite)$name %in% memberships$edge
expansion <- bipartite_projection(bipartite, multiplicity = TRUE, which = "false")
eigenvector <- eigen_centrality(expansion, weights = E(expansion)$weight)$vector
incidence <- as_biadjacency_matrix(bipartite)
walk <- incidence %*% (t(incidence) / colSums(incidence))
walk_graph <- graph_from_adjacency_matrix(walk, mode = "directed", weighted = TRUE, diag = TRUE)
result <- data.frame(
  node = V(expansion)$name,
  clique = eigenvector / sqrt(sum(eigenvector^2)),
  pagerank = page_rank(walk_graph, damping = 0.85)$vector[V(expansion)$name],
  katz = alpha_centrality(expansion, alpha = 0.05, exo = 1, weights = E(expansion)$weight) - 1)
result
