library(igraph)
memberships <- read.csv("../data/memberships.csv")
bipartite <- graph_from_data_frame(memberships, directed = FALSE)
V(bipartite)$type <- V(bipartite)$name %in% memberships$edge
E(bipartite)$weight <- 1
edges <- unique(memberships$edge)
result <- data.frame(
  edge = edges,
  size = unname(degree(bipartite, edges)),
  weight = unname(strength(bipartite, edges)),
  n_incident_edges = ego_size(bipartite, order = 2, nodes = edges, mindist = 2),
  n_neighbors = ego_size(bipartite, order = 3, nodes = edges, mindist = 3))
result
