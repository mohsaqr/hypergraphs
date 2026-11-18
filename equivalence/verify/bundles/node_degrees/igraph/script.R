library(igraph)
tribunals <- read.csv("../data/icsid_tribunals.csv", encoding = "UTF-8")
bipartite <- graph_from_data_frame(tribunals, directed = FALSE)
arbitrators <- unique(tribunals$arbitrator)
cases <- neighborhood(bipartite, order = 1, nodes = arbitrators, mindist = 1)
result <- data.frame(
  node = arbitrators,
  hyperdegree = unname(degree(bipartite, arbitrators)),
  strength = vapply(cases, function(own) sum(degree(bipartite, own)), numeric(1)),
  max_edge_size = vapply(cases, function(own) max(degree(bipartite, own)), numeric(1)),
  n_neighbors = ego_size(bipartite, order = 2, nodes = arbitrators, mindist = 2))
head(result, 10)
