library(igraph)
memberships <- read.csv("../data/memberships.csv")
bipartite <- graph_from_data_frame(memberships, directed = FALSE)
V(bipartite)$type <- V(bipartite)$name %in% memberships$edge
expansion <- bipartite_projection(bipartite, multiplicity = FALSE, which = "false")
result <- data.frame(node = V(expansion)$name,
                     projection = transitivity(expansion, type = "local"))
result
