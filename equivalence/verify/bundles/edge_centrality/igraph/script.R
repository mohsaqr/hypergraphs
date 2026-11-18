library(igraph)
memberships <- read.csv("../data/memberships.csv")
bipartite <- graph_from_data_frame(memberships, directed = FALSE)
V(bipartite)$type <- V(bipartite)$name %in% memberships$edge
line <- bipartite_projection(bipartite, multiplicity = FALSE, which = "true")
is_connected(line)
result <- data.frame(edge = V(line)$name,
                     betweenness = betweenness(line, normalized = TRUE),
                     closeness = closeness(line, normalized = TRUE))
result
