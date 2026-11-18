library(igraph)
memberships <- read.csv("../data/memberships.csv")
bipartite <- graph_from_data_frame(memberships, directed = FALSE)
V(bipartite)$type <- V(bipartite)$name %in% memberships$edge
projection <- bipartite_projection(bipartite, multiplicity = TRUE, which = "true")
projection <- delete_edges(projection, E(projection)[weight < 2])
pairs <- as_data_frame(projection, what = "edges")
result <- data.frame(pair = paste(pmin(pairs$from, pairs$to), pmax(pairs$from, pairs$to), sep = " -- "),
                     weight = pairs$weight)
result
