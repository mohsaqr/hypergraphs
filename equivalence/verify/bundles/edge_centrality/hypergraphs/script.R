library(hypergraphs)
memberships <- read.csv("../data/memberships.csv")
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
hg_edge_centrality(hg, s = 1, measure = c("betweenness", "closeness"), normalized = TRUE)
result <- transform(hg_edge_centrality(hg, s = 1, measure = c("betweenness", "closeness"), normalized = TRUE),
                    quantity = paste(measure, edge, sep = ": "), kind = measure)
result
