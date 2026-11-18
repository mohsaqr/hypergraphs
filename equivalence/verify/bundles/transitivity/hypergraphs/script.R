library(hypergraphs)
memberships <- read.csv("../data/memberships.csv")
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
result <- hg_transitivity(hg, type = "projection", what = "nodes")
result
hg_transitivity(hg, type = "projection", what = "summary")
