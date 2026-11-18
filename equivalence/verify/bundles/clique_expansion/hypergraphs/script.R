library(hypergraphs)
memberships <- read.csv("../data/memberships.csv")
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
network <- pairwise_network(hg, type = "clique")
hg_get(network, what = "edges")
result <- transform(hg_get(network, what = "edges"), pair = paste(from, to, sep = " -- ")) |>
  subset(select = c(pair, weight))
result
