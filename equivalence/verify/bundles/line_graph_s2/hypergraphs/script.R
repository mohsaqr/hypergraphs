library(hypergraphs)
memberships <- read.csv("../data/memberships.csv")
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
hg_line_graph(hg, s = 2, what = "edges")
result <- transform(hg_line_graph(hg, s = 2, what = "edges"), pair = paste(from, to, sep = " -- ")) |>
  subset(select = c(pair, weight))
result
