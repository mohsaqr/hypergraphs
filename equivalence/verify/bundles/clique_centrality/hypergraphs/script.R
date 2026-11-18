library(hypergraphs)
memberships <- read.csv("../data/memberships.csv")
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
result <- hg_centrality(hg, type = c("clique", "pagerank", "katz"), damping = 0.85, alpha = 0.05,
                        tol = 1e-12, max_iter = 10000L)
result
hg_pagerank(hg, damping = 0.85)
