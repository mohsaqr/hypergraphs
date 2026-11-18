library(hypergraphs)
tribunals <- read.csv("../data/icsid_tribunals.csv", encoding = "UTF-8")
hg <- group_hypergraph(tribunals, node = "arbitrator", hyperedge = "case")
hg
result <- hg_measures(hg, what = "nodes")
head(result, 10)
