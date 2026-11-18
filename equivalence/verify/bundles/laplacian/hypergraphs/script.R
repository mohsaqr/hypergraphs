library(hypergraphs)
memberships <- read.csv("../data/memberships.csv")
hg <- group_hypergraph(memberships, node = "node", hyperedge = "edge")
laplacian <- hg_laplacian(hg, type = "zhou")
laplacian
result <- transform(as.data.frame(as.table(laplacian)), quantity = paste(Var1, Var2, sep = ", "),
                    value = Freq, kind = "laplacian")
head(result)
