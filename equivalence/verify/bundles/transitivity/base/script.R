memberships <- read.csv("../data/memberships.csv")
incidence <- unclass(table(memberships$node, memberships$edge))
adjacency <- (tcrossprod(incidence) > 0) * 1
diag(adjacency) <- 0
neighbours <- rowSums(adjacency)
result <- data.frame(node = rownames(adjacency),
                     projection = diag(adjacency %*% adjacency %*% adjacency) / (neighbours * (neighbours - 1)))
result
