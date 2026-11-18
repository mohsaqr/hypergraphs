memberships <- read.csv("../data/memberships.csv")
incidence <- unclass(table(memberships$node, memberships$edge))
nodes <- nrow(incidence)
adjacency <- tcrossprod(incidence)
diag(adjacency) <- 0
leading <- abs(eigen(adjacency, symmetric = TRUE)$vectors[, 1])
walk <- incidence %*% (t(incidence) / colSums(incidence))
transition <- walk / rowSums(walk)
result <- data.frame(
  node = rownames(incidence),
  clique = leading / sqrt(sum(leading^2)),
  pagerank = solve(diag(nodes) - 0.85 * t(transition), rep((1 - 0.85) / nodes, nodes)),
  katz = solve(diag(nodes) - 0.05 * adjacency, rep(1, nodes)) - 1)
result
