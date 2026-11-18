memberships <- read.csv("../data/memberships.csv")
incidence <- unclass(table(memberships$node, memberships$edge))
together <- tcrossprod(incidence) > 0
diag(together) <- FALSE
result <- data.frame(
  edge = colnames(incidence),
  size = colSums(incidence),
  weight = colSums(incidence),
  n_incident_edges = rowSums(crossprod(incidence) > 0) - 1,
  n_neighbors = colSums(together %*% incidence > 0 & incidence == 0))
result
