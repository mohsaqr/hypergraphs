tribunals <- read.csv("../data/icsid_tribunals.csv", encoding = "UTF-8")
incidence <- unclass(table(tribunals$arbitrator, tribunals$case))
sizes <- colSums(incidence)
together <- tcrossprod(incidence) > 0
diag(together) <- FALSE
result <- data.frame(
  node = rownames(incidence),
  hyperdegree = rowSums(incidence),
  strength = as.vector(incidence %*% sizes),
  max_edge_size = apply(incidence, 1, function(row) max(sizes[row > 0])),
  n_neighbors = rowSums(together))
head(result, 10)
