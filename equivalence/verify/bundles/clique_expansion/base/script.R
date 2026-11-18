memberships <- read.csv("../data/memberships.csv")
incidence <- unclass(table(memberships$node, memberships$edge))
shared <- tcrossprod(incidence)
shared
pairs <- which(upper.tri(shared) & shared > 0, arr.ind = TRUE)
result <- data.frame(pair = paste(rownames(shared)[pairs[, "row"]], colnames(shared)[pairs[, "col"]], sep = " -- "),
                     weight = shared[pairs])
result
