memberships <- read.csv("../data/memberships.csv")
incidence <- unclass(table(memberships$node, memberships$edge))
shared <- crossprod(incidence)
pairs <- which(upper.tri(shared) & shared >= 2, arr.ind = TRUE)
result <- data.frame(pair = paste(rownames(shared)[pairs[, "row"]], colnames(shared)[pairs[, "col"]], sep = " -- "),
                     weight = shared[pairs])
result
