memberships <- read.csv("../data/memberships.csv")
incidence <- unclass(table(memberships$node, memberships$edge))
degree_root <- diag(1 / sqrt(rowSums(incidence)))
size_inverse <- diag(1 / colSums(incidence))
laplacian <- diag(nrow(incidence)) - degree_root %*% incidence %*% size_inverse %*% t(incidence) %*% degree_root
dimnames(laplacian) <- list(rownames(incidence), rownames(incidence))
laplacian
result <- transform(as.data.frame(as.table(laplacian)), quantity = paste(Var1, Var2, sep = ", "),
                    value = Freq, kind = "laplacian")
head(result)
