# Reference inputs and computations shared by the equivalence tests, written in base R. Nothing here calls
# hypergraphs: each function builds the input both sides read, or is one side of a comparison.

#' The small hypergraph most tests use, as a membership table
#'
#' Seven nodes and seven hyperedges of sizes 2, 3 and 4, with repeated pairs (a and b share two hyperedges) so
#' that multiplicities differ from 0/1.
#' @return A data frame `node`, `edge`, one row per membership.
ref_memberships <- function() {
  sets <- list(e1 = c("a", "b", "c"), e2 = c("a", "b", "d"), e3 = c("a", "c", "e"),
               e4 = c("b", "d", "e"), e5 = c("e", "f"), e6 = c("c", "d", "f", "g"),
               e7 = c("a", "g"))
  data.frame(node = unlist(sets, use.names = FALSE), edge = rep(names(sets), lengths(sets)),
             stringsAsFactors = FALSE)
}

#' The incidence matrix of a membership table: nodes in rows, hyperedges in columns
#'
#' @param memberships A data frame with the node and hyperedge columns.
#' @param node,edge The column names.
#' @return An integer matrix with dimnames, nodes and hyperedges in sorted order.
ref_incidence <- function(memberships, node = "node", edge = "edge") {
  unclass(table(memberships[[node]], memberships[[edge]]))
}

#' The node pairs of a symmetric matrix that pass a test, as one row per unordered pair
#'
#' @param shared A square symmetric matrix with dimnames.
#' @param keep A function of the matrix returning a logical matrix of the pairs to keep.
#' @param name The name of the value column.
#' @return A data frame `from`, `to`, `<name>`, `from` before `to` in the matrix's order.
ref_pairs_upper <- function(shared, keep = \(m) m > 0, name = "weight") {
  pairs <- which(upper.tri(shared) & keep(shared), arr.ind = TRUE)
  out <- data.frame(from = rownames(shared)[pairs[, "row"]], to = colnames(shared)[pairs[, "col"]],
                    value = shared[pairs], stringsAsFactors = FALSE)
  names(out)[[3L]] <- name
  out
}

#' Write unordered pairs with the smaller name first
#'
#' Names are compared by their bytes (radix order), not by the locale's collation, so a pair of accented or
#' mixed-case names is written the same way on every machine (`pmin()` on strings follows the locale).
#'
#' @param from,to Character vectors.
#' @return A data frame `from`, `to` with `from` before `to` in radix (C-locale) order.
ref_ordered_pairs <- function(from, to) {
  levels <- sort(unique(c(from, to)), method = "radix")
  first <- match(from, levels) <= match(to, levels)
  data.frame(from = ifelse(first, from, to), to = ifelse(first, to, from), stringsAsFactors = FALSE)
}

#' A square matrix as one row per cell
#'
#' @param cells A matrix with dimnames.
#' @param name The name of the value column.
#' @param row,column The names of the two key columns.
#' @return A data frame `<row>`, `<column>`, `<name>`.
ref_long <- function(cells, name = "value", row = "from", column = "to") {
  cells <- unclass(as.matrix(cells))
  out <- data.frame(rownames(cells)[as.vector(row(cells))], colnames(cells)[as.vector(col(cells))],
                    as.vector(cells), stringsAsFactors = FALSE)
  names(out) <- c(row, column, name)
  out
}
