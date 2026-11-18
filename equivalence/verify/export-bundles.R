# Exports hypergraphs comparison cases as bundles for Prova Verify (format:
# plan/24-verify.md in the prova repository), with Prova's R writer.
#
# Every side of a case is an R script, stored in the bundle and run there in its
# own R process: its console is the side's raw output and the table it leaves in
# `result` holds the side's values. The references are igraph (the package the
# hypergraphs tests already use for the clique expansion, PageRank, Katz and
# local transitivity) and a short base R computation written from the
# definition on the incidence matrix.
#
# Run from the project root:
#   Rscript equivalence/verify/export-bundles.R [case ...]
source(Sys.getenv("PROVA_VERIFY_R", file.path("..", "prova", "reporters", "verify.R")))

.verify_base <- file.path("equivalence", "verify")
.verify_generator <- "hypergraphs equivalence/verify/export-bundles.R"

# Tolerances, fixed before any case was run.
# Counts (degrees, sizes, intersection sizes) are integers on every side: they
# must be equal, so the tolerance is 0.
.tol_count <- 0
# Closed-form floating point quantities (shortest-path centralities, local
# transitivity, Laplacian entries) differ only by rounding in a handful of
# operations on numbers of order 1: 1e-12.
.tol_closed <- 1e-12
# Eigenvector, PageRank and Katz scores are iterative on the package's side
# (power iteration, asked here for tol = 1e-12) and on igraph's (ARPACK, PRPACK);
# the error of a power iteration stopped at 1e-12 is that divided by the
# spectral gap, so 1e-8 is the agreement that can be asked of it.
.tol_iterative <- 1e-8

#' The small hypergraph most cases use, written to a file
#'
#' Seven nodes and seven hyperedges of sizes 2, 3 and 4, with repeated pairs
#' (a and b share two hyperedges) so that multiplicities differ from 0/1.
#' @return The path of a membership table `memberships.csv` (`node`, `edge`).
.verify_memberships <- function() {
  sets <- list(e1 = c("a", "b", "c"), e2 = c("a", "b", "d"), e3 = c("a", "c", "e"),
               e4 = c("b", "d", "e"), e5 = c("e", "f"), e6 = c("c", "d", "f", "g"),
               e7 = c("a", "g"))
  memberships <- data.frame(node = unlist(sets, use.names = FALSE),
                            edge = rep(names(sets), lengths(sets)))
  file <- file.path(tempdir(), "memberships.csv")
  utils::write.csv(memberships, file, row.names = FALSE)
  file
}

#' The bundled ICSID tribunals as a membership table, written to a file
#'
#' One row per seat: 441 arbitrators (nodes) on 742 cases (hyperedges).
#' @return The path of `icsid_tribunals.csv` (`arbitrator`, `case`).
.verify_tribunals <- function() {
  environment <- new.env()
  load(file.path("data", "icsid_tribunals.rda"), envir = environment)
  tribunals <- subset(environment$icsid_tribunals, select = c(arbitrator, case))
  stopifnot("an arbitrator sits once on a case" = anyDuplicated(tribunals) == 0L)
  file <- file.path(tempdir(), "icsid_tribunals.csv")
  utils::write.csv(tribunals, file, row.names = FALSE, fileEncoding = "UTF-8")
  file
}

.read_memberships <- "memberships <- read.csv(\"../data/memberships.csv\")"
.read_tribunals <- "tribunals <- read.csv(\"../data/icsid_tribunals.csv\", encoding = \"UTF-8\")"
.build_hg <- "hg <- group_hypergraph(memberships, node = \"node\", hyperedge = \"edge\")"
.build_incidence <- "incidence <- unclass(table(memberships$node, memberships$edge))"
.build_bipartite <- c(
  "bipartite <- graph_from_data_frame(memberships, directed = FALSE)",
  "V(bipartite)$type <- V(bipartite)$name %in% memberships$edge")

.hypergraphs <- function(script, key) {
  # The version is read from the source tree that `load = "."` runs, not from an
  # installed copy (which may be absent or older).
  verify_side("hypergraphs", package = "hypergraphs", primary = TRUE, load = ".", key = key,
              version = unname(read.dcf("DESCRIPTION", fields = "Version")[1L, "Version"]),
              script = c("library(hypergraphs)", script))
}
.igraph <- function(script, key) {
  verify_side("igraph", package = "igraph", authority = "open-source", key = key,
              script = c("library(igraph)", script))
}
.base <- function(script, key) {
  verify_side("base", name = "Base R, from the definition", version = as.character(getRversion()),
              authority = "independent", key = key, script = script)
}

#' Every case: its words, data, tolerance and the script of each side
#'
#' A script's paths are relative to its own folder inside the bundle. A script
#' leaves a data frame `result`: wide with one key column, or long with
#' `quantity` and `value`.
#' @return A named list of cases.
.verify_cases <- function() {
  list(
    node_degrees = list(
      title = "Node degrees of the ICSID tribunal hypergraph",
      description = "Hyperdegree (number of hyperedges), strength (sum of the sizes of a node's hyperedges), largest hyperedge and number of distinct neighbours of the 441 arbitrators in the bundled icsid_tribunals data, where a case is a hyperedge. igraph computes them on the bipartite incidence graph; base R from the incidence matrix.",
      covers = c("group_hypergraph()", "hg_measures()"),
      data = c("ICSID tribunals (arbitrator, case)" = "tribunals"), tolerance = .tol_count,
      sides = list(
        .hypergraphs(key = "node", c(
          .read_tribunals,
          "hg <- group_hypergraph(tribunals, node = \"arbitrator\", hyperedge = \"case\")",
          "hg",
          "result <- hg_measures(hg, what = \"nodes\")",
          "head(result, 10)")),
        .igraph(key = "node", c(
          .read_tribunals,
          "bipartite <- graph_from_data_frame(tribunals, directed = FALSE)",
          "arbitrators <- unique(tribunals$arbitrator)",
          "cases <- neighborhood(bipartite, order = 1, nodes = arbitrators, mindist = 1)",
          "result <- data.frame(",
          "  node = arbitrators,",
          "  hyperdegree = unname(degree(bipartite, arbitrators)),",
          "  strength = vapply(cases, function(own) sum(degree(bipartite, own)), numeric(1)),",
          "  max_edge_size = vapply(cases, function(own) max(degree(bipartite, own)), numeric(1)),",
          "  n_neighbors = ego_size(bipartite, order = 2, nodes = arbitrators, mindist = 2))",
          "head(result, 10)")),
        .base(key = "node", c(
          .read_tribunals,
          "incidence <- unclass(table(tribunals$arbitrator, tribunals$case))",
          "sizes <- colSums(incidence)",
          "together <- tcrossprod(incidence) > 0",
          "diag(together) <- FALSE",
          "result <- data.frame(",
          "  node = rownames(incidence),",
          "  hyperdegree = rowSums(incidence),",
          "  strength = as.vector(incidence %*% sizes),",
          "  max_edge_size = apply(incidence, 1, function(row) max(sizes[row > 0])),",
          "  n_neighbors = rowSums(together))",
          "head(result, 10)")))),
    edge_sizes = list(
      title = "Hyperedge sizes and overlaps",
      description = "Size, weight (the sum of a hyperedge's incidence weights, here its size), number of hyperedges it meets and number of outside nodes adjacent to its members, for a seven-edge hypergraph with sizes 2, 3 and 4. igraph computes them on the bipartite incidence graph; base R from the incidence matrix.",
      covers = c("group_hypergraph()", "hg_edges()"),
      data = c(Memberships = "memberships"), tolerance = .tol_count,
      sides = list(
        .hypergraphs(key = "edge", c(
          .read_memberships, .build_hg,
          "hg",
          "result <- hg_edges(hg, what = \"edges\")",
          "result")),
        .igraph(key = "edge", c(
          .read_memberships, .build_bipartite,
          "E(bipartite)$weight <- 1",
          "edges <- unique(memberships$edge)",
          "result <- data.frame(",
          "  edge = edges,",
          "  size = unname(degree(bipartite, edges)),",
          "  weight = unname(strength(bipartite, edges)),",
          "  n_incident_edges = ego_size(bipartite, order = 2, nodes = edges, mindist = 2),",
          "  n_neighbors = ego_size(bipartite, order = 3, nodes = edges, mindist = 3))",
          "result")),
        .base(key = "edge", c(
          .read_memberships, .build_incidence,
          "together <- tcrossprod(incidence) > 0",
          "diag(together) <- FALSE",
          "result <- data.frame(",
          "  edge = colnames(incidence),",
          "  size = colSums(incidence),",
          "  weight = colSums(incidence),",
          "  n_incident_edges = rowSums(crossprod(incidence) > 0) - 1,",
          "  n_neighbors = colSums(together %*% incidence > 0 & incidence == 0))",
          "result")))),
    clique_expansion = list(
      title = "Clique (two-section) expansion",
      description = "The weighted pairwise network of a hypergraph: two nodes are joined when they share a hyperedge, with the number of shared hyperedges as weight. igraph gives it as the bipartite projection of the incidence graph with multiplicities; base R as the off-diagonal of H H'.",
      covers = c("pairwise_network()", "hg_get()"),
      data = c(Memberships = "memberships"), tolerance = .tol_count,
      sides = list(
        .hypergraphs(key = "pair", c(
          .read_memberships, .build_hg,
          "network <- pairwise_network(hg, type = \"clique\")",
          "hg_get(network, what = \"edges\")",
          "result <- transform(hg_get(network, what = \"edges\"), pair = paste(from, to, sep = \" -- \")) |>",
          "  subset(select = c(pair, weight))",
          "result")),
        .igraph(key = "pair", c(
          .read_memberships, .build_bipartite,
          "projection <- bipartite_projection(bipartite, multiplicity = TRUE, which = \"false\")",
          "pairs <- as_data_frame(projection, what = \"edges\")",
          "result <- data.frame(pair = paste(pmin(pairs$from, pairs$to), pmax(pairs$from, pairs$to), sep = \" -- \"),",
          "                     weight = pairs$weight)",
          "result")),
        .base(key = "pair", c(
          .read_memberships, .build_incidence,
          "shared <- tcrossprod(incidence)",
          "shared",
          "pairs <- which(upper.tri(shared) & shared > 0, arr.ind = TRUE)",
          "result <- data.frame(pair = paste(rownames(shared)[pairs[, \"row\"]], colnames(shared)[pairs[, \"col\"]], sep = \" -- \"),",
          "                     weight = shared[pairs])",
          "result")))),
    clique_centrality = list(
      title = "Clique eigenvector, PageRank and Katz centrality",
      description = "Three node centralities of a hypergraph. Clique eigenvector centrality is the leading eigenvector of the weighted clique expansion, at unit Euclidean length. PageRank (damping 0.85) is the hypergraph random walk, which for a binary incidence is the walk on W[u, v] = sum over shared hyperedges of 1 / size, self-loops included (Chitra and Raphael 2019). Katz is (I - alpha A)^-1 1 - 1 on the weighted clique expansion, alpha = 0.05.",
      covers = c("hg_centrality()", "hg_pagerank()"),
      data = c(Memberships = "memberships"), tolerance = .tol_iterative,
      sides = list(
        .hypergraphs(key = "node", c(
          .read_memberships, .build_hg,
          "result <- hg_centrality(hg, type = c(\"clique\", \"pagerank\", \"katz\"), damping = 0.85, alpha = 0.05,",
          "                        tol = 1e-12, max_iter = 10000L)",
          "result",
          "hg_pagerank(hg, damping = 0.85)")),
        .igraph(key = "node", c(
          .read_memberships, .build_bipartite,
          "expansion <- bipartite_projection(bipartite, multiplicity = TRUE, which = \"false\")",
          "eigenvector <- eigen_centrality(expansion, weights = E(expansion)$weight)$vector",
          "incidence <- as_biadjacency_matrix(bipartite)",
          "walk <- incidence %*% (t(incidence) / colSums(incidence))",
          "walk_graph <- graph_from_adjacency_matrix(walk, mode = \"directed\", weighted = TRUE, diag = TRUE)",
          "result <- data.frame(",
          "  node = V(expansion)$name,",
          "  clique = eigenvector / sqrt(sum(eigenvector^2)),",
          "  pagerank = page_rank(walk_graph, damping = 0.85)$vector[V(expansion)$name],",
          "  katz = alpha_centrality(expansion, alpha = 0.05, exo = 1, weights = E(expansion)$weight) - 1)",
          "result")),
        .base(key = "node", c(
          .read_memberships, .build_incidence,
          "nodes <- nrow(incidence)",
          "adjacency <- tcrossprod(incidence)",
          "diag(adjacency) <- 0",
          "leading <- abs(eigen(adjacency, symmetric = TRUE)$vectors[, 1])",
          "walk <- incidence %*% (t(incidence) / colSums(incidence))",
          "transition <- walk / rowSums(walk)",
          "result <- data.frame(",
          "  node = rownames(incidence),",
          "  clique = leading / sqrt(sum(leading^2)),",
          "  pagerank = solve(diag(nodes) - 0.85 * t(transition), rep((1 - 0.85) / nodes, nodes)),",
          "  katz = solve(diag(nodes) - 0.05 * adjacency, rep(1, nodes)) - 1)",
          "result")))),
    line_graph = list(
      title = "Line graph (s = 1) of a hypergraph",
      description = "Hyperedges become vertices, joined when they share at least one node, with the number of shared nodes as weight. igraph gives it as the other bipartite projection of the incidence graph; base R as the off-diagonal of H'H.",
      covers = "hg_line_graph()",
      data = c(Memberships = "memberships"), tolerance = .tol_count,
      sides = list(
        .hypergraphs(key = "pair", c(
          .read_memberships, .build_hg,
          "hg_line_graph(hg, s = 1, what = \"edges\")",
          "result <- transform(hg_line_graph(hg, s = 1, what = \"edges\"), pair = paste(from, to, sep = \" -- \")) |>",
          "  subset(select = c(pair, weight))",
          "result")),
        .igraph(key = "pair", c(
          .read_memberships, .build_bipartite,
          "projection <- bipartite_projection(bipartite, multiplicity = TRUE, which = \"true\")",
          "pairs <- as_data_frame(projection, what = \"edges\")",
          "result <- data.frame(pair = paste(pmin(pairs$from, pairs$to), pmax(pairs$from, pairs$to), sep = \" -- \"),",
          "                     weight = pairs$weight)",
          "result")),
        .base(key = "pair", c(
          .read_memberships, .build_incidence,
          "shared <- crossprod(incidence)",
          "shared",
          "pairs <- which(upper.tri(shared) & shared >= 1, arr.ind = TRUE)",
          "result <- data.frame(pair = paste(rownames(shared)[pairs[, \"row\"]], colnames(shared)[pairs[, \"col\"]], sep = \" -- \"),",
          "                     weight = shared[pairs])",
          "result")))),
    line_graph_s2 = list(
      title = "s-line graph, s = 2",
      description = "The same hypergraph with hyperedges joined only when they share at least two nodes. igraph and base R keep the pairs of the s = 1 line graph whose weight is at least 2.",
      covers = "hg_line_graph()",
      variant = list(of = "line_graph", option = "s", value = "2"),
      data = c(Memberships = "memberships"), tolerance = .tol_count,
      sides = list(
        .hypergraphs(key = "pair", c(
          .read_memberships, .build_hg,
          "hg_line_graph(hg, s = 2, what = \"edges\")",
          "result <- transform(hg_line_graph(hg, s = 2, what = \"edges\"), pair = paste(from, to, sep = \" -- \")) |>",
          "  subset(select = c(pair, weight))",
          "result")),
        .igraph(key = "pair", c(
          .read_memberships, .build_bipartite,
          "projection <- bipartite_projection(bipartite, multiplicity = TRUE, which = \"true\")",
          "projection <- delete_edges(projection, E(projection)[weight < 2])",
          "pairs <- as_data_frame(projection, what = \"edges\")",
          "result <- data.frame(pair = paste(pmin(pairs$from, pairs$to), pmax(pairs$from, pairs$to), sep = \" -- \"),",
          "                     weight = pairs$weight)",
          "result")),
        .base(key = "pair", c(
          .read_memberships, .build_incidence,
          "shared <- crossprod(incidence)",
          "pairs <- which(upper.tri(shared) & shared >= 2, arr.ind = TRUE)",
          "result <- data.frame(pair = paste(rownames(shared)[pairs[, \"row\"]], colnames(shared)[pairs[, \"col\"]], sep = \" -- \"),",
          "                     weight = shared[pairs])",
          "result")))),
    edge_centrality = list(
      title = "s-betweenness and s-closeness of hyperedges (s = 1)",
      description = "Shortest-path centralities of the hyperedges on the unweighted 1-line graph, with the normalisation of NetworkX: betweenness scaled by 2 / ((n - 1)(n - 2)) and closeness (n - 1) / sum of distances. igraph computes both on the line graph it builds itself as a bipartite projection. The line graph of this hypergraph is connected, so no correction for components enters.",
      covers = c("hg_edge_centrality()", "hg_line_graph()"),
      data = c(Memberships = "memberships"), tolerance = .tol_closed,
      sides = list(
        .hypergraphs(key = NULL, c(
          .read_memberships, .build_hg,
          "hg_edge_centrality(hg, s = 1, measure = c(\"betweenness\", \"closeness\"), normalized = TRUE)",
          "result <- transform(hg_edge_centrality(hg, s = 1, measure = c(\"betweenness\", \"closeness\"), normalized = TRUE),",
          "                    quantity = paste(measure, edge, sep = \": \"), kind = measure)",
          "result")),
        .igraph(key = "edge", c(
          .read_memberships, .build_bipartite,
          "line <- bipartite_projection(bipartite, multiplicity = FALSE, which = \"true\")",
          "is_connected(line)",
          "result <- data.frame(edge = V(line)$name,",
          "                     betweenness = betweenness(line, normalized = TRUE),",
          "                     closeness = closeness(line, normalized = TRUE))",
          "result")))),
    transitivity = list(
      title = "Local transitivity on the projection",
      description = "The local clustering coefficient of each node in the simple (unweighted) clique expansion: the share of pairs of a node's neighbours that are themselves neighbours. igraph computes it on the projection it builds itself; base R from the adjacency matrix as diag(A^3) / (k (k - 1)).",
      covers = "hg_transitivity()",
      data = c(Memberships = "memberships"), tolerance = .tol_closed,
      sides = list(
        .hypergraphs(key = "node", c(
          .read_memberships, .build_hg,
          "result <- hg_transitivity(hg, type = \"projection\", what = \"nodes\")",
          "result",
          "hg_transitivity(hg, type = \"projection\", what = \"summary\")")),
        .igraph(key = "node", c(
          .read_memberships, .build_bipartite,
          "expansion <- bipartite_projection(bipartite, multiplicity = FALSE, which = \"false\")",
          "result <- data.frame(node = V(expansion)$name,",
          "                     projection = transitivity(expansion, type = \"local\"))",
          "result")),
        .base(key = "node", c(
          .read_memberships, .build_incidence,
          "adjacency <- (tcrossprod(incidence) > 0) * 1",
          "diag(adjacency) <- 0",
          "neighbours <- rowSums(adjacency)",
          "result <- data.frame(node = rownames(adjacency),",
          "                     projection = diag(adjacency %*% adjacency %*% adjacency) / (neighbours * (neighbours - 1)))",
          "result")))),
    laplacian = list(
      title = "Normalised hypergraph Laplacian (Zhou, Huang and Scholkopf 2007)",
      description = "Every entry of I - Dv^-1/2 H W De^-1 H' Dv^-1/2 with unit hyperedge weights, where Dv holds the node degrees and De the hyperedge sizes. No package on this machine computes it, so the reference is base R written from the formula of the paper.",
      covers = "hg_laplacian()",
      data = c(Memberships = "memberships"), tolerance = .tol_closed,
      sides = list(
        .hypergraphs(key = NULL, c(
          .read_memberships, .build_hg,
          "laplacian <- hg_laplacian(hg, type = \"zhou\")",
          "laplacian",
          "result <- transform(as.data.frame(as.table(laplacian)), quantity = paste(Var1, Var2, sep = \", \"),",
          "                    value = Freq, kind = \"laplacian\")",
          "head(result)")),
        .base(key = NULL, c(
          .read_memberships, .build_incidence,
          "degree_root <- diag(1 / sqrt(rowSums(incidence)))",
          "size_inverse <- diag(1 / colSums(incidence))",
          "laplacian <- diag(nrow(incidence)) - degree_root %*% incidence %*% size_inverse %*% t(incidence) %*% degree_root",
          "dimnames(laplacian) <- list(rownames(incidence), rownames(incidence))",
          "laplacian",
          "result <- transform(as.data.frame(as.table(laplacian)), quantity = paste(Var1, Var2, sep = \", \"),",
          "                    value = Freq, kind = \"laplacian\")",
          "head(result)"))))
  )
}

#' Export one case
#'
#' @param name The case's name.
#' @return One row: the bundle, its sides and how its compared values came out.
export_case <- function(name) {
  cases <- .verify_cases()
  stopifnot("`name` must be one exported case" = length(name) == 1L && name %in% names(cases))
  case <- cases[[name]]
  files <- c(memberships = .verify_memberships(), tribunals = .verify_tribunals())
  data <- stats::setNames(unname(files[case$data]), names(case$data))
  verify_export(bundle = file.path(.verify_base, "bundles", name), title = case$title,
                description = case$description, sides = case$sides, data = data,
                tolerance = case$tolerance, covers = case$covers, variant = case$variant,
                generator = .verify_generator)
}

# igraph is the one reference package; without it nothing here is a comparison
# with an established program, so the export stops with a classed condition.
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop(errorCondition("The reference package igraph is not installed; no case was exported.",
                      class = "hypergraphs_verify_no_reference", call = NULL))
}

requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) requested <- names(.verify_cases())
print(do.call(rbind, lapply(requested, export_case)), row.names = FALSE)
