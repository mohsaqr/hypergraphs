# Sub-hypergraphs: keep chosen hyperedges, the hyperedges inside a node set,
# or the hyperedges of one source (the citation blocks of one decision).

# Rebuild the derived fields of a net_hg from a new incidence matrix
# (dense or sparse), keeping edge-level metadata aligned with the columns.
.thg_rebuild <- function(hg, incidence, keep_edges) {
  edge_names <- colnames(incidence)
  old_nodes <- hg$nodes
  hg$incidence <- incidence
  hg$nodes <- rownames(incidence)
  hg$n_nodes <- nrow(incidence)
  hg$n_hyperedges <- ncol(incidence)
  hg$hyperedges <- .thg_edge_members(incidence)
  hg$size_distribution <- .thg_size_distribution(hg$hyperedges)
  # node-aligned fields follow the kept rows by name: the planted blocks of a
  # random_hypergraph(type = "sbm") must describe the surviving nodes only
  if (!is.null(hg$blocks)) {
    hg$blocks <- stats::setNames(
      unname(hg$blocks)[match(hg$nodes, old_nodes)], hg$nodes
    )
  }
  if (!is.null(hg$edge_multiplicity)) {
    hg$edge_multiplicity <- hg$edge_multiplicity[keep_edges]
  }
  # the window counts of a window hypergraph are its hyperedge weights
  if (!is.null(hg$window_counts)) {
    hg$window_counts <- hg$window_counts[keep_edges]
  }
  if (!is.null(hg$edge_data)) {
    hg$edge_data <- hg$edge_data[
      match(edge_names, as.character(hg$edge_data$edge)), , drop = FALSE
    ]
    rownames(hg$edge_data) <- NULL
  }
  hg
}

# Size distribution of a hyperedge member list, named "size_<k>" as every
# constructor names it; integer(0) when there are no hyperedges.
.thg_size_distribution <- function(hyperedges) {
  sizes <- lengths(hyperedges)
  if (!length(sizes)) return(integer(0L))
  size_tab <- table(sizes)
  stats::setNames(as.integer(size_tab), paste0("size_", names(size_tab)))
}

# A bare net_hg from an incidence matrix (dense or sparse), with every field
# the object contract requires derived from the matrix itself, so no row or
# column (isolated vertex, empty hyperedge) can be lost on the way.
.thg_from_incidence <- function(incidence, params) {
  hyperedges <- .thg_edge_members(incidence)
  structure(
    list(
      hyperedges = hyperedges,
      incidence = incidence,
      nodes = rownames(incidence) %||% character(0L),
      n_nodes = nrow(incidence),
      n_hyperedges = ncol(incidence),
      size_distribution = .thg_size_distribution(hyperedges),
      params = params
    ),
    class = "net_hg"
  )
}

# Member index vectors of every column, without touching one column at a
# time (which is slow on a sparse matrix with tens of thousands of columns).
.thg_edge_members <- function(incidence) {
  m <- ncol(incidence)
  if (m == 0L) return(list())
  if (methods::is(incidence, "sparseMatrix")) {
    triplet <- methods::as(methods::as(incidence, "generalMatrix"),
                           "TsparseMatrix")
    nz <- triplet@x != 0
    i <- triplet@i[nz] + 1L
    j <- triplet@j[nz] + 1L
  } else {
    nz <- which(incidence != 0, arr.ind = TRUE)
    i <- nz[, 1L]
    j <- nz[, 2L]
  }
  members <- split(i, factor(j, levels = seq_len(m)))
  lapply(unname(members), function(v) sort(unique(v)))
}

#' Sub-hypergraph by hyperedges, nodes or hyperedge attributes
#'
#' Keeps part of a hypergraph and returns it as a `net_hg` with the
#' same incidence weights, edge metadata and multiplicities. Three selectors
#' combine by intersection. `edges` names the hyperedges to keep. `nodes`
#' keeps the hyperedges whose members all lie in the set, the induced
#' sub-hypergraph of HypergraphX. `where` keeps the hyperedges whose
#' attributes take given values, so `where = c(citing = "153-001")` on a
#' citation-block hypergraph is the hypergraph of one citing decision
#' (Coupette et al. 2024, Figure 3).
#'
#' @param hg A `net_hg` (from [group_hypergraph()],
#'   [hg_snapshot()], [text_hypergraph()], ...).
#' @param edges Hyperedge names to keep: a character vector, or a data.frame
#'   with an `edge` column such as the table [hg_edge_centrality()] or
#'   [hg_edges()] returns, so a ranking can be passed straight through.
#' @param nodes Character vector of node names. Only hyperedges contained in
#'   this set are kept, and the node set of the result is exactly `nodes`.
#' @param where Hyperedge attribute values to keep: a named vector or list,
#'   one element per attribute column of the edge metadata (the columns
#'   [temporal_hypergraph()] keeps because they are constant within a
#'   hyperedge), each holding the value or values to keep.
#' @param size Hyperedge sizes to keep, as a vector of member counts (distinct
#'   members): `size = 3` keeps the hyperedges with exactly three members, the
#'   3-uniform hypergraph that [hg_motifs()] needs.
#' @param component `"all"` (default) keeps every connected component;
#'   `"largest"` keeps only the largest connected component, in which two
#'   nodes are connected when a chain of shared hyperedges joins them.
#'   Spectral clustering and label spreading need a connected hypergraph
#'   (Hayashi et al. 2020), and a corpus of short texts often holds a few
#'   texts that share no word with the rest. Applied after the other
#'   filters.
#' @param drop_isolated Drop nodes that belong to no retained hyperedge?
#'   Default `TRUE`. Ignored when `nodes` is given.
#' @return A `net_hg` whose incidence matrix is the selected
#'   sub-matrix of the input, sparse if the input is sparse. Edge metadata
#'   (`edge_data`), duplicate multiplicities (`edge_multiplicity`) and the
#'   window counts of a [window_hypergraph()] are subset alongside. The
#'   nodes the subset removes are recorded in `params$subset$removed`, and
#'   the labels of removed nodes are set aside with a
#'   `hypergraphs_dropped_documents` warning by [hg_classify()],
#'   [hg_keywords()] and the other verbs that take labels, so the table that
#'   built the hypergraph can be passed back whole.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#'
#'   Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#'   random walks, Laplacians, and clustering. *Proceedings of CIKM 2020*,
#'   495-504. \doi{10.1145/3340531.3412034}
#' @examples
#' dat <- data.frame(
#'   member = c("a", "b", "c", "b", "c", "d", "d", "e"),
#'   event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3"),
#'   kind = c("x", "x", "x", "x", "x", "x", "y", "y")
#' )
#' hg <- group_hypergraph(dat, node = "member", hyperedge = "event")
#' hg_subset(hg, edges = c("e1", "e2"))
#' hg_subset(hg, nodes = c("b", "c", "d", "e"))
#' hg_subset(hg, where = c(kind = "y"))
#' hg_subset(hg, edges = c("e1", "e3"), component = "largest")
#' @export
hg_subset <- function(hg, edges = NULL, nodes = NULL, where = NULL,
                      size = NULL, component = c("all", "largest"),
                      drop_isolated = TRUE) {
  .thg_check_hg(hg)
  stopifnot(
    "`drop_isolated` must be TRUE or FALSE" =
      is.logical(drop_isolated) && length(drop_isolated) == 1L &&
      !is.na(drop_isolated)
  )
  if (!is.character(component) || !all(component %in% c("all", "largest"))) {
    .thg_bad_input("`component` must be \"all\" or \"largest\"")
  }
  component <- match.arg(component)
  no_filter <- is.null(edges) && is.null(nodes) && is.null(where) &&
    is.null(size)
  if (no_filter && identical(component, "all")) {
    .thg_bad_input(paste0("supply at least one of `edges`, `nodes`, `where`, ",
                          "`size` or `component = \"largest\"`"))
  }
  if (no_filter) return(.thg_largest_component(hg))
  if (!is.null(size) && (!is.numeric(size) || anyNA(size) || any(size < 1))) {
    .thg_bad_input("`size` must be a vector of positive member counts")
  }
  edge_names <- colnames(hg$incidence)
  keep_edge <- rep(TRUE, hg$n_hyperedges)

  if (!is.null(edges)) {
    if (is.data.frame(edges)) {
      if (!"edge" %in% names(edges)) {
        .thg_bad_input("an `edges` data.frame needs an `edge` column")
      }
      edges <- edges$edge
    }
    edges <- unique(as.character(edges))
    unknown <- setdiff(edges, edge_names)
    if (length(unknown)) {
      .thg_bad_input(sprintf("unknown hyperedge(s): %s",
                             paste(utils::head(unknown, 5L), collapse = ", ")))
    }
    keep_edge <- keep_edge & edge_names %in% edges
  }
  if (!is.null(where)) {
    where <- as.list(where)
    if (is.null(names(where)) || any(!nzchar(names(where)))) {
      .thg_bad_input("`where` must be named by hyperedge attribute")
    }
    unknown <- setdiff(names(where), setdiff(names(hg$edge_data), "edge"))
    if (length(unknown)) {
      .thg_bad_input(sprintf("no hyperedge attribute named %s",
                             paste(unknown, collapse = ", ")))
    }
    row <- match(edge_names, as.character(hg$edge_data$edge))
    for (attribute in names(where)) {
      values <- hg$edge_data[[attribute]][row]
      wanted <- where[[attribute]]
      # a value no hyperedge takes would return an empty hypergraph that
      # reads like a real result; a misspelt value is an error
      absent <- wanted[!wanted %in% values]
      if (length(absent)) .thg_bad_input(.thg_where_absent(attribute, absent, values))
      keep_edge <- keep_edge & !is.na(values) & values %in% wanted
    }
  }
  if (!is.null(size)) {
    members <- as.numeric(Matrix::colSums(.thg_binary(hg$incidence)))
    keep_edge <- keep_edge & members %in% size
  }
  if (!is.null(nodes)) {
    nodes <- as.character(nodes)
    unknown <- setdiff(nodes, hg$nodes)
    if (length(unknown)) {
      .thg_bad_input(sprintf("unknown node(s): %s",
                             paste(utils::head(unknown, 5L), collapse = ", ")))
    }
    inside <- hg$nodes %in% nodes
    b <- .thg_binary(hg$incidence)
    outside_count <- as.numeric(Matrix::colSums(b[!inside, , drop = FALSE]))
    keep_edge <- keep_edge & outside_count == 0
  }

  incidence <- hg$incidence[, keep_edge, drop = FALSE]
  keep_node <- if (!is.null(nodes)) {
    hg$nodes %in% nodes
  } else if (drop_isolated) {
    as.numeric(Matrix::rowSums(.thg_binary(incidence))) > 0
  } else {
    rep(TRUE, hg$n_nodes)
  }
  incidence <- incidence[keep_node, , drop = FALSE]
  out <- .thg_rebuild(hg, incidence, keep_edge)
  out$params$subset <- list(edges = edges, nodes = nodes, where = where,
                            size = size, drop_isolated = drop_isolated,
                            removed = c(hg$params$subset$removed,
                                        setdiff(hg$nodes, out$nodes)))
  if (is.list(hg$text)) out$text <- .thg_subset_text(hg$text, out)
  if (identical(component, "largest")) out <- .thg_largest_component(out)
  out
}

# The text layer of a subset text hypergraph: the weight table keeps the rows
# whose node and hyperedge survive, the document, sentence and vocabulary
# tables follow, and the vocabulary counts are recounted over the kept rows,
# so the layer matches the incidence matrix again.
.thg_subset_text <- function(text, out) {
  columns <- switch(
    text$construction,
    bag = if (identical(text$nodes, "word")) c(node = "word", edge = "doc")
          else c(node = "doc", edge = "word"),
    knn = c(node = "doc", edge = "edge"),
    c(node = "word", edge = "edge")
  )
  weights <- text$weights
  keep <- weights[[columns[["node"]]]] %in% out$nodes &
    weights[[columns[["edge"]]]] %in% colnames(out$incidence)
  text$weights <- weights[keep, , drop = FALSE]
  rownames(text$weights) <- NULL
  if (is.data.frame(text$sentences)) {
    text$sentences <- text$sentences[text$sentences$edge %in%
                                       colnames(out$incidence), , drop = FALSE]
  }
  kept_docs <- if ("doc" %in% names(text$weights)) {
    unique(text$weights$doc)
  } else if (is.data.frame(text$sentences)) {
    unique(text$sentences$doc)
  }
  if (!is.null(kept_docs) && is.data.frame(text$documents)) {
    text$documents <- text$documents[text$documents$doc %in% kept_docs, ,
                                     drop = FALSE]
    rownames(text$documents) <- NULL
  }
  if (is.data.frame(text$vocabulary) && nrow(text$vocabulary) > 0L &&
      "word" %in% names(text$weights)) {
    vocabulary <- text$vocabulary[text$vocabulary$word %in%
                                    text$weights$word, , drop = FALSE]
    if (all(c("count", "doc") %in% names(text$weights))) {
      vocabulary$count <- as.integer(tapply(text$weights$count,
                                            text$weights$word,
                                            sum)[vocabulary$word])
      vocabulary$doc_freq <- as.integer(tapply(text$weights$doc,
                                               text$weights$word,
                                               \(d) length(unique(d)))[
                                                 vocabulary$word])
    }
    rownames(vocabulary) <- NULL
    text$vocabulary <- vocabulary
  }
  text
}

# The largest connected component: nodes joined by a chain of shared
# hyperedges. A connected hypergraph is returned unchanged.
.thg_largest_component <- function(hg) {
  membership <- .thg_binary(hg$incidence)
  adjacency <- Matrix::tcrossprod(membership) > 0
  labels <- .thg_components(adjacency)
  largest <- which.max(tabulate(labels))
  if (all(labels == largest)) return(hg)
  hg_subset(hg, nodes = hg$nodes[labels == largest])
}

# Message for `where` values that no hyperedge takes: the closest values of
# the attribute (approximate matching, as for a misspelt name), else the
# first few in sorted order.
.thg_where_absent <- function(attribute, absent, values) {
  pool <- sort(unique(as.character(values[!is.na(values)])))
  close <- unique(unlist(lapply(as.character(absent), \(a) {
    agrep(a, pool, value = TRUE, ignore.case = TRUE, max.distance = 0.2)
  })))
  hint <- utils::head(if (length(close)) close else pool, 5L)
  sprintf("no hyperedge has %s = %s; %s %s",
          attribute,
          paste0("\"", absent, "\"", collapse = ", "),
          if (length(close)) "did you mean" else "values include",
          paste0("\"", hint, "\"", collapse = ", "))
}
