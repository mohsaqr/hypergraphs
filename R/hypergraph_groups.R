# ---- Bipartite group hypergraph (EG-6) -----------------------------------
# Direct constructor: long-format event data with member + group columns
# becomes a net_hg where each group is a hyperedge spanning all
# members that appeared in it.

# Explode a delimited membership column into one row per member. Bibliographic
# exports ship reference lists and descriptor sets as ";"-joined strings
# (EUR-Lex `citationcelex` / `eurovoc`, Scopus and WoS reference fields), so
# the split, the trim and the drop-empties belong here rather than in every
# caller's preamble.
.thg_expand_delimited <- function(data, column, separator) {
  if (!is.character(separator) || length(separator) != 1L ||
      is.na(separator) || !nzchar(separator)) {
    .thg_bad_input("`separator` must be a single non-empty string")
  }
  if (!is.character(column) || length(column) != 1L ||
      !column %in% names(data)) {
    .thg_bad_input(
      "`separator` needs `actor` to name a column of `data` to split"
    )
  }
  parts <- strsplit(as.character(data[[column]]), separator, fixed = TRUE)
  parts <- lapply(parts, \(x) {
    x <- trimws(x)
    x[nzchar(x)]
  })
  kept <- lengths(parts)
  if (sum(kept) == 0L) {
    .thg_bad_input("no member remains after splitting on `separator`")
  }
  out <- data[rep(seq_len(nrow(data)), kept), , drop = FALSE]
  out[[column]] <- unlist(parts, use.names = FALSE)
  rownames(out) <- NULL
  out
}

#' Hypergraph from co-occurrence data or an edge list
#'
#' Constructs a [net_hg][network_hypergraph] the way a network is
#' defined from data. **Co-presence data** name an `actor` and a `group`:
#' every actor sharing one value of `group` (a session, a team, a citation
#' block) belongs to one hyperedge. An **edge
#' list** names `from` and `to`, and every row is a hyperedge of size two.
#' An optional `weight` column produces a weighted incidence matrix.
#'
#' @param data Data frame in long format, one row per actor-in-group or per
#'   edge; or a clustering of sequences -- a mixture Markov
#'   fit (`net_mmm`), a distance clustering (`net_clustering`), or the
#'   per-cluster networks built from either (`netobject_group`). See
#'   "Clustered sequences" below.
#' @param actor Character. Name of the column whose values become the
#'   hypergraph's nodes (members, participants, actors).
#' @param group Character. Name of the column whose shared values bind
#'   actors into one hyperedge (groups, sessions, teams) -- Dynet's
#'   co-presence vocabulary. When neither pair of columns is named, `from`
#'   and `to` are detected case-insensitively from the alias table
#'   [temporal_hypergraph()] uses (`source`/`target`, `sender`/`receiver`,
#'   ...), and failing that `actor` and `group`.
#' @param weight Character or `NULL`. If supplied, the column is summed per
#'   `(actor, hyperedge)` pair to produce a weighted incidence matrix. Default
#'   `NULL` produces a 0/1 binary incidence matrix.
#' @param nodes Optional vector giving the complete node universe. This keeps
#'   nodes with no observed group memberships as zero-incidence rows, which is
#'   needed for representations such as citation hypergraphs where every
#'   decision is a node but some decisions are never cited.
#' @param separator Split the `actor` column on this string, one row per
#'   member, before building. Bibliographic exports ship a hyperedge's members
#'   as a single delimited cell -- EUR-Lex `citationcelex` and `eurovoc`,
#'   Scopus and Web of Science reference and keyword fields -- so
#'   `separator = ";"` replaces the caller's own split, trim and
#'   drop-empties. Members empty after trimming are dropped, and a row left
#'   with no member contributes no hyperedge.
#' @param sparse Logical. Store incidence as a sparse `Matrix`? Use this for
#'   large, sparse event data such as the full GFCC citation-block corpus.
#' @param from,to Column names of a pairwise edge list, as an alternative to
#'   `actor` and `cooccur_by`.
#' @param member,cooccur_by Deprecated names of `actor` and `group`; using
#'   them warns with a `hypergraphs_deprecated` condition.
#' @param top The number of most frequent sets kept as hyperedges, in each
#'   value of `by` or over all groups; `Inf` keeps every set. For a data
#'   frame, giving `top` or `by` counts sets (see "Frequent sets within
#'   groups"), and `NULL` (default) with no `by` builds one hyperedge per
#'   group. Clustered sequences and topic models keep `8` by default.
#' @param states For counted sets: the states to keep, such as the
#'   events of interest. Every other state is removed from each sequence's
#'   set before counting, and a sequence left with no state is not counted.
#'   `NULL` (default) keeps every state. For a topic model, the topics to
#'   keep.
#' @param threshold Topic model only: the share at which a topic counts as
#'   present in a document, one number in (0, 1]. Required for a topic
#'   model and refused for any other input.
#' @param by Count frequent sets instead of building one hyperedge per
#'   group: the name of a column constant within each `group` (an outcome, a
#'   cluster, a year). Each value of `group` is then one set, the `actor`
#'   values it holds, and the `top` most frequent sets within each value of
#'   `by` become hyperedges (see "Frequent sets within groups"). For a topic
#'   model, a column of the documents table of the fitted hypergraph.
#' @param min_size For counted sets (`by` or a topic model): the smallest set
#'   counted (default `1`); `2` keeps only the sets of two or more. The share
#'   of a set is still taken over all sets of its group.
#'
#' @return A `net_hg` object with the same structure produced by
#'   [network_hypergraph()] (`hyperedges`, `incidence`, `nodes`, `n_nodes`,
#'   `n_hyperedges`, `size_distribution`, `params`), plus `edge_data` when
#'   `data` carries hyperedge attributes (see Details). The `params` list
#'   records `source = "group_hypergraph"` and the original column names.
#'   For clustered sequences the object also has `group_sizes` and
#'   `state_counts` (read them with `hg_get()`) and `params` records
#'   `source = "clustered_sequences"`, `top`, `states` and
#'   `unit = "sequences"`.
#'
#' @details
#' The bipartite representation preserves the full group structure without
#' projecting to a pairwise network. A group of three members A, B, C
#' produces a single 3-hyperedge containing all three, not three pairwise
#' edges AB, AC, BC. This avoids information loss when group interactions are
#' the primary unit of analysis (Perc et al. 2013).
#'
#' Unlike [network_hypergraph()] (which derives hyperedges from a network's
#' clique structure), `group_hypergraph()` takes group memberships
#' directly. The two functions are complementary:
#' \itemize{
#'   \item `group_hypergraph()` - when group membership is observed
#'     (sessions, transactions, co-authorships).
#'   \item `network_hypergraph()` - when only pairwise interactions are
#'     observed and triadic structure must be inferred from triangles.
#' }
#'
#' Rows with `NA` in the actor, hyperedge or weight column are dropped
#' silently.
#'
#' @section Frequent sets within groups:
#' With `by` or `top`, a data frame is read as transactions. Each value of `group`
#' (a trial, a session, a basket) is one set, the distinct `actor` values it
#' holds, and the sets are counted within each value of `by`, which must be
#' constant within a `group`. The `top` most frequent sets of each value of
#' `by` become hyperedges, counted exactly as for clustered sequences below,
#' and `hg_get(hg, what = "sets")` reads them with their `count` and `share`
#' of the value's sets. Without `by`, the sets are counted over all groups
#' together, under one value named after `group` (`"All sessions"`).
#' `plot(hg, group = )` draws the sets of one value, every node sized by the
#' sets that contain it. `states` keeps only
#' the listed `actor` values, and `min_size` the sets of at least that many.
#'
#' @section Topic combinations:
#' Given a mixed-membership topic model fitted by [hg_topics()], every
#' document is reduced to the set of topics whose share in it is at least
#' `threshold`, the thresholded topic presence used to build topic
#' co-occurrence networks (Abuhay et al. 2017; Cassi et al. 2017). A
#' document on love, hate and romance gives the set of those three topics.
#' The `top` most frequent sets become hyperedges, counted exactly as for
#' clustered sequences below, with each document as one transaction, and the
#' topics are the nodes. `top = Inf` keeps every set, and `min_size = 2`
#' keeps only the sets that combine two or more topics. `by` names a column of the documents table of
#' the fitted hypergraph, such as a publication year or an author, and the
#' sets are then counted within each of its values; `NULL` (default) counts
#' them over all documents. A document with no topic at the threshold is not
#' counted. Read the sets with `hg_get(hg, what = "sets")`; `count` is the
#' number of documents with exactly that set and `share` its proportion of
#' the group's documents.
#'
#' @section Clustered sequences:
#' Given a clustering of sequences, every sequence of every group is reduced to
#' the set of its distinct states (order and repetition dropped; `NA` and
#' empty cells ignored), and the `top` most frequent sets of each group
#' become hyperedges -- frequent-itemset support counting with each sequence
#' as one transaction (Agrawal & Srikant 1994), restricted to the sets that
#' occur exactly. Sets of equal count are ranked by their name (states
#' sorted and joined by `" + "`). Groups are named `"Cluster 1"`,
#' `"Cluster 2"`, ... for a `net_mmm` or `net_clustering`, and by the list
#' names of a `netobject_group`, so renamed groups carry through.
#' Integer-coded states are decoded to their labels. The objects are read
#' by their structure; the package that fitted them is not needed.
#'
#' The hyperedges are named `"<group>: <set>"` and carry `group`, `set` and
#' `count` (sequences with exactly that set) as hyperedge attributes, so
#' `plot()` colours and titles them by `count` in `"sequences"`, and
#' `plot(hg, group = "Cluster 1")` draws one group's sets with each node
#' sized by the sequences of that group containing the state. Read the
#' tables with `hg_get(hg, what = "sets")` (one row per hyperedge:
#' its group, set, size, count and share of the group's sequences) and
#' `hg_get(hg, what = "state_counts")` (one row per group and
#' state). A malformed clustering (no `$data`, assignments that do not match
#' it, unnamed networks) raises `hypergraphs_bad_input`.
#'
#' Every other column of `data` that is constant within a hyperedge (a
#' session's date, a team's department) is kept as a hyperedge attribute in
#' the result's `edge_data` table, one row per hyperedge, where
#' [hg_subset()]'s `where` argument and [hg_project()]'s `edge_source` can
#' use it. Columns that vary within a hyperedge describe memberships, not
#' hyperedges, and are left out.
#'
#' @seealso [network_hypergraph()] for the clique-based constructor,
#'   [temporal_hypergraph()] for the same inputs with a clock.
#'
#' @examples
#' df <- data.frame(
#'   person = c("Alice", "Bob", "Carol", "Alice", "Bob",
#'              "Dave", "Carol", "Dave", "Eve"),
#'   session = c("S1", "S1", "S1", "S2", "S2",
#'               "S3", "S3", "S3", "S3")
#' )
#' hg <- group_hypergraph(df, actor = "person", group = "session")
#' print(hg)
#' summary(hg)
#'
#' contacts <- data.frame(from = c("a", "b"), to = c("b", "c"))
#' group_hypergraph(contacts, from = "from", to = "to")
#'
#' # a clustering of sequences, shaped as a mixture Markov fit (net_mmm)
#' fit <- structure(list(
#'   data = data.frame(V1 = c("a", "a", "b", "a", "c"),
#'                     V2 = c("b", "b", "c", "c", "a"),
#'                     V3 = c("a", NA, "a", "b", "b")),
#'   assignments = c(1L, 1L, 2L, 1L, 2L), k = 2L
#' ), class = "net_mmm")
#' sets <- group_hypergraph(fit, top = 3)
#' hg_get(sets, what = "sets")
#' plot(sets, group = "Cluster 1")
#'
#' @references
#' Agrawal, R., & Srikant, R. (1994). Fast algorithms for mining association
#' rules in large databases. In \emph{Proceedings of the 20th International
#' Conference on Very Large Data Bases (VLDB)} (pp. 487-499). Morgan
#' Kaufmann.
#'
#' Abuhay, T. M., Kovalchuk, S. V., Bochenina, K., Kampis, G.,
#' Krzhizhanovskaya, V. V., & Lees, M. H. (2017). Analysis of computational
#' science papers from ICCS 2001-2016 using topic modeling and graph theory.
#' \emph{Procedia Computer Science}, 108, 7-17.
#' \doi{10.1016/j.procs.2017.05.183}
#'
#' Cassi, L., Lahatte, A., Rafols, I., Sautier, P., & de Turckheim, E.
#' (2017). Improving fitness: Mapping research priorities against societal
#' needs on obesity. \emph{Journal of Informetrics}, 11(4), 1095-1113.
#' \doi{10.1016/j.joi.2017.09.010}
#'
#' Perc, M., Gomez-Gardenes, J., Szolnoki, A., Floria, L. M., & Moreno, Y.
#' (2013). Evolutionary dynamics of group interactions on structured
#' populations: a review. \emph{Journal of the Royal Society Interface}
#' 10(80), 20120997. \doi{10.1098/rsif.2012.0997}
#'
#' @note Dense and sparse paths are tested for exact equality. The sparse
#'   path is additionally exercised by the full GFCC reproduction from the
#'   Legal Hypergraphs Zenodo archive (3,618 nodes, 46,165 hyperedges and
#'   77,187 nonzero incidences).
#'
#'   A dense incidence with more than `.Machine$integer.max` cells cannot be
#'   addressed by the flat cell index, and would exhaust memory well before
#'   that. It raises the classed error `hypergraphs_dense_too_large` rather than
#'   attempting the allocation; pass `sparse = TRUE` for data at that scale.
#'
#' @export
group_hypergraph <- function(data, actor = NULL, group = NULL, weight = NULL,
                             nodes = NULL, sparse = FALSE, separator = NULL,
                             from = NULL, to = NULL,
                             member = NULL, cooccur_by = NULL,
                             top = NULL, states = NULL, threshold = NULL,
                             min_size = 1L, by = NULL) {
  if (inherits(data, "net_hg_topics")) {
    return(.thg_topic_set_hypergraph(data, threshold = threshold,
                                     by = by, top = top %||% 8L,
                                     states = states, min_size = min_size))
  }
  if (!is.null(threshold)) {
    .thg_bad_input("`threshold` applies to a topic model fitted by hg_topics()")
  }
  if (inherits(data, c("net_mmm", "net_clustering", "netobject_group"))) {
    if (!is.null(by) || !identical(min_size, 1L)) {
      .thg_bad_input(paste0("`by` and `min_size` apply to a data.frame or a ",
                            "topic model; clustered sequences are counted ",
                            "within their clusters"))
    }
    return(.thg_sequence_set_hypergraph(data, top = top %||% 8L,
                                        states = states))
  }
  counted <- !is.null(by) || !is.null(top)
  if (!counted && (!is.null(states) || !identical(min_size, 1L))) {
    .thg_bad_input(paste0("`states` and `min_size` apply when sets are ",
                          "counted: give `by` or `top`, a topic model, or ",
                          "clustered sequences"))
  }
  stopifnot(is.data.frame(data))
  if (!is.null(separator)) {
    data <- .thg_expand_delimited(data, actor %||% member, separator)
  }
  if (!is.null(member)) {
    .thg_deprecated("member", "actor", "group_hypergraph")
    if (is.null(actor)) actor <- member
  }
  if (!is.null(cooccur_by)) {
    .thg_deprecated("cooccur_by", "group", "group_hypergraph")
    if (is.null(group)) group <- cooccur_by
  }
  member <- actor
  if (is.null(member) && is.null(group) && is.null(from) && is.null(to)) {
    # nothing named: detect an edge list, then co-presence columns
    from <- .thg_match_column(data, "from")
    to <- .thg_match_column(data, "to", exclude = from)
    if (is.null(from) || is.null(to)) {
      from <- NULL
      to <- NULL
      member <- .thg_match_column(data, "actor")
      group <- .thg_match_column(data, "group", exclude = member)
    }
  }
  if (!is.null(from) || !is.null(to)) {
    stopifnot(
      "name either `from` and `to` or `actor` and `group`, not both" =
        is.null(member) && is.null(group),
      "`from` and `to` must both name columns of `data`" =
        is.character(from) && length(from) == 1L && from %in% names(data) &&
        is.character(to) && length(to) == 1L && to %in% names(data)
    )
    edge_id <- paste0("e", seq_len(nrow(data)))
    long <- data.frame(
      actor = c(as.character(data[[from]]), as.character(data[[to]])),
      edge = c(edge_id, edge_id), stringsAsFactors = FALSE
    )
    if (!is.null(weight)) {
      stopifnot(is.character(weight), length(weight) == 1L, weight %in% names(data))
      long$weight <- c(data[[weight]], data[[weight]])
      weight <- "weight"
    }
    # the remaining columns ride along as candidate edge attributes
    for (column in setdiff(names(data), c(from, to, weight))) {
      long[[column]] <- c(data[[column]], data[[column]])
    }
    data <- long
    member <- "actor"
    group <- "edge"
  }
  stopifnot(
    "`actor` must name one column of `data`" =
      is.character(member) && length(member) == 1L && member %in% names(data),
    "`group` must name one column of `data`" =
      is.character(group) && length(group) == 1L && group %in% names(data)
  )
  if (counted) {
    return(.thg_frame_set_hypergraph(data, member, group, by, top %||% 8L,
                                     states, min_size))
  }
  stopifnot(
    is.null(weight) ||
      (is.character(weight) && length(weight) == 1L && weight %in% names(data)),
    is.null(nodes) || is.atomic(nodes),
    is.logical(sparse), length(sparse) == 1L, !is.na(sparse)
  )

  cols <- c(member, group, weight)
  candidates <- setdiff(names(data), cols)
  d <- data[, c(cols, candidates), drop = FALSE]
  d <- d[stats::complete.cases(d[, cols, drop = FALSE]), , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No complete observations after dropping NAs.", call. = FALSE)
  }

  d[[member]] <- as.character(d[[member]])
  d[[group]]  <- as.character(d[[group]])

  if (!is.null(nodes)) {
    nodes <- as.character(nodes)
    if (anyNA(nodes) || any(!nzchar(nodes)) || anyDuplicated(nodes)) {
      stop("`nodes` must contain unique, non-missing node names.", call. = FALSE)
    }
    missing_members <- setdiff(unique(d[[member]]), nodes)
    if (length(missing_members)) {
      stop("Every observed member must occur in `nodes`.", call. = FALSE)
    }
  }
  member_levels <- sort(if (is.null(nodes)) unique(d[[member]]) else nodes)
  group_levels  <- sort(unique(d[[group]]))
  n_members <- length(member_levels)
  n_groups  <- length(group_levels)

  # Map values to row/col indices, then accumulate into the flat cell index.
  # A (member, group) pair may repeat across rows, so duplicated cells must be
  # summed BEFORE assignment -- assigning by an index vector keeps the last
  # write, not the total.
  mi <- match(d[[member]], member_levels)
  gj <- match(d[[group]],  group_levels)
  if (sparse) {
    incidence <- Matrix::sparseMatrix(
      i = mi, j = gj,
      x = if (is.null(weight)) rep.int(1, length(mi)) else as.numeric(d[[weight]]),
      dims = c(n_members, n_groups),
      dimnames = list(member_levels, group_levels)
    )
    if (is.null(weight) && length(incidence@x)) incidence@x[] <- 1
    incidence <- Matrix::drop0(incidence)
  } else {
    # The flat cell index is integer arithmetic, so a dense incidence with more
    # than .Machine$integer.max cells cannot be addressed at all -- and would
    # need >16 Gb before it got that far. Refuse with a pointer to the sparse
    # path rather than overflow to NA indices or exhaust memory.
    if (as.double(n_members) * n_groups > .Machine$integer.max) {
      stop(errorCondition(
        sprintf(
          "a dense incidence of %d members x %d groups (%.3g cells) cannot be built; use `sparse = TRUE`",
          n_members, n_groups, as.double(n_members) * n_groups
        ),
        class = "hypergraphs_dense_too_large", call = NULL
      ))
    }
    cell <- (gj - 1L) * n_members + mi
    if (is.null(weight)) {
      counts <- tabulate(cell, nbins = n_members * n_groups)
      incidence <- matrix(as.integer(counts > 0L), n_members, n_groups,
                          dimnames = list(member_levels, group_levels))
    } else {
      incidence <- matrix(0, n_members, n_groups,
                          dimnames = list(member_levels, group_levels))
      acc <- rowsum(as.numeric(d[[weight]]), cell, reorder = FALSE)
      incidence[as.integer(rownames(acc))] <- acc[, 1L]
    }
  }

  # Drop hyperedges that ended up empty (e.g. all-zero weight)
  he_sizes_pre <- if (sparse) Matrix::colSums(incidence > 0) else
    colSums(incidence > 0)
  keep <- he_sizes_pre > 0
  incidence <- incidence[, keep, drop = FALSE]
  group_levels <- group_levels[keep]
  n_groups <- length(group_levels)

  # Every other column that is constant within a group describes that
  # hyperedge (as in temporal_hypergraph()) and is kept as an attribute;
  # a column that varies within a group describes memberships and is left
  # out. Without such columns the object keeps its original layout.
  attributes <- Filter(function(column) {
    by_group <- split(d[[column]], d[[group]])
    !any(vapply(by_group, function(x) length(unique(x[!is.na(x)])) > 1L,
                logical(1L)))
  }, candidates)
  edge_data <- NULL
  if (length(attributes)) {
    edge_data <- d[match(group_levels, d[[group]]), c(group, attributes),
                   drop = FALSE]
    names(edge_data)[1L] <- "edge"
    rownames(edge_data) <- NULL
  }

  # Member indices of every column from the non-zero cells at once; reading
  # a sparse matrix one column at a time is slow with tens of thousands of
  # hyperedges.
  hyperedges <- .thg_edge_members(incidence)

  he_sizes <- vapply(hyperedges, length, integer(1L))
  size_dist <- if (length(he_sizes)) {
    tab <- table(he_sizes)
    out <- as.integer(tab)
    names(out) <- paste0("size_", names(tab))
    out
  } else {
    integer(0L)
  }

  out <- list(
    hyperedges        = hyperedges,
    incidence         = incidence,
    nodes             = member_levels,
    n_nodes           = n_members,
    n_hyperedges      = n_groups,
    size_distribution = size_dist,
    params = list(
      source         = "group_hypergraph",
      member         = member,
      group          = group,
      weight         = weight,
      nodes          = nodes,
      sparse         = sparse,
      n_observations = nrow(d)
    )
  )
  if (!is.null(edge_data)) out$edge_data <- edge_data
  structure(out, class = "net_hg")
}

# ---- Frequent state sets of clustered sequences --------------------------
# Each sequence of each group reduces to the set of its distinct states; the
# `top` most frequent sets of a group become hyperedges carrying the group and
# the number of sequences with exactly that set. This is support counting of
# itemsets (Agrawal & Srikant 1994) with each sequence as one transaction,
# restricted to the observed transactions themselves.
.thg_sequence_set_hypergraph <- function(x, top = 8L, states = NULL) {
  if (!is.numeric(top) || length(top) != 1L || !is.finite(top) || top < 1 ||
      top != round(top)) {
    .thg_bad_input("`top` must be one whole number of at least 1")
  }
  if (!is.null(states) && (!is.character(states) || !length(states) ||
                           anyNA(states))) {
    .thg_bad_input("`states` must be a character vector of state names to keep")
  }
  out <- .thg_set_hypergraph(.coerce_grouped_sequences(x), top, states,
                             item_order = sort)
  out$params <- c(out$params, list(
    source = "clustered_sequences",
    input = class(x)[1L],
    top = as.integer(top),
    states = states,
    unit = "sequences"
  ))
  out$params <- out$params[!duplicated(names(out$params), fromLast = TRUE)]
  out
}

# The shared core: `grouped` is a named list of groups, each a list of
# transactions (character vectors of items). Each transaction reduces to the
# set of its distinct items, ordered by `item_order`, and the `top` most
# frequent sets of each group become hyperedges.
.thg_set_hypergraph <- function(grouped, top, states, item_order,
                                min_size = 1L, prefix = TRUE) {
  labels <- names(grouped)
  per_group <- lapply(labels, \(g) {
    sets <- lapply(grouped[[g]], \(v) {
      v <- unique(v)
      if (!is.null(states)) v <- v[v %in% states]
      item_order(v)
    })
    sets <- sets[lengths(sets) >= min_size]
    keys <- vapply(sets, paste, character(1L), collapse = " + ")
    # table() orders the sets by name; the stable order() then keeps that
    # name order among sets of equal count
    tab <- table(keys)
    ranked <- order(-as.vector(tab))
    kept <- names(tab)[utils::head(ranked, min(top, length(ranked)))]
    counts <- as.integer(tab[kept])
    members <- sets[match(kept, keys)]
    in_state <- table(unlist(sets, use.names = FALSE))
    list(
      members = if (length(kept)) data.frame(
        state = unlist(members, use.names = FALSE),
        edge = rep(if (prefix) paste0(g, ": ", kept) else kept,
                   lengths(members)),
        group = g,
        set = rep(kept, lengths(members)),
        count = rep(counts, lengths(members)),
        stringsAsFactors = FALSE
      ),
      sizes = data.frame(group = g, sequences = length(grouped[[g]]),
                         with_states = length(sets), sets = length(tab),
                         stringsAsFactors = FALSE),
      nodes = data.frame(group = rep(g, length(in_state)),
                         node = names(in_state),
                         count = as.integer(in_state),
                         stringsAsFactors = FALSE)
    )
  })
  members <- do.call(rbind, lapply(per_group, `[[`, "members"))
  if (is.null(members) || !nrow(members)) {
    .thg_bad_input(paste0("no sequence or document keeps a state or topic: ",
                          "nothing to build a hyperedge from"))
  }
  out <- group_hypergraph(members, actor = "state", group = "edge")
  out$group_sizes <- do.call(rbind, lapply(per_group, `[[`, "sizes"))
  out$state_counts <- do.call(rbind, lapply(per_group, `[[`, "nodes"))
  rownames(out$state_counts) <- NULL
  out
}

# ---- Frequent sets of a data frame, counted within groups -----------------
# Each value of `group` (a trial, a session, a basket) is one transaction, the
# set of `actor` values it holds; the transactions are counted within the
# values of `by` (an outcome, a cluster, a year), as for clustered sequences.
.thg_frame_set_hypergraph <- function(data, actor, group, by, top, states,
                                      min_size) {
  if (!is.null(by) &&
      (!is.character(by) || length(by) != 1L || !by %in% names(data))) {
    .thg_bad_input("`by` must name one column of `data`")
  }
  if (!is.numeric(top) || length(top) != 1L || is.na(top) || top < 1 ||
      (is.finite(top) && top != round(top))) {
    .thg_bad_input("`top` must be one whole number of at least 1, or Inf")
  }
  if (!is.numeric(min_size) || length(min_size) != 1L ||
      !is.finite(min_size) || min_size < 1 || min_size != round(min_size)) {
    .thg_bad_input("`min_size` must be one whole number of at least 1")
  }
  if (!is.null(states) && (!is.character(states) || !length(states) ||
                           anyNA(states))) {
    .thg_bad_input("`states` must be a character vector of values to keep")
  }
  d <- data[stats::complete.cases(data[, c(actor, group, by)]),
            c(actor, group, by), drop = FALSE]
  if (!nrow(d)) {
    .thg_bad_input("no complete row of `actor`, `group` and `by`")
  }
  items <- as.character(d[[actor]])
  transaction <- as.character(d[[group]])
  # without `by`, the sets are counted over every group together
  label <- if (is.null(by)) {
    rep(sprintf("All %ss", group), nrow(d))
  } else {
    as.character(d[[by]])
  }
  spread <- tapply(label, transaction, \(v) length(unique(v)))
  if (any(spread > 1L)) {
    .thg_bad_input(sprintf(paste0(
      "`%s` must be constant within each `%s`; %d of them have more than one ",
      "value"), by, group, sum(spread > 1L)))
  }
  transactions <- split(items, factor(transaction, levels = unique(transaction)))
  of_transaction <- label[match(names(transactions), transaction)]
  levels <- .thg_kw_natural(unique(of_transaction))
  grouped <- lapply(stats::setNames(levels, levels),
                    \(g) unname(transactions[of_transaction == g]))
  out <- .thg_set_hypergraph(grouped, top, states, item_order = sort,
                             min_size = as.integer(min_size),
                             prefix = !is.null(by))
  out$params <- c(out$params, list(
    source = "frame_sets",
    member = actor,
    group = group,
    by = by,
    top = top,
    min_size = as.integer(min_size),
    states = states,
    unit = paste0(group, "s"),
    item = actor
  ))
  out$params <- out$params[!duplicated(names(out$params), fromLast = TRUE)]
  out
}

# ---- Topic combinations of a mixed-membership topic model -----------------
# Each document is a transaction whose items are the topics with a share of
# at least `threshold` (the thresholded topic presence of Abuhay et al. 2017
# and Cassi et al. 2017); the frequent sets are counted as for sequences.
.thg_topic_set_hypergraph <- function(x, threshold, by, top, states,
                                      min_size) {
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0 || threshold > 1) {
    .thg_bad_input(paste0("a topic model needs `threshold`, the share at ",
                          "which a topic counts as present in a document, ",
                          "one number in (0, 1]"))
  }
  if (!is.numeric(top) || length(top) != 1L || is.na(top) || top < 1 ||
      (is.finite(top) && top != round(top))) {
    .thg_bad_input("`top` must be one whole number of at least 1, or Inf")
  }
  if (!is.numeric(min_size) || length(min_size) != 1L ||
      !is.finite(min_size) || min_size < 1 || min_size != round(min_size)) {
    .thg_bad_input("`min_size` must be one whole number of at least 1")
  }
  if (!is.null(states) && (!is.character(states) || !length(states) ||
                           anyNA(states))) {
    .thg_bad_input("`states` must be a character vector of topics to keep")
  }
  shares <- x$shares
  present <- shares[shares$share >= threshold, , drop = FALSE]
  documents <- unique(shares$node)
  topics_of <- split(present$topic, factor(present$node, levels = documents))
  labels <- if (is.null(by)) {
    rep("All documents", length(documents))
  } else {
    meta <- x$documents_table
    if (!is.character(by) || length(by) != 1L || is.null(meta) ||
        !by %in% names(meta)) {
      .thg_bad_input(sprintf(paste0(
        "`by` must name a column of the documents table of the fitted ",
        "hypergraph (%s)"),
        if (is.null(meta)) "none was kept" else
          paste(setdiff(names(meta), "doc"), collapse = ", ")))
    }
    value <- meta[[by]][match(documents, meta$doc)]
    if (anyNA(value)) {
      .thg_bad_input(sprintf("`%s` is missing for some documents", by))
    }
    as.character(value)
  }
  levels <- .thg_kw_natural(unique(labels))
  grouped <- lapply(stats::setNames(levels, levels),
                    \(g) unname(topics_of[labels == g]))
  natural <- \(v) v[order(as.integer(sub("^Topic ", "", v)))]
  # ungrouped sets are named by their topics alone
  out <- .thg_set_hypergraph(grouped, top, states, item_order = natural,
                             min_size = as.integer(min_size),
                             prefix = !is.null(by))
  out$params <- c(out$params, list(
    source = "topic_sets",
    input = "net_hg_topics",
    threshold = threshold,
    by = by,
    top = top,
    min_size = as.integer(min_size),
    states = states,
    unit = "documents",
    item = "topic"
  ))
  out$params <- out$params[!duplicated(names(out$params), fromLast = TRUE)]
  out
}
