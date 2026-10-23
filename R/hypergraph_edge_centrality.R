# Walk-based hyperedge centralities from the paper. hypergraphs constructs the
# s-line graph; cograph supplies the ordinary graph centrality kernel.

#' s-betweenness and s-closeness centrality of hyperedges
#'
#' Hyperedges become vertices in the unweighted s-line graph: two are
#' adjacent when they share at least `s` nodes. Betweenness and closeness are
#' then ordinary shortest-path vertex centralities on that graph. At `s = 1`
#' these are the paper's hyperedge 1-betweenness and 1-closeness.
#'
#' @param hg A static `net_hg` or a [temporal_hypergraph()].
#' @param s One or more positive integer intersection thresholds.
#' @param measure Any of `"betweenness"` and `"closeness"`.
#' @param normalized Use NetworkX/HypergraphX normalization: undirected
#'   betweenness is scaled by `2 / ((n - 1) * (n - 2))`, and closeness uses
#'   the Wasserman--Faust correction for disconnected line graphs.
#' @param top Optional number of highest-scoring hyperedges to retain per
#'   `(time, s, measure)` group.
#' @param start,end,step,window,at Measurement grid passed to
#'   [hg_snapshots()] when `hg` is temporal.
#' @param snapshot_mode,multiedges Snapshot `mode` (`"active"` or
#'   `"cumulative"`) and multi-edge handling passed to
#'   [hg_snapshots()] when `hg` is temporal.
#' @return A tidy data frame with `edge`, `s`, `measure`, and `value`; temporal
#'   input adds a leading `time` column, the snapshot time: a date for a
#'   calendar hypergraph, the number on its clock otherwise.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#' hypergraphs. *Philosophical Transactions of the Royal Society A*,
#' 382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' h <- group_hypergraph(
#'   data.frame(member = c("a", "b", "c", "b", "c", "d"),
#'              event = rep(c("e1", "e2"), each = 3)),
#'   node = "member", hyperedge = "event"
#' )
#' hg_edge_centrality(h, s = 1)
#' @export
hg_edge_centrality <- function(hg, s = 1L,
                               measure = c("betweenness", "closeness"),
                               normalized = TRUE, top = NULL, start = NULL,
                               end = NULL, step = NULL, window = NULL, at = NULL,
                               snapshot_mode = c("active", "cumulative"),
                               multiedges = TRUE) {
  measure <- match.arg(measure, several.ok = TRUE)
  snapshot_mode <- .thg_check_mode(snapshot_mode, "hg_edge_centrality", "snapshot_mode")
  if (!is.numeric(s) || length(s) < 1L || any(!is.finite(s)) ||
      any(s < 1) || any(s > .Machine$integer.max) ||
      any(abs(s - round(s)) > sqrt(.Machine$double.eps))) {
    .thg_bad_input("`s` must contain positive whole numbers")
  }
  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    .thg_bad_input("`normalized` must be TRUE or FALSE")
  }
  if (!is.null(top)) top <- .ho_check_count(top, "top", min = 1)

  if (inherits(hg, "net_temporal_hypergraph")) {
    snaps <- hg_snapshots(hg, start = start, end = end, step = step,
                                  window = window, at = at, mode = snapshot_mode,
                                  multiedges = multiedges)
    rows <- lapply(seq_along(snaps), function(i) {
      ans <- hg_edge_centrality(snaps[[i]], s = s, measure = measure,
                                normalized = normalized, top = top)
      # an empty snapshot gives zero rows; the time is repeated to match
      data.frame(time = rep(.thg_snapshot_time(snaps, i), nrow(ans)), ans,
                 row.names = NULL, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  .thg_check_hg(hg)
  if (hg$n_hyperedges == 0L) {
    return(data.frame(edge = character(), s = integer(), measure = character(),
                      value = numeric(), stringsAsFactors = FALSE))
  }
  rows <- list()
  cursor <- 1L
  for (ss in as.integer(s)) {
    line <- hg_line_graph(hg, s = ss, what = "matrix")
    line <- (line != 0) * 1
    n_line <- nrow(line)
    # HypergraphX delegates these measures to NetworkX.  In particular,
    # NetworkX normalizes undirected betweenness by 2 / ((n - 1)(n - 2))
    # and applies the Wasserman--Faust disconnected-graph correction to
    # closeness.  cograph deliberately uses a general max-normalization for
    # most measures, so request raw graph scores here and apply the paper's
    # conventions explicitly.
    # The line graph goes to cograph as an edge list, which a sparse
    # incidence yields without densifying; a hyperedge with no s-neighbour
    # is absent from it and scores 0 below, as in NetworkX.
    triplet <- .thg_general_triplet(line)
    upper <- triplet@i < triplet@j
    line_edges <- data.frame(from = rownames(line)[triplet@i[upper] + 1L],
                             to = rownames(line)[triplet@j[upper] + 1L],
                             weight = rep(1, sum(upper)),
                             stringsAsFactors = FALSE)
    component_size <- NULL
    if (normalized && "closeness" %in% measure && n_line > 0L) {
      membership <- .thg_edge_list_components(triplet@i[upper] + 1L,
                                              triplet@j[upper] + 1L, n_line)
      component_size <- tabulate(membership)[membership]
      names(component_size) <- rownames(line)
    }
    for (metric in measure) {
      value <- if (nrow(line_edges) == 0L) {
        stats::setNames(numeric(n_line), rownames(line))
      } else {
        switch(metric,
          betweenness = cograph::centrality_betweenness(
            line_edges, weighted = FALSE, directed = FALSE, normalized = FALSE
          ),
          closeness = cograph::centrality_closeness(
            line_edges, weighted = FALSE, directed = FALSE,
            normalized = normalized
          )
        )
      }
      value <- stats::setNames(as.numeric(value[rownames(line)]),
                               rownames(line))
      if (normalized && identical(metric, "betweenness")) {
        scale <- if (n_line > 2L) 2 / ((n_line - 1) * (n_line - 2)) else 1
        value <- value * scale
      }
      if (normalized && identical(metric, "closeness") && n_line > 1L) {
        value <- value * (component_size[names(value)] - 1) / (n_line - 1)
      }
      value[!is.finite(value)] <- 0
      tab <- data.frame(
        edge = colnames(line), s = ss, measure = metric,
        value = as.numeric(value[colnames(line)]), stringsAsFactors = FALSE
      )
      tab <- tab[.ho_order(-tab$value, tab$edge), , drop = FALSE]
      if (!is.null(top)) tab <- utils::head(tab, as.integer(top))
      rows[[cursor]] <- tab
      cursor <- cursor + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Connected components of an undirected graph on `n` vertices given as an
# edge list (`from`, `to`: vertex indices), by min-label propagation: every
# vertex starts with its own index and repeatedly takes the smallest label
# among itself and its neighbours. Each pass costs O(edges), the number of
# passes is bounded by the graph diameter plus one, and an isolated vertex
# costs nothing, so a sparse or mostly disconnected line graph stays cheap.
# Returns component ids 1, 2, ... in order of first appearance.
.thg_edge_list_components <- function(from, to, n) {
  label <- seq_len(n)
  ends <- c(from, to)
  other <- c(to, from)
  # bounded by n passes; each pass either lowers a label or stops the loop
  for (pass in seq_len(n)) {
    candidate <- label[other]
    ord <- .ho_order(ends, candidate)
    first <- !duplicated(ends[ord])
    target <- ends[ord][first]
    updated <- label
    updated[target] <- pmin(label[target], candidate[ord][first])
    if (identical(updated, label)) break
    label <- updated
  }
  match(label, unique(label))
}
