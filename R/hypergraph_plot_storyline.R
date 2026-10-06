# The storyline view of a temporal hypergraph (Tanahashi & Ma 2012): every
# node is a line through time, and a hyperedge gathers the lines of its
# members at its time. plot.net_temporal_hypergraph(type = "storyline")
# routes here. Columns are the hyperedges in order of their start; the
# horizontal spacing is by order, not by elapsed time.

# Line colours: the Okabe-Ito palette without its yellow, which is too faint
# for a thin line on white. Each line's points also take a shape, so no node
# is told by colour alone, and the legend names them. Eight colours and nine
# shapes cycle with different periods, so the first 72 lines all differ in
# their pair of colour and shape.
.thg_storyline_colours <- function(n) {
  rep_len(.thg_okabe_ito[seq_len(8L)], n)
}
.thg_storyline_shapes <- function(n) {
  rep_len(c(16, 17, 15, 18, 8, 4, 3, 1, 0), n)
}

# Rows of the lines, one per (node, column) between the node's first and
# last hyperedge. One sweep carries an ordering of the nodes through the
# columns; at each column the members of its hyperedge move to the median
# of their current positions and keep their relative order (the barycentre
# rule of Sugiyama et al. 1981). A node between them moves aside. The row
# of a node is its rank among the nodes alive at that column. The sweep is
# started from the order of first appearance and then, `sweeps` times,
# from the mean row of each node in the previous layout; the layout with
# the fewest crossings is kept (the earliest on ties), so the result is
# deterministic.
.thg_storyline_layout <- function(node, col, n_col, sweeps = 4L) {
  nodes <- sort(unique(node))
  members <- split(node, factor(col, levels = seq_len(n_col)))
  first <- tapply(col, node, min)
  last <- tapply(col, node, max)
  alive <- lapply(seq_len(n_col), \(j) nodes[first[nodes] <= j & last[nodes] >= j])

  gather <- function(pos, j) {
    mem <- members[[j]]
    if (length(mem) < 2L) return(pos)
    centre <- stats::median(pos[mem])
    offsets <- rank(pos[mem], ties.method = "first") - (length(mem) + 1) / 2
    key <- ifelse(pos == centre, pos + 0.5, pos)
    key[mem] <- centre + offsets * 1e-3
    stats::setNames(rank(key, ties.method = "first"), names(pos))
  }
  sweep <- function(start) {
    orders <- Reduce(gather, seq_len(n_col), accumulate = TRUE,
                     init = start)[-1L]
    rows <- do.call(rbind, lapply(seq_len(n_col), \(j) {
      data.frame(node = alive[[j]], col = j,
                 y = as.numeric(rank(orders[[j]][alive[[j]]])))
    }))
    rows[order(rows$node, rows$col), , drop = FALSE]
  }
  crossings <- function(rows) {
    pairs <- lapply(seq_len(n_col - 1L), \(j) {
      here <- rows[rows$col == j, , drop = FALSE]
      there <- rows[rows$col == j + 1L, , drop = FALSE]
      both <- intersect(here$node, there$node)
      if (length(both) < 2L) return(0)
      a <- here$y[match(both, here$node)]
      b <- there$y[match(both, there$node)]
      flips <- sign(outer(a, a, "-")) != sign(outer(b, b, "-"))
      sum(flips) / 2
    })
    sum(unlist(pairs))
  }
  next_start <- function(rows) {
    mean_row <- tapply(rows$y, rows$node, mean)
    key <- mean_row[nodes] + first[nodes] * 1e-6
    stats::setNames(rank(key, ties.method = "first"), nodes)
  }

  appearance <- stats::setNames(
    rank(first[nodes] + seq_along(nodes) * 1e-9, ties.method = "first"), nodes
  )
  first_layout <- sweep(appearance)
  # Reduce() over no sweeps returns its init itself, not a list holding it
  layouts <- if (sweeps < 1L) list(first_layout) else
    Reduce(\(previous, i) sweep(next_start(previous)), seq_len(sweeps),
           accumulate = TRUE, init = first_layout)
  counts <- vapply(layouts, crossings, numeric(1))
  best <- layouts[[which.min(counts)]]
  attr(best, "crossings") <- min(counts)
  best
}

# Re-space the lines of each column so that the gap between two neighbouring
# lines shrinks linearly with the number of plotted hyperedges their nodes
# share: a full row for no shared hyperedge, `floor` for the most shared by
# any pair in the plot, 1 - (1 - floor) * shared / max_shared in between.
# A linear scale keeps every level of strength visible; a gap such as
# 1 / (1 + shared) cut at a floor merges all strong ties into one. The order
# of the lines, and so the crossings, are those of the layout.
.thg_storyline_spacing <- function(rows, node, edge, floor = 0.35) {
  nodes <- sort(unique(node))
  member <- unclass(table(factor(node, levels = nodes), edge)) > 0
  shared <- tcrossprod(member * 1)
  diag(shared) <- 0
  strongest <- max(shared)
  respaced <- lapply(split(rows, rows$col), \(r) {
    r <- r[order(r$y), , drop = FALSE]
    if (nrow(r) > 1L) {
      pairs <- cbind(utils::head(r$node, -1L), utils::tail(r$node, -1L))
      gaps <- if (strongest > 0) 1 - (1 - floor) * shared[pairs] / strongest else
        rep(1, nrow(pairs))
      r$y <- 1 + c(0, cumsum(gaps))
    }
    r
  })
  out <- do.call(rbind, respaced)
  out <- out[order(out$node, out$col), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "crossings") <- attr(rows, "crossings")
  out
}

# Checks shared by the storyline of a temporal hypergraph and of a static
# hypergraph with ordered hyperedges.
.thg_storyline_checks <- function(top, edge_labels, point_size, width_by) {
  if (!is.null(width_by) && !identical(width_by, "degree")) {
    .thg_bad_input("`width_by` must be NULL or \"degree\" for a storyline")
  }
  if (!is.null(top) && (!is.numeric(top) || length(top) != 1L ||
                        !is.finite(top) || top < 1 ||
                        abs(top - round(top)) > sqrt(.Machine$double.eps))) {
    .thg_bad_input("`top` must be NULL or one positive whole number")
  }
  if (!is.logical(edge_labels) || length(edge_labels) != 1L || is.na(edge_labels)) {
    .thg_bad_input("`edge_labels` must be TRUE or FALSE for a storyline")
  }
  if (!is.numeric(point_size) || length(point_size) != 1L ||
      !is.finite(point_size) || point_size <= 0) {
    .thg_bad_input("`point_size` must be one positive number")
  }
  invisible(TRUE)
}

# The `top` nodes with the most hyperedges in `m` (node, edge), ties by name,
# with the counts that rank them.
.thg_storyline_nodes <- function(m, top) {
  counts <- table(m$node)
  ranked <- names(counts)[order(-as.integer(counts), names(counts))]
  list(counts = counts,
       chosen = if (is.null(top)) ranked else utils::head(ranked, top))
}

# Storyline of a temporal hypergraph: columns are hyperedges in order of
# their start, labelled with the hyperedge and its start time.
.thg_plot_storyline <- function(x, top = 8L, start = NULL, end = NULL,
                                edge_labels = TRUE, point_size = 2.5,
                                spacing = c("even", "strength"),
                                width_by = NULL) {
  spacing <- match.arg(spacing)
  .thg_storyline_checks(top, edge_labels, point_size, width_by)
  lo <- .thg_as_time(start, x, "start") %||% -Inf
  hi <- .thg_as_time(end, x, "end") %||% Inf

  m <- unique(x$memberships[c("node", "edge", "start")])
  edge_start <- tapply(m$start, m$edge, min)
  in_window <- names(edge_start)[edge_start >= lo & edge_start <= hi]
  m <- unique(m[m$edge %in% in_window, c("node", "edge"), drop = FALSE])
  if (nrow(m) == 0L) {
    .thg_bad_input("no hyperedge begins between `start` and `end`")
  }

  picked <- .thg_storyline_nodes(m, top)
  when <- .thg_calendar(edge_start, x$origin, x$time_unit)
  edge_text <- stats::setNames(
    if (edge_labels) sprintf("%s  %s", names(edge_start), format(when)) else
      rep("", length(edge_start)),
    names(edge_start))
  .thg_storyline_draw(m, picked$chosen, picked$counts, edge_start, edge_text,
                      spacing, width_by, point_size,
                      node_title = x$params$node %||% "node",
                      edge_title = x$params$hyperedge %||% "hyperedge")
}

# Storyline of a static hypergraph whose hyperedges have an order: the
# stored order, or `sort_by` ascending (a position, a time, any value), as
# the ordered windows of window_hypergraph(collapse = FALSE).
.thg_plot_storyline_hg <- function(x, sort_by = NULL, edge_labels = TRUE,
                                   top = 8L, point_size = 2.5,
                                   spacing = c("even", "strength"),
                                   width_by = NULL) {
  spacing <- match.arg(spacing)
  .thg_storyline_checks(top, edge_labels, point_size, width_by)
  edge_names <- colnames(x$incidence) %||% paste0("h", seq_len(x$n_hyperedges))
  m <- unique(data.frame(
    node = x$nodes[unlist(x$hyperedges, use.names = FALSE)],
    edge = rep(edge_names, lengths(x$hyperedges))
  ))
  if (nrow(m) == 0L) .thg_bad_input("the hypergraph has no memberships to plot")
  picked <- .thg_storyline_nodes(m, top)
  key <- if (is.null(sort_by)) seq_along(edge_names) else
    .thg_edge_aesthetic(x, sort_by, "sort_by")
  edge_key <- stats::setNames(key, edge_names)
  edge_text <- stats::setNames(if (edge_labels) edge_names else
    rep("", length(edge_names)), edge_names)
  .thg_storyline_draw(m, picked$chosen, picked$counts, edge_key, edge_text,
                      spacing, width_by, point_size,
                      node_title = x$params$node %||% "node",
                      edge_title = x$params$hyperedge %||% "hyperedge")
}

# Lay out and draw a storyline: `m` holds the memberships (node, edge),
# `chosen` the nodes to draw, `counts` their numbers of hyperedges,
# `edge_key` the value that orders the hyperedges (ties by name) and
# `edge_text` the label of each column, both named by hyperedge.
.thg_storyline_draw <- function(m, chosen, counts, edge_key, edge_text,
                                spacing, width_by, point_size, node_title,
                                edge_title) {
  m <- m[m$node %in% chosen, , drop = FALSE]
  edges <- unique(m$edge)
  edge_levels <- edges[order(edge_key[edges], edges)]
  m$col <- match(m$edge, edge_levels)
  n_col <- length(edge_levels)
  rows <- .thg_storyline_layout(m$node, m$col, n_col)
  if (identical(spacing, "strength")) {
    rows <- .thg_storyline_spacing(rows, m$node, m$edge)
  }

  path <- rbind(transform(rows, x = col - 0.3), transform(rows, x = col + 0.3))
  path <- path[order(path$node, path$x), , drop = FALSE]
  hits <- merge(m[c("node", "col")], rows, by = c("node", "col"))
  span <- function(f) as.numeric(tapply(hits$y, factor(hits$col, levels = seq_len(n_col)), f))
  bars <- data.frame(col = seq_len(n_col), lo = span(min), hi = span(max))
  bars <- bars[!is.na(bars$lo) & bars$hi > bars$lo, , drop = FALSE]

  colours <- stats::setNames(.thg_storyline_colours(length(chosen)), chosen)
  shapes <- stats::setNames(.thg_storyline_shapes(length(chosen)), chosen)
  legend_title <- node_title
  x_labels <- unname(edge_text[edge_levels])
  ink <- "#2B2B2B"

  # a line's width follows its node's number of hyperedges in the period,
  # the count `top` ranks by; the widest line (2) stays narrower than the
  # points, so a hyperedge's point never disappears into its line
  lines <- if (is.null(width_by)) {
    ggplot2::geom_path(
      data = path,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, group = .data$node,
                             colour = .data$node),
      linewidth = 1
    )
  } else {
    path$degree <- as.integer(counts[path$node])
    list(
      ggplot2::geom_path(
        data = path,
        mapping = ggplot2::aes(x = .data$x, y = .data$y, group = .data$node,
                               colour = .data$node, linewidth = .data$degree)
      ),
      ggplot2::scale_linewidth_continuous(
        range = c(0.7, 2),
        breaks = unique(round(pretty(range(path$degree)))),
        name = paste0(edge_title, "s")
      )
    )
  }

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = bars,
      mapping = ggplot2::aes(x = .data$col, xend = .data$col,
                             y = .data$lo - 0.25, yend = .data$hi + 0.25),
      linewidth = 6, lineend = "round", colour = "grey85"
    ) +
    lines +
    ggplot2::geom_point(
      data = hits,
      mapping = ggplot2::aes(x = .data$col, y = .data$y, colour = .data$node,
                             shape = .data$node),
      size = point_size
    ) +
    ggplot2::scale_colour_manual(values = colours, breaks = chosen,
                                 name = legend_title) +
    ggplot2::scale_shape_manual(values = shapes, breaks = chosen,
                                name = legend_title) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2L, byrow = TRUE,
                                                   order = 1L),
                    shape = ggplot2::guide_legend(nrow = 2L, byrow = TRUE,
                                                  order = 1L),
                    linewidth = ggplot2::guide_legend(order = 2L)) +
    ggplot2::scale_y_reverse(breaks = NULL,
                             expand = ggplot2::expansion(add = c(0.6, 0.9))) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_col), labels = x_labels,
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7, colour = ink, angle = 90,
                                          hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      # the width legend on a row of its own below the names
      legend.box = if (is.null(width_by)) NULL else "vertical"
    )
}
