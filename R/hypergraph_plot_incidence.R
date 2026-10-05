# The incidence view of a hypergraph, after UpSet (Lex et al. 2014): one row
# per node, one column per hyperedge, and a vertical bar joining the members
# of each hyperedge. plot.net_hg(type = "incidence") routes here. The figure
# grows linearly in the number of hyperedges, so it stays legible where hulls
# merge into one shape.

# Shapes paired with the discrete colours, so a group is never told by colour
# alone; the first five take a fill.
.thg_incidence_shapes <- c(21, 24, 22, 23, 25, 16, 17, 15, 18)

# Above this many cells the grey grid of empty cells is left out: it would
# be hundreds of thousands of points that carry no information.
.thg_incidence_max_grid <- 20000L

.thg_plot_incidence <- function(x, color_by = NULL, sort_by = NULL,
                                labels = TRUE, edge_labels = FALSE,
                                label_size = NULL, edge_label_size = NULL,
                                legend_title = NULL) {
  sizes_ok <- vapply(list(label_size, edge_label_size), \(s) {
    is.null(s) || (is.numeric(s) && length(s) == 1L && is.finite(s) && s > 0)
  }, logical(1))
  if (!all(sizes_ok)) {
    .thg_bad_input("`label_size` and `edge_label_size` must be single positive numbers")
  }
  edge_names <- colnames(x$incidence) %||% paste0("h", seq_len(x$n_hyperedges))
  members <- data.frame(
    edge = rep(edge_names, lengths(x$hyperedges)),
    node = x$nodes[unlist(x$hyperedges, use.names = FALSE)]
  )
  if (nrow(members) == 0L) {
    .thg_bad_input("the hypergraph has no memberships to plot")
  }

  # rows: nodes by decreasing hyperdegree, then name
  degree <- tabulate(match(members$node, x$nodes), nbins = x$n_nodes)
  node_levels <- x$nodes[order(-degree, x$nodes)]
  members$row <- match(members$node, node_levels)
  row_degree <- degree[match(node_levels, x$nodes)]

  # columns: by `sort_by` (counts largest first, as hg_get() sorts; times --
  # dates, or a clock column such as a snapshot's numeric `start` -- and text
  # ascending), then by the highest member row, then by name, so hyperedges
  # that share the busiest nodes stand together
  by_edge <- split(members$row, factor(members$edge, levels = edge_names))
  first_row <- vapply(by_edge, \(r) if (length(r)) min(r) else Inf, numeric(1))
  key <- .thg_edge_aesthetic(x, sort_by, "sort_by")
  col_order <- if (is.null(key)) {
    order(first_row, edge_names)
  } else if (is.numeric(key) && !(is.character(sort_by) &&
                                   length(sort_by) == 1L &&
                                   .thg_norm_name(sort_by) %in% .thg_clock_names())) {
    order(-key, first_row, edge_names)
  } else {
    order(key, first_row, edge_names)
  }
  edge_levels <- edge_names[col_order]
  members$col <- match(members$edge, edge_levels)
  n_row <- length(node_levels)
  n_col <- length(edge_levels)

  by_col <- split(members$row, factor(members$col, levels = seq_len(n_col)))
  spans <- data.frame(
    col = seq_len(n_col),
    lo = vapply(by_col, \(r) if (length(r)) min(r) else NA_real_, numeric(1)),
    hi = vapply(by_col, \(r) if (length(r)) max(r) else NA_real_, numeric(1))
  )
  spans <- spans[!is.na(spans$lo) & spans$hi > spans$lo, , drop = FALSE]

  # colour (and, when discrete, shape) per hyperedge
  values <- .thg_edge_aesthetic(x, color_by, "color_by")
  legend_name <- legend_title %||%
    (if (is.character(color_by) && length(color_by) == 1L) color_by else "group")
  if (!is.null(values)) {
    by_name <- stats::setNames(values, edge_names)
    members$value <- unname(by_name[members$edge])
    spans$value <- unname(by_name[edge_levels[spans$col]])
    if (!is.numeric(values)) {
      fill <- .thg_fill_colours(values)
      group_levels <- names(fill$palette)
      members$value <- factor(members$value, levels = group_levels)
      spans$value <- factor(spans$value, levels = group_levels)
      shapes <- stats::setNames(
        rep_len(.thg_incidence_shapes, length(group_levels)), group_levels
      )
    }
  }

  ink <- "#2B2B2B"
  p <- ggplot2::ggplot()
  if (n_row * n_col <= .thg_incidence_max_grid) {
    p <- p + ggplot2::geom_point(
      data = expand.grid(row = seq_len(n_row), col = seq_len(n_col)),
      mapping = ggplot2::aes(x = .data$col, y = .data$row),
      colour = "grey88", size = 1.4
    )
  }
  if (is.null(values)) {
    p <- p +
      ggplot2::geom_segment(
        data = spans,
        mapping = ggplot2::aes(x = .data$col, xend = .data$col,
                               y = .data$lo, yend = .data$hi),
        colour = ink, linewidth = 0.9
      ) +
      ggplot2::geom_point(
        data = members, mapping = ggplot2::aes(x = .data$col, y = .data$row),
        colour = ink, size = 2.2
      )
  } else {
    p <- p +
      ggplot2::geom_segment(
        data = spans,
        mapping = ggplot2::aes(x = .data$col, xend = .data$col, y = .data$lo,
                               yend = .data$hi, colour = .data$value),
        linewidth = 0.9
      )
    if (is.numeric(values)) {
      p <- p +
        ggplot2::geom_point(
          data = members,
          mapping = ggplot2::aes(x = .data$col, y = .data$row,
                                 colour = .data$value),
          size = 2.2
        ) +
        ggplot2::scale_colour_gradientn(colours = .thg_okabe_ito_ramp(),
                                        name = legend_name)
    } else {
      p <- p +
        ggplot2::geom_point(
          data = members,
          mapping = ggplot2::aes(x = .data$col, y = .data$row,
                                 colour = .data$value, fill = .data$value,
                                 shape = .data$value),
          size = 2.3
        ) +
        ggplot2::scale_colour_manual(values = fill$palette, name = legend_name,
                                     drop = FALSE) +
        ggplot2::scale_fill_manual(values = fill$palette, name = legend_name,
                                   drop = FALSE) +
        ggplot2::scale_shape_manual(values = shapes, name = legend_name,
                                    drop = FALSE) +
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2L),
                        fill = ggplot2::guide_legend(ncol = 2L),
                        shape = ggplot2::guide_legend(ncol = 2L))
    }
  }

  # hyperdegree bars beside the rows, numbered
  bar_left <- n_col + 1
  bar_width <- max(2, 0.12 * n_col)
  bars <- data.frame(row = seq_len(n_row), degree = row_degree,
                     right = bar_left + bar_width * row_degree / max(row_degree))
  p <- p +
    ggplot2::geom_rect(
      data = bars,
      mapping = ggplot2::aes(xmin = bar_left, xmax = .data$right,
                             ymin = .data$row - 0.35, ymax = .data$row + 0.35),
      fill = "grey65", inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = bars,
      mapping = ggplot2::aes(x = .data$right + 0.2, y = .data$row,
                             label = .data$degree),
      hjust = 0, size = 2.6, colour = ink
    ) +
    ggplot2::annotate("text", x = bar_left, y = 0, label = "Hyperdegree",
                      hjust = 0, vjust = 0, size = 3, colour = ink)

  # hyperedge-size bars above the columns, when sizes differ
  edge_size <- lengths(x$hyperedges)[match(edge_levels, edge_names)]
  top_space <- 0
  if (length(unique(edge_size)) > 1L) {
    top_space <- max(2, 0.12 * n_row)
    sized <- data.frame(col = seq_len(n_col), size = edge_size,
                        top = 0.4 - top_space * edge_size / max(edge_size))
    p <- p +
      ggplot2::geom_rect(
        data = sized,
        mapping = ggplot2::aes(xmin = .data$col - 0.35, xmax = .data$col + 0.35,
                               ymin = .data$top, ymax = 0.4),
        fill = "grey65", inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = sized,
        mapping = ggplot2::aes(x = .data$col, y = .data$top - 0.15,
                               label = .data$size),
        vjust = 0, size = 2.6, colour = ink
      ) +
      ggplot2::annotate("text", x = bar_left, y = 0.4 - top_space,
                        label = "Size", hjust = 0, size = 3, colour = ink)
  }

  row_text <- .thg_incidence_text(labels, node_levels, "labels", "node")
  col_text <- if (isFALSE(edge_labels)) rep("", n_col) else
    .thg_incidence_text(edge_labels, edge_levels, "edge_labels", "hyperedge")
  row_pt <- if (is.null(label_size)) 8 else label_size * ggplot2::.pt
  col_pt <- if (is.null(edge_label_size)) 7 else edge_label_size * ggplot2::.pt
  p +
    ggplot2::scale_y_reverse(
      breaks = seq_len(n_row), labels = row_text,
      expand = ggplot2::expansion(add = c(0.7, 0.7 + top_space + 1))
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n_col), labels = col_text,
      expand = ggplot2::expansion(add = c(0.7, bar_width + 1.2))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = row_pt, colour = ink),
      axis.text.x = ggplot2::element_text(size = col_pt, colour = ink,
                                          angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    )
}

# Axis text for the incidence view: TRUE writes the names, FALSE none, and a
# vector named by node (or hyperedge) replaces the names it covers.
.thg_incidence_text <- function(spec, keys, arg, what) {
  if (isTRUE(spec)) return(keys)
  if (isFALSE(spec)) return(rep("", length(keys)))
  if (is.null(names(spec))) {
    .thg_bad_input(sprintf("`%s` must be TRUE, FALSE or a vector named by %s",
                           arg, what))
  }
  replacement <- unname(spec[keys])
  ifelse(is.na(replacement), keys, as.character(replacement))
}
