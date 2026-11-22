# ---------------------------------------------------------------------------
# hg_compare() on a hypergraph: permutation and cluster-bootstrap comparison
# of node centralities between two groups of hyperedges
# ---------------------------------------------------------------------------

#' Compare node centralities between two groups of hyperedges
#'
#' Tests whether the nodes of a hypergraph hold a different position in two
#' groups of its hyperedges, such as the messages of high and low performers.
#' The hyperedges of each group form a hypergraph of their own; a centrality
#' is computed for every node in each, and the statistic is the difference,
#' first group minus second. Its null distribution comes from permuting the
#' group labels, and its interval from a bootstrap.
#'
#' **Measures.** `"degree"` is the share of the group's hyperedges that
#' contain the node. `"clique"` is the clique-motif eigenvector centrality,
#' scaled to unit length, and `"pagerank"` the EDVW hypergraph PageRank,
#' summing to one, both as in [hg_centrality()] (Benson 2019; Chitra and
#' Raphael 2019). All three are normalised within a group, so a group with
#' more hyperedges does not score higher for that reason; the unnormalised
#' centralities of [hg_centrality()] are not offered for that reason. A node
#' that occurs in no hyperedge of a group scores 0 there.
#'
#' **Nesting.** Hyperedges often nest in larger units: messages in prompts,
#' prompts in students. `unit` names the attribute that carries one group
#' label (by default the hyperedge itself) and `cluster` the attribute whose
#' observations are not exchangeable with one another (a student). When the
#' group varies within clusters, the labels of the units are shuffled among
#' the units of the same cluster, a restricted permutation that keeps every
#' cluster's composition (Good 2005, ch. 5); a cluster whose units all share
#' one group contributes no variation. When the group is constant within
#' every cluster, the clusters themselves carry the labels and are shuffled.
#' Without `cluster` the units are shuffled freely. The p-value is
#' `(1 + extreme) / (n_perm + 1)` (Phipson and Smyth 2010), and the
#' p-values of all nodes and measures are adjusted together by the method of
#' Benjamini and Hochberg (1995). The scheme is the permutation comparison of
#' network centralities of van Borkulo et al. (2023) applied to hypergraphs.
#'
#' **Intervals.** The interval of each difference is a percentile interval
#' from a cluster bootstrap: clusters (or units, without `cluster`) are drawn
#' with replacement and keep all of their hyperedges (Davison and Hinkley
#' 1997, sec. 3.8). A draw that leaves one group without hyperedges has no
#' difference; such draws are left out, counted in `n_boot_used`, and
#' reported with the warning `hypergraphs_undefined_statistic`.
#'
#' @param x A hypergraph (`net_hg`) whose hyperedges carry the attributes
#'   named by `group`, `cluster` and `unit`, such as a [hypergraph()] of a
#'   membership table whose further columns are constant within each
#'   hyperedge (see [hg_get()] with `what = "edge_data"`).
#' @param group The hyperedge attribute that defines the groups.
#' @param groups The two groups to compare, as a character vector of length
#'   2; the first is the one the differences start from. May be omitted when
#'   `group` takes exactly two values.
#' @param cluster `NULL` (default) or the hyperedge attribute of the
#'   independent clusters, such as a student.
#' @param unit `NULL` (default; each hyperedge) or the hyperedge attribute of
#'   the unit that carries one group label, such as a prompt holding several
#'   messages. The group must be constant within each unit, and each unit
#'   must lie in one cluster.
#' @param type Measures to compare: any of `"degree"`, `"clique"` and
#'   `"pagerank"`. Default all three.
#' @param n_perm,n_boot Integers >= 2. Label permutations and bootstrap draws.
#'   Default `1000` each.
#' @param level Confidence level of the bootstrap intervals. Default `0.95`.
#' @param alpha Significance level for the `significant` flag on the adjusted
#'   p-values. Default `0.05`.
#' @param damping PageRank damping factor, as in [hg_centrality()]. Default
#'   `0.85`.
#' @param seed `NULL` or one whole number. The caller's random stream is
#'   restored on exit.
#' @param ... Must be empty: an argument that only the memory-network method
#'   takes raises `hypergraphs_bad_input`.
#' @return An object of class `hypergraphs_hg_comparison` with `print`,
#'   `summary` and `plot` methods. [hg_get()] returns its table, one row per
#'   node and measure, with `node`, `measure`, the node's value in each group
#'   (two columns named after the groups), `diff` (first group minus
#'   second), `lower` and `upper` (bootstrap interval), `p_value`, `p_adj`
#'   (Benjamini-Hochberg), `significant` and `higher` (the group with the
#'   larger value, or `"equal"`).
#'
#'   Raises `hypergraphs_bad_input` when an attribute is missing, `groups`
#'   does not name two groups of `group`, the group varies within a unit or
#'   a unit spans clusters. Hyperedges without a value of `group` are left
#'   out with the warning `hypergraphs_dropped_rows`.
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#' rate. \emph{Journal of the Royal Statistical Society B}, 57(1), 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Benson, A. R. (2019). Three hypergraph eigenvector centralities.
#' \emph{SIAM Journal on Mathematics of Data Science} 1(2), 293-312.
#'
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. \emph{Proceedings of the 36th
#' International Conference on Machine Learning}, PMLR 97, 1172-1181.
#'
#' Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and
#' their Application}. Cambridge University Press.
#'
#' Good, P. (2005). \emph{Permutation, Parametric and Bootstrap Tests of
#' Hypotheses} (3rd ed.). Springer.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero. \emph{Statistical Applications in Genetics and Molecular Biology}
#' 9(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio,
#' P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023). Comparing
#' network structures on three aspects: A permutation test.
#' \emph{Psychological Methods} 28(6), 1273-1285. \doi{10.1037/met0000476}
#' @examples
#' moves <- data.frame(
#'   message = rep(paste0("m", 1:16), each = 2),
#'   code = rep(c("plan", "ask", "plan", "check", "ask", "check",
#'                "plan", "ask"), 4),
#'   student = rep(paste0("s", 1:4), each = 8),
#'   level = rep(c("high", "high", "low", "low"), 8)
#' )
#' messages <- hypergraph(moves, node = "code", hyperedge = "message")
#' comparison <- hg_compare(messages, group = "level", cluster = "student",
#'                          n_perm = 99, n_boot = 99, seed = 1)
#' comparison
#' hg_get(comparison, measure = "degree")
#' @seealso [hg_centrality()], [hg_compare()]
#' @export
hg_compare.net_hg <- function(x, group, groups = NULL, cluster = NULL,
                              unit = NULL,
                              type = c("degree", "clique", "pagerank"),
                              n_perm = 1000L, n_boot = 1000L, level = 0.95,
                              alpha = 0.05, damping = 0.85, seed = NULL,
                              ...) {
  .ho_no_dots(..., .for = "a hypergraph")
  .thg_check_hg(x)
  if (missing(group)) {
    .thg_bad_input("`group` must name the hyperedge attribute that defines the groups")
  }
  measures <- c("degree", "clique", "pagerank")
  if (!is.character(type) || !length(type) || anyNA(type) ||
      !all(type %in% measures) || anyDuplicated(type) > 0L) {
    .thg_bad_input(sprintf("`type` must name each measure once, among: %s",
                           paste(measures, collapse = ", ")))
  }
  n_perm <- .ho_check_count(n_perm, "n_perm", min = 2)
  n_boot <- .ho_check_count(n_boot, "n_boot", min = 2)
  level <- .ho_check_number(level, "level", min = 0, max = 1)
  alpha <- .ho_check_number(alpha, "alpha", min = 0, max = 1)
  if (level <= 0 || level >= 1) .thg_bad_input("`level` must lie strictly between 0 and 1")
  if (alpha <= 0 || alpha >= 1) .thg_bad_input("`alpha` must lie strictly between 0 and 1")
  damping <- .ho_check_number(damping, "damping", min = 0, max = 1)
  .ho_check_seed(seed)

  design <- .hg_cmp_design(x, group, groups, cluster, unit)
  groups <- design$groups
  reserved <- c("node", "measure", "diff", "lower", "upper", "p_value",
                "p_adj", "significant", "higher")
  if (any(groups %in% reserved)) {
    .thg_bad_input(sprintf("a group may not be called %s",
                           paste(intersect(groups, reserved), collapse = ", ")))
  }
  edges <- design$edges
  b <- .thg_binary(x$incidence[, edges, drop = FALSE])
  nodes <- .ho_sort(x$nodes[as.numeric(Matrix::rowSums(b)) > 0])
  score <- function(idx) .hg_cmp_scores(x, idx, nodes, type, damping)
  difference <- function(labels) {
    first <- edges[labels == groups[1L]]
    second <- edges[labels == groups[2L]]
    if (!length(first) || !length(second)) return(NULL)
    score(first) - score(second)
  }

  .ho_local_seed(seed)
  value_1 <- score(edges[design$labels == groups[1L]])
  value_2 <- score(edges[design$labels == groups[2L]])
  observed <- value_1 - value_2

  permuted <- vapply(seq_len(n_perm), \(r) {
    difference(.hg_cmp_shuffle(design))
  }, observed)
  tolerance <- sqrt(.Machine$double.eps)
  extreme <- rowSums(abs(permuted) >= c(abs(observed)) - tolerance, dims = 2L)
  p_value <- (1 + extreme) / (n_perm + 1)

  draws <- lapply(seq_len(n_boot), \(r) .hg_cmp_resample(design))
  boot <- lapply(draws, \(d) {
    first <- d$edges[d$labels == groups[1L]]
    second <- d$edges[d$labels == groups[2L]]
    if (!length(first) || !length(second)) return(NULL)
    score(first) - score(second)
  })
  used <- !vapply(boot, is.null, logical(1L))
  if (sum(used) < 2L) {
    .thg_bad_input(paste0("fewer than two bootstrap draws hold both groups: ",
                          "the groups share too few clusters to resample"))
  }
  if (!all(used)) {
    warning(warningCondition(sprintf(
      "%d of %d bootstrap draws hold only one group and are left out of the intervals",
      sum(!used), n_boot), class = "hypergraphs_undefined_statistic", call = NULL))
  }
  boot <- simplify2array(boot[used])
  tail_prob <- (1 - level) / 2
  lower <- apply(boot, c(1L, 2L), stats::quantile, probs = tail_prob,
                 names = FALSE)
  upper <- apply(boot, c(1L, 2L), stats::quantile, probs = 1 - tail_prob,
                 names = FALSE)

  diff <- as.vector(observed)
  table <- data.frame(
    node = rep(nodes, length(type)),
    measure = rep(type, each = length(nodes)),
    first = as.vector(value_1), second = as.vector(value_2),
    diff = diff, lower = as.vector(lower), upper = as.vector(upper),
    p_value = as.vector(p_value), stringsAsFactors = FALSE
  )
  names(table)[3:4] <- groups
  table$p_adj <- stats::p.adjust(table$p_value, method = "BH")
  table$significant <- table$p_adj < alpha
  table$higher <- ifelse(abs(diff) < tolerance, "equal",
                         ifelse(diff > 0, groups[1L], groups[2L]))

  structure(list(
    table = table, groups = groups, group = group, cluster = cluster,
    unit = unit, permuted = design$permuted, type = type,
    n_hyperedges = stats::setNames(
      c(sum(design$labels == groups[1L]), sum(design$labels == groups[2L])),
      groups),
    n_units = length(design$unit_ids), n_clusters = design$n_clusters,
    n_perm = n_perm, n_boot = n_boot, n_boot_used = sum(used),
    level = level, alpha = alpha
  ), class = "hypergraphs_hg_comparison")
}

# The design of a comparison: the hyperedges of the two groups (column
# indices of the incidence), their group labels, the unit and cluster each
# lies in, and how labels are permuted. Validates the attributes.
.hg_cmp_design <- function(x, group, groups, cluster, unit) {
  ed <- x$edge_data
  attributes <- setdiff(names(ed), "edge")
  check <- function(value, arg, null_ok) {
    if (is.null(value) && null_ok) return(invisible(NULL))
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !value %in% attributes) {
      .thg_bad_input(sprintf(
        "`%s` must name one hyperedge attribute; the hypergraph has %s", arg,
        if (length(attributes)) paste(attributes, collapse = ", ") else "none"))
    }
  }
  check(group, "group", FALSE)
  check(cluster, "cluster", TRUE)
  check(unit, "unit", TRUE)
  row <- match(colnames(x$incidence), as.character(ed$edge))
  labels <- as.character(ed[[group]][row])
  available <- .ho_sort(unique(labels[!is.na(labels)]))
  if (is.null(groups)) {
    if (length(available) != 2L) {
      .thg_bad_input(sprintf(
        "`%s` has %d values (%s); choose two with `groups =`", group,
        length(available), paste(available, collapse = ", ")))
    }
    groups <- available
  }
  groups <- as.character(groups)
  if (length(groups) != 2L || anyNA(groups) || anyDuplicated(groups) ||
      !all(groups %in% available)) {
    .thg_bad_input(sprintf("`groups` must be two of: %s",
                           paste(available, collapse = ", ")))
  }
  if (anyNA(labels)) {
    warning(warningCondition(sprintf(
      "%d of %d hyperedges have no value of `%s` and are left out",
      sum(is.na(labels)), length(labels), group),
      class = "hypergraphs_dropped_rows", call = NULL))
  }
  edges <- which(!is.na(labels) & labels %in% groups)
  labels <- labels[edges]
  units <- if (is.null(unit)) colnames(x$incidence)[edges] else
    as.character(ed[[unit]][row][edges])
  clusters <- if (is.null(cluster)) rep("all", length(edges)) else
    as.character(ed[[cluster]][row][edges])
  if (anyNA(units) || anyNA(clusters)) {
    .thg_bad_input("`unit` and `cluster` must have a value on every hyperedge compared")
  }
  spread <- function(values, by) {
    any(tapply(values, by, \(v) length(unique(v))) > 1L)
  }
  if (spread(labels, units)) {
    .thg_bad_input(sprintf("`%s` must be constant within each `%s`", group, unit))
  }
  if (spread(clusters, units)) {
    .thg_bad_input(sprintf("each `%s` must lie in one `%s`", unit, cluster))
  }
  unit_ids <- unique(units)
  first <- match(unit_ids, units)
  unit_group <- labels[first]
  unit_cluster <- clusters[first]
  # the group constant within every cluster: the clusters carry the labels
  between <- !is.null(cluster) && !spread(unit_group, unit_cluster)
  list(groups = groups, edges = edges, labels = labels,
       unit_of_edge = match(units, unit_ids), unit_ids = unit_ids,
       unit_group = unit_group, unit_cluster = unit_cluster,
       cluster_of_edge = match(clusters, unique(clusters)),
       n_clusters = length(unique(clusters)),
       permuted = if (between) "clusters" else if (is.null(cluster)) "units"
                  else "units within clusters")
}

# One permutation of the group labels of the compared hyperedges.
.hg_cmp_shuffle <- function(design) {
  unit_group <- design$unit_group
  shuffled <- if (identical(design$permuted, "clusters")) {
    cluster_ids <- unique(design$unit_cluster)
    cluster_group <- unit_group[match(cluster_ids, design$unit_cluster)]
    cluster_group[sample.int(length(cluster_group))][
      match(design$unit_cluster, cluster_ids)]
  } else {
    strata <- split(seq_along(unit_group), design$unit_cluster)
    out <- unit_group
    positions <- unlist(strata, use.names = FALSE)
    drawn <- unlist(lapply(strata, \(s) s[sample.int(length(s))]),
                    use.names = FALSE)
    out[positions] <- unit_group[drawn]
    out
  }
  shuffled[design$unit_of_edge]
}

# One cluster-bootstrap draw: clusters drawn with replacement, each keeping
# all of its compared hyperedges with their own labels.
.hg_cmp_resample <- function(design) {
  by_cluster <- split(seq_along(design$edges), design$cluster_of_edge)
  drawn <- sample.int(length(by_cluster), replace = TRUE)
  rows <- unlist(by_cluster[drawn], use.names = FALSE)
  list(edges = design$edges[rows], labels = design$labels[rows])
}

# The measures of every node in the hypergraph made of the hyperedges `idx`
# (column indices, repeats allowed): a node x measure matrix, 0 for a node in
# no hyperedge of `idx`.
.hg_cmp_scores <- function(x, idx, nodes, type, damping) {
  incidence <- x$incidence[, idx, drop = FALSE]
  colnames(incidence) <- paste0("e", seq_along(idx))
  counts <- as.numeric(Matrix::rowSums(.thg_binary(incidence)))
  out <- matrix(0, length(nodes), length(type),
                dimnames = list(nodes, type))
  at <- match(x$nodes, nodes)
  inside <- !is.na(at)
  if ("degree" %in% type) {
    out[at[inside], "degree"] <- counts[inside] / length(idx)
  }
  centralities <- intersect(type, c("clique", "pagerank"))
  if (length(centralities)) {
    present <- counts > 0
    part <- x
    part$edge_data <- NULL
    part <- .thg_rebuild(part, incidence[present, , drop = FALSE], idx)
    cen <- .hg_centrality_fit(part, type = centralities, damping = damping)
    out[match(cen$node, nodes), centralities] <- as.matrix(cen[centralities])
  }
  out
}

# ---------------------------------------------------------------------------
# S3 methods: hypergraphs_hg_comparison
# ---------------------------------------------------------------------------

#' Table of a hypergraph comparison
#'
#' @param x A `hypergraphs_hg_comparison` from [hg_compare()] on a
#'   hypergraph.
#' @param what `NULL` or `"nodes"` (the only table): one row per node and
#'   measure.
#' @param ... Must be empty: an argument this method does not take raises
#'   `hypergraphs_bad_input`.
#' @param measure `NULL` (default, every measure) or the measures to keep.
#' @param significant Logical. `TRUE` keeps the rows whose adjusted p-value
#'   falls below the object's `alpha`. Default `FALSE`.
#' @param sort_by `NULL` (measure, then node; default), `"abs_diff"`
#'   (largest absolute difference first) or `"p_adj"` (smallest first),
#'   ties broken by measure and node.
#' @param top Integer or `NULL`. The first `top` rows, after the filters
#'   and `sort_by`.
#' @return A data.frame, one row per node and measure; see
#'   [hg_compare.net_hg()] for the columns.
#' @examples
#' moves <- data.frame(
#'   message = rep(paste0("m", 1:16), each = 2),
#'   code = rep(c("plan", "ask", "plan", "check", "ask", "check",
#'                "plan", "ask"), 4),
#'   student = rep(paste0("s", 1:4), each = 8),
#'   level = rep(c("high", "high", "low", "low"), 8)
#' )
#' messages <- hypergraph(moves, node = "code", hyperedge = "message")
#' comparison <- hg_compare(messages, group = "level", cluster = "student",
#'                          n_perm = 99, n_boot = 99, seed = 1)
#' hg_get(comparison, sort_by = "abs_diff", top = 3)
#' @export
hg_get.hypergraphs_hg_comparison <- function(x, what = NULL, ...,
                                             measure = NULL,
                                             significant = FALSE,
                                             sort_by = NULL, top = NULL) {
  .ho_no_dots(..., .for = "hg_get() on a `hypergraphs_hg_comparison`")
  .ho_match_what(what, "nodes")
  if (!is.logical(significant) || length(significant) != 1L ||
      is.na(significant)) {
    .ho_bad_input("`significant` must be TRUE or FALSE")
  }
  out <- x$table
  if (!is.null(measure)) {
    if (!is.character(measure) || !length(measure) ||
        !all(measure %in% x$type)) {
      .ho_bad_input(sprintf("`measure` must be among: %s",
                            paste(x$type, collapse = ", ")))
    }
    out <- out[out$measure %in% measure, , drop = FALSE]
  }
  if (isTRUE(significant)) out <- out[out$significant, , drop = FALSE]
  if (!is.null(sort_by)) {
    sort_by <- .ho_match_what(sort_by, c("abs_diff", "p_adj"))
    key <- if (identical(sort_by, "abs_diff")) -abs(out$diff) else out$p_adj
    out <- out[.ho_order(key, match(out$measure, x$type), out$node), ,
               drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Print a hypergraph comparison
#'
#' @param x A `hypergraphs_hg_comparison` object.
#' @param n Number of rows of the table to print. Default `10`.
#' @param ... Additional arguments (ignored).
#' @return The input `x`, invisibly.
#' @examples
#' moves <- data.frame(
#'   message = rep(paste0("m", 1:16), each = 2),
#'   code = rep(c("plan", "ask", "plan", "check", "ask", "check",
#'                "plan", "ask"), 4),
#'   student = rep(paste0("s", 1:4), each = 8),
#'   level = rep(c("high", "high", "low", "low"), 8)
#' )
#' messages <- hypergraph(moves, node = "code", hyperedge = "message")
#' print(hg_compare(messages, group = "level", n_perm = 99, n_boot = 99,
#'                  seed = 1), n = 5)
#' @export
print.hypergraphs_hg_comparison <- function(x, n = 10L, ...) {
  cat(sprintf("Hypergraph comparison of `%s`: %s (%d hyperedges) vs %s (%d hyperedges)\n",
              x$group, x$groups[1L], x$n_hyperedges[[1L]], x$groups[2L],
              x$n_hyperedges[[2L]]))
  cat(sprintf("  %d permutations of %s; %d of %d bootstrap draws (%g%% intervals)\n",
              x$n_perm, x$permuted, x$n_boot_used, x$n_boot, 100 * x$level))
  cat(sprintf("  Significant (BH, alpha = %.2f): %d of %d node-measure pairs\n",
              x$alpha, sum(x$table$significant), nrow(x$table)))
  .ho_print_table(x, n)
  invisible(x)
}

#' @rdname result-summary
#' @export
summary.hypergraphs_hg_comparison <- function(object, ...) {
  tab <- object$table
  by_measure <- do.call(rbind, lapply(object$type, \(m) {
    d <- tab[tab$measure == m, , drop = FALSE]
    data.frame(measure = m, n_nodes = nrow(d),
               n_significant = sum(d$significant),
               max_abs_diff = max(abs(d$diff)), stringsAsFactors = FALSE)
  }))
  .ho_summary(object, list(by_measure = by_measure))
}

#' Plot a hypergraph comparison
#'
#' The difference of every node, first group minus second, with its
#' bootstrap interval, one panel per measure. Significance (BH-adjusted) is
#' shown by both colour and point shape.
#'
#' @param x A `hypergraphs_hg_comparison` object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object, invisibly.
#' @examples
#' moves <- data.frame(
#'   message = rep(paste0("m", 1:16), each = 2),
#'   code = rep(c("plan", "ask", "plan", "check", "ask", "check",
#'                "plan", "ask"), 4),
#'   student = rep(paste0("s", 1:4), each = 8),
#'   level = rep(c("high", "high", "low", "low"), 8)
#' )
#' messages <- hypergraph(moves, node = "code", hyperedge = "message")
#' plot(hg_compare(messages, group = "level", n_perm = 99, n_boot = 99,
#'                 seed = 1))
#' @export
plot.hypergraphs_hg_comparison <- function(x, ...) {
  d <- x$table
  d$node <- factor(d$node, levels = rev(.ho_sort(unique(d$node))))
  d$measure <- factor(d$measure, levels = x$type)
  d$sig_f <- factor(ifelse(d$significant, "significant", "not significant"),
                    levels = c("significant", "not significant"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$diff, y = .data$node)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "#999999") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                   colour = .data$sig_f),
      width = 0.25, orientation = "y") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$sig_f,
                                     shape = .data$sig_f), size = 2.4) +
    ggplot2::scale_colour_manual(
      values = c(significant = "#D55E00", `not significant` = "#0072B2"),
      name = NULL, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(significant = 17, `not significant` = 16),
      name = NULL, drop = FALSE) +
    ggplot2::facet_wrap(~ measure, scales = "free_x") +
    ggplot2::labs(
      x = sprintf("Difference, %s minus %s (%g%% interval)",
                  x$groups[1L], x$groups[2L], 100 * x$level),
      y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  print(p)
  invisible(p)
}
