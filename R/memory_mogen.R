# ---- hg_get() for net_mogen ----
#
# mogen() returns the multi-order generative model unchanged; its tables are
# read here. "paths" replaces a raw-data path counter: the k-step path counts
# are the counts the fitted order-(k - 1) layer already holds, so they are
# read off the model rather than recomputed from the sequences.

#' Tables of a multi-order generative model
#'
#' @param x A `net_mogen` object from [mogen()].
#' @param what `"orders"` (default) for the per-order model-selection table,
#'   `"transitions"` for the fitted transition probabilities of one layer,
#'   `"paths"` for the frequency of every observed path of `k` states, or
#'   `"pathways"` for the layer's transitions written as pathways
#'   (`"a b -> c"`), the form [simplicial()] turns into simplices.
#' @param ... For `what = "pathways"`, `min_prob` (drop transitions below
#'   this probability). Any other argument raises `hypergraphs_bad_input`.
#' @param order Integer or `NULL`. For `what = "transitions"` and
#'   `"pathways"`: which layer to read, a whole number from `0` to the
#'   model's highest order. Default `NULL` uses the selected order. The
#'   order-0 layer has no memory and no transitions, so at order 0 both
#'   tables have no rows. `order`, `k` and `min_count` given with a table
#'   they do not apply to raise `hypergraphs_bad_input`.
#' @param k Integer vector, for `what = "paths"`: the path length or
#'   lengths in states, each between `2` and the model's highest order plus
#'   one, e.g. `k = 3` or `k = 2:4`. Default `2` (single transitions).
#' @param min_count Integer, for `what = "transitions"` and `"pathways"`:
#'   drop transitions observed fewer times. Default `1`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter. Default `NULL` returns every row.
#' @return A data.frame. For `what = "orders"`, one row per candidate order:
#'   `order`, `log_likelihood`, `aic`, `bic`, `df`, `layer_df`, `p_value`
#'   and `optimal` (logical, `TRUE` for the selected order). `layer_df` counts
#'   the parameters of the order's own layer: one per path of `order`
#'   transitions the observed first-order network allows, less one per row of
#'   the layer's transition matrix (Scholtes 2017, Eq. 7 and its footnote:
#'   on a complete network this is the Markov chain count
#'   \eqn{|V|^k (|V| - 1)}{|V|^k (|V| - 1)}). `df` sums `layer_df` over
#'   orders 0 to `order` (Eq. 8). `p_value` is the
#'   likelihood-ratio test of the model of highest order `order` against the
#'   one of highest order `order - 1`, on `layer_df` degrees of freedom
#'   (Scholtes 2017, Eq. 6); it is `NA` at orders 0 and 1, where the models
#'   are not nested, and where the statistic or `layer_df` is not positive.
#'   With `criterion = "lrt"` the selected order is the highest whose
#'   `p_value` is below `lrt_alpha`, and 1 if there is none. For
#'   `what = "transitions"`, one row per transition: `path`, `count`,
#'   `probability`, `from`, `to`. For `what = "paths"`, one row per distinct
#'   path of each length in `k` (`top` applies within each length), grouped
#'   by length and most frequent first (ties alphabetical): `k`, `path`,
#'   `count`, `proportion` (share of all `k`-state paths, rounded to 4
#'   decimals). For `what = "pathways"`, one row per transition with the
#'   column `pathway`.
#' @references
#' Scholtes, I. (2017). When is a network a network? Multi-order graphical
#' model selection in pathways and temporal networks. \emph{Proceedings of
#' KDD 2017}, 1037--1046. \doi{10.1145/3097983.3098145}
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' mg <- mogen(seqs, max_order = 2)
#' hg_get(mg)
#' hg_get(mg, what = "transitions", top = 5)
#' hg_get(mg, what = "paths", k = 3)
#' @export
hg_get.net_mogen <- function(x, what = c("orders", "transitions", "paths",
                                         "pathways"), ...,
                             order = NULL, k = 2L, min_count = 1L,
                             top = NULL) {
  what <- .ho_match_what(what)
  applies <- switch(what, orders = character(0L), paths = "k",
                    c("order", "min_count"))
  given <- c(order = !is.null(order), k = !missing(k),
             min_count = !missing(min_count))
  .hon_refuse_args(given[setdiff(names(given), applies)], what, "net_mogen")
  if (what == "pathways") {
    .hon_allowed_dots(list(...), "min_prob",
                      "hg_get(what = \"pathways\") on a `net_mogen`")
  } else {
    .ho_no_dots(..., .for = sprintf("hg_get(what = \"%s\") on a `net_mogen`", what))
  }
  if (what %in% c("transitions", "pathways")) {
    order <- .hmg_layer(x, order)
    min_count <- .ho_check_count(min_count, "min_count")
  }
  if (what == "transitions") {
    if (order == 0L) {
      # the order-0 layer is the state distribution: it has no transitions
      return(.ho_top(data.frame(path = character(0L), count = integer(0L),
                                probability = numeric(0L),
                                from = character(0L), to = character(0L),
                                stringsAsFactors = FALSE), top))
    }
    return(.ho_top(Nestimate::mogen_transitions(
      x, order = order, min_count = min_count), top))
  }
  if (what == "pathways") {
    return(.ho_top(.ho_pathway_table(Nestimate::pathways(
      x, order = order, min_count = min_count, ...)), top))
  }
  if (what == "paths") {
    if (!is.numeric(k) || !length(k) || anyNA(k)) {
      .ho_bad_input("`k` must be a non-empty numeric vector of path lengths")
    }
    tables <- lapply(unique(k), \(kk) {
      counts <- .ho_top(.hmg_path_counts(x, kk), top)
      data.frame(k = rep(as.integer(kk), nrow(counts)), counts,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, tables)
    rownames(out) <- NULL
    return(out)
  }
  out <- data.frame(
    order          = as.integer(x$orders),
    log_likelihood = as.numeric(x$log_likelihood),
    aic            = as.numeric(x$aic),
    bic            = as.numeric(x$bic),
    df             = as.integer(x$dof),
    layer_df       = as.integer(x$layer_dof),
    p_value        = as.numeric(x$lrt_p),
    stringsAsFactors = FALSE
  )
  out$optimal <- out$order == x$optimal_order
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' The layer of a multi-order model to read: the selected order, or a
#' validated whole number from 0 to the highest fitted order
#' @noRd
.hmg_layer <- function(x, order) {
  if (is.null(order)) return(as.integer(x$optimal_order))
  max_order <- max(x$orders)
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
      abs(order - round(order)) > sqrt(.Machine$double.eps) || order < 0 ||
      order > max_order) {
    .ho_bad_input(sprintf(paste0(
      "`order` must be one whole number from 0 to %d (the orders the ",
      "model fitted)"), max_order))
  }
  as.integer(order)
}

#' Path counts of length k read off a fitted multi-order model
#'
#' The order-(k - 1) layer of a MOGen model counts every transition between
#' overlapping (k - 1)-grams, i.e. every observed path of k states. The table
#' is laid out as a raw-data k-gram count: most frequent first, ties in
#' sorted path order, proportion of all k-state paths rounded to 4 decimals.
#'
#' @param x A `net_mogen`.
#' @param k Integer path length in states.
#' @return data.frame: path, count, proportion.
#' @noRd
.hmg_path_counts <- function(x, k) {
  max_k <- max(x$orders) + 1L
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != round(k) ||
      k < 2L || k > max_k) {
    stop(errorCondition(
      sprintf(paste0("`k` must be a single whole number from 2 to %d (the ",
                     "model's highest order plus one); refit mogen() with a ",
                     "larger `max_order` for longer paths."), max_k),
      class = "hypergraphs_bad_input", call = NULL))
  }
  tr <- Nestimate::mogen_transitions(x, order = as.integer(k) - 1L,
                                     min_count = 1L)
  if (!nrow(tr)) {
    return(data.frame(path = character(0), count = integer(0),
                      proportion = numeric(0), stringsAsFactors = FALSE))
  }
  counts <- tapply(as.integer(tr$count), tr$path, sum)
  counts <- counts[.ho_sort(names(counts))]
  counts <- counts[.ho_order(as.integer(counts), decreasing = TRUE)]
  out <- data.frame(path = names(counts), count = as.integer(counts),
                    proportion = round(as.numeric(counts) / sum(counts), 4),
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
