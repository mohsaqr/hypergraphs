# ---- hypergraph(): the main constructor ------------------------------------
# One verb for every input. Each method hands its input to the constructor of
# that kind of data and returns that constructor's result unchanged.

#' Build a hypergraph from any input
#'
#' `hypergraph()` is the main constructor of the package. It reads the
#' input and builds the hypergraph the input describes, through the
#' constructor of that kind of data:
#'
#' * A data frame in long format with `actor` and `group` (or `from` and
#'   `to`) describes observed groups, and every group becomes a hyperedge
#'   ([group_hypergraph()]). With `by` or `top` the groups are read as sets
#'   and the most frequent sets become hyperedges. With `time`, `start` or
#'   `end` the hyperedges carry a clock ([temporal_hypergraph()]). With
#'   `window`, `step` or `action` the data are sequences, and every window of
#'   consecutive actions becomes a hyperedge ([window_hypergraph()]).
#' * A list of sequences is read the same way ([window_hypergraph()]).
#' * A network, as a weight matrix, a sparse matrix, a `netobject` or a
#'   `cograph_network`, has its cliques promoted to hyperedges
#'   ([network_hypergraph()]); with `window`, a model object built from
#'   sequences is read as sequences instead.
#' * A topic model fitted by [hg_topics()], or a clustering of sequences,
#'   becomes the hypergraph of its frequent sets ([group_hypergraph()]).
#'
#' Every argument in `...` is passed to that constructor, whose
#' documentation describes it, and the result is the constructor's own.
#'
#' @param data The input: a data frame, a list of sequences, a network, a
#'   topic model or a clustering of sequences.
#' @param ... Arguments of the constructor that the input selects.
#' @return A `net_hg` object (a `net_temporal_hypergraph` for data with a
#'   clock), identical to the result of the constructor called directly.
#' @section Conditions:
#' `hypergraphs_bad_input` for an input of a class no constructor reads,
#' naming the class, and the conditions of the constructor called.
#' @references
#' Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M., Patania,
#' A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
#' interactions: Structure and dynamics. \emph{Physics Reports}, 874, 1-92.
#' \doi{10.1016/j.physrep.2020.05.004}
#' @seealso [group_hypergraph()], [window_hypergraph()],
#'   [temporal_hypergraph()], [network_hypergraph()], [text_hypergraph()]
#'   for a corpus.
#' @examples
#' meetings <- data.frame(
#'   person = c("Alice", "Bob", "Carol", "Alice", "Bob", "Dave", "Eve"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m3", "m3"))
#' hypergraph(meetings, actor = "person", group = "meeting")
#'
#' sessions <- list(c("a", "b", "c", "a"), c("b", "c", "d"))
#' hypergraph(sessions, window = 2)
#'
#' weights <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3,
#'                   dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
#' hypergraph(weights)
#' @export
hypergraph <- function(data, ...) UseMethod("hypergraph")

# arguments that select a route for a data frame or a model object
.hgm_window_args <- c("window", "step", "action", "session",
                      "time_threshold", "min_weight")
.hgm_clock_args <- c("time", "start", "end", "observation_start",
                     "observation_end", "time_unit")

#' @export
hypergraph.data.frame <- function(data, ...) {
  given <- names(list(...))
  if (any(given %in% .hgm_window_args)) {
    return(window_hypergraph(data, ...))
  }
  if (any(given %in% .hgm_clock_args)) {
    return(temporal_hypergraph(data, ...))
  }
  group_hypergraph(data, ...)
}

#' @export
hypergraph.list <- function(data, ...) window_hypergraph(data, ...)

.hgm_network <- function(data, ...) {
  if (any(names(list(...)) %in% .hgm_window_args)) {
    return(window_hypergraph(data, ...))
  }
  network_hypergraph(data, ...)
}

#' @export
hypergraph.matrix <- function(data, ...) .hgm_network(data, ...)

# clique finding works on a dense adjacency; a sparse matrix is converted
#' @export
hypergraph.Matrix <- function(data, ...) .hgm_network(as.matrix(data), ...)

#' @export
hypergraph.netobject <- function(data, ...) .hgm_network(data, ...)

#' @export
hypergraph.cograph_network <- function(data, ...) .hgm_network(data, ...)

#' @export
hypergraph.net_hg_topics <- function(data, ...) group_hypergraph(data, ...)

#' @export
hypergraph.net_mmm <- function(data, ...) group_hypergraph(data, ...)

#' @export
hypergraph.net_clustering <- function(data, ...) group_hypergraph(data, ...)

#' @export
hypergraph.netobject_group <- function(data, ...) group_hypergraph(data, ...)

#' @export
hypergraph.default <- function(data, ...) {
  .thg_bad_input(sprintf(paste0(
    "hypergraph() reads a data frame, a list of sequences, a network, a ",
    "topic model or a clustering of sequences, not an object of class %s"),
    paste(class(data), collapse = "/")))
}
