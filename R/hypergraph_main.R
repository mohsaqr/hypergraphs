# ---- hypergraph(): the main constructor ------------------------------------
# One verb for every input. Each method hands its input to the constructor of
# that kind of data and returns that constructor's result unchanged.

#' Build a hypergraph from any input
#'
#' `hypergraph()` is the main constructor of the package. It reads the
#' input and builds the hypergraph the input describes, through the
#' constructor of that kind of data. A data frame is read in one of three
#' formats, each with its own argument names, and the format is recognised
#' from the arguments given:
#'
#' * **Membership data** name a `node` and a `hyperedge` (or `from` and `to`
#'   for an edge list): every value of `hyperedge` becomes a hyperedge of the
#'   nodes it holds ([group_hypergraph()]). With `time`, `start` or `end` the
#'   hyperedges carry a clock ([temporal_hypergraph()]).
#' * **Event data** name an `action`, the column of what happened, with the
#'   `session` and `actor` it belongs to, in the vocabulary of the memory
#'   family. Without `window`, the actions of each session become one
#'   hyperedge; a session is read within its actor when both are given, and
#'   with `actor` alone each actor's actions become one hyperedge. With
#'   `window`, every window of consecutive actions within an actor becomes a
#'   hyperedge ([window_hypergraph()]), ordered by `time`.
#' * With `group`, `top` or `min_share`, the sets of either format are
#'   counted and the most frequent become hyperedges; `group` names the
#'   comparison variable within whose values the sets are counted.
#'
#' A list of sequences is read as event data with `window`. A network, as a
#' weight matrix, a sparse matrix, a `netobject` or a `cograph_network`, has
#' its cliques promoted to hyperedges ([network_hypergraph()]); with
#' `window`, a model object built from sequences is read as sequences
#' instead. A topic model fitted by [hg_topics()], or a clustering of
#' sequences, becomes the hypergraph of its frequent sets
#' ([group_hypergraph()]).
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
#' naming the class; for event data given a clock without `window`, since the
#' set of a session has no order; and the conditions of the constructor
#' called.
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
#' hypergraph(meetings, node = "person", hyperedge = "meeting")
#'
#' visits <- data.frame(
#'   user = c("u1", "u1", "u1", "u1", "u2", "u2"),
#'   visit = c(1, 1, 2, 2, 1, 1),
#'   page = c("home", "cart", "home", "help", "home", "cart"))
#' hypergraph(visits, action = "page", actor = "user", session = "visit")
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
.hgm_window_args <- c("window", "step", "time_threshold", "min_weight")
.hgm_clock_args <- c("time", "start", "end", "observation_start",
                     "observation_end", "time_unit")

#' @export
hypergraph.data.frame <- function(data, ...) {
  given <- names(list(...))
  if (any(given %in% .hgm_window_args)) {
    return(window_hypergraph(data, ...))
  }
  if ("action" %in% given) {
    return(.hgm_event_sets(data, ...))
  }
  if (any(given %in% .hgm_clock_args)) {
    return(temporal_hypergraph(data, ...))
  }
  group_hypergraph(data, ...)
}

# Event data without a window: the actions of each session (within its actor,
# when both are given), or of each actor, form one hyperedge, built as
# membership data with the action as the node.
.hgm_event_sets <- function(data, action, actor = NULL, session = NULL, ...) {
  dots <- names(list(...))
  ordered <- intersect(dots, .hgm_clock_args)
  if (length(ordered)) {
    .thg_bad_input(sprintf(paste0(
      "the actions of a session form an unordered set, so %s does not ",
      "apply; give `window` for windows of consecutive actions"),
      paste0("`", ordered, "`", collapse = ", ")))
  }
  for (column in c(action, actor, session)) {
    if (!is.character(column) || length(column) != 1L ||
        !column %in% names(data)) {
      .thg_bad_input(paste0("`action`, `actor` and `session` must each name ",
                            "one column of `data`"))
    }
  }
  if (is.null(actor) && is.null(session)) {
    .thg_bad_input(paste0("event data need a `session` or an `actor` whose ",
                          "actions form a hyperedge, or `window` for windows"))
  }
  unit <- if (!is.null(actor) && !is.null(session)) {
    paste(data[[actor]], data[[session]], sep = ".")
  } else {
    as.character(data[[session %||% actor]])
  }
  hyperedge <- session %||% actor
  data[[hyperedge]] <- unit
  group_hypergraph(data, node = action, hyperedge = hyperedge, ...)
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
