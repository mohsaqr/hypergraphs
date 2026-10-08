# ---- Shared internal helpers ----
#
# PRIVATE copies of Nestimate internals that hypergraphs' own verbs need
# (`:::` is not allowed and Nestimate does not export them):
# .coerce_sequence_input / .as_netobject / .extract_edges_from_matrix
# (sequence input of hg_bootstrap(), hg_compare(), temporal_hypergraph())
# and .validate_mcml_matrix / .wrap_netobject (pairwise_network()).
# Each is a verbatim copy of Nestimate's same-named internal; identity is
# asserted in local_testing_and_equivalence/test-identity-nestimate-memory.R.
# .coerce_grouped_sequences() is hypergraphs' own (group_hypergraph()).

#' Coerce tna or netobject to labeled sequence data.frame
#'
#' When \code{data} is a \code{tna} or \code{netobject}, extracts the
#' sequence data and converts numeric state IDs to label names. This
#' allows \code{hg_bootstrap()}, \code{hg_compare()} and
#' \code{temporal_hypergraph()} to accept model objects directly.
#'
#' @param data Input: data.frame, list, tna, or netobject.
#' @return A data.frame or list suitable for \code{.hon_parse_input()}.
#' @noRd
.coerce_sequence_input <- function(data) {
  if (inherits(data, "tna")) {
    if (is.null(data$data)) { # nocov start
      stop("tna object has no sequence data ($data). ",
           "Build the tna from sequence data, not a raw matrix.",
           call. = FALSE)
    } # nocov end
    df <- as.data.frame(data$data, stringsAsFactors = FALSE) # nocov start
    lbl <- attr(data$data, "labels") %||% data$labels
    if (!is.null(lbl) && length(lbl) > 0L &&
        (is.integer(df[[1]]) || is.numeric(df[[1]]))) {
      df[] <- lapply(df, function(col) {
        idx <- as.integer(col)
        ifelse(is.na(idx) | idx < 1L | idx > length(lbl),
               NA_character_, lbl[idx])
      })
    }
    return(df) # nocov end
  }
  if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    data <- .as_netobject(data)
  }
  if (inherits(data, "netobject")) {
    if (is.null(data$data)) { # nocov start
      stop("netobject has no sequence data ($data). ",
           "Build the network from sequence data.",
           call. = FALSE) # nocov end
    }
    df <- as.data.frame(data$data, stringsAsFactors = FALSE)
    lbl <- rownames(data$weights)
    if (!is.null(lbl) && length(lbl) > 0L &&
        (is.integer(df[[1]]) || is.numeric(df[[1]]))) {
      df[] <- lapply(df, function(col) { # nocov start
        idx <- as.integer(col)
        ifelse(is.na(idx) | idx < 1L | idx > length(lbl),
               NA_character_, lbl[idx])
      }) # nocov end
    }
    return(df)
  }
  ## Bare sequence matrix (character / logical) -> wide data.frame.
  if (is.matrix(data) && !is.numeric(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }
  data
}

#' Convert pure cograph_network to dual-class netobject/cograph_network
#'
#' Internal converter so that the builders can accept either
#' \code{netobject} or \code{cograph_network} inputs transparently.
#' Objects that already have the \code{"netobject"} class are returned
#' unchanged.
#'
#' @param x A \code{netobject} (returned unchanged) or \code{cograph_network}.
#' @return A dual-class \code{c("netobject", "cograph_network")} object.
#' @noRd
.as_netobject <- function(x) {
  if (inherits(x, "netobject")) return(x)
  if (!inherits(x, "cograph_network")) {
    stop("Expected a netobject or cograph_network.", call. = FALSE)
  }

  mat <- x$weights
  if (is.null(mat)) {
    stop("cograph_network has no $weights matrix.", call. = FALSE)
  }
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (!is.numeric(mat)) storage.mode(mat) <- "double"
  nodes_df <- x$nodes
  states <- nodes_df$label
  raw_data <- x$data
  directed <- x$directed %||% TRUE

  # Infer method from tna metadata or matrix symmetry
  tna_meta <- x$meta$tna
  method <- if (!is.null(tna_meta$method)) {
    tna_meta$method
  } else if (is.matrix(mat) && isSymmetric(mat)) {
    "co_occurrence"
  } else {
    "relative"
  }

  is_sequence_method <- method %in% c(
    "relative", "frequency", "co_occurrence", "attention"
  )

  # Decode integer-encoded tna data -> character labels
  # Only for sequence methods; association methods keep numeric data as-is
  if (!is.null(raw_data)) {
    raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE)
    if (is_sequence_method &&
        (is.integer(raw_data[[1]]) || is.numeric(raw_data[[1]]))) {
      raw_data[] <- lapply(raw_data, function(col) {
        idx <- as.integer(col)
        ifelse(is.na(idx) | idx < 1L | idx > length(states),
               NA_character_, states[idx])
      })
    }
  }

  edges <- .extract_edges_from_matrix(mat, directed = directed)

  structure(list(
    data = raw_data, weights = mat, nodes = nodes_df,
    edges = edges, directed = directed, method = method,
    params = list(), scaling = NULL, threshold = 0,
    n_nodes = length(states), n_edges = nrow(edges),
    level = NULL,
    meta = x$meta %||% list(source = "cograph", layout = NULL,
                            tna = list(method = method)),
    node_groups = x$node_groups
  ), class = c("netobject", "cograph_network"))
}

#' Extract a tidy edge list from a weight matrix
#'
#' @param mat Numeric weight matrix.
#' @param directed Logical. Directed network?
#' @return data.frame with integer \code{from}/\code{to} and \code{weight}.
#' @noRd
.extract_edges_from_matrix <- function(mat, directed = FALSE) {
  if (directed) {
    # Keep self-loops too: every non-zero entry is a real edge.
    idx <- which(mat != 0, arr.ind = TRUE)
  } else {
    # row <= col keeps the upper triangle PLUS the diagonal (one row
    # per self-loop, no double-count for undirected networks).
    idx <- which(mat != 0 & row(mat) <= col(mat), arr.ind = TRUE)
  }

  if (nrow(idx) == 0) {
    return(data.frame(
      from = integer(0), to = integer(0),
      weight = numeric(0), stringsAsFactors = FALSE
    ))
  }

  data.frame(
    from   = as.integer(idx[, 1]),
    to     = as.integer(idx[, 2]),
    weight = mat[idx],
    stringsAsFactors = FALSE
  )
}

# =========================================================================
# Matrix validation + minimal netobject wrapper (hypergraph family)
# =========================================================================

.validate_mcml_matrix <- function(mat) {
  if (any(is.na(mat) | !is.finite(mat))) {
    stop("x matrix must contain finite non-missing weights.", call. = FALSE)
  }
  row_names <- rownames(mat)
  col_names <- colnames(mat)
  has_row_names <- !is.null(row_names)
  has_col_names <- !is.null(col_names)
  if (has_row_names && (any(is.na(row_names)) || any(!nzchar(row_names)))) {
    stop("x matrix row names must not contain missing or empty values.",
         call. = FALSE)
  }
  if (has_col_names && (any(is.na(col_names)) || any(!nzchar(col_names)))) {
    stop("x matrix column names must not contain missing or empty values.",
         call. = FALSE)
  }
  if (has_row_names && anyDuplicated(row_names)) {
    stop("x matrix row names must be unique.", call. = FALSE)
  }
  if (has_col_names && anyDuplicated(col_names)) {
    stop("x matrix column names must be unique.", call. = FALSE)
  }
  if (has_row_names && has_col_names && !identical(row_names, col_names)) {
    stop("x matrix row and column names must be identical and in the same order.",
         call. = FALSE)
  }
  invisible(TRUE)
}



#' Wrap a weight matrix + optional data into a minimal netobject
#' @noRd
.wrap_netobject <- function(mat, data = NULL, method = "relative",
                            directed = TRUE, inits = NULL) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("'mat' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("'mat' must be a square matrix.", call. = FALSE)
  }
  .validate_mcml_matrix(mat)
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)) {
    stop("'method' must be a single non-missing character value.",
         call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }
  states <- rownames(mat)
  if (is.null(states)) {
    states <- as.character(seq_len(nrow(mat)))
    dimnames(mat) <- list(states, states)
  }
  edges <- .extract_edges_from_matrix(mat, directed = directed)
  nodes_df <- data.frame(
    id = seq_along(states), label = states, name = states,
    x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
  )
  if (!is.null(inits)) {
    if (!is.numeric(inits) || length(inits) != length(states) ||
        any(is.na(inits) | !is.finite(inits))) {
      stop("'inits' must be a finite numeric vector with one value per state.",
           call. = FALSE)
    }
    if (!is.null(names(inits))) {
      missing_inits <- setdiff(states, names(inits))
      extra_inits <- setdiff(names(inits), states)
      if (length(missing_inits) > 0L || length(extra_inits) > 0L) {
        stop("'inits' names must match matrix state names.", call. = FALSE)
      }
      inits <- inits[states]
    } else {
      names(inits) <- states
    }
  }

  structure(
    list(
      data       = data,
      weights    = mat,
      inits      = inits,
      nodes      = nodes_df,
      edges      = edges,
      directed   = directed,
      method     = method,
      params     = list(),
      scaling    = NULL,
      threshold  = 0,
      n_nodes    = length(states),
      n_edges    = nrow(edges),
      level      = NULL,
      meta       = list(source = "nestimate", layout = NULL,
                        tna = list(method = method)),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
}

#' Sequences of a Nestimate clustering, split by group
#'
#' Reads the three shapes a Nestimate clustering takes without calling
#' Nestimate: a mixture Markov fit (`net_mmm`) or a distance clustering
#' (`net_clustering`), each with one row of `$data` per sequence and an
#' integer `$assignments`, and the per-group networks built from either
#' (`netobject_group`, a named list of netobjects each with its own `$data`).
#' Group labels are `"Cluster k"` for the first two, as Nestimate's
#' `build_network()` names them, and the list names for a
#' `netobject_group` (so `rename_models()` names carry through).
#'
#' @param x A `net_mmm`, `net_clustering` or `netobject_group`.
#' @return A named list, one element per group in group order, each a list
#'   of character vectors (one per sequence, NA cells dropped).
#' @noRd
.coerce_grouped_sequences <- function(x) {
  rows_of <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (nrow(df) == 0L) return(list())
    cells <- do.call(cbind, lapply(df, as.character))
    keep <- !is.na(cells) & nzchar(cells)
    # column-major order keeps each row's cells in sequence order
    unname(split(cells[keep], factor(row(cells)[keep], levels = seq_len(nrow(df)))))
  }
  if (inherits(x, "netobject_group")) {
    labels <- names(x)
    if (!length(x) || is.null(labels) || anyNA(labels) || any(!nzchar(labels)) ||
        anyDuplicated(labels)) {
      .thg_bad_input("a netobject_group needs unique, non-empty group names")
    }
    ok <- vapply(x, \(net) inherits(net, c("netobject", "cograph_network")) &&
                   !is.null(net$data), logical(1L))
    if (!all(ok)) {
      .thg_bad_input(sprintf(
        "every network of the netobject_group needs its sequence data ($data); missing in: %s",
        paste(labels[!ok], collapse = ", ")))
    }
    return(stats::setNames(
      lapply(unclass(x), \(net) rows_of(.coerce_sequence_input(net))), labels))
  }
  if (!inherits(x, c("net_mmm", "net_clustering"))) {
    .thg_bad_input("expected a net_mmm, net_clustering or netobject_group")
  }
  if (is.null(x$data) || is.null(x$assignments)) {
    .thg_bad_input(sprintf("this %s carries no $data or no $assignments",
                           class(x)[1L]))
  }
  data <- as.data.frame(x$data, stringsAsFactors = FALSE)
  assignments <- x$assignments
  k <- x$k %||% max(assignments)
  if (length(assignments) != nrow(data) || anyNA(assignments) ||
      any(assignments < 1 | assignments > k | assignments != round(assignments))) {
    .thg_bad_input(sprintf(
      "`assignments` must be one whole number in 1..%d per sequence of $data (%d sequences)",
      k, nrow(data)))
  }
  # a mixture fit keeps its component networks, whose weight dimnames decode
  # integer-coded state cells as .coerce_sequence_input() does for a netobject
  labels <- if (inherits(x, "net_mmm") && length(x$models)) {
    rownames(x$models[[1L]]$weights)
  }
  if (length(labels) && ncol(data) &&
      all(vapply(data, is.numeric, logical(1L)))) {
    data[] <- lapply(data, \(col) {
      idx <- as.integer(col)
      ifelse(is.na(idx) | idx < 1L | idx > length(labels), NA_character_,
             labels[idx])
    })
  }
  sequences <- rows_of(data)
  groups <- lapply(seq_len(k), \(g) sequences[assignments == g])
  stats::setNames(groups, paste("Cluster", seq_len(k)))
}

# =========================================================================
# One sequence-input helper for every sequence-taking verb
# =========================================================================
#
# hon(), mogen(), hypa(), markov_order(), memory(),
# hg_markov_stability(), hg_bootstrap(), hg_compare(), simplicial(type =
# "window") and window_hypergraph() all read sequences through
# .ho_sequence_input(). It accepts
#   * a long event table, with `action` (the state column) and optionally
#     `actor`, `session`, `time`, `time_threshold` and `timezone`;
#   * a wide data.frame or character matrix, one sequence per row;
#   * a list of vectors, one sequence per element;
#   * a model object carrying its sequences (netobject, netobject_group,
#     tna, cograph_network).
# The long route is Nestimate's: build_network()'s detection of the columns
# named action, time and session / session_id when those arguments are
# NULL, then Nestimate::prepare() with the same arguments, whose wide
# $sequence_data is returned unchanged. The sequences are therefore the
# ones build_network(method = "relative", ...) builds from the same call.
# `session = FALSE` switches session detection off. A long table passed
# without `action` and without an `action` column is refused
# (hypergraphs_long_format).

#' Coerce any sequence input to what a sequence verb consumes
#'
#' @param data Sequence input (see the block comment above).
#' @param action,actor,time,session Long-format column names or `NULL`;
#'   `session = FALSE` switches session detection off.
#' @param time_threshold,timezone Passed to [Nestimate::prepare()].
#' @param lists `"keep"` returns a list of sequences as a list; `"wide"`
#'   pads it into a wide data.frame (one row per sequence, `NA`-padded),
#'   for consumers that read only frames.
#' @param models `"keep"` passes model objects through untouched (the
#'   consumer reads them itself); `"decode"` turns a `tna` or a bare
#'   `cograph_network` into its wide, label-decoded sequence frame, for
#'   consumers that only read a netobject's `$data`.
#' @return The wide sequence data.frame of [Nestimate::prepare()] (long
#'   input), a wide data.frame, a list, or the model object.
#' @noRd
.ho_sequence_input <- function(data, action = NULL, actor = NULL,
                               time = NULL, session = NULL,
                               time_threshold = 900, timezone = "UTC",
                               lists = c("keep", "wide"),
                               models = c("keep", "decode")) {
  lists <- match.arg(lists)
  models <- match.arg(models)
  # build_network()'s column detection: a column named action, time,
  # session or session_id (any case) is used when its argument is NULL
  # (time and session matter only for a long table, so they are looked
  # for once there is an action column)
  if (is.data.frame(data)) {
    detected <- .ho_detect_columns(data)
    action <- action %||% detected$action
    if (!is.null(action)) {
      time <- time %||% detected$time
      if (is.null(session)) session <- detected$session
    }
  }
  if (isFALSE(session)) session <- NULL
  if (is.null(action)) {
    stray <- c(actor = !is.null(actor), time = !is.null(time),
               session = !is.null(session))
    if (any(stray)) {
      .ho_bad_input(sprintf(
        "`%s` requires `action` (the state column of a long event table)",
        names(stray)[stray][1L]))
    }
    .hon_guard_long_format(data)
    if (models == "decode" &&
        (inherits(data, "tna") ||
         (inherits(data, "cograph_network") && !inherits(data, "netobject")))) {
      data <- .coerce_sequence_input(data)
    }
    if (is.matrix(data) && !is.numeric(data)) {
      data <- as.data.frame(data, stringsAsFactors = FALSE)
    }
    if (lists == "wide" && is.list(data) && !is.data.frame(data) &&
        !inherits(data, c("netobject", "netobject_group", "tna",
                          "cograph_network"))) {
      data <- .ho_wide_sequences(data)
    }
    return(data)
  }
  .ho_check_long_columns(data, action = action, actor = actor, time = time,
                         session = session)
  stopifnot(
    "`time_threshold` must be a positive number, Inf, or FALSE" =
      isFALSE(time_threshold) ||
      (is.numeric(time_threshold) && length(time_threshold) == 1L &&
         !is.na(time_threshold) && time_threshold > 0)
  )
  keys <- c(actor, session)
  if (length(keys) > 0L && anyNA(data[keys])) {
    .ho_bad_input(sprintf(paste0(
      "missing values in %s: every event needs an identifier; drop or ",
      "relabel these rows first"), paste0("`", keys, "`", collapse = ", ")))
  }
  if (is.null(actor)) {
    # build_network()'s notice: without an actor every event is one sequence
    notice <- simpleMessage(paste0(
      "A network with one long sequence is not recommended and can't be ",
      "validated using bootstrap and other confirmatory testings.\n"))
    class(notice) <- c("hypergraphs_single_sequence", class(notice))
    message(notice)
  }
  # only the columns that define the sequences, so prepare() has no other
  # columns to aggregate into session metadata
  columns <- unique(c(action, actor, time, session))
  prepared <- Nestimate::prepare(data[columns], actor = actor,
                                 action = action, time = time,
                                 session = session,
                                 time_threshold = time_threshold,
                                 timezone = timezone)
  prepared$sequence_data
}

#' Columns found by name, as build_network() finds them
#'
#' @param data A data.frame.
#' @return A list with `action`, `time` and `session`: the column whose
#'   lower-cased name is `"action"`, `"time"`, and `"session"` (else
#'   `"session_id"`), each `NULL` when there is no single such column.
#' @noRd
.ho_detect_columns <- function(data) {
  lower <- tolower(names(data))
  match1 <- \(name) {
    hit <- which(lower == name)
    if (length(hit) == 1L) names(data)[hit]
  }
  list(action = match1("action"), time = match1("time"),
       session = match1("session") %||% match1("session_id"))
}

#' Validate the long-format column arguments against a data.frame
#'
#' @param data The long table.
#' @param action,actor,time,session Column names or `NULL`; `actor` and
#'   `session` may name several columns.
#' @return `NULL`, invisibly; raises `hypergraphs_bad_input`.
#' @noRd
.ho_check_long_columns <- function(data, action, actor, time, session) {
  if (!is.data.frame(data)) {
    .ho_bad_input("`data` must be a data.frame when `action` is given")
  }
  single <- list(action = action, time = time)
  several <- list(actor = actor, session = session)
  lapply(names(single), function(arg) {
    col <- single[[arg]]
    if (is.null(col)) return(NULL)
    if (!is.character(col) || length(col) != 1L || !col %in% names(data)) {
      .ho_bad_input(sprintf("`%s` must name a column of `data`", arg))
    }
    NULL
  })
  lapply(names(several), function(arg) {
    col <- several[[arg]]
    if (is.null(col)) return(NULL)
    if (!is.character(col) || length(col) < 1L || !all(col %in% names(data))) {
      .ho_bad_input(sprintf("`%s` must name one or more columns of `data`",
                            arg))
    }
    NULL
  })
  invisible(NULL)
}

#' Pad a list of sequences into a wide data.frame
#'
#' @param sequences List of vectors.
#' @return data.frame, one row per sequence, columns `T1`..`Tn`, `NA`-padded
#'   on the right.
#' @noRd
.ho_wide_sequences <- function(sequences) {
  if (!length(sequences)) {
    .ho_bad_input("`data` holds no sequence: the list of sequences is empty")
  }
  atomic <- vapply(sequences, \(s) is.atomic(s) || is.factor(s), logical(1L))
  if (!all(atomic)) {
    .ho_bad_input(sprintf(paste0(
      "every sequence must be a vector of states; element %d is a %s"),
      which(!atomic)[1L], class(sequences[[which(!atomic)[1L]]])[1L]))
  }
  sequences <- lapply(sequences, as.character)
  if (all(lengths(sequences) == 0L)) {
    .ho_bad_input("`data` holds no state: every sequence is empty")
  }
  width <- max(lengths(sequences))
  rows <- lapply(sequences, function(s) {
    c(s, rep(NA_character_, width - length(s)))
  })
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- paste0("T", seq_len(width))
  rownames(out) <- NULL
  out
}

#' Raise a classed bad-input error
#' @param msg Message.
#' @noRd
.ho_bad_input <- function(msg) {
  stop(errorCondition(msg, class = "hypergraphs_bad_input", call = NULL))
}

# Apply `fn` to every element of `x`, serially or with parallel::mclapply
# (forking, so never on Windows, where it runs serially). Work that seeds
# itself per element returns the same result either way, so `fn` must not
# return NULL and must not draw from the caller's stream. A worker whose R
# code fails comes back as a "try-error" value and is re-raised as
# `hypergraphs_parallel_failed`. A worker the operating system killed comes
# back as NULL -- on macOS this happens when a forked child calls the
# Accelerate BLAS after the parent has (it is not fork-safe) -- and is rerun
# here, serially, with one `hypergraphs_parallel_fallback` warning: the work
# is deterministic, so the rerun gives the value the worker would have.
.ho_apply <- function(x, fn, parallel = FALSE, n_cores = 2L) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    .ho_input_error("`parallel` must be TRUE or FALSE")
  }
  n_cores <- .ho_check_count(n_cores, "n_cores")
  if (!parallel || .Platform$OS.type == "windows") return(lapply(x, fn))
  # mclapply also warns "scheduled core(s) ... encountered errors" or "did
  # not deliver results"; both outcomes are handled below, so only those
  # warnings are muffled here.
  out <- withCallingHandlers(
    parallel::mclapply(x, fn, mc.cores = n_cores),
    warning = function(w) {
      if (grepl("scheduled cores?", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  errored <- vapply(out, inherits, logical(1L), "try-error")
  if (any(errored)) {
    stop(errorCondition(
      sprintf("%d of %d parallel workers failed; first error: %s",
              sum(errored), length(out),
              conditionMessage(attr(out[[which(errored)[1L]]], "condition"))),
      class = "hypergraphs_parallel_failed", call = NULL
    ))
  }
  killed <- vapply(out, is.null, logical(1L))
  if (any(killed)) {
    warning(warningCondition(sprintf(paste0(
      "%d of %d parallel workers were killed by the operating system and ",
      "were rerun serially (on macOS the Accelerate BLAS is not fork-safe)"),
      sum(killed), length(out)),
      class = "hypergraphs_parallel_fallback", call = NULL))
    out[killed] <- lapply(x[killed], fn)
  }
  out
}

# k-means with `nstart` random starts, drawn exactly as stats::kmeans(x, k,
# nstart = nstart) draws them (one start: k rows of `x`, redrawn from the
# unique rows if two coincide; several: every start from the unique rows,
# in order, from one stream), keeping the first start with the smallest
# total within-cluster sum. The difference: a Hartigan-Wong start that
# stops before it has finished (ifault 4, "Quick-TRANSfer stage steps
# exceeded", or ifault 2, `iter.max` reached) is finished instead of
# entering the comparison unfinished: Hartigan-Wong is continued from the
# centres it reached while that lowers the within-cluster sum (at most
# `max_continue` times), and a start that still cycles is finished by
# Lloyd's algorithm from those centres (Lloyd 1982), whose steps never
# raise the same objective. Starts that finish give exactly the
# stats::kmeans() result. A start Lloyd cannot finish either raises one
# `hypergraphs_no_converge` warning.
.ho_kmeans <- function(x, k, nstart = 1L, iter.max = 100L,
                       max_continue = 20L, lloyd_iter = 1000L) {
  first <- if (nstart == 1L) x[sample.int(nrow(x), k), , drop = FALSE]
  need_unique <- is.null(first) || anyDuplicated(first) > 0L
  unique_rows <- if (need_unique) unique(x)
  draw_unique <- \() unique_rows[sample.int(nrow(unique_rows), k), , drop = FALSE]
  if (need_unique) first <- draw_unique()
  starts <- c(list(first), lapply(seq_len(nstart - 1L), \(s) draw_unique()))
  runs <- lapply(starts, \(centres) {
    run <- .ho_kmeans_once(x, centres, iter.max)
    continued <- 0L
    # each continuation restarts Hartigan-Wong from the centres reached,
    # while it still lowers the within-cluster sum; bounded by `max_continue`
    while (!run$finished && continued < max_continue) {
      nxt <- .ho_kmeans_once(x, run$fit$centers, iter.max)
      improved <- nxt$fit$tot.withinss < run$fit$tot.withinss
      run <- nxt
      continued <- continued + 1L
      if (!improved) break
    }
    if (!run$finished) {
      run <- .ho_kmeans_once(x, run$fit$centers, iter.max = lloyd_iter,
                             algorithm = "Lloyd")
    }
    run
  })
  unfinished <- !vapply(runs, \(r) r$finished, logical(1L))
  if (any(unfinished)) {
    warning(warningCondition(sprintf(
      "k-means: %d of %d starts did not converge, even by Lloyd's algorithm",
      sum(unfinished), nstart),
      class = "hypergraphs_no_converge", call = NULL))
  }
  # which.min() keeps the first of tied minima, as stats::kmeans() does
  totals <- vapply(runs, \(r) r$fit$tot.withinss, numeric(1L))
  runs[[which.min(totals)]]$fit
}

# One k-means run from given centres, as list(fit, finished). The two
# "stopped early" warnings (Hartigan-Wong's ifault 4 "Quick-TRANSfer stage
# steps exceeded", and "did not converge in ... iterations" from either
# algorithm) are what .ho_kmeans() acts on, so they are recorded in
# `finished` and muffled; any other warning passes through.
.ho_kmeans_once <- function(x, centres, iter.max,
                            algorithm = "Hartigan-Wong") {
  stopped_early <- FALSE
  fit <- withCallingHandlers(
    stats::kmeans(x, centers = centres, iter.max = iter.max,
                  algorithm = algorithm),
    warning = function(w) {
      if (grepl("^(Quick-TRANSfer stage steps exceeded|did not converge in)",
                conditionMessage(w))) {
        stopped_early <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  list(fit = fit, finished = !stopped_early)
}
