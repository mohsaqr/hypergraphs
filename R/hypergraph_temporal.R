# Temporal hypergraphs as hyperedge spells. This is deliberately an hypergraphs
# object, not a Dynet adapter: hyperedges remain first-class throughout. The
# vocabulary is Dynet's (`../temporal`): column aliases, time parsing,
# observation bounds and the step/window measurement grid follow its
# definitions so a co-presence log reads the same way in both packages.

.thg_bad_input <- function(message, class = character()) {
  stop(errorCondition(message, class = c(class, "hypergraphs_bad_input"), call = NULL))
}

.thg_deprecated <- function(old, new, fn) {
  warning(warningCondition(
    sprintf("`%s` is deprecated in %s(); use `%s` instead", old, fn, new),
    class = c("hypergraphs_deprecated", "deprecatedWarning"), call = NULL
  ))
}

# ---------------------------------------------------------------------------
# Column detection: Dynet's alias table (R/detect.R), checked case-
# insensitively after stripping non-alphanumerics; order is priority order.
# ---------------------------------------------------------------------------

.thg_aliases <- list(
  from     = c("from", "source", "sender", "tail", "ego", "actor1", "node1",
               "sourceid", "fromid", "i"),
  to       = c("to", "target", "receiver", "head", "alter", "actor2", "node2",
               "targetid", "toid", "j"),
  start    = c("start", "onset", "begin", "starttime", "startdate", "tstart",
               "fromtime"),
  end      = c("end", "terminus", "finish", "endtime", "enddate", "tend",
               "totime", "stop"),
  duration = c("duration", "dur", "elapsed", "length", "lasted"),
  time     = c("time", "timestamp", "datetime", "date", "when", "t", "occurred",
               "createdat", "posttime"),
  thread   = c("thread", "threadid", "discussion", "discussiontitle", "topic",
               "conversation", "parent", "postid", "root"),
  session  = c("session", "period", "wave", "phase", "cohort", "course"),
  group    = c("group", "event", "context", "room", "class", "meeting",
               "venue", "team", "channel", "hyperedge"),
  actor    = c("actor", "person", "member", "student", "participant", "user",
               "id", "name", "node"),
  weight   = c("weight", "weights", "strength")
)

.thg_norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

.thg_match_column <- function(data, role, exclude = character()) {
  aliases <- .thg_aliases[[role]]
  cand <- setdiff(names(data), exclude)
  hit <- match(aliases, .thg_norm_name(cand))
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0L) NULL else cand[hit[1L]]
}

# One column: an explicit name must exist as given; otherwise the first alias
# match, or NULL.
.thg_resolve_column <- function(data, given, role, exclude = character(),
                                arg = role) {
  if (!is.null(given)) {
    if (!is.character(given) || length(given) != 1L || is.na(given)) {
      .thg_bad_input(sprintf("`%s` must be a single column name", arg))
    }
    if (!given %in% names(data)) {
      .thg_bad_input(
        sprintf("column `%s` (given as `%s`) is not in the data; available: %s",
                given, arg, paste(names(data), collapse = ", ")),
        class = "hypergraphs_missing_column"
      )
    }
    return(given)
  }
  .thg_match_column(data, role, exclude)
}

# ---------------------------------------------------------------------------
# Time parsing: Dynet's rule. Numeric times pass through with unit "step";
# Date, POSIXct and character date-times become elapsed time since an origin
# in a unit chosen for the span (or requested).
# ---------------------------------------------------------------------------

.thg_time_formats <- c(
  "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M",
  "%Y-%m-%dT%H:%M", "%Y-%m-%d",
  "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M", "%Y/%m/%d",
  "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y",
  "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M", "%m/%d/%Y",
  "%d-%m-%Y %H:%M:%S", "%d-%m-%Y %H:%M", "%d-%m-%Y",
  "%d %b %Y %H:%M", "%d %B %Y %H:%M", "%d %b %Y", "%d %B %Y"
)

# The whole-string pattern of a format: strptime() reads a prefix and
# ignores what follows, so a string is accepted only when it matches this
# pattern from end to end.
.thg_format_pattern <- function(fmt) {
  tokens <- c("%Y" = "[0-9]{4}", "%m" = "[0-9]{1,2}", "%d" = "[0-9]{1,2}",
              "%H" = "[0-9]{1,2}", "%M" = "[0-9]{1,2}", "%S" = "[0-9]{1,2}",
              "%b" = "[[:alpha:]]+", "%B" = "[[:alpha:]]+")
  pattern <- Reduce(function(p, token) gsub(token, tokens[[token]], p, fixed = TRUE),
                    names(tokens), fmt)
  paste0("^", pattern, "$")
}

# The date part of a format: formats that share it read one way of writing a
# date, with or without a time of day.
.thg_format_family <- function(fmt) sub("[ T]%H.*$", "", fmt)

# An ISO 8601 ending of a date-time: fractional seconds and a UTC designator
# (`Z`) or offset (`+0200`, `+02:00`, `-05`... in hours and minutes). Returns
# the string without them and the seconds to add to its reading as UTC.
.thg_iso_suffix <- function(x) {
  parts <- regmatches(x, regexec(
    "^(.*?[0-9]:[0-9]{2}(?::[0-9]{2})?)([.][0-9]+)?(Z|[+-][0-9]{2}:?[0-9]{2})?$",
    x, perl = TRUE))
  shift <- vapply(parts, function(p) {
    if (!length(p)) return(0)
    fraction <- if (nzchar(p[[3L]])) as.numeric(p[[3L]]) else 0
    zone <- p[[4L]]
    offset <- if (!nzchar(zone) || identical(zone, "Z")) 0 else {
      digits <- gsub(":", "", substring(zone, 2L), fixed = TRUE)
      (if (startsWith(zone, "-")) -1 else 1) *
        (as.numeric(substr(digits, 1L, 2L)) * 3600 +
           as.numeric(substr(digits, 3L, 4L)) * 60)
    }
    fraction - offset
  }, numeric(1L))
  # a fraction belongs to the seconds: after hours and minutes alone it is
  # left in place, and the string then fails to parse
  bare <- vapply(parts, function(p) {
    if (!length(p)) return(NA_character_)
    if (nzchar(p[[3L]]) && !grepl(":[0-9]{2}:[0-9]{2}$", p[[2L]])) {
      return(NA_character_)
    }
    p[[2L]]
  }, character(1L))
  matched <- !is.na(bare)
  list(base = ifelse(matched, bare, x), shift = ifelse(matched, shift, 0))
}

# Character date-times on one reading. Every non-missing string must be read
# completely, by the formats of one date family (the family that reads the
# most strings, the earlier one on a tie): a string read only in part, or by
# no format of that family, is an error, never a missing time. Strings
# without an offset are UTC; an ISO offset or `Z` is honoured.
.thg_parse_datetime_strings <- function(x) {
  present <- !is.na(x)
  text <- trimws(x)
  iso <- .thg_iso_suffix(text)
  families <- .thg_format_family(.thg_time_formats)
  read <- lapply(.thg_time_formats, function(fmt) {
    whole <- present & grepl(.thg_format_pattern(fmt), iso$base)
    out <- rep(as.POSIXct(NA_real_, tz = "UTC"), length(x))
    if (any(whole)) {
      out[whole] <- as.POSIXct(iso$base[whole], format = fmt, tz = "UTC")
    }
    out
  })
  family_names <- unique(families)
  n_ok <- vapply(family_names, function(f) {
    ok <- Reduce(`|`, lapply(read[families == f], \(p) !is.na(p)),
                 rep(FALSE, length(x)))
    sum(ok)
  }, integer(1L))
  best <- family_names[which.max(n_ok)]
  # within the family, the first format that reads a string reads it
  out <- Reduce(function(acc, p) {
    fill <- is.na(acc) & !is.na(p)
    acc[fill] <- p[fill]
    acc
  }, read[families == best], rep(as.POSIXct(NA_real_, tz = "UTC"), length(x)))
  failed <- present & is.na(out)
  if (any(failed)) {
    .thg_bad_input(
      sprintf(paste0("could not read %d time string(s) such as '%s' as a date ",
                     "or date-time; supply numeric, Date or POSIXct times"),
              sum(failed), x[failed][1L]),
      class = "hypergraphs_unparsed_time"
    )
  }
  out + iso$shift
}

.thg_unit_seconds <- function(time_unit) {
  switch(time_unit, seconds = 1, minutes = 60, hours = 3600, days = 86400,
         weeks = 604800)
}

.thg_auto_unit <- function(span_seconds) {
  if (!is.finite(span_seconds) || span_seconds <= 0) return("seconds")
  if (span_seconds < 2 * 60) "seconds"
  else if (span_seconds < 3 * 3600) "minutes"
  else if (span_seconds < 3 * 86400) "hours"
  else "days"
}

# Parse one or more time vectors on a shared clock. `columns` is a named
# list; NA entries are allowed (an open end, a node without an entry time)
# and stay NA. Returns the parsed columns plus `unit` and `origin`.
.thg_parse_clock <- function(columns, time_unit = "auto") {
  if (!is.character(time_unit) || length(time_unit) != 1L ||
      !time_unit %in% c("auto", "seconds", "minutes", "hours", "days", "weeks")) {
    .thg_bad_input("`time_unit` must be \"auto\", \"seconds\", \"minutes\", \"hours\", \"days\" or \"weeks\"")
  }
  columns <- lapply(columns, function(x) if (is.factor(x)) as.character(x) else x)
  # a column with no value at all (an `end` column of open intervals, read
  # by data.frame() as logical NA) takes the clock of the others
  empty <- vapply(columns, function(x) all(is.na(x)), logical(1L))
  if (all(empty)) .thg_bad_input("no time value could be read")
  if (any(empty)) {
    read <- .thg_parse_clock(columns[!empty], time_unit)
    blank <- lapply(columns[empty], function(x) rep(NA_real_, length(x)))
    return(c(read[names(columns)[!empty]], blank,
             read[c("unit", "origin")])[c(names(columns), "unit", "origin")])
  }
  kinds <- vapply(columns, function(x) {
    if (is.numeric(x)) "numeric"
    else if (inherits(x, "POSIXct") || inherits(x, "Date")) "calendar"
    else if (is.character(x)) "character"
    else "unsupported"
  }, character(1L))
  if (any(kinds == "unsupported")) {
    bad <- names(columns)[kinds == "unsupported"][1L]
    .thg_bad_input(sprintf(
      "time column `%s` has unsupported type <%s>; supply numeric, Date, POSIXct or character times",
      bad, paste(class(columns[[bad]]), collapse = "/")))
  }
  if (all(kinds == "numeric")) {
    # NA is a missing time (an open end); NaN and infinite values are not
    # times at all
    bad <- vapply(columns, function(x) any(is.nan(x) | is.infinite(x)), logical(1L))
    if (any(bad)) {
      .thg_bad_input(sprintf(
        "the %s times hold NaN or infinite values; times must be finite",
        sub("_", " ", names(columns)[bad][1L], fixed = TRUE)))
    }
    return(c(lapply(columns, as.numeric), list(unit = "step", origin = 0)))
  }
  if (any(kinds == "numeric")) {
    .thg_bad_input("time columns must all be dates or all be numbers; they cannot be mixed")
  }
  parsed <- lapply(columns, function(x) {
    if (is.character(x)) .thg_parse_datetime_strings(x) else as.POSIXct(x, tz = "UTC")
  })
  all_times <- do.call(c, lapply(parsed, function(p) p[!is.na(p)]))
  if (length(all_times) == 0L) .thg_bad_input("no time value could be read")
  span <- as.numeric(difftime(max(all_times), min(all_times), units = "secs"))
  unit <- if (identical(time_unit, "auto")) .thg_auto_unit(span) else time_unit
  origin <- min(all_times)
  values <- lapply(parsed, function(p) {
    as.numeric(difftime(p, origin, units = "secs")) / .thg_unit_seconds(unit)
  })
  c(values, list(unit = unit, origin = origin))
}

# A user-supplied time (`at`, `start`, `observation_end`, ...) on the stored
# clock: numbers pass through, dates are converted for a calendar hypergraph
# and refused for a numeric one. Vectorised.
.thg_as_time <- function(v, x, arg) {
  if (is.null(v)) return(NULL)
  if (is.factor(v)) v <- as.character(v)
  calendar <- inherits(x$origin, "POSIXt")
  if (inherits(v, "Date") || inherits(v, "POSIXt") || is.character(v)) {
    if (!calendar) {
      .thg_bad_input(sprintf(
        "`%s` was given as a date, but this hypergraph's times are plain numbers",
        arg))
    }
    parsed <- if (is.character(v)) .thg_parse_datetime_strings(v) else
      as.POSIXct(v, tz = "UTC")
    if (anyNA(parsed)) .thg_bad_input(sprintf("`%s` contains missing times", arg))
    return(as.numeric(difftime(parsed, x$origin, units = "secs")) /
             .thg_unit_seconds(x$time_unit))
  }
  if (!is.numeric(v) || anyNA(v) || any(!is.finite(v))) {
    .thg_bad_input(sprintf("`%s` must be finite numbers or dates", arg))
  }
  as.numeric(v)
}

# The calendar value of stored times, for axes and labels; numeric clocks
# return the numbers unchanged.
.thg_calendar <- function(t, origin, time_unit) {
  if (!inherits(origin, "POSIXt")) return(t)
  out <- origin + t * .thg_unit_seconds(time_unit)
  if (identical(time_unit, "days") || identical(time_unit, "weeks")) as.Date(out) else out
}

# The time of snapshot `i` of a hg_snapshots() list as a reader gives it: a
# date (or date-time) for a calendar hypergraph, the number on the clock
# otherwise.
.thg_snapshot_time <- function(snaps, i) {
  .thg_calendar(snaps[[i]]$params$at, attr(snaps, "origin"),
                attr(snaps, "time_unit"))
}

.thg_clock_label <- function(origin, time_unit) {
  if (!inherits(origin, "POSIXt")) return("time (numeric steps)")
  sprintf("%s since %s", time_unit, format(origin, if (identical(time_unit, "days") ||
                                                          identical(time_unit, "weeks"))
    "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"))
}

#' Temporal hypergraph from an edge list or co-occurrence data
#'
#' Builds a temporal hypergraph the way a temporal network is defined from
#' data in this ecosystem (Dynet's `dynet()`): one row per relation, with the
#' columns that name its ends and its clock. Two input shapes are accepted.
#' An **edge list** names `from` and `to`, and every row is a hyperedge of
#' size two, so an ordinary temporal network is the same object.
#' **Membership data** name a `node` and a `hyperedge`: every node sharing
#' one value of `hyperedge` (a case, a citation block, a seminar) belongs to one
#' hyperedge. Two clocks are understood, and the one you name selects the
#' format:
#'
#' \describe{
#'   \item{interval}{`start` and `end`: the hyperedge is active on the closed
#'     interval between them. A missing `end` means it stays active through
#'     the end of observation.}
#'   \item{contact}{`time`: the hyperedge is an instantaneous event, as in a
#'     contact log -- a citation, a message, a meeting on one day. It is
#'     present at that instant only. The cumulative view in which every
#'     hyperedge stays once it has appeared (the point-aggregation model of
#'     Coupette et al. 2024) is a snapshot `mode`, not a property of the
#'     data: ask for it with `mode = "cumulative"` in [hg_snapshot()],
#'     [hg_growth()] and [hg_edges()].}
#' }
#'
#' A third shape is a **sequence table**: one row per session and one column
#' per position holding the state at that step (the wide format of tna and
#' TraMineR), a list of character vectors, or a `tna` / `netobject` model
#' built from one. It is recognised when no relational column is named or
#' detected and every column is categorical. Each session is one hyperedge
#' whose members are the states it contains, and a state's membership is a
#' contact at its position, `1` to the session's length, on a `"step"`
#' clock: the simple co-occurrence reading in which time is order.
#'
#' In any shape a membership may carry its own time, as when a log has one
#' row per attendance rather than one time per group. The hyperedge then
#' spans from its first to its last membership, each membership is present
#' on its own spell only, and a snapshot keeps the memberships present in
#' its window: `mode = "cumulative"` at step `t` is what each session had
#' shown by `t`, and `window = 3` at `t` is what it showed on steps `t` to
#' `t + 2`. When every membership carries its hyperedge's time, as a
#' tribunal or a citation block does, nothing changes.
#'
#' Column names are resolved case-insensitively from the same alias table
#' Dynet uses, so `Sender`/`Receiver`, `source`/`target`, `onset`/`terminus`
#' and `timestamp` are understood without being spelled out; a name you give
#' explicitly must exist as written. Times may be numeric, `Date`, `POSIXct`
#' or character date-time strings. Numeric times are kept as they are and
#' reported in `"step"` units; calendar input is converted to elapsed time
#' since the earliest time in the data, in a unit chosen for the span
#' (`"days"` beyond three days) or given as `time_unit`, and the unit and
#' origin are stored in the object and shown by `print()`. Every time you
#' pass later -- `at`, `start`, `end`, the observation bounds -- may be a
#' date for a calendar hypergraph or a number on that clock, and every
#' `time` column in a result is on that clock.
#'
#' `observation_start` and `observation_end` declare the study window when it
#' is known independently of the log, with Dynet's meaning: they bound the
#' snapshot times and the measurement grid, and an open-ended hyperedge is
#' active through `observation_end`, but the stored memberships are never
#' rewritten and `hg_get()` returns the original spells. Without them
#' the window is the span of the data.
#'
#' Every other column that is constant within a hyperedge is kept as a
#' hyperedge attribute in the edge metadata, where `plot()` can colour by it
#' and where [pairwise_network()] finds the source a citation block belongs to.
#' Columns that vary within a hyperedge, such as the seat an arbitrator held,
#' are not attributes of the hyperedge and are left out.
#'
#' @param data A data frame with one row per relation, or a sequence table
#'   (a wide data frame of states, a list of character vectors, or a `tna` /
#'   `netobject` model); see Details.
#' @param from,to Column names of a pairwise edge list. Detected from the
#'   alias table when neither is given and `node`/`hyperedge` are not named.
#' @param node,hyperedge Column names of membership data: the node, and the
#'   grouping whose shared values bind nodes into one hyperedge. Naming
#'   either selects the membership format; the other is then detected by
#'   alias if not given. When no column is named and no edge list is
#'   detected, both are detected by alias (`node`, `actor`, `person`, ...
#'   and `hyperedge`, `group`, `event`, ...), as [group_hypergraph()] does.
#' @param time Column with the instant of a contact hyperedge.
#' @param start,end Columns with the interval on which a hyperedge is active.
#'   `start` without an `end` column is read as a contact clock.
#' @param weight Optional membership-weight column.
#' @param nodes Optional node universe: a character vector of node names, or
#'   a data frame whose first column holds the names and whose `start`,
#'   `time` or `date` column, if present, gives the time each node enters (on
#'   the same clock as `data`). Nodes that never appear in a hyperedge are
#'   then kept as zero-degree nodes of every snapshot, and [hg_growth()]
#'   counts a node from its own start. Every observed node must be in the
#'   universe.
#' @param time_unit Unit for calendar times: `"auto"` (default),
#'   `"seconds"`, `"minutes"`, `"hours"`, `"days"` or `"weeks"`. Numeric
#'   times are left alone and reported as `"step"`.
#' @param observation_start,observation_end Optional bounds of the
#'   observation window, as numbers on the stored clock or as dates for a
#'   calendar hypergraph. Either may be omitted; the corresponding limit of
#'   the data is then used.
#' @param separator Split the `node` column on this string, one row per
#'   member, before building. Bibliographic exports ship a hyperedge's members
#'   as a single delimited cell -- EUR-Lex `citationcelex` and `eurovoc`,
#'   Scopus and Web of Science reference and keyword fields -- so
#'   `separator = ";"` replaces the caller's own split, trim and
#'   drop-empties. Members empty after trimming are dropped, and a row left
#'   with no member contributes no hyperedge.
#' @param sparse Store every snapshot's incidence as a sparse `Matrix`?
#'   Default `FALSE`.
#' @param cooccur_by Deprecated name of `group`; using it warns with a
#'   `hypergraphs_deprecated` condition.
#' @return A `net_temporal_hypergraph` holding the membership table (`node`,
#'   `edge`, `start`, `end`, `weight`), the edge metadata (`edge`, `start`,
#'   `end` and the hyperedge attributes), the node universe with entry
#'   times, the sorted event times, `format` (`"interval"` or `"contact"`),
#'   `time_unit`, `origin` and the `observation` bounds. `hg_get(x,
#'   what = "memberships" | "edges" | "nodes")` returns the three tables;
#'   `summary()` the one-row description.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' seats <- data.frame(
#'   case = c("A", "A", "A", "B", "B", "B"),
#'   arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
#'   constituted = c(1, 1, 1, 2, 2, 2), concluded = c(4, 4, 4, 5, 5, 5),
#'   sector = c("oil", "oil", "oil", "gas", "gas", "gas")
#' )
#' thg <- temporal_hypergraph(seats, node = "arbitrator", hyperedge = "case",
#'                            start = "constituted", end = "concluded")
#' thg
#' hg_snapshot(thg, at = 3)
#'
#' # a contact log on a calendar: instants at each date, cumulative on request
#' contacts <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
#'                        date = as.Date(c("2024-01-01", "2024-01-05", "2024-01-09")))
#' calls <- temporal_hypergraph(contacts, from = "from", to = "to", time = "date")
#' calls
#' hg_snapshot(calls, at = as.Date("2024-01-05"))
#' hg_snapshot(calls, at = as.Date("2024-01-05"), mode = "cumulative")
#' @export
temporal_hypergraph <- function(data, from = NULL, to = NULL, node = NULL,
                                hyperedge = NULL, time = NULL, start = NULL,
                                end = NULL, weight = NULL, nodes = NULL,
                                time_unit = "auto",
                                observation_start = NULL, observation_end = NULL,
                                sparse = FALSE, separator = NULL,
                                cooccur_by = NULL) {
  if (.thg_is_sequence_input(data, from, to, node, hyperedge, time, start, end)) {
    data <- .thg_sequence_memberships(data)
    node <- "state"
    hyperedge <- "sequence"
    time <- "position"
  }
  if (!is.data.frame(data) || nrow(data) == 0L) {
    .thg_bad_input("`data` must be a non-empty data.frame")
  }
  if (!is.null(separator)) {
    data <- .thg_expand_delimited(data, node %||% to, separator)
  }
  if (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse)) {
    .thg_bad_input("`sparse` must be TRUE or FALSE")
  }
  if (!is.null(cooccur_by)) {
    .thg_deprecated("cooccur_by", "hyperedge", "temporal_hypergraph")
    if (is.null(hyperedge)) hyperedge <- cooccur_by
  }
  roles <- .thg_resolve_roles(data, from, to, node, hyperedge, time, start, end,
                              weight)
  from <- roles$from; to <- roles$to; node <- roles$actor; hyperedge <- roles$group
  clock <- roles$clock; end <- roles$end; weight <- roles$weight
  edge_list <- !is.null(from)
  format <- roles$format

  # One row per membership. An edge list is unpivoted to its two ends; the
  # remaining columns ride along as candidate hyperedge attributes, kept in a
  # table of their own so that a column named `node`, `start` or `weight`
  # never replaces the structural column of that name.
  used <- c(from, to, node, hyperedge, clock, end, weight)
  candidates <- setdiff(names(data), used)
  if (edge_list) {
    edge_id <- paste0("e", seq_len(nrow(data)))
    memberships <- data.frame(
      node = c(as.character(data[[from]]), as.character(data[[to]])),
      edge = c(edge_id, edge_id),
      stringsAsFactors = FALSE
    )
    rows <- c(seq_len(nrow(data)), seq_len(nrow(data)))
  } else {
    memberships <- data.frame(
      node = as.character(data[[node]]),
      edge = as.character(data[[hyperedge]]),
      stringsAsFactors = FALSE
    )
    rows <- seq_len(nrow(data))
  }
  if (!is.null(weight) && (!is.numeric(data[[weight]]) || is.factor(data[[weight]]))) {
    .thg_bad_input("membership weights must be numbers (a factor's codes are not weights)")
  }
  if (is.numeric(data[[clock]]) && any(is.nan(data[[clock]]))) {
    .thg_bad_input(sprintf("time column `%s` holds NaN; times must be finite", clock))
  }
  memberships$start <- data[[clock]][rows]
  memberships$end <- if (is.null(end)) rep(NA, length(rows)) else data[[end]][rows]
  memberships$weight <- if (is.null(weight)) rep(1, length(rows)) else
    as.numeric(data[[weight]][rows])
  metadata <- data[rows, candidates, drop = FALSE]

  keep <- !is.na(memberships$node) & nzchar(memberships$node) &
    !is.na(memberships$edge) & nzchar(memberships$edge) &
    !is.na(memberships$start) & !is.na(memberships$weight)
  memberships <- memberships[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  rownames(memberships) <- NULL
  rownames(metadata) <- NULL
  if (nrow(memberships) == 0L) {
    .thg_bad_input("no complete memberships remain after dropping missing rows")
  }
  if (any(!is.finite(memberships$weight)) || any(memberships$weight < 0)) {
    .thg_bad_input("membership weights must be finite and non-negative")
  }

  # One clock for every time the object holds: hyperedge starts and ends and
  # the node entry times share the origin and the unit.
  observed_nodes <- sort(unique(memberships$node))
  node_data <- .thg_node_universe(nodes, observed_nodes)
  clock_columns <- list(start = memberships$start)
  if (!is.null(end)) clock_columns$end <- memberships$end
  if (!all(is.na(node_data$start))) clock_columns$node_start <- node_data$start
  parsed <- .thg_parse_clock(clock_columns, time_unit)
  timed <- unlist(parsed[setdiff(names(clock_columns), "unit")], use.names = FALSE)
  if (any(is.infinite(timed))) {
    .thg_bad_input("times must be finite; an infinite time is not on the clock")
  }
  memberships$start <- parsed$start
  memberships$end <- if (is.null(end)) rep(NA_real_, nrow(memberships)) else parsed$end
  node_data$start <- if (is.null(parsed$node_start)) rep(NA_real_, nrow(node_data)) else
    parsed$node_start
  # every membership spell on its own, before any hyperedge envelope can hide
  # a reversed one behind a valid one
  if (any(!is.na(memberships$end) & memberships$end < memberships$start)) {
    .thg_bad_input("every interval must satisfy `end >= start`")
  }

  # Times must not depend on which membership row carried them; attribute
  # columns that vary within a hyperedge are not hyperedge attributes.
  constant_within_edge <- function(values, allow_na = FALSE) {
    by_edge <- split(values, memberships$edge)
    !any(vapply(by_edge, function(x) {
      x <- if (allow_na) x[!is.na(x)] else x
      length(unique(x)) > 1L
    }, logical(1L)))
  }
  # A hyperedge's time is the hull of its memberships. When every membership
  # carries the hyperedge's time (a tribunal, a citation block) the hull is
  # that time and nothing below changes. When memberships carry their own
  # times (a session read as a sequence, a log with one row per attendance)
  # the hyperedge spans its first to its last membership and each membership
  # is present on its own spell only; a contact is present at its instant.
  # An open end (NA) is an end of its own: a membership still active is not
  # one that ended, so open and finite ends within a hyperedge differ.
  membership_times <- !constant_within_edge(memberships$start) ||
    (!is.null(end) && !constant_within_edge(memberships$end))
  if (membership_times && identical(format, "contact")) {
    memberships$end <- memberships$start
    format <- "interval"
  }
  attributes <- Filter(function(column) {
    constant_within_edge(metadata[[column]], allow_na = TRUE)
  }, candidates)

  edge_names <- sort(unique(memberships$edge))
  first <- match(edge_names, memberships$edge)
  edge_data <- .thg_edge_attributes(metadata[attributes], memberships$edge,
                                    edge_names,
                                    reserved = c("edge", "start", "end"))
  edge_data$start <- memberships$start[first]
  edge_data$end <- memberships$end[first]
  edge_data <- edge_data[, c("edge", "start", "end",
                             setdiff(names(edge_data), c("edge", "start", "end"))),
                         drop = FALSE]
  attributes <- setdiff(names(edge_data), c("edge", "start", "end"))
  if (membership_times) {
    by_edge <- factor(memberships$edge, levels = edge_names)
    edge_data$start <- as.numeric(tapply(memberships$start, by_edge, min))
    edge_data$end <- as.numeric(tapply(memberships$end, by_edge, function(e) {
      if (anyNA(e)) NA_real_ else max(e)
    }))
  }
  if (identical(format, "interval")) {
    bad_interval <- !is.na(edge_data$end) & edge_data$end < edge_data$start
    if (any(bad_interval)) .thg_bad_input("every interval must satisfy `end >= start`")
  }
  memberships <- memberships[, c("node", "edge", "start", "end", "weight"),
                             drop = FALSE]

  times <- sort(unique(c(edge_data$start, edge_data$end[!is.na(edge_data$end)],
                         memberships$start, memberships$end[!is.na(memberships$end)])))
  span <- range(c(times, node_data$start[!is.na(node_data$start)]))

  x <- structure(
    list(
      memberships = memberships,
      edge_data = edge_data,
      node_data = node_data,
      nodes = node_data$node,
      edges = edge_names,
      times = times,
      format = format,
      time_unit = parsed$unit,
      origin = parsed$origin,
      observation = c(start = span[[1L]], end = span[[2L]]),
      params = list(from = from, to = to, node = node, hyperedge = hyperedge,
                    time = clock, end = end, weight = weight,
                    attributes = attributes, sparse = sparse,
                    membership_times = membership_times,
                    nodes_given = !is.null(nodes),
                    observation_explicit = !is.null(observation_start) ||
                      !is.null(observation_end))
    ),
    class = "net_temporal_hypergraph"
  )
  lower <- .thg_as_time(observation_start, x, "observation_start")
  upper <- .thg_as_time(observation_end, x, "observation_end")
  if (!is.null(lower) && length(lower) != 1L) .thg_bad_input("`observation_start` must be one time")
  if (!is.null(upper) && length(upper) != 1L) .thg_bad_input("`observation_end` must be one time")
  x$observation <- c(start = lower %||% span[[1L]], end = upper %||% span[[2L]])
  if (x$observation[["end"]] < x$observation[["start"]]) {
    .thg_bad_input("`observation_end` is earlier than `observation_start`")
  }
  x
}

# A sequence table is states only: no relational column is named or detected,
# no clock column is detected, and every column is categorical. Anything
# else takes the relational route, where a table without ends still fails
# with the message that names them.
.thg_is_sequence_input <- function(data, from, to, actor, group, time, start, end) {
  named <- !is.null(from) || !is.null(to) || !is.null(actor) || !is.null(group) ||
    !is.null(time) || !is.null(start) || !is.null(end)
  if (named) return(FALSE)
  if (inherits(data, c("tna", "netobject", "cograph_network"))) return(TRUE)
  if (is.matrix(data)) return(!is.numeric(data))
  if (!is.data.frame(data)) return(is.list(data))
  if (ncol(data) == 0L) return(FALSE)
  categorical <- vapply(data, function(column) {
    is.character(column) || is.factor(column) || all(is.na(column))
  }, logical(1L))
  if (!all(categorical)) return(FALSE)
  detected <- function(role) !is.null(.thg_match_column(data, role))
  relational <- (detected("from") && detected("to")) ||
    (detected("actor") && detected("group"))
  !relational && !detected("time") && !detected("start")
}

# One row per state occurrence: the session is the hyperedge, the state the
# member, and the position 1..n its contact time. Sessions are named as
# window_hypergraph() names them, so the two readings of one table agree.
.thg_sequence_memberships <- function(data) {
  if (inherits(data, c("tna", "netobject", "cograph_network"))) {
    data <- .coerce_sequence_input(data)
  }
  trajectories <- .wh_parse_input(data, action = NULL, actor = NULL, time = NULL)
  trajectories <- trajectories[lengths(trajectories) > 0L]
  if (!length(trajectories)) .thg_bad_input("`data` holds no non-empty sequence")
  n <- lengths(trajectories)
  data.frame(
    state = unlist(trajectories, use.names = FALSE),
    sequence = rep(names(trajectories), n),
    position = unlist(lapply(n, seq_len), use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

# Which columns play which role, and therefore the input shape and the clock.
.thg_resolve_roles <- function(data, from, to, actor, group, time, start, end,
                               weight) {
  copresence <- !is.null(actor) || !is.null(group)
  if (copresence && (!is.null(from) || !is.null(to))) {
    .thg_bad_input("name either `from` and `to` or `node` and `hyperedge`, not both")
  }
  claimed <- character()
  if (copresence) {
    actor <- .thg_resolve_column(data, actor, "actor", arg = "node")
    group <- .thg_resolve_column(data, group, "group", exclude = actor, arg = "hyperedge")
    if (is.null(actor) || is.null(group)) {
      .thg_bad_input("membership data need both `node` and `hyperedge`; name the one that could not be detected")
    }
    claimed <- c(actor, group)
  } else {
    named <- !is.null(from) || !is.null(to)
    from <- .thg_resolve_column(data, from, "from", arg = "from")
    to <- .thg_resolve_column(data, to, "to", exclude = from, arg = "to")
    if (!named && (is.null(from) || is.null(to))) {
      # nothing named and no edge list detected: membership columns, by the
      # alias table, as group_hypergraph() detects them
      actor <- .thg_match_column(data, "actor")
      group <- .thg_match_column(data, "group", exclude = actor)
      if (is.null(actor) || is.null(group)) {
        actor <- NULL
        group <- NULL
      } else {
        from <- NULL
        to <- NULL
        claimed <- c(actor, group)
      }
    }
    if (is.null(actor) && (is.null(from) || is.null(to))) {
      .thg_bad_input("name `from` and `to` for an edge list, or `node` and `hyperedge` for membership data")
    }
    if (is.null(actor)) claimed <- c(from, to)
  }
  if (!is.null(start) && !is.null(time)) {
    .thg_bad_input("supply only one of `start` and `time`")
  }
  weight <- .thg_resolve_column(data, weight, "weight", exclude = claimed, arg = "weight")
  claimed <- c(claimed, weight)
  if (!is.null(time)) {
    clock <- .thg_resolve_column(data, time, "time", arg = "time")
    if (!is.null(end)) .thg_bad_input("`end` belongs to the interval format; name `start` with it, not `time`")
    format <- "contact"
    end <- NULL
  } else if (!is.null(start)) {
    clock <- .thg_resolve_column(data, start, "start", arg = "start")
    end <- .thg_resolve_column(data, end, "end", exclude = c(claimed, clock), arg = "end")
    format <- if (is.null(end)) "contact" else "interval"
  } else {
    detected_start <- .thg_match_column(data, "start", claimed)
    detected_end <- .thg_match_column(data, "end", c(claimed, detected_start))
    if (!is.null(end)) {
      end <- .thg_resolve_column(data, end, "end", arg = "end")
      if (is.null(detected_start)) .thg_bad_input("`end` needs a `start` column; name it")
      clock <- detected_start
      format <- "interval"
    } else if (!is.null(detected_start) && !is.null(detected_end)) {
      clock <- detected_start
      end <- detected_end
      format <- "interval"
    } else {
      clock <- .thg_match_column(data, "time", claimed) %||% detected_start
      if (is.null(clock)) {
        .thg_bad_input("no clock found: name `time` for a contact log or `start` and `end` for intervals")
      }
      end <- NULL
      format <- "contact"
    }
  }
  list(from = from, to = to, actor = actor, group = group, clock = clock,
       end = end, weight = weight, format = format)
}

# The node universe as a sorted `node` / `start` table. `start` is NA when the
# caller gave no entry times.
.thg_node_universe <- function(nodes, observed_nodes) {
  if (is.null(nodes)) {
    return(data.frame(node = observed_nodes, start = rep(NA, length(observed_nodes)),
                      stringsAsFactors = FALSE))
  }
  if (is.data.frame(nodes)) {
    if (ncol(nodes) == 0L) .thg_bad_input("a `nodes` data.frame needs a column of node names")
    node <- as.character(nodes[[1L]])
    clock <- intersect(c("start", "time", "date"), names(nodes))
    start <- if (length(clock)) nodes[[clock[[1L]]]] else rep(NA, length(node))
  } else if (is.atomic(nodes)) {
    node <- as.character(nodes)
    start <- rep(NA, length(node))
  } else {
    .thg_bad_input("`nodes` must be a character vector or a data.frame")
  }
  if (anyNA(node) || any(!nzchar(node)) || anyDuplicated(node)) {
    .thg_bad_input("`nodes` must contain unique, non-missing node names")
  }
  missing_members <- setdiff(observed_nodes, node)
  if (length(missing_members)) {
    .thg_bad_input(sprintf("%d observed node(s) are not in `nodes`",
                           length(missing_members)))
  }
  ord <- order(node)
  data.frame(node = node[ord], start = start[ord], stringsAsFactors = FALSE)
}

.thg_empty_hypergraph <- function(nodes = character(), sparse = FALSE) {
  # a snapshot without a node universe has no nodes: character(0), not NULL
  nodes <- as.character(nodes %||% character())
  n <- length(nodes)
  incidence <- if (sparse) {
    Matrix::Matrix(0, n, 0L, sparse = TRUE, dimnames = list(nodes, NULL))
  } else {
    matrix(0, n, 0L, dimnames = list(nodes, NULL))
  }
  structure(list(
    hyperedges = list(), incidence = incidence, nodes = nodes,
    n_nodes = n, n_hyperedges = 0L, size_distribution = integer(),
    params = list(source = "group_hypergraph", member = "member",
                  group = "edge", weight = NULL)
  ), class = "net_hg")
}

.thg_collapse_duplicate_edges <- function(hg) {
  if (hg$n_hyperedges < 2L) {
    hg$edge_multiplicity <- rep.int(1L, hg$n_hyperedges)
    return(hg)
  }
  members <- .thg_edge_members(hg$incidence)
  signatures <- vapply(members, paste, collapse = "\r", character(1L))
  unique_sig <- unique(signatures)
  first <- match(unique_sig, signatures)
  multiplicity <- tabulate(match(signatures, unique_sig), length(unique_sig))
  incidence <- .thg_binary(hg$incidence[, first, drop = FALSE])
  keep <- seq_len(hg$n_hyperedges) %in% first
  hg <- .thg_rebuild(hg, incidence, keep)
  hg$edge_multiplicity <- as.integer(multiplicity)
  hg
}

# ---------------------------------------------------------------------------
# The measurement grid (Dynet's `.window_spec`): `start`/`end` bound the
# measured period, `step` is how often to look, `window` how much time each
# look covers (0 = a point; "all" = one window over the whole period), and
# `at` names instants directly. Without `step` and `at`, the grid is the
# event times of the hypergraph.
# ---------------------------------------------------------------------------

.thg_check_temporal <- function(x) {
  if (!inherits(x, "net_temporal_hypergraph")) {
    .thg_bad_input("`x` must come from temporal_hypergraph()")
  }
}

.thg_grid <- function(x, start = NULL, end = NULL, step = NULL, window = NULL,
                      at = NULL) {
  whole <- is.character(window) && length(window) == 1L && !is.na(window) &&
    identical(tolower(window), "all")
  if (!whole && !is.null(window) &&
      (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 0)) {
    .thg_bad_input("`window` must be a single non-negative number, or the string \"all\"")
  }
  if (!is.null(step) && (!is.numeric(step) || length(step) != 1L ||
                         !is.finite(step) || step <= 0)) {
    .thg_bad_input("`step` must be a single positive number")
  }
  if (whole && !is.null(step)) {
    .thg_bad_input("`step` has no meaning with `window = \"all\"`, which measures one window")
  }
  if (!is.null(at) && !is.null(step)) {
    .thg_bad_input("supply either `at` (instants) or `step` (a grid), not both")
  }
  obs <- x$observation
  explicit <- isTRUE(x$params$observation_explicit)
  start <- .thg_as_time(start, x, "start")
  end <- .thg_as_time(end, x, "end")
  if (!is.null(start) && length(start) != 1L) .thg_bad_input("`start` must be one time")
  if (!is.null(end) && length(end) != 1L) .thg_bad_input("`end` must be one time")
  if (!is.null(start) && !is.null(end) && end < start) {
    .thg_bad_input(sprintf("`end` (%s) is earlier than `start` (%s)", end, start))
  }
  if (explicit) {
    query_start <- start %||% obs[["start"]]
    query_end <- end %||% obs[["end"]]
    if (query_end < obs[["start"]] || query_start > obs[["end"]]) {
      .thg_bad_input("the requested measurement range does not intersect the observation window",
                     class = "hypergraphs_outside_observation")
    }
    if (!is.null(start)) start <- max(start, obs[["start"]])
    if (!is.null(end)) end <- min(end, obs[["end"]])
  }
  start <- start %||% obs[["start"]]
  end <- end %||% obs[["end"]]

  if (!is.null(at)) {
    times <- .thg_as_time(at, x, "at")
    if (length(times) == 0L) .thg_bad_input("`at` contains no valid times")
    if (explicit && any(times < obs[["start"]] | times > obs[["end"]])) {
      .thg_bad_input("`at` lies outside the observation window",
                     class = "hypergraphs_outside_observation")
    }
    times <- sort(unique(times))
    window <- if (whole) end - start else as.numeric(window %||% 0)
  } else if (whole) {
    times <- start
    window <- end - start
  } else if (!is.null(step)) {
    times <- seq(start, end, by = step)
    window <- as.numeric(window %||% step)
  } else {
    times <- x$times[x$times >= start & x$times <= end]
    if (length(times) == 0L) times <- end
    window <- as.numeric(window %||% 0)
  }
  list(times = times, window = window, closed = whole)
}

# Has a spell begun by the end of the window from `t`? A point (window 0) and
# the closed whole-period window include their end; a positive window
# [t, t + window) excludes it. The one boundary rule for hyperedges,
# memberships and node entries alike, so growth counts and snapshots agree.
.thg_begun <- function(start, t, window, closed = FALSE) {
  upper <- t + window
  if (window > 0 && !closed) start < upper else start <= upper
}

# Logical over the edge metadata: which hyperedges a window from `t`
# measures. A point (window 0) is closed: an interval ending exactly at `t`
# is still active, a contact at `t` is present. A positive window covers
# [t, t + window), closed on the right only for the whole-period window.
# Which spells [start, end] (end NA = open) are in the window [t, t + window)
# -- closed on the right for a point or a closed grid -- or, cumulatively,
# have begun by its end.
.thg_spells_in_window <- function(start, end, t, window, mode, closed = FALSE) {
  begun <- .thg_begun(start, t, window, closed)
  if (identical(mode, "cumulative")) return(begun)
  begun & (is.na(end) | end >= t)
}

.thg_edges_in_window <- function(x, t, window, mode, closed = FALSE) {
  ed <- x$edge_data
  if (identical(x$format, "contact") && !identical(mode, "cumulative")) {
    return(.thg_begun(ed$start, t, window, closed) & ed$start >= t)
  }
  .thg_spells_in_window(ed$start, ed$end, t, window, mode, closed)
}

# The memberships of the active hyperedges; when memberships carry their own
# times, only those present in the window.
.thg_memberships_in_window <- function(x, t, window, mode, closed = FALSE) {
  keep <- .thg_edges_in_window(x, t, window, mode, closed)
  mem <- x$memberships
  present <- mem$edge %in% x$edge_data$edge[keep]
  if (isTRUE(x$params$membership_times)) {
    present <- present & .thg_spells_in_window(mem$start, mem$end, t, window, mode, closed)
  }
  present
}

# Node universe of a snapshot. Without a universe, the nodes are the members
# of the active hyperedges. With one, every node is kept; when nodes carry
# entry times, only nodes entered by the end of the window plus any active
# member are kept.
.thg_snapshot_nodes <- function(x, t, window, closed, active_members) {
  if (!isTRUE(x$params$nodes_given)) return(NULL)
  node_data <- x$node_data
  if (all(is.na(node_data$start))) return(node_data$node)
  entered <- !is.na(node_data$start) &
    .thg_begun(node_data$start, t, window, closed)
  sort(union(node_data$node[entered], unique(active_members)))
}

.thg_snapshot_at <- function(x, t, window, mode, multiedges, closed = FALSE) {
  ed <- x$edge_data
  d <- x$memberships[.thg_memberships_in_window(x, t, window, mode, closed), , drop = FALSE]
  sparse <- isTRUE(x$params$sparse)
  universe <- .thg_snapshot_nodes(x, t, window, closed, d$node)
  if (nrow(d) == 0L) {
    hg <- .thg_empty_hypergraph(universe, sparse)
  } else {
    hg <- group_hypergraph(d, node = "node", hyperedge = "edge", weight = "weight",
                           nodes = universe, sparse = sparse)
    hg$edge_data <- ed[match(colnames(hg$incidence), ed$edge), , drop = FALSE]
    rownames(hg$edge_data) <- NULL
  }
  # the data's names for nodes and hyperedges (arbitrator, case), not the
  # internal membership columns, so printed and plotted snapshots use them
  hg$params$node <- x$params$node %||% hg$params$node
  hg$params$hyperedge <- x$params$hyperedge %||% hg$params$hyperedge
  hg$params$temporal_mode <- mode
  hg$params$at <- t
  hg$params$window <- window
  hg$params$format <- x$format
  hg$params$time_unit <- x$time_unit
  hg$params$origin <- x$origin
  if (!multiedges) hg <- .thg_collapse_duplicate_edges(hg)
  hg
}

.thg_check_mode <- function(mode, fn, arg = "mode") {
  if (identical(mode, "all")) {
    .thg_deprecated(sprintf("%s = \"all\"", arg),
                    sprintf("%s = \"cumulative\" (with `at` unset for the whole period)", arg),
                    fn)
    return("cumulative")
  }
  match.arg(mode, c("active", "cumulative"))
}

.thg_check_multiedges <- function(multiedges) {
  if (!is.logical(multiedges) || length(multiedges) != 1L || is.na(multiedges)) {
    .thg_bad_input("`multiedges` must be TRUE or FALSE")
  }
}

#' Extract one snapshot from a temporal hypergraph
#'
#' A snapshot is a static `net_hg` holding the hyperedges a window
#' measures. `mode = "active"` follows the format: an interval hyperedge is
#' active on its closed interval, a contact hyperedge at its instant.
#' `mode = "cumulative"` keeps every hyperedge begun by the end of the
#' window, the growing view of a contact log; with `at` unset it is the
#' static aggregate of everything observed.
#'
#' @param x A [temporal_hypergraph()].
#' @param at One time value, as a number on the hypergraph's clock or a date
#'   for a calendar hypergraph. `NULL` (default) is the end of observation.
#' @param window How much time the snapshot covers, in the hypergraph's time
#'   unit: `0` (default) evaluates the instant `at`; a positive width covers
#'   `[at, at + window)`; `"all"` covers the whole observation period.
#' @param mode `"active"` (default) or `"cumulative"`, see Description. The
#'   former `"all"` is deprecated: it is `"cumulative"` with `at` unset.
#' @param multiedges Keep distinct edge identities with identical member sets?
#'   `TRUE` matches a multi-hypergraph; `FALSE` collapses them to one edge and
#'   records their counts in `edge_multiplicity`.
#' @return A static `net_hg` usable by every hypergraphs hypergraph verb;
#'   its `params` record `at`, `window`, `temporal_mode` and the clock.
#' @export
hg_snapshot <- function(x, at = NULL, window = 0,
                                mode = c("active", "cumulative"),
                                multiedges = TRUE) {
  .thg_check_temporal(x)
  mode <- .thg_check_mode(mode, "hg_snapshot")
  .thg_check_multiedges(multiedges)
  if (!is.null(at) && length(at) != 1L) {
    .thg_bad_input("`at` must be one time value; use hg_snapshots() for several")
  }
  whole <- is.character(window) && length(window) == 1L && !is.na(window) &&
    identical(tolower(window), "all")
  if (whole && !is.null(at)) {
    .thg_bad_input("`at` has no meaning with `window = \"all\"`, which covers the whole period")
  }
  grid <- .thg_grid(x, window = window,
                    at = if (whole) NULL else at %||% x$observation[["end"]])
  .thg_snapshot_at(x, grid$times, grid$window, mode, multiedges, grid$closed)
}

#' Extract a sequence of temporal-hypergraph snapshots
#'
#' Snapshots on a measurement grid with Dynet's meaning: `start` and `end`
#' bound the measured period (default: the observation window), `step` is
#' how often to look, `window` how much time each look covers, and `at`
#' names instants directly. Without `step` and `at`, the grid is the event
#' times of the hypergraph, each looked at as a point.
#'
#' @inheritParams hg_snapshot
#' @param start,end Bounds of the measured period, as numbers on the
#'   hypergraph's clock or dates for a calendar hypergraph. Default: the
#'   observation window.
#' @param step How often to measure, in the hypergraph's time unit. `NULL`
#'   (default) measures at every event time.
#' @param window How much time each measurement covers. Defaults to `step`
#'   when a step is given (a partition of the period) and to `0` otherwise
#'   (a point sample); larger than `step` gives a rolling window; `"all"`
#'   measures the whole period as one window.
#' @param at Instants to measure instead of a grid; may be a vector.
#' @return A named list of `net_hg` objects with class
#'   `net_hg_snapshots`, named by the window starts on the
#'   hypergraph's clock.
#' @examples
#' seats <- data.frame(
#'   case = rep(c("A", "B", "C"), each = 3),
#'   arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
#'   constituted = rep(c(1, 2, 4), each = 3), concluded = rep(c(4, 3, 6), each = 3)
#' )
#' thg <- temporal_hypergraph(seats, node = "arbitrator", hyperedge = "case",
#'                            start = "constituted", end = "concluded")
#' yearly <- hg_snapshots(thg, step = 2)
#' names(yearly)
#' @export
hg_snapshots <- function(x, start = NULL, end = NULL, step = NULL,
                                 window = NULL, mode = c("active", "cumulative"),
                                 at = NULL, multiedges = TRUE) {
  .thg_check_temporal(x)
  mode <- .thg_check_mode(mode, "hg_snapshots")
  .thg_check_multiedges(multiedges)
  grid <- .thg_grid(x, start, end, step, window, at)
  out <- lapply(grid$times, function(t) {
    .thg_snapshot_at(x, t, grid$window, mode, multiedges, grid$closed)
  })
  names(out) <- make.unique(as.character(grid$times))
  structure(out, class = c("net_hg_snapshots", "list"),
            time_unit = x$time_unit, origin = x$origin, window = grid$window)
}

#' @export
print.net_temporal_hypergraph <- function(x, n = 10L, ...) {
  cat(sprintf("Temporal hypergraph: %d nodes, %d hyperedges, %d event times\n",
              length(x$nodes), length(x$edges), length(x$times)))
  cat(if (isTRUE(x$params$membership_times)) {
    paste0("Format: interval; memberships carry their own times ",
           "(a hyperedge spans its first to its last membership)\n")
  } else if (identical(x$format, "interval")) {
    "Format: interval (a hyperedge is active from its start to its end)\n"
  } else {
    paste0("Format: contact (a hyperedge is an instantaneous event; ",
           "mode = \"cumulative\" keeps every hyperedge once it appears)\n")
  })
  cat(sprintf("Time: %s; observed from %s to %s%s\n",
              .thg_clock_label(x$origin, x$time_unit),
              format(x$observation[["start"]]), format(x$observation[["end"]]),
              if (isTRUE(x$params$observation_explicit)) " (declared)" else ""))
  .ho_print_table(x, n)
  invisible(x)
}

#' @export
summary.net_temporal_hypergraph <- function(object, ...) {
  # an edge's size counts its distinct members: a member recorded twice (a
  # repeated contact, a second spell) is one member
  distinct <- unique(object$memberships[c("edge", "node")])
  memberships_per_edge <- table(distinct$edge)
  duration <- if (all(is.na(object$edge_data$end))) NA_real_ else
    as.numeric(object$edge_data$end - object$edge_data$start)
  data.frame(
    n_nodes = length(object$nodes),
    n_hyperedges = length(object$edges),
    n_event_times = length(object$times),
    first_time = min(object$times),
    last_time = max(object$times),
    n_memberships = nrow(object$memberships),
    mean_edge_size = mean(as.numeric(memberships_per_edge)),
    median_edge_size = stats::median(as.numeric(memberships_per_edge)),
    mean_duration = if (all(is.na(duration))) NA_real_ else
      mean(as.numeric(duration), na.rm = TRUE),
    format = object$format,
    time_unit = object$time_unit,
    row.names = NULL
  )
}

#' Tidy tables of a temporal hypergraph
#'
#' @param x A [temporal_hypergraph()].
#' @param what `"memberships"` (default) for one row per node-in-hyperedge
#'   spell (`member`, `edge`, `start`, `end`, `weight`), `"edges"` for one
#'   row per hyperedge with its clock and attributes, `"nodes"` for the node
#'   universe with entry times. Times are on the hypergraph's clock.
#' @param ... Ignored.
#' @return A base `data.frame`.
#' @export
hg_get.net_temporal_hypergraph <- function(x, what = c("memberships",
                                                       "edges", "nodes"),
                                           ...) {
  what <- .ho_match_what(what)
  switch(what, memberships = x$memberships, edges = x$edge_data, nodes = x$node_data)
}

#' Plot a temporal hypergraph
#'
#' Plots the snapshot that [hg_snapshot()] takes at `at` as a hypergraph,
#' with [plot.net_hg()]: the hyperedges active then, as hulls by default or
#' as the incidence matrix with `type = "incidence"`. The snapshot keeps the
#' hyperedge attributes and the data's names for nodes and hyperedges, so
#' `color_by` can name an attribute and the legends use the data's words.
#'
#' `type = "storyline"` plots the whole history instead (Tanahashi and Ma
#' 2012). Each of the `top` nodes with the most hyperedges is a line, from
#' its first hyperedge to its last, and each hyperedge is a column, in order
#' of its start, where a grey bar gathers the lines of its members among the
#' plotted nodes. The columns are spaced by order, not by elapsed time, and
#' are labelled with the hyperedge and its start. Lines are ordered by the
#' barycentre rule of Sugiyama, Tagawa and Toda (1981): along the columns the
#' members of each hyperedge move to the median of their current rows, and
#' the ordering is refined over repeated sweeps, keeping the one with the
#' fewest line crossings. Each line takes an Okabe-Ito colour and its points
#' a shape, and the legend names the nodes in decreasing order of their
#' number of hyperedges. Eight colours and nine shapes cycle, so up to 72
#' lines differ in their pair of colour and shape.
#'
#' @param x A [temporal_hypergraph()].
#' @param at Snapshot time; defaults to the end of observation. Not used by
#'   the storyline.
#' @param mode Snapshot mode passed to [hg_snapshot()]. Not used by the
#'   storyline.
#' @param type `"hulls"` (default) or `"incidence"` plot the snapshot at
#'   `at` (see [plot.net_hg()]); `"storyline"` plots the whole history.
#' @param top For `type = "storyline"`: the number of nodes with the most
#'   hyperedges to draw as lines (default `8`, one Okabe-Ito colour each;
#'   ties broken by name), or
#'   `NULL` for every node.
#' @param start,end For `type = "storyline"`: draw only the hyperedges that
#'   begin in this period, as dates for a calendar hypergraph or numbers on
#'   its clock. `NULL` (default) leaves the period open.
#' @param method Deprecated. The snapshot used to be projected to a pairwise
#'   network; `plot(pairwise_network(hg_snapshot(x, at), type = method))`
#'   plots that network.
#' @param ... For `"hulls"` and `"incidence"`, arguments passed to
#'   [plot.net_hg()], such as `color_by` or `labels`. For `"storyline"`,
#'   `edge_labels` (`TRUE` writes the hyperedge and its start under each
#'   column), `point_size` (size of the points, default `2.5`) and
#'   `spacing`: `"even"` (default) puts neighbouring lines one row apart,
#'   `"strength"` brings two neighbouring lines closer the more plotted
#'   hyperedges their nodes share: a full row apart for none, 0.35 of a row
#'   for the largest number shared by any pair in the plot, and linearly in
#'   between, so lines that often meet run together. Spacing never changes
#'   the order of the lines. `width_by = "degree"` draws each line with a
#'   width that grows linearly with its node's number of hyperedges in the
#'   period, the count `top` ranks by, and adds a width legend; `NULL`
#'   (default) draws every line alike.
#' @return A ggplot object.
#' @section Conditions:
#' `hypergraphs_bad_input` when the snapshot has no active nodes, when no
#' hyperedge begins between `start` and `end`, for an invalid `top`, for
#' `at`, `mode` or an argument of [plot.net_hg()] with `type = "storyline"`,
#' for `top`, `start` or `end` with another type, and as [plot.net_hg()]
#' raises it; `hypergraphs_deprecated` (a warning) for `method`.
#' @references
#' Tanahashi, Y., & Ma, K.-L. (2012). Design considerations for optimizing
#' storyline visualizations. *IEEE Transactions on Visualization and
#' Computer Graphics*, 18(12), 2679-2688. \doi{10.1109/TVCG.2012.212}
#'
#' Sugiyama, K., Tagawa, S., & Toda, M. (1981). Methods for visual
#' understanding of hierarchical system structures. *IEEE Transactions on
#' Systems, Man, and Cybernetics*, 11(2), 109-125.
#' \doi{10.1109/TSMC.1981.4308636}
#' @examples
#' seats <- data.frame(
#'   case = rep(c("A", "B", "C"), each = 3),
#'   arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
#'   constituted = rep(c(1, 2, 4), each = 3),
#'   concluded = rep(c(4, 3, 6), each = 3)
#' )
#' thg <- temporal_hypergraph(seats, node = "arbitrator", hyperedge = "case",
#'                            start = "constituted", end = "concluded")
#' plot(thg, at = 2.5)
#' plot(thg, at = 2.5, type = "incidence", edge_labels = TRUE)
#' plot(thg, type = "storyline")
#' @export
plot.net_temporal_hypergraph <- function(x, at = NULL,
                                         mode = c("active", "cumulative"),
                                         type = c("hulls", "incidence", "storyline"),
                                         top = 8L, start = NULL, end = NULL,
                                         method = NULL, ...) {
  type <- match.arg(type)
  given <- names(match.call())[-1L]
  if (!is.null(method)) {
    .thg_deprecated("method", "plot(pairwise_network(hg_snapshot(x, at), type = method))",
                    "plot")
  }
  if (identical(type, "storyline")) {
    snapshot_only <- intersect(given, c("at", "mode"))
    dots <- list(...)
    dot_names <- names(dots) %||% rep("", length(dots))
    refused <- c(paste0("`", snapshot_only, "`"),
                 paste0("`", setdiff(dot_names[nzchar(dot_names)],
                                     c("edge_labels", "point_size", "spacing", "width_by")), "`"),
                 if (any(!nzchar(dot_names))) "unnamed arguments")
    refused <- refused[refused != "``"]
    if (length(refused)) {
      .thg_bad_input(sprintf("%s: not used by `type = \"storyline\"`",
                             paste(refused, collapse = ", ")))
    }
    return(do.call(.thg_plot_storyline,
                   c(list(x, top = top, start = start, end = end), dots)))
  }
  storyline_only <- intersect(given, c("top", "start", "end"))
  if (length(storyline_only)) {
    .thg_bad_input(sprintf("%s applies only to `type = \"storyline\"`",
                           paste0("`", storyline_only, "`", collapse = ", ")))
  }
  mode <- .thg_check_mode(mode, "plot")
  hg <- hg_snapshot(x, at = at, mode = mode)
  if (hg$n_nodes == 0L) .thg_bad_input("the selected snapshot has no active nodes")
  plot(hg, type = type, ...)
}
