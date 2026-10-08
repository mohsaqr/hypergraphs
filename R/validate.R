# Shared argument validators. Every check raises the classed
# `hypergraphs_bad_input` condition, so a caller can catch a broken contract
# the same way whichever verb raised it.

.ho_input_error <- function(message) {
  stop(errorCondition(message, class = "hypergraphs_bad_input", call = NULL))
}

# A single whole-number count in [min, max]. `allow_inf = TRUE` admits +Inf
# (an unlimited row count such as `n = Inf`), which is returned unchanged;
# every finite value is returned as an integer. Fractions, NA, -Inf, vectors
# and values beyond the integer range are refused rather than truncated.
.ho_check_count <- function(x, arg, min = 1, max = .Machine$integer.max,
                            allow_inf = FALSE) {
  ok <- is.numeric(x) && !is.factor(x) && length(x) == 1L && !is.na(x)
  if (ok && allow_inf && identical(as.numeric(x), Inf)) return(Inf)
  ok <- ok && is.finite(x) && abs(x - round(x)) < sqrt(.Machine$double.eps) &&
    x >= min && x <= max
  if (!ok) {
    upper <- if (allow_inf) "or Inf" else
      if (max < .Machine$integer.max) sprintf("<= %s", format(max)) else ""
    .ho_input_error(trimws(sprintf(
      "`%s` must be one whole number >= %s %s", arg, format(min), upper
    )))
  }
  as.integer(round(x))
}

# A single finite number in [min, max] (a continuous control: a tolerance,
# a probability, a threshold).
.ho_check_number <- function(x, arg, min = -Inf, max = Inf) {
  ok <- is.numeric(x) && !is.factor(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x >= min && x <= max
  if (!ok) {
    range <- if (is.finite(min) && is.finite(max)) {
      sprintf(" in [%s, %s]", format(min), format(max))
    } else if (is.finite(min)) {
      sprintf(" >= %s", format(min))
    } else if (is.finite(max)) {
      sprintf(" <= %s", format(max))
    } else {
      ""
    }
    .ho_input_error(sprintf("`%s` must be one finite number%s", arg, range))
  }
  as.numeric(x)
}

# Identifiers (node, edge, document names) that must be present, non-empty
# and unique, because they are used as keys. Returns them as character.
.ho_check_ids <- function(ids, arg) {
  ids_chr <- as.character(ids)
  if (anyNA(ids_chr) || any(!nzchar(ids_chr))) {
    .ho_input_error(sprintf("`%s` must not contain missing or empty names", arg))
  }
  dup <- unique(ids_chr[duplicated(ids_chr)])
  if (length(dup)) {
    .ho_input_error(sprintf(
      "`%s` names must be unique; repeated: %s", arg,
      paste(utils::head(dup, 5L), collapse = ", ")
    ))
  }
  ids_chr
}

# Non-negative finite numeric weights (factors refused: their codes are not
# weights). Returns them as double.
.ho_check_weights <- function(w, arg) {
  if (!is.numeric(w) || is.factor(w) || anyNA(w) || any(!is.finite(w)) ||
      any(w < 0)) {
    .ho_input_error(sprintf(
      "`%s` must be finite, non-negative numbers with no missing values", arg
    ))
  }
  as.numeric(w)
}
