# Sequence input shared by the sequence verbs

Sequence input shared by the sequence verbs

## Sequence input

Every verb that reads sequences reads them the same way:

- Long event table:

  a data.frame with one event per row. `action` names the state column,
  `actor` the column (or columns) that identify whose events they are,
  `session` the column (or columns) that split an actor's events into
  sessions, and `time` the column that orders them: timestamps in any
  common date-time format, Unix times, or positions; row order when
  `NULL`. A column named `action`, `time`, `session` or `session_id` (in
  any case) is used when its argument is `NULL`; `session = FALSE`
  switches the session detection off. When `time` is given, a gap of
  more than `time_threshold` seconds between consecutive events starts a
  new sequence (default `900`, fifteen minutes);
  `time_threshold = FALSE` keeps each actor or session in one sequence
  however long the gaps. `timezone` is the time zone of timestamps that
  carry none (default `"UTC"`). Events with a missing actor or session
  raise `hypergraphs_bad_input`; a missing action stays in its sequence
  as a gap. Without `actor` all events form one sequence, announced by
  the message `hypergraphs_single_sequence`.

- Wide data.frame or character matrix:

  one sequence per row; trailing `NA`s end a sequence.

- List:

  one character vector per sequence.

- Model object:

  a `netobject`, `netobject_group`, `tna` or `cograph_network` that
  carries its sequence data.

These are the conventions of the tna family of packages, and the same
call builds the same sequences there.

A missing state anywhere but at the end of a wide row is a gap. The
memory verbs
([`hon()`](https://pak.dynasite.org/hypergraphs/reference/hon.md),
[`mogen()`](https://pak.dynasite.org/hypergraphs/reference/mogen.md),
[`hypa()`](https://pak.dynasite.org/hypergraphs/reference/hypa.md),
[`markov_order()`](https://pak.dynasite.org/hypergraphs/reference/markov_order.md),
[`memory()`](https://pak.dynasite.org/hypergraphs/reference/memory.md),
[`hg_bootstrap()`](https://pak.dynasite.org/hypergraphs/reference/hg_bootstrap.md),
[`hg_compare()`](https://pak.dynasite.org/hypergraphs/reference/hg_compare.md))
split a sequence at every gap into its runs of observed states: no
transition is counted across or into a gap and no state `"NA"` is
created (a state spelled `"NA"` is an ordinary state).
[`hg_bootstrap()`](https://pak.dynasite.org/hypergraphs/reference/hg_bootstrap.md)
and
[`hg_compare()`](https://pak.dynasite.org/hypergraphs/reference/hg_compare.md)
still resample the original sequences, each with all its runs. Repeats
are collapsed within a run. State labels containing `" -> "` (the
notation of memory nodes) are refused with `hypergraphs_bad_input` by
the verbs whose results name paths.

A data.frame with columns named like an event table (`code`, `state`,
`user`, `timestamp`, ...) that is passed without `action =` and has no
`action` column raises `hypergraphs_long_format` (a
`hypergraphs_bad_input`): read as wide, its actor ids and times would
silently become states. `actor`, `time` or `session` without an action
column raises `hypergraphs_bad_input`.
