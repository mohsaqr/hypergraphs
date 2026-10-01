# Sessions with a coding assistant

The events of 4,000 sessions in which a user runs code with the help of
an assistant, one row per event. A session is one task; it ends `Done`
when a run passes, or `Abort` when the user exits after failed runs. A
run is a stretch of a session between two verdicts: a new run starts at
the event after each `Fail`, so the assistant's moves after a failed run
belong to the rerun they lead to. The sessions are real event sequences
of another domain with every event renamed and learner identifiers
removed; the sessions are a random 4,000 of 13,309. Rebuilt by
`data-raw/debug_events.R`.

## Usage

``` r
debug_events
```

## Format

A data frame with 25,008 rows and 5 columns:

- session:

  Session number (integer).

- run:

  Run identifier, unique over the data (`"<session>.<n>"`).

- position:

  Position of the event within its session (integer).

- event:

  The event: the user's `Open`, `Run`, `Rerun` and `Exit`, the task's
  `Spec`, the verdicts `Pass` and `Fail`, the session's end `Done` or
  `Abort`, and the assistant's `Tip`, `Query`, `Fix`, `Note`, `Counter`,
  `Reflect`, `Confirm` and `Reassure`.

- group:

  `"quick"` for a session done on the first run or after one rerun,
  `"slow"` for one that needs several reruns or is aborted.

## Source

Renamed and sampled by `data-raw/debug_events.R` (seed 20261001).

## Examples

``` r
session_sets <- group_hypergraph(debug_events, actor = "event",
                                 group = "session", by = "group")
session_sets
#> Hypergraph: 15 nodes, 16 hyperedges (sizes 4: 1, 5: 1, 7: 2, 8: 5, 9: 5, 10: 2)
#> Source: sets of event per session, counted within group (top 8 per group)
#>                                                                          hyperedge
#>                       quick: Done + Fail + Fix + Pass + Query + Rerun + Run + Spec
#>                               quick: Done + Fail + Fix + Pass + Rerun + Run + Spec
#>               quick: Done + Fail + Note + Open + Pass + Query + Rerun + Run + Spec
#>                      quick: Done + Fail + Note + Pass + Query + Rerun + Run + Spec
#>                        quick: Done + Fail + Open + Pass + Rerun + Run + Spec + Tip
#>                               quick: Done + Fail + Pass + Rerun + Run + Spec + Tip
#>                                             quick: Done + Open + Pass + Run + Spec
#>                                                    quick: Done + Pass + Run + Spec
#>  slow: Abort + Confirm + Exit + Fail + Query + Reassure + Rerun + Run + Spec + Tip
#>            slow: Abort + Exit + Fail + Fix + Query + Reassure + Rerun + Run + Spec
#>  size                                                            members weight
#>     8                     Done, Fail, Fix, Pass, Query, Rerun, Run, Spec     NA
#>     7                            Done, Fail, Fix, Pass, Rerun, Run, Spec     NA
#>     9              Done, Fail, Note, Open, Pass, Query, Rerun, Run, Spec     NA
#>     8                    Done, Fail, Note, Pass, Query, Rerun, Run, Spec     NA
#>     8                      Done, Fail, Open, Pass, Rerun, Run, Spec, Tip     NA
#>     7                            Done, Fail, Pass, Rerun, Run, Spec, Tip     NA
#>     5                                        Done, Open, Pass, Run, Spec     NA
#>     4                                              Done, Pass, Run, Spec     NA
#>    10 Abort, Confirm, Exit, Fail, Query, Reassure, Rerun, Run, Spec, Tip     NA
#>     9          Abort, Exit, Fail, Fix, Query, Reassure, Rerun, Run, Spec     NA
#> ... 6 more rows
```
