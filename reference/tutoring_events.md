# Problem steps of learners with a tutor

The events of 13,309 problem steps that learners worked through with an
AI tutor, one row per event. A step begins with the task, the learner
attempts an answer, and the answer is correct or incorrect. After an
incorrect answer the tutor responds with guidance, a question, an order,
a refutation, a prompt to reflect or comfort, and the learner reattempts
or gives up. A step ends `Completed` when an answer is correct, or
`Stopped` when it ends without a correct answer. A trial is the part of
a step between two answers: a new trial starts at the event after each
`Incorrect`, so the tutor's response to an incorrect answer belongs to
the reattempt it leads to. The events are real; each is renamed to a
word of the same meaning, similar events are merged, and learner and
skill identifiers are removed. Rebuilt by `data-raw/tutoring_events.R`.

## Usage

``` r
tutoring_events
```

## Format

A data frame with 84,356 rows and 5 columns:

- step:

  Problem step number (integer).

- trial:

  Trial identifier, unique over the data (`"<step>.<n>"`).

- position:

  Position of the event within its step (integer).

- event:

  The event: the step's `Begin` and `Task`, the learner's `Attempt`,
  `Reattempt` and `GiveUp`, the answers `Correct` and `Incorrect`, the
  step's end `Completed` or `Stopped`, and the tutor's `Guidance`,
  `Question`, `Order`, `Refute`, `Reflect` and `Comfort`.

- outcome:

  `"completed"` for a step that ends with a correct answer, `"stopped"`
  for one that ends without.

## Source

Renamed by `data-raw/tutoring_events.R`.

## Examples

``` r
step_sets <- group_hypergraph(tutoring_events, node = "event",
                              hyperedge = "trial", group = "outcome", top = 4)
step_sets
#> Hypergraph: 12 nodes, 8 hyperedges (sizes 3: 5, 4: 2, 5: 1)
#> Source: sets of event per trial, counted within each outcome (top 4 per group)
#>                                                hyperedge size
#>  completed: Attempt + Begin + Completed + Correct + Task    5
#>          completed: Attempt + Completed + Correct + Task    4
#>                    completed: Attempt + Incorrect + Task    3
#>    completed: Completed + Correct + Guidance + Reattempt    4
#>                      stopped: Attempt + Incorrect + Task    3
#>                      stopped: Comfort + GiveUp + Stopped    3
#>                stopped: Guidance + Incorrect + Reattempt    3
#>                stopped: Incorrect + Question + Reattempt    3
#>                                   members
#>  Attempt, Begin, Completed, Correct, Task
#>         Attempt, Completed, Correct, Task
#>                  Attempt, Incorrect, Task
#>   Completed, Correct, Guidance, Reattempt
#>                  Attempt, Incorrect, Task
#>                  Comfort, GiveUp, Stopped
#>            Guidance, Incorrect, Reattempt
#>            Incorrect, Question, Reattempt
```
