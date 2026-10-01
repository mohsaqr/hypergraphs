# Build the bundled `debug_events` dataset: the event sequences of 4,000
# sessions, disguised from a corpus of real tutoring steps (kept outside the
# package): every event is renamed into the vocabulary of a debugging session
# with a coding assistant, two pairs of events are merged, learner and skill
# identifiers are dropped, and 4,000 of the 13,309 steps are drawn at random.
# A run is a stretch of a session between two verdicts: a new run starts at
# the event after each Fail.
# Run from the package root: Rscript data-raw/debug_events.R

load(file.path("local_testing_and_equivalence", "real-data",
               "tutor_events_real.rda"))
real <- tutor_events

rename <- c(
  Start = "Open", Instructions = "Spec", Try = "Run", Retry = "Rerun",
  Right = "Pass", Wrong = "Fail", Solved = "Done", Hint = "Tip",
  Question = "Query", Directive = "Fix", Correction = "Fix", State = "Note",
  Disprove = "Counter", Think = "Reflect", Validation = "Confirm",
  Assurance = "Reassure", Quit = "Exit", Skip = "Exit", Terminated = "Abort")
stopifnot("every event is renamed" = all(real$event %in% names(rename)))

set.seed(20261001)
kept <- sort(sample(unique(real$step), 4000L))
real <- real[real$step %in% kept, ]
session <- match(real$step, kept)

quick <- real$outcome %in% c("Solved on the first try",
                             "Solved after one retry")
debug_events <- data.frame(
  session = session,
  run = sprintf("%04d.%d", session, real$attempt),
  position = real$position,
  event = unname(rename[real$event]),
  group = ifelse(quick, "quick", "slow"),
  stringsAsFactors = FALSE
)
debug_events <- debug_events[order(debug_events$session,
                                   debug_events$position), ]
rownames(debug_events) <- NULL

stopifnot(
  "17 events" = length(unique(debug_events$event)) == 17L,
  "no missing values" = !anyNA(debug_events),
  "one group per session" = all(tapply(debug_events$group,
                                       debug_events$session,
                                       \(v) length(unique(v))) == 1L)
)

dir.create("data", showWarnings = FALSE)
save(debug_events, file = file.path("data", "debug_events.rda"),
     compress = "xz")
cat("rows:", nrow(debug_events), "| sessions:", max(debug_events$session),
    "| runs:", length(unique(debug_events$run)),
    "| size:", file.size(file.path("data", "debug_events.rda")), "bytes\n")
print(table(tapply(debug_events$group, debug_events$session, `[`, 1L)))
