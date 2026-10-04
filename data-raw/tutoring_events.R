# Build the bundled `tutoring_events` dataset: the event sequences of all
# 13,309 problem steps of a corpus of learners working with an AI tutor (the
# real data are kept outside the package). Every event is renamed to a word of
# the same meaning, similar events are merged (19 events become 15), and
# learner and skill identifiers are dropped.
# A trial is a stretch of a step between two answers: a new trial starts at
# the event after each Incorrect, so the tutor's response to an incorrect
# answer belongs to the reattempt it leads to.
# Run from the package root: Rscript data-raw/tutoring_events.R

load(file.path("local_testing_and_equivalence", "real-data",
               "tutor_events_real.rda"))
real <- tutor_events

rename <- c(
  Start = "Begin", Instructions = "Task", Try = "Attempt",
  Retry = "Reattempt", Right = "Correct", Validation = "Correct",
  Wrong = "Incorrect", Solved = "Completed", Terminated = "Stopped",
  Hint = "Guidance", State = "Guidance", Question = "Question",
  Directive = "Order", Correction = "Order", Disprove = "Refute",
  Think = "Reflect", Assurance = "Comfort", Quit = "GiveUp", Skip = "GiveUp")
stopifnot("every event is renamed" = all(real$event %in% names(rename)))

# the step's outcome: completed with a correct answer, or stopped without one
outcome <- ifelse(real$outcome == "Terminated", "stopped", "completed")
tutoring_events <- data.frame(
  step = real$step,
  trial = real$trial,
  position = real$position,
  event = unname(rename[real$event]),
  outcome = outcome,
  stringsAsFactors = FALSE
)
tutoring_events <- tutoring_events[order(tutoring_events$step,
                                         tutoring_events$position), ]
rownames(tutoring_events) <- NULL

stopifnot(
  "15 events" = length(unique(tutoring_events$event)) == 15L,
  "no missing values" = !anyNA(tutoring_events),
  "one outcome per step" = all(tapply(tutoring_events$outcome,
                                      tutoring_events$step,
                                      \(v) length(unique(v))) == 1L),
  "a step is stopped exactly when it holds Stopped" =
    identical(tapply(tutoring_events$outcome, tutoring_events$step, `[`, 1L) ==
                "stopped",
              tapply(tutoring_events$event, tutoring_events$step,
                     \(v) "Stopped" %in% v))
)

dir.create("data", showWarnings = FALSE)
save(tutoring_events, file = file.path("data", "tutoring_events.rda"),
     compress = "xz")
cat("rows:", nrow(tutoring_events), "| steps:", max(tutoring_events$step),
    "| trials:", length(unique(tutoring_events$trial)),
    "| size:", file.size(file.path("data", "tutoring_events.rda")), "bytes\n")
print(table(tapply(tutoring_events$outcome, tutoring_events$step, `[`, 1L)))
print(sort(table(tutoring_events$event), decreasing = TRUE))
