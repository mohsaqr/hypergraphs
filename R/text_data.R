#' COVID-19 education research abstracts
#'
#' A deterministic stratified sample of 165 abstracts (up to 40 per year,
#' 2020--2024) from a Scopus export of COVID-19 education research, the same
#' source corpus used by the `sbert` package's topic-modeling articles.
#' Rebuilt by `data-raw/covid_abstracts.R`.
#'
#' @format A data frame with 165 rows and 4 columns:
#' \describe{
#'   \item{doc}{Scopus EID, the unique document identifier.}
#'   \item{title}{Article title.}
#'   \item{abstract}{Abstract text (each at least 400 characters).}
#'   \item{year}{Publication year (integer, 2020--2024).}
#' }
#' @source Scopus export of COVID-19 education research, 2020--2024.
#' @examples
#' hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc")
#' hg
"covid_abstracts"

#' A random sample of 1,000 COVID-19 education research abstracts
#'
#' A simple random sample of 1,000 abstracts from the Scopus export of
#' COVID-19 education research that the `sbert` package also uses
#' (4,187 records, 2020--2024). Records without an abstract of at least 400
#' characters are excluded before sampling, and an abstract that appears
#' twice is kept once. The years are represented in the proportions of the
#' source. Rebuilt by `data-raw/covid_sample.R`.
#'
#' @format A data frame with 1,000 rows and 4 columns:
#' \describe{
#'   \item{doc}{Scopus EID, the unique document identifier.}
#'   \item{title}{Article title.}
#'   \item{abstract}{Abstract text (each at least 400 characters).}
#'   \item{year}{Publication year (integer, 2020--2024).}
#' }
#' @source Scopus export of COVID-19 education research, 2020--2024.
#' @examples
#' hg <- text_hypergraph(covid_sample, column = "abstract", id = "doc",
#'                       stop_words = stop_words_en(), min_count = 5)
#' hg
"covid_sample"

#' Sentence embeddings of the COVID-19 abstracts
#'
#' Sentence embeddings of [covid_abstracts]' abstract texts, computed with
#' the `sbert` package's pinned `all-MiniLM-L6-v2` model (L2-normalized
#' rows). Bundled so that `text_hypergraph(construction = "knn")` runs
#' offline; rebuilt by `data-raw/covid_embeddings.R`.
#'
#' @format A numeric matrix with 165 rows (rownames = `covid_abstracts$doc`)
#'   and 384 columns.
#' @source Computed from [covid_abstracts] with `sbert::encode()`.
#' @examples
#' hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
#'                       construction = "knn", k = 10,
#'                       embeddings = covid_embeddings)
#' hg
"covid_embeddings"

#' Sessions with a coding assistant
#'
#' The events of 4,000 sessions in which a user runs code with the help of
#' an assistant, one row per event. A session is one task; it ends `Done`
#' when a run passes, or `Abort` when the user exits after failed runs. A run
#' is a stretch of a session between two verdicts: a new run starts at the
#' event after each `Fail`, so the assistant's moves after a failed run belong
#' to the rerun they lead to. The sessions are real event sequences of another
#' domain with every event renamed and learner identifiers removed; the
#' sessions are a random 4,000 of 13,309. Rebuilt by
#' `data-raw/debug_events.R`.
#'
#' @format A data frame with 25,008 rows and 5 columns:
#' \describe{
#'   \item{session}{Session number (integer).}
#'   \item{run}{Run identifier, unique over the data (`"<session>.<n>"`).}
#'   \item{position}{Position of the event within its session (integer).}
#'   \item{event}{The event: the user's `Open`, `Run`, `Rerun` and `Exit`,
#'     the task's `Spec`, the verdicts `Pass` and `Fail`, the session's end
#'     `Done` or `Abort`, and the assistant's `Tip`, `Query`, `Fix`, `Note`,
#'     `Counter`, `Reflect`, `Confirm` and `Reassure`.}
#'   \item{group}{`"quick"` for a session done on the first run or after one
#'     rerun, `"slow"` for one that needs several reruns or is aborted.}
#' }
#' @source Renamed and sampled by `data-raw/debug_events.R` (seed 20261001).
#' @examples
#' session_sets <- group_hypergraph(debug_events, actor = "event",
#'                                  group = "session", by = "group")
#' session_sets
"debug_events"
