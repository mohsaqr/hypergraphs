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

#' Problem steps of learners with a tutor
#'
#' The events of 13,309 problem steps that learners worked through with an AI
#' tutor, one row per event. A step begins with the task, the learner attempts
#' an answer, and the answer is correct or incorrect. After an incorrect answer
#' the tutor responds with guidance, a question, an order, a refutation, a
#' prompt to reflect or comfort, and the learner reattempts or gives up. A
#' step ends `Completed` when an answer is correct, or `Stopped` when it ends
#' without a correct answer.
#' A trial is the part of a step between two answers: a new trial starts at the
#' event after each `Incorrect`, so the tutor's response to an incorrect answer
#' belongs to the reattempt it leads to. The events are real; each is renamed
#' to a word of the same meaning, similar events are merged, and learner and
#' skill identifiers are removed. Rebuilt by `data-raw/tutoring_events.R`.
#'
#' @format A data frame with 84,356 rows and 5 columns:
#' \describe{
#'   \item{step}{Problem step number (integer).}
#'   \item{trial}{Trial identifier, unique over the data (`"<step>.<n>"`).}
#'   \item{position}{Position of the event within its step (integer).}
#'   \item{event}{The event: the step's `Begin` and `Task`, the learner's
#'     `Attempt`, `Reattempt` and `GiveUp`, the answers `Correct` and
#'     `Incorrect`, the step's end `Completed` or `Stopped`, and the tutor's
#'     `Guidance`, `Question`, `Order`, `Refute`, `Reflect` and `Comfort`.}
#'   \item{outcome}{`"completed"` for a step that ends with a correct
#'     answer, `"stopped"` for one that ends without.}
#' }
#' @source Renamed by `data-raw/tutoring_events.R`.
#' @examples
#' step_sets <- group_hypergraph(tutoring_events, node = "event",
#'                               hyperedge = "trial", group = "outcome", top = 4)
#' step_sets
"tutoring_events"
