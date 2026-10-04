#' English function-word stop list
#'
#' A fixed list of English function words for the `stop_words` argument of
#' [text_hypergraph()] and [clean_text()]. `type = "minimal"` (the default)
#' is a small list of articles, prepositions, conjunctions and auxiliaries,
#' deliberately minimal and versioned with the package. `type = "snowball"`
#' is the English stop list of the Snowball project (Porter 2001), 174 words
#' that add the personal pronouns, the forms of "do" and the contractions
#' ("i", "you", "don't", "i'm"), which chat messages and other informal text
#' are full of. Corpus-specific boilerplate (e.g. "study", "results" in an
#' abstract corpus) is added by the caller, as in
#' `c(stop_words_en(), "study", "results")`.
#'
#' @param type `"minimal"` (default) or `"snowball"`.
#' @return A sorted character vector of lowercase English function words.
#'   Raises `hypergraphs_bad_input` for an unknown `type`.
#' @references
#' Manning, C. D., Raghavan, P., & Schütze, H. (2008). *Introduction to
#' Information Retrieval*. Cambridge University Press.
#' \doi{10.1017/CBO9780511809071}
#'
#' Porter, M. F. (2001). Snowball: A language for stemming algorithms.
#' \url{https://snowballstem.org/texts/introduction.html}
#' @examples
#' hg <- text_hypergraph(
#'   c(a = "the salt and the soup", b = "the soup and the stars"),
#'   stop_words = stop_words_en()
#' )
#' hg_get(hg, what = "vocabulary")
#' length(stop_words_en(type = "snowball"))
#' @export
stop_words_en <- function(type = c("minimal", "snowball")) {
  type <- tryCatch(match.arg(type), error = function(e) {
    .thg_bad_input("`type` must be \"minimal\" or \"snowball\"")
  })
  if (identical(type, "snowball")) {
    # the English stop list of snowballstem.org/algorithms/english/stop.txt
    return(sort(c(
    "i", "me", "my", "myself", "we", "our", "ours", "ourselves", "you",
    "your", "yours", "yourself", "yourselves", "he", "him", "his", "himself",
    "she", "her", "hers", "herself", "it", "its", "itself", "they", "them",
    "their", "theirs", "themselves", "what", "which", "who", "whom", "this",
    "that", "these", "those", "am", "is", "are", "was", "were", "be", "been",
    "being", "have", "has", "had", "having", "do", "does", "did", "doing",
    "would", "should", "could", "ought", "i'm", "you're", "he's", "she's",
    "it's", "we're", "they're", "i've", "you've", "we've", "they've", "i'd",
    "you'd", "he'd", "she'd", "we'd", "they'd", "i'll", "you'll", "he'll",
    "she'll", "we'll", "they'll", "isn't", "aren't", "wasn't", "weren't",
    "hasn't", "haven't", "hadn't", "doesn't", "don't", "didn't", "won't",
    "wouldn't", "shan't", "shouldn't", "can't", "cannot", "couldn't",
    "mustn't", "let's", "that's", "who's", "what's", "here's", "there's",
    "when's", "where's", "why's", "how's", "a", "an", "the", "and", "but",
    "if", "or", "because", "as", "until", "while", "of", "at", "by", "for",
    "with", "about", "against", "between", "into", "through", "during",
    "before", "after", "above", "below", "to", "from", "up", "down", "in",
    "out", "on", "off", "over", "under", "again", "further", "then", "once",
    "here", "there", "when", "where", "why", "how", "all", "any", "both",
    "each", "few", "more", "most", "other", "some", "such", "no", "nor",
    "not", "only", "own", "same", "so", "than", "too", "very"
    )))
  }
  sort(c(
    "the", "a", "an", "and", "or", "of", "to", "in", "on", "for", "with",
    "by", "from", "as", "at", "this", "that", "these", "those", "is",
    "are", "was", "were", "be", "been", "being", "it", "its", "we", "our",
    "their", "they", "has", "have", "had", "not", "no", "but", "which",
    "who", "during", "into", "through", "between", "among", "also", "can",
    "could", "may", "will", "would", "than", "then", "there", "here",
    "such", "more", "most", "other", "both", "each", "all", "some", "any",
    "how", "what", "when", "where", "while", "because", "about"
  ))
}
