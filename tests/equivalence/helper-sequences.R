# Sequences with memory and k-gram counts, shared by the memory-network tests. Base R only; nothing here calls
# hypergraphs.

#' Sequences of a second-order chain on A, B, C
#'
#' After A then B comes C, after C then B comes A; every other step draws A, B, C with probabilities 0.3, 0.4, 0.3.
#' @param n_sequences,length Number and length of the sequences.
#' @param seed The seed.
#' @return A list of character vectors.
ref_second_order <- function(n_sequences = 30L, length = 15L, seed = 7L) {
  set.seed(seed)
  lapply(seq_len(n_sequences), \(i) {
    Reduce(\(so_far, t) {
      last_two <- utils::tail(so_far, 2L)
      following <- if (identical(last_two, c("A", "B"))) "C" else if (identical(last_two, c("C", "B"))) "A" else
        sample(c("A", "B", "C"), 1L, prob = c(0.3, 0.4, 0.3))
      c(so_far, following)
    }, seq_len(length - 2L), sample(c("A", "B", "C"), 2L, replace = TRUE))
  })
}

#' Random first-order sequences
#'
#' @param n_sequences,length,states,seed As named.
#' @return A list of character vectors.
ref_random_sequences <- function(n_sequences, length, states, seed) {
  set.seed(seed)
  lapply(seq_len(n_sequences), \(i) sample(states, length, replace = TRUE))
}

#' Every k-gram of a list of sequences, counted
#'
#' @param sequences A list of character vectors.
#' @param k The number of states in a k-gram.
#' @param separator What joins the states of a k-gram.
#' @return A named numeric vector of counts.
ref_kgrams <- function(sequences, k, separator = " -> ") {
  grams <- unlist(lapply(Filter(\(s) length(s) >= k, sequences), \(s) {
    starts <- seq_len(length(s) - k + 1L)
    vapply(starts, \(start) paste(s[start:(start + k - 1L)], collapse = separator), character(1))
  }))
  counted <- table(grams)
  stats::setNames(as.numeric(counted), names(counted))
}
