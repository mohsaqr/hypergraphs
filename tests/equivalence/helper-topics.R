# Reference code for the topic-model tests (test-hg_topics.R, test-hg_topic_search.R). Base R only; nothing here
# calls hypergraphs. The functions take the package's tables as plain data frames.

#' A corpus with planted, separable topics
#'
#' `k` topics over `k * anchors` anchor words (each used by one topic only, 70% of the topic's mass, spread evenly)
#' and `shared` words every topic uses (30% of its mass, exponential weights). Half the documents are pure (one
#' topic each, in turn); the other half mix the topics with weights uniform on the simplex. Anchor words and pure
#' documents make the factorisation identifiable up to the order of the topics (separability), so the maximum-
#' likelihood topics converge to the planted ones as the documents grow.
#'
#' @param seed Seed of the draws.
#' @param n_docs,len Number of documents and tokens per document.
#' @param k,anchors,shared Topics, anchor words per topic, shared words.
#' @return A list: `X` (documents x words counts), `H` (planted P(w | z), k x words), `theta` (documents x k),
#'   `words`, `memberships` (long table `doc`, `word`, `count` of the non-zero cells).
ref_topic_corpus <- function(seed, n_docs = 60L, len = 1000L, k = 3L, anchors = 6L, shared = 12L) {
  stopifnot("`n_docs` must be even" = n_docs %% 2L == 0L)
  set.seed(seed)
  words <- c(sprintf("t%dw%02d", rep(seq_len(k), each = anchors), rep(seq_len(anchors), k)), sprintf("sw%02d", seq_len(shared)))
  v <- length(words)
  H <- t(vapply(seq_len(k), \(z) {
    p <- numeric(v)
    p[(z - 1L) * anchors + seq_len(anchors)] <- 0.7 / anchors
    common <- k * anchors + seq_len(shared)
    weights <- stats::rexp(shared)
    p[common] <- 0.3 * weights / sum(weights)
    p
  }, numeric(v)))
  dimnames(H) <- list(seq_len(k), words)
  dominant <- rep(seq_len(k), length.out = n_docs)
  theta <- t(vapply(seq_len(n_docs), \(d) {
    if (d <= n_docs / 2L) return(as.numeric(seq_len(k) == dominant[[d]]))
    g <- stats::rexp(k)
    g / sum(g)
  }, numeric(k)))
  probabilities <- theta %*% H
  X <- t(vapply(seq_len(n_docs), \(d) as.numeric(stats::rmultinom(1L, len, probabilities[d, ])), numeric(v)))
  dimnames(X) <- list(sprintf("d%03d", seq_len(n_docs)), words)
  cells <- which(X > 0, arr.ind = TRUE)
  list(X = X, H = H, theta = theta, words = words,
       memberships = data.frame(doc = rownames(X)[cells[, 1L]], word = words[cells[, 2L]], count = X[cells]))
}

#' A wide matrix from a long table
#'
#' @param table A data frame. @param row,column,value Its column names.
#' @param rows,columns The row and column names of the result, in order; absent cells are 0.
#' @return A numeric matrix.
ref_wide <- function(table, row, column, value, rows, columns) {
  out <- matrix(0, length(rows), length(columns), dimnames = list(rows, columns))
  out[cbind(match(table[[row]], rows), match(table[[column]], columns))] <- table[[value]]
  out
}

#' Generalised Kullback-Leibler divergence D(X || X_hat) (Lee and Seung 2001)
#'
#' sum over cells of x log(x / x_hat) - x + x_hat, with 0 log 0 = 0.
#'
#' @param X,fitted Matrices of the same shape.
#' @return One number.
ref_kl <- function(X, fitted) {
  positive <- X > 0
  sum(X[positive] * log(X[positive] / fitted[positive])) - sum(X) + sum(fitted)
}

#' The fitted counts of probabilistic latent semantic analysis from the topic shares and word probabilities
#'
#' X_hat[d, w] = n_d sum_z P(z | d) P(w | z), n_d the document's token count (Hofmann 1999; Gaussier and Goutte
#' 2005: the KL factorisation is the PLSA maximum-likelihood fit, whose row totals are the documents' counts at
#' the fixed point of the multiplicative updates).
#'
#' @param X The counts. @param shares Long table `node`, `topic`, `share`. @param words Long table `topic`,
#'   `word`, `probability`.
#' @return The fitted matrix, shaped as `X`.
ref_plsa_fitted <- function(X, shares, words) {
  topics <- sort(unique(shares$topic))
  theta <- ref_wide(shares, "node", "topic", "share", rownames(X), topics)
  phi <- ref_wide(words, "topic", "word", "probability", topics, colnames(X))
  rowSums(X) * (theta %*% phi)
}

#' UMass coherence of one topic (Mimno et al. 2011)
#'
#' sum over m = 2..M and l = 1..m-1 of log((D(w_m, w_l) + 0.01) / (D(w_l) + 0.01)), D the number of documents that
#' contain the word(s), the words in the topic's rank order; the smoothing 0.01 is stm's.
#'
#' @param present Documents x words presence (0/1). @param top The topic's top words, in rank order.
#' @return One number.
ref_umass <- function(present, top) {
  co <- crossprod(present[, top, drop = FALSE])
  pairs <- utils::combn(length(top), 2L)
  sum(log(co[cbind(pairs[2L, ], pairs[1L, ])] + 0.01) - log(diag(co)[pairs[1L, ]] + 0.01))
}

#' The Pareto frontier of two criteria, both to be maximised
#'
#' A point is on the frontier when no other point is larger on both (Roberts et al. 2014).
#'
#' @param a,b Numeric vectors of equal length.
#' @return A logical vector.
ref_frontier <- function(a, b) {
  vapply(seq_along(a), \(i) !any(a > a[[i]] & b > b[[i]]), logical(1))
}
