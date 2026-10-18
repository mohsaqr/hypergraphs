# ---- Non-negative matrix-factorization clustering ----------------------

# Internal symmetric NMF solver for Hayashi et al.'s RDC-Sym objective
#   min_{F >= 0} ||S - F F'||_F^2.
# The multiplicative update follows the positive/negative gradient split;
# a square-root step accounts for updating both tied factors simultaneously.
.hl_symnmf <- function(S, k, nstart = 1L, seed = NULL, max_iter = 500L,
                       tol = 1e-6, epsilon = 1e-12, parallel = FALSE,
                       n_cores = 2L) {
  n <- nrow(S)
  objective <- function(F) sum((S - tcrossprod(F))^2)
  # Every start's initial factor is drawn first, in start order (from its
  # own seed, or one after the other from the caller's stream); the updates
  # draw nothing, so the starts can then run in parallel with the result of
  # the serial run.
  scale <- sqrt(max(mean(S), epsilon) / k)
  initial <- lapply(seq_len(nstart), \(restart) {
    if (!is.null(seed)) set.seed(as.integer(seed) + restart - 1L)
    matrix(stats::runif(n * k, 0.5, 1.5) * scale, nrow = n)
  })
  fits <- .ho_apply(seq_len(nstart), \(restart) {
    F <- initial[[restart]]
    history <- objective(F)
    converged <- FALSE
    # Multiplicative updates: each step starts from the previous factor, so
    # the iterations cannot be vectorised; bounded by `max_iter`.
    for (iteration in seq_len(max_iter)) {
      numerator <- S %*% F
      denominator <- F %*% crossprod(F)
      candidate <- F * sqrt(pmax(numerator, 0) /
                              pmax(denominator, epsilon))
      value <- objective(candidate)
      # Numerical guard: backtrack to a monotone non-negative step.
      if (value > history[length(history)] * (1 + 1e-12)) {
        proposal <- candidate
        step <- 0.5
        repeat {
          candidate <- (1 - step) * F + step * proposal
          value <- objective(candidate)
          if (value <= history[length(history)] * (1 + 1e-12) ||
              step < 2^-20) break
          step <- step / 2
        }
      }
      history <- c(history, value)
      change <- abs(history[length(history) - 1L] - value) /
        max(history[length(history) - 1L], epsilon)
      F <- pmax(candidate, epsilon)
      if (change <= tol) {
        converged <- TRUE
        break
      }
    }
    list(factor = F, objective = utils::tail(history, 1L),
         objective_history = history, iterations = iteration,
         converged = converged, restart = restart)
  }, parallel = parallel, n_cores = n_cores)
  # the first start with the lowest objective, as the serial loop kept it
  fits[[which.min(vapply(fits, \(f) f$objective, numeric(1L)))]]
}
