/*
 * Expected counts of the KL-NMF topic model (hg_topics()) at the non-zero
 * cells of the document x word matrix: for every cell (i, j),
 *   max(sum_t Wt[t, i] * Htt[t, j], epsilon).
 * This is exactly what the R expression
 *   pmax(colSums(Wt[, i] * Htt[, j]), epsilon)
 * computes, without gathering two topics x cells matrices: each product is
 * rounded to double on its own (as R's `*` does), and the products are
 * added in the accumulator colSums() uses on this build of R -- long double
 * where R has it (`long_double = TRUE`), double otherwise -- in topic order,
 * so the values are bit-identical to the R expression. Four cells are
 * summed side by side, each in its own accumulator and its own order, so
 * their dependent additions overlap without reordering any of them.
 */
#include <R.h>
#include <Rinternals.h>

/* Never fuse a multiply and an add into one rounding (FMA): R rounds every
   product to double before colSums() adds it. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

#define HG_FITTED_KERNEL(ACC)                                                \
  {                                                                          \
    R_xlen_t c = 0;                                                          \
    for (; c + 3 < n; c += 4) {                                              \
      const double *w0 = w + (R_xlen_t) k * (row[c] - 1);                    \
      const double *w1 = w + (R_xlen_t) k * (row[c + 1] - 1);                \
      const double *w2 = w + (R_xlen_t) k * (row[c + 2] - 1);                \
      const double *w3 = w + (R_xlen_t) k * (row[c + 3] - 1);                \
      const double *h0 = h + (R_xlen_t) k * (col[c] - 1);                    \
      const double *h1 = h + (R_xlen_t) k * (col[c + 1] - 1);                \
      const double *h2 = h + (R_xlen_t) k * (col[c + 2] - 1);                \
      const double *h3 = h + (R_xlen_t) k * (col[c + 3] - 1);                \
      ACC s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;                            \
      for (int t = 0; t < k; t++) {                                          \
        double p0 = w0[t] * h0[t];                                           \
        double p1 = w1[t] * h1[t];                                           \
        double p2 = w2[t] * h2[t];                                           \
        double p3 = w3[t] * h3[t];                                           \
        s0 += p0;                                                            \
        s1 += p1;                                                            \
        s2 += p2;                                                            \
        s3 += p3;                                                            \
      }                                                                      \
      value[c] = (double) s0;                                                \
      value[c + 1] = (double) s1;                                            \
      value[c + 2] = (double) s2;                                            \
      value[c + 3] = (double) s3;                                            \
    }                                                                        \
    for (; c < n; c++) {                                                     \
      const double *wc = w + (R_xlen_t) k * (row[c] - 1);                    \
      const double *hc = h + (R_xlen_t) k * (col[c] - 1);                    \
      ACC s = 0.0;                                                           \
      for (int t = 0; t < k; t++) {                                          \
        double p = wc[t] * hc[t];                                            \
        s += p;                                                              \
      }                                                                      \
      value[c] = (double) s;                                                 \
    }                                                                        \
  }

SEXP hg_tm_fitted(SEXP Wt, SEXP Htt, SEXP i, SEXP j, SEXP epsilon,
                  SEXP long_double) {
  const int k = Rf_nrows(Wt);
  const R_xlen_t n = XLENGTH(i);
  const double *w = REAL(Wt);
  const double *h = REAL(Htt);
  const int *row = INTEGER(i);
  const int *col = INTEGER(j);
  const double floor_value = Rf_asReal(epsilon);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *value = REAL(out);
  if (Rf_asLogical(long_double) == TRUE) {
    HG_FITTED_KERNEL(long double)
  } else {
    HG_FITTED_KERNEL(double)
  }
  for (R_xlen_t c = 0; c < n; c++) {
    if (value[c] < floor_value) value[c] = floor_value;
  }
  UNPROTECT(1);
  return out;
}
