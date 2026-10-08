#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP hg_tm_fitted(SEXP Wt, SEXP Htt, SEXP i, SEXP j, SEXP epsilon,
                  SEXP long_double);

static const R_CallMethodDef call_methods[] = {
  {"hg_tm_fitted", (DL_FUNC) &hg_tm_fitted, 6},
  {NULL, NULL, 0}
};

void R_init_hypergraphs(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
