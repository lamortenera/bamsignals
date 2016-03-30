#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP bamsignals_checkList(SEXP lSEXP, SEXP ssSEXP);
SEXP bamsignals_fastWidth(SEXP lSEXP, SEXP ssSEXP);
SEXP bamsignals_pileup_core(SEXP bampathSEXP, SEXP grSEXP, SEXP tlen_filterSEXP, SEXP mapqualSEXP, SEXP binsizeSEXP, SEXP shiftSEXP, SEXP ssSEXP, SEXP requiredFSEXP, SEXP filteredFSEXP, SEXP pe_midSEXP, SEXP maxgapSEXP);
SEXP bamsignals_coverage_core(SEXP bampathSEXP, SEXP grSEXP, SEXP tlen_filterSEXP, SEXP mapqualSEXP, SEXP requiredFSEXP, SEXP filteredFSEXP, SEXP tspanSEXP, SEXP maxgapSEXP);
SEXP bamsignals_writeSamAsBamAndIndex(SEXP sampathSEXP, SEXP bampathSEXP);

R_CallMethodDef callMethods[]  = {
  {"bamsignals_checkList", (DL_FUNC) &bamsignals_checkList, 2},
  {"bamsignals_fastWidth", (DL_FUNC) &bamsignals_fastWidth, 2},
  {"bamsignals_pileup_core", (DL_FUNC) &bamsignals_pileup_core, 11},
  {"bamsignals_coverage_core", (DL_FUNC) &bamsignals_coverage_core, 8},
  {"bamsignals_writeSamAsBamAndIndex", (DL_FUNC) &bamsignals_writeSamAsBamAndIndex, 2},
  {NULL, NULL, 0}
};

void R_init_bamsignals(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

