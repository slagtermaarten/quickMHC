// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// match_seq_ext_Rcpp
Rcpp::LogicalVector match_seq_ext_Rcpp(const Rcpp::CharacterVector& query_peps, const Rcpp::CharacterVector& ref_list, const Rcpp::NumericMatrix& scorematrix, const bool debug);
RcppExport SEXP _quickMHC_match_seq_ext_Rcpp(SEXP query_pepsSEXP, SEXP ref_listSEXP, SEXP scorematrixSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type query_peps(query_pepsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ref_list(ref_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type scorematrix(scorematrixSEXP);
    Rcpp::traits::input_parameter< const bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(match_seq_ext_Rcpp(query_peps, ref_list, scorematrix, debug));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_quickMHC_match_seq_ext_Rcpp", (DL_FUNC) &_quickMHC_match_seq_ext_Rcpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_quickMHC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
