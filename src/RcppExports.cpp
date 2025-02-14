// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// C_lasmerge_neighbors
NumericVector C_lasmerge_neighbors(S4 las, S4 las_from, const unsigned int k, const std::string func, const std::string var);
RcppExport SEXP _lasRANN_C_lasmerge_neighbors(SEXP lasSEXP, SEXP las_fromSEXP, SEXP kSEXP, SEXP funcSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type las(lasSEXP);
    Rcpp::traits::input_parameter< S4 >::type las_from(las_fromSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const std::string >::type func(funcSEXP);
    Rcpp::traits::input_parameter< const std::string >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(C_lasmerge_neighbors(las, las_from, k, func, var));
    return rcpp_result_gen;
END_RCPP
}
// C_lasmerge_neighbors_Rfunc
NumericVector C_lasmerge_neighbors_Rfunc(S4 las, S4 las_from, const unsigned int k, Rcpp::Function func, const std::string var);
RcppExport SEXP _lasRANN_C_lasmerge_neighbors_Rfunc(SEXP lasSEXP, SEXP las_fromSEXP, SEXP kSEXP, SEXP funcSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type las(lasSEXP);
    Rcpp::traits::input_parameter< S4 >::type las_from(las_fromSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type func(funcSEXP);
    Rcpp::traits::input_parameter< const std::string >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(C_lasmerge_neighbors_Rfunc(las, las_from, k, func, var));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_class_WANN();

static const R_CallMethodDef CallEntries[] = {
    {"_lasRANN_C_lasmerge_neighbors", (DL_FUNC) &_lasRANN_C_lasmerge_neighbors, 5},
    {"_lasRANN_C_lasmerge_neighbors_Rfunc", (DL_FUNC) &_lasRANN_C_lasmerge_neighbors_Rfunc, 5},
    {"_rcpp_module_boot_class_WANN", (DL_FUNC) &_rcpp_module_boot_class_WANN, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_lasRANN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
