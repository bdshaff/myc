// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// myc_proj
NumericVector myc_proj(NumericVector u, NumericVector a);
RcppExport SEXP _myc_myc_proj(SEXP uSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(myc_proj(u, a));
    return rcpp_result_gen;
END_RCPP
}
// myc_matmult
NumericMatrix myc_matmult(NumericMatrix A, NumericMatrix B);
RcppExport SEXP _myc_myc_matmult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(myc_matmult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// myc_diag
NumericVector myc_diag(NumericMatrix A);
RcppExport SEXP _myc_myc_diag(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(myc_diag(A));
    return rcpp_result_gen;
END_RCPP
}
// myc_qr
List myc_qr(NumericMatrix A);
RcppExport SEXP _myc_myc_qr(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(myc_qr(A));
    return rcpp_result_gen;
END_RCPP
}
// myc_eigen
List myc_eigen(NumericMatrix A, double margin);
RcppExport SEXP _myc_myc_eigen(SEXP ASEXP, SEXP marginSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type margin(marginSEXP);
    rcpp_result_gen = Rcpp::wrap(myc_eigen(A, margin));
    return rcpp_result_gen;
END_RCPP
}
// myc_dist
List myc_dist(NumericMatrix A);
RcppExport SEXP _myc_myc_dist(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(myc_dist(A));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _myc_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_myc_myc_proj", (DL_FUNC) &_myc_myc_proj, 2},
    {"_myc_myc_matmult", (DL_FUNC) &_myc_myc_matmult, 2},
    {"_myc_myc_diag", (DL_FUNC) &_myc_myc_diag, 1},
    {"_myc_myc_qr", (DL_FUNC) &_myc_myc_qr, 1},
    {"_myc_myc_eigen", (DL_FUNC) &_myc_myc_eigen, 2},
    {"_myc_myc_dist", (DL_FUNC) &_myc_myc_dist, 1},
    {"_myc_rcpp_hello_world", (DL_FUNC) &_myc_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_myc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
