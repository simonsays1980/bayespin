// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// simulateEKOP_cc
arma::imat simulateEKOP_cc(const unsigned int nobs, const double alpha, const double epsilon, const double delta, const double mu, const double T);
RcppExport SEXP _bayespin_simulateEKOP_cc(SEXP nobsSEXP, SEXP alphaSEXP, SEXP epsilonSEXP, SEXP deltaSEXP, SEXP muSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateEKOP_cc(nobs, alpha, epsilon, delta, mu, T));
    return rcpp_result_gen;
END_RCPP
}
// simulateEKOPMis_cc
arma::imat simulateEKOPMis_cc(const unsigned int nobs, const double alpha, const double epsilon, const double delta, const double mu, const double mis, const double T);
RcppExport SEXP _bayespin_simulateEKOPMis_cc(SEXP nobsSEXP, SEXP alphaSEXP, SEXP epsilonSEXP, SEXP deltaSEXP, SEXP muSEXP, SEXP misSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type mis(misSEXP);
    Rcpp::traits::input_parameter< const double >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateEKOPMis_cc(nobs, alpha, epsilon, delta, mu, mis, T));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayespin_simulateEKOP_cc", (DL_FUNC) &_bayespin_simulateEKOP_cc, 6},
    {"_bayespin_simulateEKOPMis_cc", (DL_FUNC) &_bayespin_simulateEKOPMis_cc, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayespin(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}