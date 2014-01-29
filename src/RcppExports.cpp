// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rstrauss_BD
List rstrauss_BD(double beta, double gamma, double R, NumericVector win, int toroidal, int iter, int dbg);
RcppExport SEXP rstrauss_rstrauss_BD(SEXP betaSEXP, SEXP gammaSEXP, SEXP RSEXP, SEXP winSEXP, SEXP toroidalSEXP, SEXP iterSEXP, SEXP dbgSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< double >::type R(RSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type win(winSEXP );
        Rcpp::traits::input_parameter< int >::type toroidal(toroidalSEXP );
        Rcpp::traits::input_parameter< int >::type iter(iterSEXP );
        Rcpp::traits::input_parameter< int >::type dbg(dbgSEXP );
        List __result = rstrauss_BD(beta, gamma, R, win, toroidal, iter, dbg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rstrauss_DCFTP
List rstrauss_DCFTP(double beta, double gamma, double R, NumericVector win, int toroidal, int T0, int dbg, int maxtry);
RcppExport SEXP rstrauss_rstrauss_DCFTP(SEXP betaSEXP, SEXP gammaSEXP, SEXP RSEXP, SEXP winSEXP, SEXP toroidalSEXP, SEXP T0SEXP, SEXP dbgSEXP, SEXP maxtrySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< double >::type R(RSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type win(winSEXP );
        Rcpp::traits::input_parameter< int >::type toroidal(toroidalSEXP );
        Rcpp::traits::input_parameter< int >::type T0(T0SEXP );
        Rcpp::traits::input_parameter< int >::type dbg(dbgSEXP );
        Rcpp::traits::input_parameter< int >::type maxtry(maxtrySEXP );
        List __result = rstrauss_DCFTP(beta, gamma, R, win, toroidal, T0, dbg, maxtry);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rstrauss_MH
List rstrauss_MH(int n, double gamma, double R, NumericVector win, int toroidal, int iter, int dbg);
RcppExport SEXP rstrauss_rstrauss_MH(SEXP nSEXP, SEXP gammaSEXP, SEXP RSEXP, SEXP winSEXP, SEXP toroidalSEXP, SEXP iterSEXP, SEXP dbgSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< double >::type R(RSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type win(winSEXP );
        Rcpp::traits::input_parameter< int >::type toroidal(toroidalSEXP );
        Rcpp::traits::input_parameter< int >::type iter(iterSEXP );
        Rcpp::traits::input_parameter< int >::type dbg(dbgSEXP );
        List __result = rstrauss_MH(n, gamma, R, win, toroidal, iter, dbg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
