// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_EZ_cpp
arma::cube get_EZ_cpp(Rcpp::List Xlist, arma::vec beta, arma::mat ab, arma::mat U, arma::mat V);
RcppExport SEXP amen_get_EZ_cpp(SEXP XlistSEXP, SEXP betaSEXP, SEXP abSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Xlist(XlistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ab(abSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(get_EZ_cpp(Xlist, beta, ab, U, V));
    return rcpp_result_gen;
END_RCPP
}
// rbeta_ab_rep_fc_cpp
List rbeta_ab_rep_fc_cpp(arma::cube zCube, arma::cube XrCube, arma::cube XcCube, arma::cube mXCube, arma::cube mXtCube, arma::cube xxCube, arma::cube xxTCube, arma::mat iSe2, arma::mat Sabs, int k, arma::mat G, arma::mat e, arma::vec colE);
RcppExport SEXP amen_rbeta_ab_rep_fc_cpp(SEXP zCubeSEXP, SEXP XrCubeSEXP, SEXP XcCubeSEXP, SEXP mXCubeSEXP, SEXP mXtCubeSEXP, SEXP xxCubeSEXP, SEXP xxTCubeSEXP, SEXP iSe2SEXP, SEXP SabsSEXP, SEXP kSEXP, SEXP GSEXP, SEXP eSEXP, SEXP colESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type zCube(zCubeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type XrCube(XrCubeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type XcCube(XcCubeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type mXCube(mXCubeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type mXtCube(mXtCubeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type xxCube(xxCubeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type xxTCube(xxTCubeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type iSe2(iSe2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sabs(SabsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type colE(colESEXP);
    rcpp_result_gen = Rcpp::wrap(rbeta_ab_rep_fc_cpp(zCube, XrCube, XcCube, mXCube, mXtCube, xxCube, xxTCube, iSe2, Sabs, k, G, e, colE));
    return rcpp_result_gen;
END_RCPP
}
// rrho_mh_rep_cpp
double rrho_mh_rep_cpp(arma::cube ET, double rho, double s2);
RcppExport SEXP amen_rrho_mh_rep_cpp(SEXP ETSEXP, SEXP rhoSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type ET(ETSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(rrho_mh_rep_cpp(ET, rho, s2));
    return rcpp_result_gen;
END_RCPP
}
// rs2_rep_fc_cpp
double rs2_rep_fc_cpp(arma::cube ET, arma::mat rhoMat);
RcppExport SEXP amen_rs2_rep_fc_cpp(SEXP ETSEXP, SEXP rhoMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type ET(ETSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type rhoMat(rhoMatSEXP);
    rcpp_result_gen = Rcpp::wrap(rs2_rep_fc_cpp(ET, rhoMat));
    return rcpp_result_gen;
END_RCPP
}
// rUV_sym_fc_cpp
arma::mat rUV_sym_fc_cpp(arma::mat E, arma::mat U, arma::mat V, double s2, bool shrink);
RcppExport SEXP amen_rUV_sym_fc_cpp(SEXP ESEXP, SEXP USEXP, SEXP VSEXP, SEXP s2SEXP, SEXP shrinkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type E(ESEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< bool >::type shrink(shrinkSEXP);
    rcpp_result_gen = Rcpp::wrap(rUV_sym_fc_cpp(E, U, V, s2, shrink));
    return rcpp_result_gen;
END_RCPP
}
// Xbeta_cpp
arma::mat Xbeta_cpp(arma::cube X, arma::vec beta);
RcppExport SEXP amen_Xbeta_cpp(SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(Xbeta_cpp(X, beta));
    return rcpp_result_gen;
END_RCPP
}
