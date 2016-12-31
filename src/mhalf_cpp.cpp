//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Symmetric square root of a matrix
//' 
//' Computes the symmetric square root of a positive definite matrix
//' 
//' 
//' @usage mhalf_cpp(M)
//' @param M a positive definite matrix
//' @return a matrix \code{H} such that \code{H^2} equals \code{M}
//' @author Peter Hoff, Shahryar Minhas
//' @export mhalf_cpp
// [[Rcpp::export]]

arma::mat mhalf_cpp(
	arma::mat M
	) {

	arma::vec eigVal;
	arma::mat eigVec;
	arma::eig_sym(eigVal, eigVec, M);

	arma::mat eigValDiagMat = pow(diagmat(eigVal), .5);
	arma::mat tmp = eigVec * eigValDiagMat * eigVec.t();

	return(tmp);
}