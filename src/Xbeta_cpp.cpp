//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Linear combinations of submatrices of an array
//' 
//' Computes a matrix of expected values based on an array X of predictors and a
//' vector beta of regression coefficients.
//' 
//' 
//' @usage Xbeta_cpp(X, beta)
//' @param X an n by n by p array
//' @param beta a p by 1 vector
//' @return An n by n matrix
//' @author Peter Hoff, Shahryar Minhas
//' @export Xbeta_cpp
// [[Rcpp::export]]

arma::mat Xbeta_cpp(
	arma::cube X, arma::vec beta
	) {

	int n = X.n_rows;
	int p = beta.n_elem;
	arma::mat XB = arma::zeros(n, n);
	for (int k=0 ; k < p ; ++k){
		XB = XB + beta(k) * X.slice(k);
	}
  
	return( XB );
}
