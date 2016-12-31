//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Simulation from a multivariate normal distribution
//' 
//' Simulates a matrix where the rows are i.i.d. samples from a multivariate
//' normal distribution
//' 
//' 
//' @usage rmvnorm_cpp(n, mu, Sigma)
//' @param n sample size
//' @param mu multivariate mean vector
//' @param Sigma covariance matrix
//' @return a matrix with \code{n} rows
//' @author Peter Hoff, Shahryar Minhas
//' @export rmvnorm_cpp
// [[Rcpp::export]]

arma::mat rmvnorm_cpp(
	int n, arma::vec mu, arma::mat Sigma
	) {

	arma::mat E = rnorm( n * mu.size() ) ; E.reshape(n, mu.size());
	arma::mat tmp = ( E * chol(Sigma) ).t();
	for(int i=0 ; i < tmp.n_rows ; i++){
		tmp.row(i) = tmp.row(i) + mu[i];
	}

	return( tmp );
}