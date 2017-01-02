//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {

    // clone a into b to leave a alone
    Rcpp::NumericVector b = Rcpp::clone(a);

    std::random_shuffle(b.begin(), b.end(), randWrapper);

    return b;
}

//' Gibbs sampling of U and V
//' 
//' A Gibbs sampler for updating the multiplicative effect matrices U and V
//' in the symmetric case. In this case \code{U\%*\%t(V)} is symmetric, so
//' this is parameterized as \code{V=U\%*\%L} where \code{L} is the 
//' diagonal matrix of eigenvalues of \code{U\%*\%t(V)}. 
//' 
//' @usage rUV_sym_fc_cpp(E, U, V, s2 = 1, shrink=TRUE)
//' @param E square residual relational matrix
//' @param U current value of U
//' @param V current value of V
//' @param s2 dyadic variance
//' @param shrink adaptively shrink the factors with a hierarchical prior
//' @return \item{U}{a new value of U} \item{V}{a new value of V}
//' @author Peter Hoff, Shahryar Minhas
//' @examples
//' 
//' U0<-matrix(rnorm(30,2),30,2) ; V0<-U0%*%diag(c(3,-2)) 
//' E<- U0%*%t(V0) + matrix(rnorm(30^2),30,30) 
//' rUV_sym_fc(E,U0,V0) 
//' 
//' @export rUV_sym_fc_cpp
// [[Rcpp::export]]

arma::mat rUV_sym_fc_cpp(
	arma::mat E, arma::mat U, arma::mat V, 
	double s2, bool shrink) {

	int R = U.n_cols; int n = U.n_rows;
	arma::mat L = diagmat( V.row(0)/U.row(0)  );
	L.replace(datum::nan, 0);

	arma::vec ivU = arma::zeros(R);
	arma::mat ivDiagMat;
	arma::rowvec scale;
	double shape;

	if(shrink){
		shape = (2+n)/2;
		scale = (1+sum(pow(U,2),0))/2;
		for(int r=0 ; r<R ; r++){
			ivU[r] = R::rgamma( shape, 1/scale[r] );
		}
		ivDiagMat = diagmat(ivU);
	}

	double ugh = n; arma::vec tmp(R) ; tmp.fill( 1/ugh );
	if(!shrink){
		ivDiagMat = diagmat(tmp);
	}
	
	return( ivDiagMat );	
}
