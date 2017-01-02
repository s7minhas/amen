//Includes/namespaces
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::export]]
IntegerVector sample_num(IntegerVector x, int size, bool replace, 
	NumericVector prob = NumericVector::create() ) {
  IntegerVector ret = Rcpp::RcppArmadillo::sample(x, size, replace, prob) ;
  return ret ;
}

// [[Rcpp::export]]
bool matMultVec(arma::mat x, arma::vec y){
	int nRows = x.n_rows;
	int nCols = x.n_cols;

	arma::mat xy = arma::zeros(nRows,nCols);
	if(nRows <= nCols){
		for(int r = 0 ; r<nRows ; r++){
		  xy.row(r) = x.row(r) % y.t();
		}		
	}
	if(nCols < nRows){
		for(int c = 0 ; c<nCols ; c++){
		  xy.col(c) = x.col(c) % y;
		}		
	}		

	return(xy);
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

	Rcpp::IntegerVector loopIDs = rep(
		sample_num( seq_len(U.n_rows), U.n_rows, FALSE  ), 4);

	// int i = loopIDs[2];
	arma::mat eui = arma::zeros(n,R);
	arma::vec erow = E.row(9).t();
	for(int r = 0 ; r<R ; r++){
	  eui.col(r) = U.col(r) % erow;
	}

	arma::rowvec euicolsum = sum(eui,0);

	arma::mat l = L % trans((euicolsum - U.row(9) * E(9,9))/s2);

	return( L );	
}
