//Includes/namespaces
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

IntegerVector sample_num(IntegerVector x, int size, bool replace, 
	NumericVector prob = NumericVector::create() ) {
  IntegerVector ret = Rcpp::RcppArmadillo::sample(x, size, replace, prob) ;
  return ret ;
}

arma::mat matMultVec(arma::mat x, arma::vec y){
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
//' rUV_sym_fc_cpp(E,U0,V0) 
//' 
//' @export rUV_sym_fc_cpp
// [[Rcpp::export]]

List rUV_sym_fc_cpp(
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
			ivU[r] = rgamma( shape, 1/scale[r] )[0];
		}
		ivDiagMat = diagmat(ivU);
	}

	double ugh = n; arma::vec tmp(R) ; tmp.fill( 1/ugh );
	if(!shrink){
		ivDiagMat = diagmat(tmp);
	}

	// Rcpp::IntegerVector loopIDs = rep(
	// 	sample_num( seq_len(U.n_rows), U.n_rows, FALSE  ), 4);

	// need to add in loops ids to simulate for(i in rep(sample(1:n),4))
	// just index loopIDs in line 90 before running fn code
	for(int i = 0; i<n ; i++){
	// for(int s = 0; s<loopIDs ; s++){
	// 	int i = loopIDs[s];
		arma::vec erow = E.row(i).t();
		arma::mat eui = matMultVec(U, erow);
	
		arma::rowvec euicolsum = sum(eui,0);
	
		arma::mat l = L * trans((euicolsum - U.row(i) * E(i,i))/s2);
		arma::mat iQ = inv( ivDiagMat + L * ( (U.t() * U) - (U.row(i).t() * U.row(i)) ) * L/s2 );
		arma::vec randNormDraw = rnorm(R);
		U.row(i) = trans(iQ * l + trans(chol(iQ)) * randNormDraw);
	}
	
	arma::mat tmponesmat = trimatl(arma::ones(E.n_rows, E.n_cols));
	arma::uvec tmpindex = find(tmponesmat==0);	

	// arma::vec l_vec; arma::vec iq_vec;
	// for(int r = 0 ; r<R ; r++){
	int r = 0;
		arma::mat Usmall = U; Usmall.shed_col(r);
		arma::mat Lsmall = L; Lsmall.shed_col(r); Lsmall.shed_row(r);
		arma::mat Er = E - Usmall * Lsmall * Usmall.t();
		arma::mat lMat = Er % ( U.col(r) * U.col(r).t() );
		double l = accu(lMat.elem( tmpindex ))/s2;
		arma::mat uut2 = pow(U.col(r) * U.col(r).t(), 2);
		double iq = 1/(1+accu(uut2.elem(tmpindex))/s2);
		// l_vec[r]=l; iq_vec[r]=iq;
		NumericVector rand = rnorm(iq*l, pow(iq,.5) );
		L(r,r) = rand[0];
	// }

	return(
		Rcpp::List::create(
			Rcpp::Named("U")=U,
			Rcpp::Named("V")=U*L,
			Rcpp::Named("Usmall")=Usmall,
			Rcpp::Named("Lsmall")=Lsmall,
			Rcpp::Named("Er")=Er,
			Rcpp::Named("lMat")=lMat,
			Rcpp::Named("l")=l,
			Rcpp::Named("iq")=iq,
			Rcpp::Named("rand")=rand,
			Rcpp::Named("L")=L
			// Rcpp::Named("l_vec")=l_vec,
			// Rcpp::Named("iq_vec")=iq_vec
			// ,Rcpp::Named("loopIDs")=loopIDs
			)		
		);	
}
