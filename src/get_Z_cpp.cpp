//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Simulate Z based on a probit model
//' 
//' Simulates a random latent matrix Z given its expectation, dyadic correlation
//' and a binary relational matrix Y
//' 
//' 
//' @usage rZ_bin_fc_cpp(Z, EZ, rho, Y)
//' @param Z a square matrix, the current value of Z
//' @param EZ expected value of Z
//' @param rho dyadic correlation
//' @param Y square binary relational matrix
//' @return a square matrix , the new value of Z
//' @author Peter Hoff, Shahryar Minhas
//' @export rZ_bin_fc_cpp
// [[Rcpp::export]]

Rcpp::List rZ_bin_fc_cpp(
	arma::mat Z, arma::mat EZ, double rho, arma::mat Y
	) {

	double sz = pow( ( 1 - pow(rho,2) ), .5);
	arma::uvec ut = find(trimatl(arma::ones(EZ.n_rows, EZ.n_cols))==0);
	arma::uvec lt = find(trimatu(arma::ones(EZ.n_rows, EZ.n_cols))==0);
	Y.replace(datum::nan, -1);	

	int y=-1;
	double lb = -1*datum::inf; double ub = datum::inf;
	IntegerVector Yyindex = as<IntegerVector>(wrap( find(Y==y) ));
	IntegerVector utindex = as<IntegerVector>(wrap( ut ));	
	arma::uvec up = as<arma::uvec>(wrap( intersect(Yyindex, utindex) ));
	arma::mat Zt = Z.t(); arma::mat EZt = EZ.t();
	arma::vec ez = EZ.elem(up) + rho * ( Zt.elem(up) - EZt.elem(up) );

	// arma::vec ezUp = EZ.elem( ut ); arma::vec yy = Y.elem( ut );
	// arma::mat EZt = EZ.t(); arma::vec eztUp = EZt.elem( ut );
	// arma::vec ezUpY = ezUp.elem(find(yy==y)); arma::vec eztUpY = eztUp.elem(find(yy==y));
	// arma::mat Zt = Z.t(); arma::vec ztUp = Zt.elem( ut ); arma::vec ztUpY = ztUp.elem(find(yy==y));
	// arma::vec ez = ezUpY + rho * ( ztUpY - eztUpY );

	NumericVector x1; x1 = (lb-ez)/sz; NumericVector x2; x2 = (ub-ez)/sz;
	NumericVector rUnifLo = pnorm(x1,0.0,1.0); NumericVector rUnifHi = pnorm(x2,0.0,1.0);
	NumericVector rUnifDraw(up.n_elem);
	for(int ii=0;ii<up.n_elem;++ii){ rUnifDraw[ii]=runif(1,rUnifLo[ii],rUnifHi[ii])[0]; }
	NumericVector qnormDraw = qnorm(rUnifDraw,0.0,1.0); arma::vec szQnorm = sz * qnormDraw;
	Z.elem(up) = ez + szQnorm;

	return(
		Rcpp::List::create(
			Rcpp::Named("sz") = sz,
			Rcpp::Named("ez") = ez,
			Rcpp::Named("rUnifLo") = rUnifLo,
			Rcpp::Named("rUnifHi") = rUnifHi,
			Rcpp::Named("rUnifDraw") = rUnifDraw,
			Rcpp::Named("qnormDraw") = qnormDraw,
			// Rcpp::Named("newZ") = newZ,
			// Rcpp::Named("test") = test,
			// Rcpp::Named("ut") = ut,
			Rcpp::Named("Z") = Z
			)
		);
}