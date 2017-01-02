//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Metropolis update for dyadic correlation with independent replicate data
//' 
//' Metropolis update for dyadic correlation with independent replicate data. 
//' 
//' 
//' @usage rrho_mh_rep_cpp(E.T, rho, s2 = 1)
//' @param E.T Array of square residual relational matrix series. The third
//' dimension of the array is for different replicates. Each slice of the array
//' according to the third dimension is a square residual relational matrix. 
//' @param rho current value of rho
//' @param s2 current value of s2
//' @return a new value of rho
//' @author Peter Hoff, Yanjun He, Shahryar Minhas
//' @export rrho_mh_rep_cpp
// [[Rcpp::export]]

double rrho_mh_rep_cpp(
	arma::cube ET, double rho, double s2
	) {

	int N = ET.n_slices;
	arma::mat tmp = trimatl(arma::ones(ET.n_rows, ET.n_cols));
	arma::uvec tmpindex = find(tmp==0);
	arma::mat EM;
	
	for(int t=0 ; t<N ; ++t){
		arma::mat E = ET.slice(t); arma::mat Et = E.t();
		arma::vec eUpperTri = E.elem( tmpindex );
		arma::vec etUpperTri = Et.elem( tmpindex );
		EM = join_cols(EM, join_rows(eUpperTri, etUpperTri) / pow(s2, .5));
	}

	double emcp = accu( EM.col(0) % EM.col(1) );
	double emss = accu( pow(EM, 2) );

	int m = EM.n_rows;
	arma::mat emCor = arma::cor(EM);
	double sr = 2*( 1 - pow(emCor(0,1), 2) )/pow(m, .5);

	double rho1 = R::pnorm(1, (1-rho)/sr, 1, 0, FALSE );	

	return(rho1);
}