//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Gibbs update for dyadic variance with independent replicate relational data
//' 
//' Gibbs update for dyadic variance with independent replicate relational data
//' 
//' 
//' @usage rs2_rep_fc_cpp(E.T, rho)
//' @param E.T Array of square residual relational matrix series. The third
//' dimension of the array is for different replicates. Each slice of the array
//' according to the third dimension is a square residual relational matrix
//' @param rho s2 matrix
//' @return a new value of s2
//' @author Peter Hoff, Yanjun He, Shahryar Minhas
//' @export rs2_rep_fc_cpp
// [[Rcpp::export]]

double rs2_rep_fc_cpp(
	arma::cube ET, arma::mat rhoMat
	) {
	
	int N = ET.n_slices;

	arma::vec eigVal; arma::mat eigVec;
	arma::eig_sym(eigVal, eigVec, rhoMat);
	arma::mat eigValDiagMat = pow(diagmat(eigVal), .5);
	arma::mat H = eigVec * eigValDiagMat * eigVec.t();

	arma::mat tmp = trimatl(arma::ones(ET.n_rows, ET.n_cols));
	arma::uvec tmpindex = find(tmp==0);
	arma::mat EM;
	for(int t=0 ; t<N ; ++t){
		arma::mat E = ET.slice(t); arma::mat Et = E.t();
		arma::vec eUpperTri = E.elem( tmpindex );
		arma::vec etUpperTri = Et.elem( tmpindex );
		EM = join_cols(EM, join_rows(eUpperTri, etUpperTri) * H);
	}

	double shape = EM.n_elem + 1 ;
	double shape2 = shape/2;
	double scale = accu( pow(EM, 2) )/2;
	double s2 = 1/R::rgamma(shape2, 1/scale);
	return( s2 );
}