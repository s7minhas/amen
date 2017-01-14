//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Generates parameters to aid in Gibbs sampling of regression coefficient
//' with independent replicate relational data
//' 
//' @param Z.T n x n x T array, with the third dimension for replicates.
//' Each slice of the array is a (latent) normal relational matrix, with
//' multiplicative effects subtracted out
//' @param XrCube n x p x T row covariate array generated within ame_repL fn
//' @param XcCube n x p x T column covariate array generated within ame_repL fn
//' @param mXCube n^2 x p x T design array generated within ame_repL fn
//' @param mXtCube n^2 x p x T transposed design array generated within ame_repL fn
//' @param xxCube p x p x T regression sum of squares array generated within ame_repl fn
//' @param xxTCube p x p x T transposed regression sum of squares array generated 
//' within ame_repl fn
//' @param iSe2 variance matrix
//' @return list of parameters to feed into gibbs sampling of covariate,
//' a random effects, and b random effects
//' @author Shahryar Minhas
//' @export rbeta_rep_cpp
// [[Rcpp::export]]

List rbeta_rep_cpp(
	arma::cube zCube, arma::cube XrCube, arma::cube XcCube, 
	arma::cube mXCube, arma::cube mXtCube, 
	arma::cube xxCube, arma::cube xxTCube,
	double td, double to
	) {

  int N = zCube.n_slices;  
  int n = zCube.n_rows;
  int p = XrCube.n_cols;

  arma::vec lb = arma::zeros(p);
  arma::mat Qb = arma::zeros(p,p);
  arma::vec ZrT = arma::zeros(n);
  arma::vec ZcT = arma::zeros(n);
  arma::mat XrT = arma::zeros(n,p);
  arma::mat XcT = arma::zeros(n,p);  

  for(int t=0 ; t < N ; ++t){
    arma::mat mXs = td * mXCube.slice(t) + to*mXtCube.slice(t);
    arma::mat XXs = (pow(to, 2)+pow(td, 2))*xxCube.slice(t) + 2*to*td*xxTCube.slice(t);
    arma::mat Zs = td*zCube.slice(t) + to*zCube.slice(t).t();

    arma::vec zr = sum(zCube.slice(t),1);
    arma::vec zc = sum(zCube.slice(t),0).t();

    if(p > 0){
      lb = lb + (mXs.t() * vectorise(Zs));
      Qb = Qb + XXs + ( xxCube.slice(t)/mXs.n_rows )/zCube.n_slices;
    }

    arma::mat Xsr = td * XrCube.slice(t) + to*XcCube.slice(t);
    arma::mat Xsc = td * XcCube.slice(t) + to*XrCube.slice(t);

    ZrT = ZrT + zr;
    ZcT = ZcT + zc;
    XrT = XrT + Xsr;
    XcT = XcT + Xsc;
  }

return(
  Rcpp::List::create(
    Rcpp::Named("lb")=lb,
    Rcpp::Named("Qb")=Qb,
    Rcpp::Named("ZrT")=ZrT,
    Rcpp::Named("ZcT")=ZcT,
    Rcpp::Named("XrT")=XrT,
    Rcpp::Named("XcT")=XcT
    )
  );
}
