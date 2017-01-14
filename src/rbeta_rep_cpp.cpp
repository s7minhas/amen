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
  arma::cube ZT, double to, double td, 
  arma::cube Xr, arma::cube Xc, arma::cube mX,
  arma::cube mXt, arma::cube XX, arma::cube XXt
  ){

  const int N = ZT.n_slices;  
  const int n = ZT.n_rows;
  const int p = Xr.n_cols;

  arma::vec lb = arma::zeros(p);
  arma::mat Qb = arma::zeros(p,p);
  arma::vec ZrT = arma::zeros(n);
  arma::vec ZcT = arma::zeros(n);
  arma::mat XrT = arma::zeros(n,p);
  arma::mat XcT = arma::zeros(n,p);

  for(int t=0 ; t < N ; ++t){
    arma::mat Z = ZT.slice(t);
    arma::mat mXs = td * mX.slice(t) + to * mXt.slice(t);
    arma::mat XXs = (pow(to,2)+pow(td,2))*XX.slice(t) + 2*to*td*XXt.slice(t);
    arma::mat Zs = td*Z + to*trans(Z);
    arma::vec zr = sum(Zs,1); arma::vec zc=trans(sum(Zs,0));

    if(p>0){
      lb = lb + (trans(mXs) * vectorise(Zs));
      Qb = Qb + XXs + (XX.slice(t)/mXs.n_rows)/ZT.n_slices;
    }

    arma::mat Xsr = td*Xr.slice(t) + to*Xc.slice(t);
    arma::mat Xsc = td*Xc.slice(t) + to*Xr.slice(t);

    ZrT = ZrT+zr;
    ZcT = ZcT+zc;
    XrT = XrT+Xsr;
    XcT = XcT+Xsc;
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