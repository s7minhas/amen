//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Gibbs sampling of additive row and column effects and regression coefficient
//' with independent replicate relational data
//' 
//' Simulates from the joint full conditional distribution of (a,b,beta),
//' assuming same additive row and column effects and regression coefficient
//' across replicates. 
//' 
//' 
//' @usage rbeta_ab_rep_cpp(Z.T,Sab,rho,X.T,s2=1)
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
//' @param Sabs row and column covariance
//' @param k dimensions for row and column random effects
//' @param G eigenvalue calcs from Sab
//' @param e n x k gaussian error matrix
//' @param colE column sums of e
//' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
//' \item{b}{additive column effects}
//' @author Peter Hoff, Yanjun He, Shahryar Minhas
//' @export rbeta_ab_rep_fc_cpp
// [[Rcpp::export]]

List rbeta_ab_rep_fc_cpp(
	arma::cube zCube, arma::cube XrCube, arma::cube XcCube, 
	arma::cube mXCube, arma::cube mXtCube, 
	arma::cube xxCube, arma::cube xxTCube,
	arma::mat iSe2, arma::mat Sabs, int k, arma::mat G, 
	arma::mat e, arma::vec colE
	) {

  double td = iSe2(0,0);
  double to = iSe2(0,1);	

  int N = zCube.n_slices;  
  int n = zCube.n_rows;
  int p = XrCube.n_cols;

  arma::vec lb = arma::zeros(p);
  arma::mat Qb = arma::zeros(p, p);
  arma::vec ZrT = arma::zeros(n);
  arma::vec ZcT = arma::zeros(n);
  arma::mat XrT = arma::zeros(n,p);
  arma::mat XcT = arma::zeros(n,p);  

  for(int t=0 ; t < N ; ++t){
    arma::mat mXs = td * mXCube.slice(t) + to*mXtCube.slice(t);
    arma::mat XXs = (pow(to, 2)+pow(td, 2))*xxCube.slice(t) + 2*to*td*xxTCube.slice(t);
    arma::mat Zs = td*zCube.slice(t) + to*zCube.slice(t).t();

    arma::vec zr(n);
    for (int i=0 ; i < n ; ++i){
      zr[i] = arma::sum(zCube.slice(t).row(i));
    }

    arma::vec zc(n);
    for (int j=0 ; j < n ; ++j){
      zc[j] = arma::sum(zCube.slice(t).col(j));
    }

    if ( p > 0 ) {
      lb = lb + mXs.t() * vectorise(Zs);
      Qb = Qb + XXs + ( xxCube.slice(t)/mXs.n_rows )/N;
    }

    arma::mat Xsr = td * XrCube.slice(t) + to*XcCube.slice(t);
    arma::mat Xsc = td * XcCube.slice(t) + to*XrCube.slice(t);

    ZrT = ZrT + zr;
    ZcT = ZcT + zc;
    XrT = XrT + Xsr;
    XcT = XcT + Xsc;
  }

  arma::mat ab = arma::zeros(n,2);

  // if ( k > 0 ){
  	
  	arma::mat K = arma::mat("0 1; 1 0");
  	arma::mat idmatk = eye<mat>(k,k);

  	arma::mat A = N*n*G.t()*G + idmatk;
  	arma::mat B = N*G.t()*K*G;
  	arma::mat iA0 = inv(A);
  	arma::mat C0 = -inv(A + n*B)*B*iA0;

  	arma::mat iA = G * iA0 * G.t();
  	arma::mat C = G * C0 * G.t();

	arma::mat idmatn = eye<mat>(n,n);
	arma::mat onesmatn = ones<mat>(n,n);
	arma::mat H = arma::kron(iA, idmatn) + arma::kron(C,onesmatn);
	arma::mat Hrr = H.submat(0,0,n-1,n-1);
	arma::mat Hrc = H.submat(0,n,n-1,2*n-1);
	arma::mat Hcr = H.submat(n,0,2*n-1,n-1);
	arma::mat Hcc = H.submat(n,n,2*n-1,2*n-1);
	Qb = Qb-XrT.t()*Hrr*XrT-XcT.t()*Hcr*XrT-XrT.t()*Hrc*XcT-XcT.t()*Hcc*XcT;
	lb = lb-XrT.t()*Hrr*ZrT-XcT.t()*Hcr*ZrT-XrT.t()*Hrc*ZcT-XcT.t()*Hcc*ZcT;
	// }

  // if ( p > 0 ) {

  	arma::mat Vb = inv(Qb);
  	arma::mat Mb = Vb * lb;
  	arma::vec beta = as<arma::vec>(amen::rmvnorm(1,Mb,Vb));
  // }

  // if ( p == 0 ) {

  // 	arma::vec beta = arma::zeros(1);
  // }

  arma::mat RrT = ZrT - XrT * beta; 
  arma::mat RcT = ZcT - XcT * beta; 

  arma::mat RTcrossiA0G = join_rows(RrT, RcT) * (iA0 * G.t()).t();
  arma:vec RrTC0G = sum(accu(RrT) * C0 * G.t(),1);
  arma::mat m = arma::zeros(n,RTcrossiA0G.n_cols);
  for( int r=0 ; r < RTcrossiA0G.n_cols ; r++ ) {
    m.col(r) = RTcrossiA0G.col(r) + RrTC0G[r];
  }

  arma::mat hiA0 = as<arma::mat>(amen::mhalf(iA0));
  arma::mat ehiA0 = (e * hiA0).t();
  arma::mat iA0nCo = as<arma::mat>(amen::mhalf(iA0+n*C0));
  arma::mat hiA0nCo = (hiA0-iA0nCo)/n;
  arma::vec ugh = hiA0nCo * colE;
  arma::mat w = arma::zeros(n,RTcrossiA0G.n_cols);
  for( int r=0 ; r < RTcrossiA0G.n_cols ; r++ ) {
    w.col(r) = m.col(r) + (ehiA0.row(r) - ugh[r]).t();
  }

  arma::mat abVec = w * G.t() * inv(iSe2);
  arma::vec a = abVec.col(0);
  arma::vec b = abVec.col(1);

  return(Rcpp::List::create(beta, a, b));

}