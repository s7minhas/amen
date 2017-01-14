#' Gibbs sampling of additive row and column effects and regression coefficient
#' with independent replicate relational data
#' 
#' Simulates from the joint full conditional distribution of (a,b,beta),
#' assuming same additive row and column effects and regression coefficient
#' across replicates. 
#' 
#' 
#' @usage rbetaCPP_ab_rep_fc(Z.T,Sab,rho,s2=1,...)
#' @param Z.T n x n x T array, with the third dimension for replicates.
#' Each slice of the array is a (latent) normal relational matrix, with
#' multiplicative effects subtracted out
#' @param XrLong n x p x T row covariate array generated within ame_repL fn
#' @param XcLong n x p x T column covariate array generated within ame_repL fn
#' @param mXLong n^2 x p x T design array generated within ame_repL fn
#' @param mXtLong n^2 x p x T transposed design array generated within ame_repL fn
#' @param xxLong p x p x T regression sum of squares array generated within ame_repl fn
#' @param xxTLong p x p x T transposed regression sum of squares array generated 
#' within ame_repl fn
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author Peter Hoff, Yanjun He, Shahryar Minhas
#' @export rbetaCPP_ab_rep_fc
rbetaCPP_ab_rep_fc <-
  function(Z.T,Sab,rho,s2=1,
    XrLong,XcLong,mXLong,mXtLong,xxLong,xxTLong) 
  {
    ###
    N<-dim(Z.T)[3]
    p<-dim(XrLong)[2]
    Se<-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2<-mhalf(solve(Se))
    td<-iSe2[1,1] ; to<-iSe2[1,2]
    Sabs<-iSe2%*%Sab%*%iSe2
    tmp<-eigen(Sabs)
    k<-sum(zapsmall(tmp$val)>0 )
    ###
    
    betaEstimates<-rbeta_rep_cpp(ZT=Z.T, to=to, td=td, Xr=XrLong,
      Xc=XcLong, mX=mXLong, mXt=mXtLong, XX=xxLong, 
      XXt=xxTLong)
    lb = c(betaEstimates$lb)
    Qb = c(betaEstimates$Qb)
    Zr.T = c(betaEstimates$ZrT)
    Zc.T = c(betaEstimates$ZcT)    
    Xr.T = betaEstimates$XrT
    Xc.T = betaEstimates$XcT
    rm(betaEstimates)

    ## dyadic and prior contributions  

    ## row and column reduction
    ab<-matrix(0,nrow(Z.T),2)
    if(k>0)
    {
      n<-nrow(Z.T[,,1])
      G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
      K<-matrix(c(0,1,1,0),2,2)
      
      A<-N*n*t(G)%*%G + diag(k)
      B<-N*t(G)%*%K%*%G
      iA0<-solve(A)
      C0<- -solve(A+ n*B)%*%B%*%iA0
      
      iA<-G%*%iA0%*%t(G)
      C<-G%*%C0%*%t(G)
      
      #BigMatrix<-rbind(cbind(n*diag(n),matrix(1,n,n)),cbind(matrix(1,n,n),n*diag(n)))
      #Gn<-G%x%diag(n)
      #V.inv<-N*crossprod(crossprod(BigMatrix,Gn),Gn)+diag(k*n)
      #V<-solve(V.inv)
      #H<-tcrossprod(tcrossprod(Gn,V),Gn)
      H<-iA%x%diag(n)+C%x%matrix(1,nrow=n,ncol=n)
      Hrr<-H[1:n,1:n]
      Hrc<-H[1:n,(n+1):(2*n)]
      Hcr<-H[(n+1):(2*n),1:n]
      Hcc<-H[(n+1):(2*n),(n+1):(2*n)]
      Qb<-Qb-t(Xr.T)%*%Hrr%*%Xr.T-t(Xc.T)%*%Hcr%*%Xr.T-t(Xr.T)%*%Hrc%*%Xc.T-t(Xc.T)%*%Hcc%*%Xc.T
      lb<-lb-t(Xr.T)%*%Hrr%*%Zr.T-t(Xc.T)%*%Hcr%*%Zr.T-t(Xr.T)%*%Hrc%*%Zc.T-t(Xc.T)%*%Hcc%*%Zc.T
    } 

    ##
    if(p>0) 
    {
      V.b<-solve(Qb)
      M.b<-V.b%*%lb
      beta<-c(rmvnorm(1,M.b,V.b))
    }
    if(p==0){ beta<-numeric(0) } 

    
    Rr.T<-Zr.T-Xr.T%*%matrix(beta,ncol=1)
    Rc.T<-Zc.T-Xc.T%*%matrix(beta,ncol=1)
    #M.ab<-V%*%(t(G)%x%diag(n))%*%matrix(c(Rr.T,Rc.T),ncol=1)
    #W<-rmvnorm(1,mu=M.ab,V)
    #ab.vec.t<-(G%x%diag(n))%*%matrix(W,ncol=1)
    #ab.vec<-solve(iSe2)%*%matrix(ab.vec.t,nrow=2,byrow=T)
    
    m<-t(t(crossprod(rbind(c(Rr.T),c(Rc.T)),t(iA0%*%t(G)))) + rowSums(sum(Rr.T)*C0%*%t(G)) )
    hiA0<-mhalf(iA0)
    e<-matrix(rnorm(n*k),n,k) 
    w<-m+ t( t(e%*%hiA0) - c(((hiA0-mhalf(iA0+n*C0))/n)%*% colSums(e) ) )
    ab.vec<- t(w%*%t(G)%*%solve(iSe2))
    
    return( list(beta=beta,a=ab.vec[1,],b=ab.vec[2,] ) )
  }
