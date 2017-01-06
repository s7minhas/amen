#' AME model fitting routine for replicated relational data
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to replicated relational data of
#' various types. 
#' 
#' This command provides posterior inference for parameters in AME models of
#' independent replicated relational data, assuming one of six possible data
#' types/models:
#' 
#' "nrm": A normal AME model.
#' 
#' "bin": A binary probit AME model.
#' 
#' "ord": An ordinal probit AME model. An intercept is not identifiable in this
#' model.
#' 
#' "cbin": An AME model for censored binary data.  The value of 'odmax'
#' specifies the maximum number of links each row may have.
#' 
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#' 
#' "rrl": An AME model based on the row ranks. This is appropriate if the
#' relationships across rows are not directly comparable in terms of scale. An
#' intercept, row random effects and row regression effects are not estimable
#' for this model.
#' 
#' @usage ame_repL(Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(model=="rrl")
#' , cvar = TRUE, dcor = !symmetric, nvar=TRUE,  R = 0, model="nrm",
#' intercept=!is.element(model,c("rrl","ord")),
#' symmetric=FALSE,
#' odmax=NULL, seed = 1,
#' nscan = 10000, burn = 500, odens = 25, plot=TRUE, print = TRUE, gof=TRUE)
#' @param Y a T length list of n x n relational matrices, where T corresponds to the number of replicates (over time, for example). See
#' model below for various data types.
#' @param Xdyad a T length list of n x n x pd arrays of covariates
#' @param Xrow a T length list of n x pr matrices of nodal row covariates
#' @param Xcol a T length list of n x pc matrices of nodal column covariates
#' @param rvar logical: fit row random effects (asymmetric case)?
#' @param cvar logical: fit column random effects (asymmetric case)? 
#' @param dcor logical: fit a dyadic correlation (asymmetric case)?
#' @param nvar logical: fit nodal random effects (symmetric case)? 
#' @param R integer: dimension of the multiplicative effects (can be zero)
#' @param model character: one of "nrm","bin","ord","cbin","frn","rrl" - see
#' the details below
#' @param intercept logical: fit model with an intercept?
#' @param symmetric logical: Is the sociomatrix symmetric by design?
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' models
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param plot logical: plot results while running?
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @return \item{BETA}{posterior samples of regression coefficients}
#' \item{VC}{posterior samples of the variance parameters}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u} 
#' \item{V}{posterior mean of multiplicative column effects v (asymmetric case)}
#' \item{UVPM}{posterior mean of UV}
#' \item{ULUPM}{posterior mean of ULU (symmetric case)} 
#' \item{L}{posterior mean of L (symmetric case)} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix} \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' @author Peter Hoff, Yanjun He, Shahryar Minhas
#' @examples
#' 
#' data(YX_bin_list) 
#' fit<-ame_repL(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,model="bin")
#' # you should run the Markov chain much longer than this
#' 
#' @export ame_repL
ame_repL<-function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
           rvar = !(model=="rrl") , cvar = TRUE, dcor = !symmetric, 
           nvar=TRUE, 
           R = 0,
           model="nrm",
           intercept=!is.element(model,c("rrl","ord")), 
           symmetric=FALSE, 
           odmax=NULL,
           seed = 1, nscan = 10000, burn = 500, odens = 25,
           plot=TRUE, print = TRUE, gof=TRUE)
{ 

  # set random seed 
  set.seed(seed)

  # get actor info
  actorByYr <- lapply(Y, rownames)
  actorSet <- sort(unique(unlist( actorByYr ))) ; n <- length(actorSet)

  # reset odmax param
  odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), n )

  # check formatting of input objects
  checkFormat(Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol)

  # set diag to NA 
  N<-length(Y) ; pdLabs <- names(Y) ; Y<-lapply(Y, function(y){diag(y)=NA; return(y)})

  # convert into large array format
  tmp <- array(NA, dim=c(n,n,N),
    dimnames=list( actorSet, actorSet, names(Y) ) )
  for(t in 1:N){ tmp[rownames(Y[[t]]),rownames(Y[[t]]),t] <- Y[[t]] }
  Y <- tmp ; rm(tmp)

  if(!is.null(Xdyad)){
    tmp <- array(NA, dim=c(n,n,dim(Xdyad[[1]])[3],N),
      dimnames=list( actorSet, actorSet, dimnames(Xdyad[[1]])[[3]], pdLabs ) )
    for(t in 1:N){
      tmp[rownames(Xdyad[[t]]),rownames(Xdyad[[t]]),,t] <- Xdyad[[t]] }
    Xdyad <- tmp ; rm(tmp) }

  if(!is.null(Xrow)){
    tmp <- array(NA, dim=c(n, dim(Xrow[[1]])[2], N),
      dimnames=list( actorSet, colnames(Xrow[[1]]), pdLabs) )
    for(t in 1:N){
      tmp[rownames(Xrow[[t]]),,t] <- Xrow[[t]] }
    Xrow <- tmp ; rm(tmp) }

  if(!is.null(Xcol)){
    tmp <- array(NA, dim=c(n, dim(Xcol[[1]])[2], N),
      dimnames=list( actorSet, colnames(Xcol[[1]]), pdLabs) )
    for(t in 1:N){
      tmp[rownames(Xcol[[t]]),,t] <- Xcol[[t]] }
    Xcol <- tmp ; rm(tmp) }

  # force binary if binary model specified 
  if(is.element(model,c("bin","cbin"))) { Y<-1*(Y>0) } 
   
  # observed and max outdegrees 
  if(is.element(model,c("cbin","frn","rrl")) ){
    odobs<-apply(Y>0,c(1,3),sum,na.rm=TRUE) 
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y[,,1])) } 
  }

  # some settings for symmetric case
  if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }

  # construct design matrix    
  pr<-length(Xrow[,,1])/n
  pc<-length(Xcol[,,1])/n
  pd<-length(Xdyad[,,,1])/n^2

  X<-array(dim=c(n,n,pr+pc+pd+intercept,N)) 
  for (t in 1:N ){ 
    Xt<-design_array_listwisedel(Xrow[,,t],Xcol[,,t],Xdyad[,,,t],intercept,n) 

    # re-add intercept if it was removed
    if(dim(Xt)[3]<dim(X)[3] ){
      tmp<-array( dim=dim(Xt)+c(0,0,1),
        dimnames=list(NULL,NULL,c('intercept',dimnames(Xt)[[3]])) ) 
      tmp[,,1]<-1 ; tmp[,,-1]<-Xt   
      Xt<-tmp
    }

    X[,,,t]<-Xt
  } 
  if( (pr+pc+pd+intercept)>0 ){ dimnames(X)[[3]]<-dimnames(Xt)[[3]] }
  dimnames(X)[[4]]<-dimnames(Y)[[3]]

  # add NAs from design array to Y
  if( (pr + pc + pd) > 0){
    for(t in 1:N){
      for(p in 1:dim(X)[3]){
        # find NAs
        naMat <- X[,,p,t] ; naMat[!is.na(naMat)] <- 1 ; naMat[is.nan(naMat)] <- NA
        # add NAs to Y
        Y[,,t] <- Y[,,t] * naMat
      }
    } ; rm(naMat)
  }

  # turn NAs in design array to zero
  X[is.na(X)]<-0  

  # useful design array transformations to do up front
  ## these are used in calculation of rbeta
  Xlist <- lapply(1:N, function(t){ array(X[,,,t], dim=dim(X)[1:3], dimnames=dimnames(X)[1:3]) })
  XrLong <- apply(X, c(1,3,4), sum)                 # row sum
  XcLong <- apply(X, c(2,3,4), sum)                 # col sum
  mXLong <- apply(X, c(3,4), c)                     # design matrix
  mXtLong <- apply(aperm(X, c(2,1,3,4)), c(3,4), c) # dyad-transposed design matrix
  # regression sums of squares
  xxLong <- array(apply(mXLong, 3, function(x){ t(x)%*%x }), dim=c(dim(X)[3],dim(X)[3],N))
  xxTLong <- array(unlist(lapply(1:N, function(t){ t(mXLong[,,t]) %*% mXtLong[,,t] })), dim=dim(xxLong))
  
  # design matrix warning for rrl
  if(model=='rrl'){
    if( any(apply(apply(X,c(1,3),var),2,sum)==0) & !any( apply(X,c(3),function(x){var(c(x))})==0) ){
      cat("WARNING: row effects are not estimable using this procedure ","\n")  
    }
  }

  # design matrix warning for rrl and ord
  if( is.element(model,c("ord","rrl")) ){
    if( any( apply(X,c(3),function(x){var(c(x))})==0 ) ){
      cat("WARNING: an intercept is not estimable using this procedure ","\n")
    }
  }
    
  # construct matrix of ranked nominations for frn, rrl   
  if(is.element(model,c("frn","rrl")) ){
    ymx<-max(apply(1*(Y>0),c(1,3),sum,na.rm=TRUE))
    YL<-list()
    for (t in 1:N ){
      YL.t<-NULL
      warn<-FALSE
      for(i in 1:nrow(Y[,,1]) ){
        yi<-Y[i,,t] ; rnkd<-which( !is.na(yi)&yi>0 )
        if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
        yi[rnkd]<-rank(yi[rnkd],ties.method="random")
        Y[i,,t]<-yi
        YL.t<-rbind(YL.t, match(1:ymx,yi))
      }
      YL[[t]]<-YL.t
      if(warn){cat("WARNING: Random reordering used to break ties in ranks\n")}
    }
  }
    
  # starting Z values
  Z<-array(dim=dim(Y))
  for (t in 1:N ){
    if(model=="nrm"){Z[,,t]<-Y[,,t] }
    if(model=="ord"){Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),ncol(Y[,,t]))} 
    if(model=="rrl" ){  
      Z[,,t]<-matrix(t(apply(Y[,,t],1,zscores)),nrow(Y[,,t]),ncol(Y[,,t])) 
    }  
    if(model=="bin" ){ 
      Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),nrow(Y[,,t])) 
      # zyMax <- max(Z[,,t][Y[,,t]==0],na.rm=TRUE)
      zyMax <- ifelse( sum(Y[,,t]==0, na.rm=TRUE)!=0, max(Z[,,t][Y[,,t]==0],na.rm=TRUE), 0)
      # zyMin <- min(Z[,,t][Y[,,t]==1],na.rm=TRUE)
      zyMin <- ifelse( sum(Y[,,t]==1, na.rm=TRUE)!=0, max(Z[,,t][Y[,,t]==1],na.rm=TRUE), 0)
      z01<-.5*(zyMax+zyMin ) 
      Z[,,t]<-Z[,,t] - z01
    } 
      
    if(is.element(model,c("cbin","frn")) ){
      Z[,,t]<-Y[,,t]
      for(i in 1:nrow(Y[,,t]) ){
        yi<-Y[i,,t]
        zi<-zscores(yi)
        rnkd<-which( !is.na(yi) & yi>0 ) 
        if(length(rnkd)>0 && min(zi[rnkd])<0 ){ 
          zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 
        }
          
        if(length(rnkd)<odmax[i] ){
          urnkd<-which( !is.na(yi) & yi==0 ) 
          if(max(zi[urnkd])>0) { zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 }
        }
          
        Z[i,,t]<-zi
      } 
    }
  }

  # starting values for missing entries  
  ZA<-Z
  for (t in 1:N){ 
    mu<-mean(Z[,,t],na.rm=TRUE) 
    a<-rowMeans(Z[,,t],na.rm=TRUE) ; b<-colMeans(Z[,,t],na.rm=TRUE)
    # a[is.na(a)] <- mean(a, na.rm=TRUE) ; b[is.na(b)] <- mean(b, na.rm=TRUE)
    a[is.na(a)] <- 0 ; b[is.na(b)] <- 0
    ZA[,,t]<-mu + outer(a,b,"+")
  }
  Z[is.na(Z)]<-ZA[is.na(Z)]
    
  # other starting values
  beta<-rep(0,dim(X)[3]) 
  s2<-1 
  rho<-0
  Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
  U<-V<-matrix(0, nrow(Y[,,1]), R) 
    
  #  output items
  BETA <- matrix(nrow = nscan/odens, ncol = dim(X)[3] - pr*symmetric)
  VC<-matrix(nrow=nscan/odens,ncol=5-3*symmetric)   
  UVPS <- U %*% t(V) * 0 
  APS<-BPS<- rep(0,nrow(Y[,,1]))  
  YPS<-array(0,dim=dim(Y)) ; dimnames(YPS)<-dimnames(Y)
  GOF <- matrix(NA, nrow=(nscan/odens)+1, ncol=4,
    dimnames=list(c('obs',1:(nscan/odens)),c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")))
  GOF[1,] <- rowMeans(apply(Y,3,gofstats))
  names(APS)<-names(BPS)<-rownames(U)<-rownames(V)<-rownames(Y[,,1])
   
  # names of parameters, asymmetric case  
  if(!symmetric ){ 
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve") 
    colnames(BETA) <- dimnames(X)[[3]] 
  }   

  # names of parameters, symmetric case
  if(symmetric ){ 
    colnames(VC) <- c("va", "ve")  
    rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
    bnames<-dimnames(X)[[3]]
    bni<-bnames[1*intercept] 
    bnn<-gsub("row",bnames[rb],replacement="node") 
    bnd<-bnames[-c(1*intercept,rb,cb)]
    colnames(BETA)<-c(bni,bnn,bnd) 
  }    

  # MCMC
  have_coda<-suppressWarnings(try(requireNamespace("coda",quietly = TRUE),silent=TRUE)) 

  for (s in 1:(nscan + burn) ){ 

    # update Z
    E.nrm<-array(dim=dim(Z))
    EZt <- get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
    for(t in 1:N ){
      EZ<-EZt[,,t]
      if(model=="nrm" ){ 
        Z[,,t]<-rZ_nrm_fc(Z[,,t],EZ,rho,s2,Y[,,t]) ; E.nrm[,,t]<-Z[,,t]-EZ
      }
      if(model=="bin"){ Z[,,t]<-rZ_bin_fc(Z[,,t],EZ,rho,Y[,,t]) }
      if(model=="ord"){ Z[,,t]<-rZ_ord_fc(Z[,,t],EZ,rho,Y[,,t]) }
      if(model=="cbin"){Z[,,t]<-rZ_cbin_fc(Z[,,t],EZ,rho,Y[,,t],odmax,odobs)}
      if(model=="frn" ){ 
        Z[,,t]<-rZ_frn_fc(Z[,,t],EZ,rho,Y[,,t],YL[[t]],odmax,odobs)
      }
      if(model=="rrl"){ Z[,,t]<-rZ_rrl_fc(Z[,,t],EZ,rho,Y[,,t],YL[[t]]) } 
    }

    # update s2
    if (model=="nrm"){ s2<-rs2_rep_fc_cpp(E.nrm,solve(matrix(c(1,rho,rho,1),2,2))) }
      
    # update beta, a b
    iSe2<-mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)) ; Sabs<-iSe2%*%Sab%*%iSe2
    tmpz<-eigen(Sabs) ; k<-sum(zapsmall(tmpz$val)>0 ) ; G<-tmpz$vec[,1:k] %*% sqrt(diag(tmpz$val[1:k],nrow=k))
    tmp <- rbeta_ab_rep_fc_cpp( 
      zCube=sweep(Z,c(1,2),U%*%t(V)), XrCube=XrLong, XcCube=XcLong, 
      mXCube=mXLong, mXtCube=mXtLong, xxCube=xxLong, xxTCube=xxTLong,
      iSe2=iSe2, Sabs=Sabs, k=k, G=G )
    beta <- tmp$beta
    a <- tmp$a * rvar
    b <- tmp$b * cvar 
    if(symmetric){ a<-b<-(a+b)/2 }
 
    # update Sab - full SRM
    if(rvar & cvar & !symmetric ){
      Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+nrow(Z[,,1])))
    }

    # update Sab - rvar only
    if (rvar & !cvar & !symmetric ){
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(a^2))/2)
    }
     
    # update Sab - cvar only 
    if (!rvar & cvar & !symmetric ){
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(b^2))/2)
    }
     
    # update Sab - symmetric case
    if(symmetric & nvar ){
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]   
    }

    # update rho
    if( dcor ){
      E.T <- Z - get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
      rho<-rrho_mh_rep_cpp(E.T, rho,s2)
    }
     
    # shrink rho - symmetric case 
    if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }

    # update U,V
    if (R > 0){
      E <- Z-get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U*0, V*0 )
      shrink <- (s>.5*burn)

      if(symmetric ){ 
        EA<-apply(E,c(1,2),mean) ; EA<-.5*(EA+t(EA))
        UV<-rUV_sym_fc_cpp(EA, U, V, s2/dim(E)[3], shrink, rep(sample(1:nrow(E)),4)-1)
      }
      if(!symmetric){
        UV <- rUV_rep_fc_cpp(E2, U, V, rho, s2, 
          mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)),
          maxmargin=1e-6, shrink, sample(1:R)-1 )
      }

      U<-UV$U ; V<-UV$V
    }

    # burn-in countdown
    if(s%%odens==0&s<=burn){cat(round(100*s/burn,2)," pct burnin complete \n")}
  
    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn ){ 

      # save BETA and VC - symmetric case 
      if(symmetric ){
        br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2
        sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
        BETA[(s-burn)/odens,]<-beta
        VC[(s-burn)/odens,]<-c(Sab[1,1],s2)
      }
    
      # save BETA and VC - asymmetric case 
      if(!symmetric ){
        BETA[(s-burn)/odens,]<-beta
        VC[(s-burn)/odens,]<- c(Sab[upper.tri(Sab, diag = T)], rho,s2)        
      }

      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b 
        
      # simulate from posterior predictive 
      Ys<-array(dim=dim(Z))
      EZ <- get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
      for(t in 1:N ){
        EZt<-EZ[,,t]
        if(symmetric){ EZt<-(EZt+t(EZt))/2 }
        if(model=="bin"){ Ys[,,t]<-simY_bin(EZt,rho) }
        if(model=="cbin"){ Ys[,,t]<-1*(simY_frn(EZt,rho,odmax,YO=Y[,,t])>0)}
        if(model=="frn"){ Ys[,,t]<-simY_frn(EZt,rho,odmax,YO=Y[,,t]) }
        if(model=="rrl"){ Ys[,,t]<-simY_rrl(EZt,rho,odobs,YO=Y[,,t] ) }
        if(model=="nrm"){ Ys[,,t]<-simY_nrm(EZt,rho,s2) }
        if(model=="ord"){ Ys[,,t]<-simY_ord(EZt,rho,Y[,,t]) }
    
        if(symmetric ){  
          Yst<-Ys[,,t] ; Yst[lower.tri(Yst)]<-0 ; Ys[,,t]<-Yst+t(Yst)
        }
      } 

      # update posterior sum
      YPS<-YPS+Ys

      # save posterior predictive GOF stats
      if(gof){ Ys[is.na(Y)]<-NA ; GOF[((s-burn)/odens)+1,]<-rowMeans(apply(Ys,3,gofstats)) }
       
      # print MC progress 
      if(print){
        cat(s,round(apply(BETA,2,mean),2),":",round(apply(VC,2,mean),2),"\n")
        if (have_coda & nrow(VC) > 3 & length(beta)>0){
          cat(round(coda::effectiveSize(BETA)), "\n")
        }
      }

      # plot MC progress
      if(plot){
        # plot VC
        par(mfrow=c(1+2*gof,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
        mVC <- apply(VC, 2, median)
        matplot(VC, type = "l", lty = 1)
        abline(h = mVC, col = 1:length(mVC)) 
       
        # plot BETA
        if(length(beta)>0){
          mBETA <- apply(BETA, 2, median)
          matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
          abline(h = mBETA, col = 1:length(mBETA))
          abline(h = 0, col = "gray") 
        } 
        
        # plot GOF
        if(gof){
          for(k in 1:4){
            hist(GOF[-1,k],xlim=range(GOF[,k]),main="",prob=TRUE,
                 xlab=colnames(GOF)[k],col="lightblue",ylab="",yaxt="n")  
            abline(v=GOF[1,k],col="red") 
          }
        } 
      } 

    } # end post burnin save


  } # end MCMC  
    
  # output 

  # posterior means 
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)  
  UVPM<-UVPS/nrow(VC)
  YPM<-YPS/nrow(VC) 
  EZ<-array(dim=dim(Y)) 
  for (t in 1:N){
    EZ[,,t]<-Xbeta_cpp(array(X[,,,t],dim(X)[1:3]),apply(BETA,2,mean)) + 
      outer(APM,BPM,"+")+UVPM 
  }

  names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(Y)[[1]]
  dimnames(YPM)<-dimnames(EZ)<-dimnames(Y) 
  rownames(BETA)<-NULL  

  # asymmetric output
  if(!symmetric){
    UDV<-svd(UVPM)
    U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    rownames(U)<-rownames(V)<-dimnames(Y)[[1]]

    # reformat EZ and YPM as list objects
    EZ <- lapply(1:dim(EZ)[3], function(t){
      actorT <- actorByYr[[t]]
      ez <- EZ[actorT,actorT,t]
      diag(ez) <- NA
      return( ez ) }) ; names(EZ) <- pdLabs

    YPM <- lapply(1:dim(YPM)[3], function(t){
      actorT <- actorByYr[[t]]
      ypm <- YPM[actorT,actorT,t]
      diag(ypm) <- NA
      return( ypm ) }) ; names(YPM) <- pdLabs    

    # create fitted object
    fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
                YPM=YPM,GOF=GOF) 
  }
 
  # symmetric output
  if(symmetric){
    ULUPM<-UVPM 
    eULU<-eigen(ULUPM) 
    eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
    U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
    L<-eULU$val[eR]   
    rownames(U)<-rownames(ULUPM)<-colnames(ULUPM)<-dimnames(Y)[[1]]
    for(t in 1:N){ 
      EZ[,,t]<-.5*(EZ[,,t]+t(EZ[,,t]))
      YPM[,,t]<-.5*(YPM[,,t]+t(YPM[,,t]))
    }

    # reformat EZ and YPM as list objects
    EZ <- lapply(1:dim(EZ)[3], function(t){
      actorT <- actorByYr[[t]]
      ez <- EZ[actorT,actorT,t]
      diag(ez) <- NA
      return( ez ) }) ; names(EZ) <- pdLabs

    YPM <- lapply(1:dim(YPM)[3], function(t){
      actorT <- actorByYr[[t]]
      ypm <- YPM[actorT,actorT,t]
      diag(ypm) <- NA
      return( ypm ) }) ; names(YPM) <- pdLabs    

    # create fitted object
    fit<-list(BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
      YPM=YPM,GOF=GOF)
  } 

  class(fit) <- "ame"
  fit

}