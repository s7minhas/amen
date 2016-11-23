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
#' @author Peter Hoff, Yanjun He
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

  # create full frame of actors
  fullActorSet <- unique(unlist(lapply(Y,rownames)))

  # reset odmax param
  odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), length(fullActorSet) )

  # make sure inputted data is formatted correctly
  if( !is.list(Y) ){ stop('Y needs to be inputted as a list.') }
  # if( is.element(FALSE, unlist(lapply(Y, function(y){c(class(y)=='matrix', nrow(y) == ncol(y))}))) ){
  #   stop('List objects in Y need to be inputted as n x n matrices.') }
  if(!is.null(Xdyad)){ if( !is.list(Xdyad) ){ stop('Xdyad needs to be inputted as a list.') } }
  if(!is.null(Xrow)){ if( !is.list(Xrow) ){ stop('Xrow needs to be inputted as a list.') } }
  if(!is.null(Xcol)){ if( !is.list(Xcol) ){ stop('Xcol needs to be inputted as a list.') } }

  # set diag to NA 
  N<-length(Y) ; Y<-lapply(Y, function(y){diag(y)=NA; return(y)})

  # make sure actor labels are provided
  if( do.call('sum', lapply(Y, function(y) !is.null( rownames(y) ) ) )!=length(Y) | 
    do.call('sum', lapply(Y, function(y) !is.null( colnames(y) ) ) )!=length(Y) ){
    stop('Actor labels need to be provided for all Y list objects.') }

  if(!is.null(Xdyad)){
  if( do.call('sum', lapply(Xdyad, function(x) !is.null( dimnames(x)[[1]] ) ) )!=length(Y) | 
    do.call('sum', lapply(Xdyad, function(x) !is.null( dimnames(x)[[2]] ) ) )!=length(Y) ){
    stop('Actor labels need to be provided for all Xdyad list objects.') } }

  if(!is.null(Xrow)){
    if( do.call('sum', lapply(Xrow, function(x) !is.null(rownames(x)) )) != length(Y) ){
      stop('Actor labels need to be provided for all Xrow list objects.') } }

  if(!is.null(Xcol)){
    if( do.call('sum', lapply(Xcol, function(x) !is.null(rownames(x)) )) != length(Y) ){
      stop('Actor labels need to be provided for all Xcol list objects.') } }

  # make sure actor labels appear in same order across objects
  for(t in 1:N){
    tNames <- rownames(Y[[t]]) ; check <- identical(tNames, colnames(Y[[t]]))
    if(!is.null(Xdyad)){ check <- c(check, identical( tNames, dimnames(Xdyad[[t]])[[1]] )) }
    if(!is.null(Xdyad)){ check <- c(check, identical( tNames, dimnames(Xdyad[[t]])[[2]] )) }  
    if(!is.null(Xrow)){ check <- c(check, identical(tNames, rownames(Xrow[[t]]) )) }
    if(!is.null(Xcol)){ check <- c(check, identical(tNames, rownames(Xcol[[t]]) )) }
    if( sum(check)/length(check) != 1 ){
      stop('Actor labels are not identical across inputted data within time periods.') } } ; rm(list=c('tNames','check','t'))

  # force binary if binary model specified 
  if(is.element(model, c('bin','cbin'))){ Y<-lapply(Y,function(y){1*(y>0)})}

  # observed and max outdegrees 
  if(is.element(model,c("cbin","frn","rrl")))
  {
  odobs<-do.call('cbind',lapply(Y, function(y){ apply(y>0,1,sum,na.rm=TRUE)}))
  if(length(odmax)==1) { odmax<-rep(odmax,length(fullActorSet)) } 
  }

  # some settings for symmetric case
  if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }

  # construct design matrix    
  pr<-ifelse(is.null(Xrow),0,ncol(Xrow[[1]]))
  pc<-ifelse(is.null(Xcol),0,ncol(Xcol[[1]]))
  pd<-ifelse(is.null(Xdyad),0,dim(Xdyad[[1]])[3])

  X<-lapply(1:N, function(t){
    Xt<-design_array(Xrow[[t]], Xcol[[t]], Xdyad[[t]], intercept, nrow(Y[[t]]))
    dimnames(Xt)[[1]] <- dimnames(Xt)[[2]] <- rownames(Xdyad[[t]])

    # re-add intercept if it was removed
    if(dim(Xt)[3] < sum(pr,pc,pd,intercept) )
    {
      tmp <- array(dim=dim(Xt)+c(0,0,1))
      tmp[,,1]<-1 ; tmp[,,-1]<-Xt
      Xt<-tmp
    }
    return(Xt)
  }) ; names(X) <- names(Y)

  # design matrix warnings
  if( is.element( model, c('ord','rrl') ) & sum(pr,pc,pd,intercept)>0 ){

    # construct feature array with missing
    xArr <- array(NA, 
      dim=c(length(fullActorSet),length(fullActorSet),dim(X[[1]])[3],N),
      dimnames=list( fullActorSet, fullActorSet, dimnames(X[[1]])[[3]], names(Y) ) )
    for(t in 1:N){ for(v in 1:dim(xArr)[3]){
      x <- X[[t]][,,v] ; missing <- setdiff(fullActorSet, rownames(x))
      for(m in missing){ x<-cbind(x,NA) ; x<-rbind(x,NA) ; rownames(x)[nrow(x)] <- colnames(x)[ncol(x)] <- m }
      xArr[,,v,t]<-x } } 

    # design matrix warning for rrl
    if( model=='rrl' & any(apply(apply(xArr,c(1,3),var,na.rm=TRUE),2,sum,na.rm=TRUE)==0) 
      & !any( apply(xArr,c(3),function(x){var(c(x),na.rm=TRUE)})==0 ) ){
      cat("WARNING: row effects are not estimable using this procedure ","\n")
    }

    # design matrix warning for rrl and ord
    if( any( apply(xArr,c(3),function(x){var(c(x),na.rm=TRUE)})==0 ) ){
      cat("WARNING: an intercept is not estimable using this procedure ","\n")
    } } ; rm(list=c('t','v','x','missing','xArr'))

  # construct matrix of ranked nominations for frn, rrl   
  if(is.element(model,c("frn","rrl")))
  {
  ymx<-max(do.call('cbind', lapply( Y, function(y){ apply(1*(y>0),1,sum,na.rm=TRUE) } )))
  YL <- list()
  for(t in 1:N){
    YL.t<-NULL
    warn<-FALSE
    for(i in 1:nrow(Y[[t]]))
    {
      yi<-Y[[t]][i,] ; rnkd<-which( !is.na(yi)&yi>0 )
      if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
      yi[rnkd]<-rank(yi[rnkd],ties.method="random")
      Y[[t]][i,]<-yi
      YL.t<-rbind(YL.t, match(1:ymx,yi))
    }
    if(warn){cat("WARNING: Random reordering used to break ties in ranks\n")}
    YL[[t]] <- YL.t
    }
  }

  # starting Z values
  Z<-lapply(Y, function(y){
    if(model=="nrm"){ return(y) }
    if(model=="ord"){ return( matrix(zscores(y),nrow(y),ncol(y)) ) }   
    if(model=='rrl'){ return( matrix(t(apply(y,1,zscores)),nrow(y),ncol(y)) ) }
    if(model=='bin'){
      z<-matrix(zscores(y), nrow(y), ncol(y))
      z01<-.5*(max(z[y==0],na.rm=TRUE) + min(z[y==1],na.rm=TRUE))
      return( z - z01 ) }

    if(is.element(model,c('cbin','frn'))){
      z <- y
      for(i in 1:nrow(y)){
        yi<-y[i,]
        zi<-zscores(yi)
        rnkd<-which( !is.na(yi) & yi>0 )
        if(length(rnkd)>0 && min(zi[rnkd])<0){
          zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3
        }

        if(length(rnkd)<odmax[i]){
          urnkd<-which( !is.na(yi) & yi==0)
          if(max(zi[urnkd]>0)){ zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3}
        }
        z[i,]<-zi
      }
      return(z) }
    })

  # add dimnames to Z
  Z <- lapply(1:N, function(t){ rownames(Z[[t]])=colnames(Z[[t]])=rownames(Y[[t]]) ; return(Z[[t]]) })

  # starting values for missing entries  
  missStart <- lapply(Z, function(z){
    mu<-mean(z,na.rm=TRUE)
    a<-rowMeans(z,na.rm=TRUE) ; b<-colMeans(z,na.rm=TRUE)  
    ab<-matrix(cbind(a,b), nrow=length(a), ncol=2, dimnames=list(names(a), c('a','b')))  
    za<-mu + outer(a,b,'+')
    z[is.na(z)] <- za[is.na(z)]
    return(list(z=z,ab=ab)) })
  Z<-lapply(missStart,function(x)x$z) ; ab<-do.call('rbind',lapply(1:N,function(t)cbind(missStart[[t]]$ab,t)))

  # other starting values
  beta<-rep(0,dim(X[[1]])[3]) 
  s2<-1 
  rho<-0
  ab <- do.call('rbind',lapply(
    split(ab,rownames(ab)),function(x){x=matrix(x,ncol=3);x[order(x[,3],decreasing=TRUE),1:2][1,]}))
  a <- ab[fullActorSet,1] ; b <- ab[fullActorSet,2] ; rm(list='ab')
  Sab<-tcrossprod(c(rvar,cvar))*cov( cbind(a,b) )
  U<-V<-matrix(0, length(fullActorSet), R, dimnames=list(fullActorSet,NULL)) 

  #  output items
  BETA <- matrix(nrow = 0, ncol = dim(X[[1]])[3] - pr*symmetric)
  VC<-matrix(nrow=0,ncol=5-3*symmetric) 
  UVPS <- U %*% t(V) * 0 
  APS<-BPS<- rep(0,length(fullActorSet)) ; names(APS)<-names(BPS)<-fullActorSet
  YPS <- lapply(Y, function(y){ y*0 })
  GOF<-matrix(rowMeans( do.call('cbind', lapply(Y, function(y){ gofstats(y) }) ) ),1,4)  
  rownames(GOF)<-"obs"
  colnames(GOF)<- c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")
  names(APS)<-names(BPS)<-rownames(U)<-rownames(V)<-fullActorSet

  # names of parameters, asymmetric case  
  if(!symmetric)
  { 
  colnames(VC) <- c("va", "cab", "vb", "rho", "ve") 
  colnames(BETA) <- dimnames(X[[1]])[[3]] 
  }   

  # names of parameters, symmetric case
  if(symmetric)
  { 
  colnames(VC) <- c("va", "ve")  
  rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
  bnames<-dimnames(X[[1]])[[3]]
  bni<-bnames[1*intercept] 
  bnn<-gsub("row",bnames[rb],replacement="node") 
  bnd<-bnames[-c(1*intercept,rb,cb)]
  colnames(BETA)<-c(bni,bnn,bnd) 
  }    

  # MCMC
  have_coda<-suppressWarnings(
             try(requireNamespace("coda",quietly = TRUE),silent=TRUE)) 

  for (s in 1:(nscan + burn)) 
  { 
    
    # update Z
    EZ <- lapply(1:N, function(t){
      tActors <- rownames(Y[[t]])
      ez<-Xbeta(X[[t]], beta)+ outer(a[tActors], b[tActors],"+")+ U[tActors,]%*%t(V[tActors,])
      return(ez) })

    Z <- lapply(1:N, function(t){
      if(model=="nrm"){ z<-rZ_nrm_fc(Z[[t]],EZ[[t]],rho,s2,Y[[t]]) }
      if(model=="bin"){ z<-rZ_bin_fc(Z[[t]],EZ[[t]],rho,Y[[t]]) }
      if(model=="ord"){ z<-rZ_ord_fc(Z[[t]],EZ[[t]],rho,Y[[t]]) }
      if(model=="cbin"){ z<-rZ_cbin_fc(Z[[t]],EZ[[t]],rho,Y[[t]],odmax,odobs) }
      if(model=="frn"){ z<-rZ_frn_fc(Z[[t]],EZ[[t]],rho,Y[[t]],YL[[t]],odmax,odobs) }
      if(model=="rrl"){ z<-rZ_rrl_fc(Z[[t]],EZ[[t]],rho,Y[[t]],YL[[t]]) } 
      return(z) })

    # update s2
    if(model=="nrm"){
      E.nrm <- lapply(1:N, function(t){ return( Z[[t]]-EZ[[t]] ) })  
      s2<-rs2_repL_fc(E.nrm,rho) }
      
    # update beta, a b
    tmp <- rbeta_ab_repL_fc(
      lapply(Z, function(z){tActors<-rownames(z); return(z - U[tActors,]%*%t(V[tActors,])) }),
      Sab, rho, X, s2)
    beta <- tmp$beta
    a <- tmp$a * rvar
    b <- tmp$b * cvar 
    if(symmetric){ a<-b<-(a+b)/2 }

    # update Sab - full SRM
    if(rvar & cvar & !symmetric )
    {
      Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+length(fullActorSet) ))
    }

    # update Sab - rvar only
    if (rvar & !cvar & !symmetric) 
    {
      Sab[1, 1] <- 1/rgamma(1, (1 + length(fullActorSet))/2, (1 + sum(a^2))/2)
    }
     
    # update Sab - cvar only 
    if (!rvar & cvar & !symmetric) 
    {
      Sab[2, 2] <- 1/rgamma(1, (1 + length(fullActorSet))/2, (1 + sum(b^2))/2)
    }
     
    # update Sab - symmetric case
    if(symmetric & nvar)
    {
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+length(fullActorSet))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]   
    }

    # update rho
    if (dcor) 
    {
      E.T <- lapply(1:N, function(t){
        tActors <- rownames(Y[[t]])
        ez<-Xbeta(X[[t]], beta)+ outer(a[tActors], b[tActors],"+")+ U[tActors,]%*%t(V[tActors,])
        return( Z[[t]] - ez ) })
      rho<-rrho_mh_repL(E.T, rho,s2)
    }
     
    # shrink rho - symmetric case 
    if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }

    # update U,V
    if (R > 0)
    {
      E <- array(NA, dim=c(length(fullActorSet),length(fullActorSet),N), 
        dimnames=list(fullActorSet,fullActorSet,names(Y)))
      for(t in 1:N){
        tActors <- rownames(Y[[t]])
        E[tActors,tActors,t] <- Z[[t]] - ( Xbeta(X[[t]], beta) + outer(a[tActors], b[tActors],"+") )
      }
      EA<-apply(E,c(1,2),mean,na.rm=TRUE)
      shrink<- (s>.5*burn)

      if(symmetric)
      { 
        EA<-.5*(EA+t(EA))
        UV<-rUV_sym_fc(EA, U, V, s2/dim(E)[3],shrink) 
      }
      if(!symmetric){
        for(t in 1:N){ E[,,t][is.na(E[,,t])] <- EA[is.na(E[,,t])] }
        UV<-rUV_rep_fc(E, U, V,rho, s2,shrink)
      }

      U<-UV$U ; V<-UV$V
    }

    # burn-in countdown
    if(s%%odens==0&s<=burn){cat(round(100*s/burn,2)," pct burnin complete \n")}

    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) 
    { 
      # save BETA and VC - symmetric case 
      if(symmetric) 
      {
        br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2
        sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
        BETA<-rbind(BETA,sbeta)  

        VC<-rbind(VC,c(Sab[1,1],s2) ) 
      }

      # save BETA and VC - asymmetric case 
      if(!symmetric)
      {
        BETA<-rbind(BETA, beta)
        VC<-rbind(VC, c(Sab[upper.tri(Sab, diag = T)], rho,s2)) 
      }

      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b 
        
      # simulate from posterior predictive 
      EZ <- lapply(1:N, function(t){
        tActors <- rownames(Y[[t]])
        ez<-Xbeta(X[[t]], beta)+ outer(a[tActors], b[tActors],"+")+ U[tActors,]%*%t(V[tActors,])
        return(ez) })
      if(symmetric){ EZ <- lapply(EZ, function(ez){ (ez+t(ez))/2 }) }

      Ys <- lapply(1:N, function(t){
        if(model=="bin"){ ys<-simY_bin(EZ[[t]],rho) }
        if(model=="cbin"){ ys<-1*(simY_frn(EZ[[t]],rho,odmax,YO=Y[[t]])>0)}
        if(model=="frn"){ ys<-simY_frn(EZ[[t]],rho,odmax,YO=Y[[t]]) }
        if(model=="rrl"){ ys<-simY_rrl(EZ[[t]],rho,odobs,YO=Y[[t]] ) }
        if(model=="nrm"){ ys<-simY_nrm(EZ[[t]],rho,s2) }
        if(model=="ord"){ ys<-simY_ord(EZ[[t]],rho,Y[[t]]) }
        if(symmetric){ Yst<-ys ; Yst[lower.tri(Yst)]<-0 ; ys<-Yst+t(Yst) }
        return(ys) })

      # update posterior sum
      YPS<-lapply(1:N, function(t){ YPS[[t]]+Ys[[t]] } )

      # save posterior predictive GOF stats
      if(gof){
        Ys <- lapply(1:N, function(t){ ys<-Ys[[t]] ; ys[is.na(Y[[t]])]<-NA ; return(ys) })
        GOF <- rbind(GOF, rowMeans(do.call('cbind', lapply(Ys, gofstats))))
      }
       
      # print MC progress 
      if(print) 
      {
        cat(s,round(apply(BETA,2,mean),2),":",round(apply(VC,2,mean),2),"\n")
        if (have_coda & nrow(VC) > 3 & length(beta)>0) 
        {
          cat(round(coda::effectiveSize(BETA)), "\n")
        }
      }

      # plot MC progress
      if(plot) 
      {
        # plot VC
        par(mfrow=c(1+2*gof,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
        mVC <- apply(VC, 2, median)
        matplot(VC, type = "l", lty = 1)
        abline(h = mVC, col = 1:length(mVC)) 
       
        # plot BETA
        if(length(beta)>0) 
        {
          mBETA <- apply(BETA, 2, median)
          matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
          abline(h = mBETA, col = 1:length(mBETA))
          abline(h = 0, col = "gray") 
        } 
        
        # plot GOF
        if(gof)
        {
          for(k in 1:4)
          {
            hist(GOF[-1,k],xlim=range(GOF[,k]),main="",prob=TRUE,
                 xlab=colnames(GOF)[k],col="lightblue",ylab="",yaxt="n")  
            abline(v=GOF[1,k],col="red") 
          }
        } 
      } # end plot
    } # end param save
  } # end MCMC  

  # organize output into fit object

  # posterior means 
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)  
  UVPM<-UVPS/nrow(VC)
  YPM<-lapply(YPS, function(yps){yps/nrow(VC)})
  EZ <- lapply(1:N, function(t){
    tActors <- rownames(Y[[t]])
    ez<-Xbeta(X[[t]], apply(BETA, 2, mean))+ outer(APM[tActors], BPM[tActors],"+")+ UVPM[tActors,tActors]
    return(ez) })
  rownames(BETA)<-NULL  

  # asymmetric output
  if(!symmetric)
  {
  UDV<-svd(UVPM)
  U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
  V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
  rownames(U)<-rownames(V)<-rownames(UVPM)
  fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
              YPM=YPM,GOF=GOF) 
  }

  # symmetric output
  if(symmetric) 
  {
  ULUPM<-UVPM 
  eULU<-eigen(ULUPM) 
  eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
  U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
  L<-eULU$val[eR]   
  rownames(U)<-rownames(ULUPM)
  for(t in 1:N)
  { 
    EZ[[t]]<-.5*(EZ[[t]]+t(EZ[[t]]))
    YPM[[t]]<-.5*(YPM[[t]]+t(YPM[[t]]))
  }  
  fit<-list(BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
            YPM=YPM,GOF=GOF)
  } 

  class(fit) <- "ame"
  return( fit )

} # end function