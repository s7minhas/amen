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
#' fit<-ame_repL2(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,model="bin")
#' # you should run the Markov chain much longer than this
#' 
#' @export ame_repL2
ame_repL2<-function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
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


# rm(list=ls())
# library(amen)

# data(YX_bin_long)

# Y=YX_bin_long$Y ; Xdyad=YX_bin_long$X ; Xrow=NULL ;
# Xcol=NULL ; R=2
# burn=1 ; nscan=5 ; odens=1 ; model="bin" ; symmetric=TRUE
# rvar = !(model=="rrl"); cvar = TRUE ; dcor = !symmetric
# nvar=TRUE
# intercept=!is.element(model,c("rrl","ord"))
# intercept=FALSE
# # odmax=rep(max(apply(Y>0,c(1,3),sum,na.rm=TRUE)),nrow(Y[,,1]))
# seed=6886
# plot=FALSE ;  print = FALSE ;  gof=TRUE

# # restructure Y
# Y <- lapply(1:dim(YX_bin_long$Y)[3], function(t){YX_bin_long$Y[,,t]})
# Xdyad <- lapply(1:dim(YX_bin_long$X)[4], function(t){YX_bin_long$X[,,,t]})

# # add labels
# set.seed(6886) ; actors <- as.character( sample(300:700,size=50,replace=FALSE) )
# Y <- lapply(Y, function(y){ rownames(y) <- actors ; colnames(y) <- actors ; return(y) })
# varLabs = paste0('dyadVar',1:3)
# Xdyad <- lapply(Xdyad, function(x){ rownames(x) <- actors ; colnames(x) <- actors ; dimnames(x)[[3]] <- varLabs ; return(x) })

# YX_bin_list <- list(Y=Y, X=Xdyad)
    
# set.seed(6886) ; toRem = sample(1:length(actors), 5)
# toRem = actors[toRem]
# Y[[1]] = Y[[1]][-which(actors %in% toRem),-which(actors %in% toRem)]
# Y[[2]] = Y[[2]][-which(actors %in% toRem),-which(actors %in% toRem)]
# Xdyad[[1]] = Xdyad[[1]][-which(actors %in% toRem),-which(actors %in% toRem),]
# Xdyad[[2]] = Xdyad[[2]][-which(actors %in% toRem),-which(actors %in% toRem),]

# load buthe milner data
# load('~/Dropbox/Research/netsMatter/replications/mansfield_milner_2012/inputData/amenData.rda')
# Y=yList ; Xdyad = xDyadList ; Xrow = xNodeList ; seed = 6886
  # set random seed 
  set.seed(seed)

  # get actor info
  actorByYr <- lapply(Y, rownames)
  actorSet <- sort(unique(unlist( actorByYr ))) ; nActors <- length(actorSet)

  # reset odmax param
  odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), nActors )

  # make sure inputted data is formatted correctly
  if( !is.list(Y) ){ stop('Y needs to be inputted as a list.') }
  if(!is.null(Xdyad)){ if( !is.list(Xdyad) ){ stop('Xdyad needs to be inputted as a list.') } }
  if(!is.null(Xrow)){ if( !is.list(Xrow) ){ stop('Xrow needs to be inputted as a list.') } }
  if(!is.null(Xcol)){ if( !is.list(Xcol) ){ stop('Xcol needs to be inputted as a list.') } }

  # set diag to NA 
  N<-length(Y) ; pdLabs <- names(Y) ; Y<-lapply(Y, function(y){diag(y)=NA; return(y)})

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

  # convert into array format for easier processing
  tmp <- array(NA, dim=c(nActors,nActors,N),
    dimnames=list( actorSet, actorSet, names(Y) ) )
  for(t in 1:N){ tmp[rownames(Y[[t]]),rownames(Y[[t]]),t] <- Y[[t]] }
  Y <- tmp ; rm(tmp)

  if(!is.null(Xdyad)){
    tmp <- array(NA, dim=c(nActors,nActors,dim(Xdyad[[1]])[3],N),
      dimnames=list( actorSet, actorSet, dimnames(Xdyad[[1]])[[3]], pdLabs ) )
    for(t in 1:N){
      tmp[rownames(Xdyad[[t]]),rownames(Xdyad[[t]]),,t] <- Xdyad[[t]] }
    Xdyad <- tmp ; rm(tmp) }

  if(!is.null(Xrow)){
    tmp <- array(NA, dim=c(nActors, dim(Xrow[[1]])[2], N),
      dimnames=list( actorSet, colnames(Xrow[[1]]), pdLabs) )
    for(t in 1:N){
      tmp[rownames(Xrow[[t]]),,t] <- Xrow[[t]] }
    Xrow <- tmp ; rm(tmp) }

  if(!is.null(Xcol)){
    tmp <- array(NA, dim=c(nActors, dim(Xcol[[1]])[2], N),
      dimnames=list( actorSet, colnames(Xcol[[1]]), pdLabs) )
    for(t in 1:N){
      tmp[rownames(Xcol[[t]]),,t] <- Xcol[[t]] }
    Xcol <- tmp ; rm(tmp) }

  # force binary if binary model specified 
  if(is.element(model,c("bin","cbin"))) { Y<-1*(Y>0) } 
   
  # observed and max outdegrees 
  if(is.element(model,c("cbin","frn","rrl")))
  {
    odobs<-apply(Y>0,c(1,3),sum,na.rm=TRUE) 
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y[,,1])) } 
  }

  # some settings for symmetric case
  if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }

  # construct design matrix    
  n<-nrow(Y[,,1])
  pr<-length(Xrow[,,1])/n
  pc<-length(Xcol[,,1])/n
  pd<-length(Xdyad[,,,1])/n^2

  X<-array(dim=c(n,n,pr+pc+pd+intercept,N)) 
  for (t in 1:N)
  { 
    Xt<-design_array_listwisedel(Xrow[,,t],Xcol[,,t],Xdyad[,,,t],intercept,n) 

    # re-add intercept if it was removed
    if(dim(Xt)[3]<dim(X)[3])
    {
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

  # design matrix warning for rrl
  if( model=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0)
                   & !any( apply(X,c(3),function(x){var(c(x))})==0) )
  {
    cat("WARNING: row effects are not estimable using this procedure ","\n")
  }

  # design matrix warning for rrl and ord
  if( is.element(model,c("ord","rrl")) & 
      any( apply(X,c(3),function(x){var(c(x))})==0 ) )
  {
    cat("WARNING: an intercept is not estimable using this procedure ","\n")
  }

    
  # construct matrix of ranked nominations for frn, rrl   
  if(is.element(model,c("frn","rrl")))
  {
    ymx<-max(apply(1*(Y>0),c(1,3),sum,na.rm=TRUE))
    YL<-list()
    for (t in 1:N) 
    {
      YL.t<-NULL
      warn<-FALSE
      for(i in 1:nrow(Y[,,1]))
      {
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
  for (t in 1:N)
  {
    if(model=="nrm"){Z[,,t]<-Y[,,t] }
    if(model=="ord"){Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),ncol(Y[,,t]))} 
    if(model=="rrl")
    {  
      Z[,,t]<-matrix(t(apply(Y[,,t],1,zscores)),nrow(Y[,,t]),ncol(Y[,,t])) 
    }  
    if(model=="bin")
    { 
      Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),nrow(Y[,,t])) 
      # zyMax <- max(Z[,,t][Y[,,t]==0],na.rm=TRUE)
      zyMax <- ifelse( sum(Y[,,t]==0, na.rm=TRUE)!=0, max(Z[,,t][Y[,,t]==0],na.rm=TRUE), 0)
      # zyMin <- min(Z[,,t][Y[,,t]==1],na.rm=TRUE)
      zyMin <- ifelse( sum(Y[,,t]==1, na.rm=TRUE)!=0, max(Z[,,t][Y[,,t]==1],na.rm=TRUE), 0)
      z01<-.5*(zyMax+zyMin ) 
      Z[,,t]<-Z[,,t] - z01
    } 
      
    if(is.element(model,c("cbin","frn")))
    {
      Z[,,t]<-Y[,,t]
      for(i in 1:nrow(Y[,,t]))
      {
        yi<-Y[i,,t]
        zi<-zscores(yi)
        rnkd<-which( !is.na(yi) & yi>0 ) 
        if(length(rnkd)>0 && min(zi[rnkd])<0)
        { 
          zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 
        }
          
        if(length(rnkd)<odmax[i]) 
        {
          urnkd<-which( !is.na(yi) & yi==0 ) 
          if(max(zi[urnkd])>0) { zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 }
        }
          
        Z[i,,t]<-zi
      } 
    }
  }
    
  # starting values for missing entries  
  ZA<-Z
  for (t in 1:N)
  { 
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
  BETA <- matrix(nrow = 0, ncol = dim(X)[3] - pr*symmetric)
  VC<-matrix(nrow=0,ncol=5-3*symmetric) 
  UVPS <- U %*% t(V) * 0 
  APS<-BPS<- rep(0,nrow(Y[,,1]))  
  YPS<-array(0,dim=dim(Y)) ; dimnames(YPS)<-dimnames(Y)
  GOF<-matrix(rowMeans(apply(Y,3,gofstats)),1,4)  
  rownames(GOF)<-"obs"
  colnames(GOF)<- c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")
  names(APS)<-names(BPS)<-rownames(U)<-rownames(V)<-rownames(Y[,,1])
   
  # names of parameters, asymmetric case  
  if(!symmetric)
  { 
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve") 
    colnames(BETA) <- dimnames(X)[[3]] 
  }   

  # names of parameters, symmetric case
  if(symmetric)
  { 
    colnames(VC) <- c("va", "ve")  
    rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
    bnames<-dimnames(X)[[3]]
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
   # s=1
    # update Z
    E.nrm<-array(dim=dim(Z))
    for (t in 1:N)
    {
      EZ<-Xbeta(array(X[,,,t],dim(X)[1:3]), beta)+ outer(a, b,"+")+ U%*%t(V)
      if(model=="nrm")
      { 
        Z[,,t]<-rZ_nrm_fc(Z[,,t],EZ,rho,s2,Y[,,t]) ; E.nrm[,,t]<-Z[,,t]-EZ
      }
      if(model=="bin"){ Z[,,t]<-rZ_bin_fc(Z[,,t],EZ,rho,Y[,,t]) }
      if(model=="ord"){ Z[,,t]<-rZ_ord_fc(Z[,,t],EZ,rho,Y[,,t]) }
      if(model=="cbin"){Z[,,t]<-rZ_cbin_fc(Z[,,t],EZ,rho,Y[,,t],odmax,odobs)}
      if(model=="frn")
      { 
        Z[,,t]<-rZ_frn_fc(Z[,,t],EZ,rho,Y[,,t],YL[[t]],odmax,odobs)
      }
      if(model=="rrl"){ Z[,,t]<-rZ_rrl_fc(Z[,,t],EZ,rho,Y[,,t],YL[[t]]) } 
    }

    # update s2
    if (model=="nrm") s2<-rs2_rep_fc(E.nrm,rho) 
      
    # update beta, a b
    tmp <- rbeta_ab_rep_fc(sweep(Z,c(1,2),U%*%t(V)), Sab, rho, X, s2)
    beta <- tmp$beta
    a <- tmp$a * rvar
    b <- tmp$b * cvar 
    if(symmetric){ a<-b<-(a+b)/2 }
 
    # update Sab - full SRM
    if(rvar & cvar & !symmetric )
    {
      Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+nrow(Z[,,1])))
    }

    # update Sab - rvar only
    if (rvar & !cvar & !symmetric) 
    {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(a^2))/2)
    }
     
    # update Sab - cvar only 
    if (!rvar & cvar & !symmetric) 
    {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(b^2))/2)
    }
     
    # update Sab - symmetric case
    if(symmetric & nvar)
    {
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]   
    }
 
    # update rho
    if (dcor) 
    {
      E.T<-array(dim=dim(Z))
      for (t in 1:N)
      {
        E.T[,,t]<-Z[,,t]-(Xbeta(array(X[,,,t],dim(X)[1:3]),beta) + 
                          outer(a, b, "+") + U %*% t(V))
      }
      rho<-rrho_mh_rep(E.T, rho,s2)
    }
     
    # shrink rho - symmetric case 
    if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }



    # update U,V
    if (R > 0)
    {
      E<-array(dim=dim(Z))
      for(t in 1:N){E[,,t]<-Z[,,t]-(Xbeta(array(X[,,,t],dim(X)[1:3]),beta)+
                    outer(a, b, "+"))}
      shrink<- (s>.5*burn)

      if(symmetric)
      { 
        EA<-apply(E,c(1,2),mean) ; EA<-.5*(EA+t(EA))
        UV<-rUV_sym_fc(EA, U, V, s2/dim(E)[3],shrink) 
      }
      if(!symmetric){UV<-rUV_rep_fc(E, U, V,rho, s2,shrink) }

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
      EZ<-Ys<-array(dim=dim(Z))
      for (t in 1:N)
      {
        EZ[,,t]<-Xbeta(array(X[,,,t],dim(X)[1:3]),beta) + 
                  outer(a, b, "+") + U %*% t(V)
        if(symmetric){ EZ[,,t]<-(EZ[,,t]+t(EZ[,,t]))/2 }

        if(model=="bin"){ Ys[,,t]<-simY_bin(EZ[,,t],rho) }
        if(model=="cbin"){ Ys[,,t]<-1*(simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t])>0)}
        if(model=="frn"){ Ys[,,t]<-simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t]) }
        if(model=="rrl"){ Ys[,,t]<-simY_rrl(EZ[,,t],rho,odobs,YO=Y[,,t] ) }
        if(model=="nrm"){ Ys[,,t]<-simY_nrm(EZ[,,t],rho,s2) }
        if(model=="ord"){ Ys[,,t]<-simY_ord(EZ[,,t],rho,Y[,,t]) }
    
        if(symmetric)
        {  
          Yst<-Ys[,,t] ; Yst[lower.tri(Yst)]<-0 ; Ys[,,t]<-Yst+t(Yst)
        }

      } 

      # update posterior sum
      YPS<-YPS+Ys

      # save posterior predictive GOF stats
      if(gof){Ys[is.na(Y)]<-NA ;GOF<-rbind(GOF,rowMeans(apply(Ys,3,gofstats)))}
       
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

      } 


    } 


  } # end MCMC  
    
  # output 

  # posterior means 
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)  
  UVPM<-UVPS/nrow(VC)
  YPM<-YPS/nrow(VC) 
  EZ<-array(dim=dim(Y)) 
  for (t in 1:N)
  {
    EZ[,,t]<-Xbeta(array(X[,,,t],dim(X)[1:3]),apply(BETA,2,mean)) + 
      outer(APM,BPM,"+")+UVPM 
  }

  names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(Y)[[1]]
  dimnames(YPM)<-dimnames(EZ)<-dimnames(Y) 
  rownames(BETA)<-NULL  

  # asymmetric output
  if(!symmetric)
  {
    UDV<-svd(UVPM)
    U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    rownames(U)<-rownames(V)<-dimnames(Y)[[1]]
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
    rownames(U)<-rownames(ULUPM)<-colnames(ULUPM)<-dimnames(Y)[[1]]
    for(t in 1:N)
    { 
      EZ[,,t]<-.5*(EZ[,,t]+t(EZ[,,t]))
      YPM[,,t]<-.5*(YPM[,,t]+t(YPM[,,t]))
    }  
    fit<-list(BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
              YPM=YPM,GOF=GOF)
  } 

# summStats = function(x){ c( mu=mean(x), sd=sd(x), med=median(x), quantile(x, probs=c(0.025,0.975)) ) }
# round(t(apply(fit$BETA, 2, summStats)),2)
# library(reshape2) ; library(ggplot2)
# ugh = melt(fit$BETA)
# ggplot(ugh, aes(x=Var1, y=value)) + geom_line() + facet_wrap(~Var2, scales='free')

  class(fit) <- "ame"
  fit

}