library(amen)
Xcol=NULL ; R=2
burn=1 ; nscan=5 ; odens=1 ; model="bin" ; symmetric=TRUE
rvar = !(model=="rrl"); cvar = TRUE ; dcor = !symmetric
nvar=TRUE
intercept=!is.element(model,c("rrl","ord"))
intercept=FALSE
# odmax=rep(max(apply(Y>0,c(1,3),sum,na.rm=TRUE)),nrow(Y[,,1]))
seed=6886
plot=FALSE ;  print = FALSE ;  gof=TRUE

load('~/Dropbox/Research/netsMatter/replications/mansfield_milner_2012/inputData/amenData.rda')
Y=yList ; Xdyad = xDyadList ; Xrow = xNodeList ; seed = 6886
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
s=1
ameIterOld = function(){
  E.nrm<-array(dim=dim(Z))
  for(t in 1:N ){
    EZ<-Xbeta(array(X[,,,t],dim(X)[1:3]), beta)+ outer(a, b,"+")+ U%*%t(V) # move out of loop
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

  if (model=="nrm"){ s2<-rs2_rep_fc(E.nrm, solve(matrix(c(1,rho,rho,1),2,2)) ) }

  tmp <- rbeta_ab_rep_fc( sweep(Z,c(1,2),U%*%t(V)), Sab, rho, X, s2 )
  beta <- c(tmp$beta)
  a <- c(tmp$a) * rvar
  b <- c(tmp$b) * cvar 
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

  if( dcor ){
          E.T<-array(dim=dim(Z))
          for (t in 1:N ){
            E.T[,,t]<-Z[,,t]-(Xbeta(array(X[,,,t],dim(X)[1:3]),beta) + # move out of loop
                              outer(a, b, "+") + U %*% t(V))
          }
          rho<-rrho_mh_rep(E.T, rho,s2)
  }

  if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }

  if (R > 0){
    E<-array(dim=dim(Z))
    for(t in 1:N){E[,,t]<-Z[,,t]-(Xbeta(array(X[,,,t],dim(X)[1:3]),beta)+ # move out of loop
                  outer(a, b, "+"))}
    shrink <- (s>.5*burn)
    if(symmetric ){ 
      EA<-apply(E,c(1,2),mean) ; EA<-.5*(EA+t(EA))
      UV<-rUV_sym_fc(EA, U, V, s2/dim(E)[3], shrink) }
    if(!symmetric){  UV<-rUV_rep_fc(E2, U, V,rho, s2,shrink) }
    U<-UV$U ; V<-UV$V
  }  
}

ameIterNew = function(){
  E.nrm<-array(dim=dim(Z))
  EZt = get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
  for(t in 1:N ){
    EZ<-EZt[,,t] # move out of loop
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
  if (model=="nrm"){ s2<-rs2_rep_fc_cpp(E.nrm, solve(matrix(c(1,rho,rho,1),2,2)) ) }

  iSe2<-mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)) ; Sabs<-iSe2%*%Sab%*%iSe2
  tmpz<-eigen(Sabs) ; k<-sum(zapsmall(tmpz$val)>0 ) ; G<-tmpz$vec[,1:k] %*% sqrt(diag(tmpz$val[1:k],nrow=k))
  tmp <- rbeta_ab_rep_fc_cpp(
    zCube=sweep(Z,c(1,2),U%*%t(V)), XrCube=XrLong, XcCube=XcLong, 
    mXCube=mXLong, mXtCube=mXtLong, xxCube=xxLong, xxTCube=xxTLong,
    iSe2=iSe2, Sabs=Sabs, k=k, G=G )
  beta <- c(tmp$beta)
  a <- c(tmp$a) * rvar
  b <- c(tmp$b) * cvar 
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
    E.T = Z - get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )  
    rho<-rrho_mh_rep_cpp(E.T, rho,s2)
  }

  # shrink rho - symmetric case 
  if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }

  # update U,V
  if (R > 0){
    E = Z-get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U*0, V*0 )
    shrink <- (s>.5*burn)
    if(symmetric ){ 
      EA<-apply(E,c(1,2),mean) ; EA<-.5*(EA+t(EA))
      UV<-rUV_sym_fc_cpp(EA, U, V, s2/dim(E)[3], shrink, rep(sample(1:nrow(E)),4)-1)
    }
    if(!symmetric){
      UV = rUV_rep_fc_cpp(E2, U, V, rho, s2, 
        mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)), 
        maxmargin=1e-6, shrink, sample(1:R)-1 )
    }
    U<-UV$U ; V<-UV$V
  }
}

library(rbenchmark)
benchmark(ameIterOld(), ameIterNew(), replications = 100)
# > library(rbenchmark)
# > benchmark(ameIterOld(), ameIterNew(), replications = 100)
#           test replications  elapsed relative user.self sys.self user.child  sys.child
# 2 ameIterNew()          100  114.840    1.000    92.169   22.889          0  0
# 1 ameIterOld()          100 1752.906   15.264  1175.685  572.889          0  0