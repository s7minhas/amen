rm(list=ls())
library(devtools)
devtools::install_github('s7minhas/amen')
library(amen)

data(YX_bin_long)

Y=YX_bin_long$Y ; Xdyad=YX_bin_long$X ; Xrow=NULL ; Xcol=NULL ; R=2
burn=5 ; nscan=5 ; odens=1 ; model="bin" ; symmetric=FALSE
rvar = !(model=="rrl"); cvar = TRUE ; dcor = !symmetric
nvar=TRUE
intercept=!is.element(model,c("rrl","ord"))
odmax=rep(max(apply(Y>0,c(1,3),sum,na.rm=TRUE)),nrow(Y[,,1]))
seed=6886
plot=FALSE ;  print = FALSE ;  gof=TRUE

# restructure Y
Y <- lapply(1:dim(YX_bin_long$Y)[3], function(t){YX_bin_long$Y[,,t]})
Xdyad <- lapply(1:dim(YX_bin_long$X)[4], function(t){YX_bin_long$X[,,,t]})

# add labels
set.seed(6886) ; actors <- as.character( sample(300:700,size=50,replace=FALSE) )
Y <- lapply(Y, function(y){ dimnames(y)[[1]] <- actors ; dimnames(y)[[2]] <- actors ; return(y) })
Xdyad <- lapply(Xdyad, function(x){ dimnames(x)[[1]] <- actors ; dimnames(x)[[2]] <- actors ; return(x) })

# change actor composition


################################################################
# reset odmax param
odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), nrow(Y[[1]]) )

# set random seed 
set.seed(seed)

# make sure actor labels are provided
dvRowLabCnt <- do.call('sum', lapply(Y, function(y) !is.null( dimnames(y)[[1]] ) ) )
dvColLabCnt <- do.call('sum', lapply(Y, function(y) !is.null( dimnames(y)[[2]] ) ) )
if( dvRowLabCnt!=length(Y) | dvColLabCnt!=length(Y) ){
  stop('Actor labels are missing from Y.') }

if(!is.null(Xdyad)){
xDyadRowLabCnt <- do.call('sum', lapply(Xdyad, function(x) !is.null( dimnames(x)[[1]] ) ) )
xDyadvColLabCnt <- do.call('sum', lapply(Xdyad, function(x) !is.null( dimnames(x)[[2]] ) ) )
if( xDyadRowLabCnt!=length(Y) | xDyadvColLabCnt!=length(Y) ){
  stop('Actor labels are missing from Xdyad.') }
}

if(!is.null(Xrow)){
  xRowLabCnt <- do.call('sum', lapply(Xrow, function(x) !is.null(rownames(x)) ))
  if( xRowLabCnt != length(Y) ){
    stop('Actor labels are missing from Xrow.') }
}

if(!is.null(Xcol)){
  xColLabCnt <- do.call('sum', lapply(Xcol, function(x) !is.null(rownames(x)) ))
  if( xColLabCnt != length(Y) ){
    stop('Actor labels are missing from Xcol.') }
}

# create full frame of actors
fullActorSet <- unique(unlist(lapply(Y,rownames)))

# set diag to NA 
N<-length(Y) ; Y<-lapply(Y, function(y){diag(y)=NA; return(y)})

# force binary if binary model specified 
if(is.element(model, c('bin','cbin'))){ Y<-lapply(Y,function(y){1*(y>0)})}

# observed and max outdegrees 
if(is.element(model,c("cbin","frn","rrl")))
{
odobs<-do.call('cbind',lapply(Y, function(y){ apply(y>0,1,sum,na.rm=TRUE)}))
if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y[[1]])) } 
}

# some settings for symmetric case
if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }

# construct design matrix    
pr<-ifelse(is.null(Xrow),0,ncol(Xrow[[1]]))
pc<-ifelse(is.null(Xcol),0,ncol(Xcol[[1]]))
pd<-ifelse(is.null(Xdyad),0,dim(Xdyad[[1]])[3])

X<-lapply(1:N, function(t){
  Xt<-design_array(Xrow[[t]], Xcol[[t]], Xdyad[[t]], intercept, nrow(Y[[t]]))

  # re-add intercept if it was removed
  if(dim(Xt)[3] < sum(pr,pc,pd,intercept) )
  {
    tmp <- array(dim=dim(Xt)+c(0,0,1))
    tmp[,,1]<-1 ; tmp[,,-1]<-Xt
    Xt<-tmp
  }
  return(Xt)
})
names(X) <- dimnames(Y)[[3]]

# # design matrix warning for rrl
# if( model=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0)
#                & !any( apply(X,c(3),function(x){var(c(x))})==0) )
# {
# cat("WARNING: row effects are not estimable using this procedure ","\n")
# }

# # design matrix warning for rrl and ord
# if( is.element(model,c("ord","rrl")) & 
#   any( apply(X,c(3),function(x){var(c(x))})==0 ) )
# {
# cat("WARNING: an intercept is not estimable using this procedure ","\n")
# }


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
    return( z - z01 )
  }

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
    return(z)
  }
  })

# row and column means
#### need to pull out "last" value for every actor's a/b into n vector for prior
ab<-lapply(Z, function(z){
  a<-rowMeans(z, na.rm=TRUE) ; b<-colMeans(z, na.rm=TRUE)
  return( list(a=a,b=b) ) })
a <- ab[[length(ab)]]$a ; b <- ab[[length(ab)]]$b

# starting values for missing entries  
Z <- lapply(Z, function(z){
  mu<-mean(z,na.rm=TRUE)
  a<-rowMeans(z,na.rm=TRUE) ; b<-colMeans(z,na.rm=TRUE)  
  za<-mu + outer(a,b,'+')
  z[is.na(z)] <- za[is.na(z)]
  return(z)
  })

# other starting values
beta<-rep(0,dim(X[[1]])[3]) 
s2<-1 
rho<-0
Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
U<-V<-matrix(0, length(fullActorSet), R) 


#  output items
BETA <- matrix(nrow = 0, ncol = dim(X[[1]])[3] - pr*symmetric)
VC<-matrix(nrow=0,ncol=5-3*symmetric) 
UVPS <- U %*% t(V) * 0 
APS<-BPS<- rep(0,length(fullActorSet))  
YPS <- lapply(Y, function(y){ y*0 })
GOF<-matrix(rowMeans( do.call('cbind', lapply(Y, function(y){ gofstats(y) }) ) ),1,4)  
rownames(GOF)<-"obs"
colnames(GOF)<- c("sd.rowmean","sd.colmean","dyad.dep","triad.dep")
names(APS)<-names(BPS)<-rownames(U)<-rownames(V)<-rownames(Y[,,1])

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

class(fit) <- "ame"