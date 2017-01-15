rm(list=ls())
library(amen)
library(rbenchmark)

############################################################
# Mansfield milner example
load('~/Dropbox/Research/netsMatter/replications/mansfield_milner_2012/inputData/amenData.rda')
set.seed(6886) ; fitOrig <- ame_repL(
  Y=yList, Xdyad=xDyadList, Xrow=xNodeList,
  model='bin', symmetric=TRUE, R=2, intercept=FALSE,
  burn=0, nscan=200, odens=25, seed=6886,
  plot=TRUE, print=FALSE, gof=TRUE
  )
############################################################

############################################################
## Running with changing actor composition
# load data
data(YX_bin_long) ; data(YX_bin_list) # same as other but replicates are stored in list

# randomly delete some nodes
yL = YX_bin_list$Y
xDyad = YX_bin_list$X
actors = rownames(yL[[1]])

set.seed(6886) ; toRem = sample(1:length(actors), 5)
yL[[1]] = yL[[1]][-toRem,-toRem]
xDyad[[1]] = xDyad[[1]][-toRem,-toRem,]

# run mod
fitList<-ame_repL(yL,xDyad,R=2,
	model='bin',
	burn=10,nscan=40,odens=1,plot=FALSE, print=FALSE)
############################################################

##############################
# fn to compare mod results
applyR = function(...){round(apply(...),4)}
meanR = function(...){round(mean(...),4)}
compareResults <- function(mod1, mod2, symm){
  beta = cbind(applyR(mod1$BETA, 2, mean), applyR(mod2$BETA, 2, mean))
  vc = cbind(applyR(mod1$VC, 2, mean), applyR(mod2$VC, 2, mean))
  apm = c(meanR(mod1$APM), meanR(mod2$APM))
  if(!symm){bpm = c(meanR(mod1$BPM), meanR(mod2$BPM))}
  u = cbind(applyR(mod1$U, 2, mean), applyR(mod2$U, 2, mean))
  if(!symm){v = cbind(applyR(mod1$V, 2, mean), applyR(mod2$V, 2, mean))}
  if(symm){return(list(beta=beta,vc=vc,apm=apm,u=u))}
  if(!symm){return(list(beta=beta,vc=vc,apm=apm,bpm=bpm,u=u,v=v))}
}
##############################

##############################
# load data
data(YX_bin_long) ; data(YX_bin_list) # same as other but replicates are stored in list

Y = YX_bin_list$Y
Xdyad = YX_bin_list$X
Xrow=NULL ; Xcol=NULL ; intercept=TRUE
# get actor info
actorByYr <- lapply(Y, rownames)
actorSet <- sort(unique(unlist( actorByYr ))) ; n <- length(actorSet)
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

set.seed(6886)
benchmark(
  fitOrig <- ame_rep(
    Y, Xdyad, R=2,
    model='nrm', seed=6886,symmetric=FALSE,intercept=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
  ),
  fitList <- ame_repL(
    Y=YX_bin_list$Y, Xdyad=YX_bin_list$X, R=2,
    model='nrm', seed=6886,symmetric=FALSE,intercept=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
  ),
  replications=1
)
compareResults(fitOrig, fitList, FALSE)

set.seed(6886)
benchmark(
  fitOrig <- ame_rep(
    Y, Xdyad, R=2,
    model='nrm', seed=6886,symmetric=TRUE,intercept=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
  ),
  fitList <- ame_repL(
    Y=YX_bin_list$Y, Xdyad=YX_bin_list$X, R=2,
    model='nrm', seed=6886,symmetric=TRUE,intercept=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
  ),
  replications=1
)
compareResults(fitOrig, fitList, TRUE)
##############################

##############################
# load some other data
data(dutchcollege)

# code from hoff vignette pg 39-40
Y<-1*( dutchcollege$Y >= 2 )[,,2:7]
n<-dim(Y)[1] ; t<-dim(Y)[3]
# nodal covariates
Xnode<-dutchcollege$X[,1:2]
Xnode<-array(Xnode,dim=c(n,ncol(Xnode),t))
dimnames(Xnode)[[2]]<-c("male","smoker")
# dyadic covariates
Xdyad<-array(dim=c(n,n,5,t))
Xdyad[,,1,]<-1*( dutchcollege$Y >= 2 )[,,1:6]
Xdyad[,,2,]<-array(apply(Xdyad[,,1,],3,t),dim=c(n,n,t))
Xdyad[,,3,]<-tcrossprod(Xnode[,1,1])
Xdyad[,,4,]<-tcrossprod(Xnode[,2,1])
Xdyad[,,5,]<-outer( dutchcollege$X[,3],dutchcollege$X[,3],"==")
dimnames(Xdyad)[[3]]<-c("Ylag","tYlag","bothmale","bothsmoke","sameprog")
# create list version
# note that list version of fn requires actor labels
set.seed(6886) ; actors=as.character(sample(2:600,n))
yL <- lapply(1:t, function(x){ y=Y[,,x] ; rownames(y)=colnames(y)=actors ; y})
xNodeL <- lapply(1:t, function(x){ xnode=Xnode[,,x] ; rownames(xnode)=actors ; xnode })
xDyadL <- lapply(1:t, function(x){ xdyad=Xdyad[,,,x] ; rownames(xdyad)=colnames(xdyad)=actors ; xdyad})

benchmark(
  fitOrig<-ame_rep(
    Y,Xdyad,Xnode,Xnode,R=0, model='bin',
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  fitList<-ame_repL(
    yL,xDyadL,xNodeL,xNodeL,R=0, model='bin',
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  replications=1
)
compareResults(fitOrig, fitList, FALSE)

benchmark(
  fitOrig<-ame_rep(
    Y,Xdyad,Xnode,Xnode,R=2, model='bin',
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  fitList<-ame_repL(
    yL,xDyadL,xNodeL,xNodeL,R=2, model='bin',
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  replications=1
)
compareResults(fitOrig, fitList, FALSE)

benchmark(
  fitOrig<-ame_rep(
    Y,Xdyad,Xnode,Xnode,R=2,model='bin', symmetric=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  fitList<-ame_repL(
    yL,xDyadL,xNodeL,xNodeL,R=2,model='bin', symmetric=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  replications=1
)
compareResults(fitOrig, fitList, TRUE)
##############################

##############################
# load some other data
data(coldwar)

# code modified from hoff vignette pg 43-44
Y=sign( coldwar$cc ) ; n=nrow(Y) ; t=dim(Y)[3]
# nodal covariates
Xn <- array(NA,dim=c(n,2,t),dimnames=list(rownames(Y),c('lgdp','polity'),NULL))
for(x in 1:t){ Xn[,'lgdp',x] <- log(coldwar$gdp[,x]) ; Xn[,'polity',x] <- coldwar$polity[,x] }
# dyadic covariates
Xd<-array(dim=c(nrow(Y),nrow(Y),1,t),dimnames=list(rownames(Y),rownames(Y),'ldist',NULL))
Xd[,,'ldist',] <- log(coldwar$distance)
# create list version
yL <- lapply(1:t, function(x){ Y[,,x] })
xNodeL <- lapply(1:t, function(x){ Xn[,,x] })
xDyadL <- lapply(1:t, function(x){ array(Xd[,,,x],dim=c(n,n,1),dimnames=list(rownames(Y),rownames(Y),'ldist')) })

# ordinal ame_rep with symmetric=TRUE, R=1
benchmark(
  fitOrig<-ame_rep(
    Y,Xd,Xn,R=1, model="ord",symmetric=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  fitList<-ame_repL(
    yL,xDyadL,xNodeL,R=1, model="ord",symmetric=TRUE,
    burn=1000,nscan=2000,odens=25,plot=FALSE, print=FALSE
    ),
  replications=1
)
compareResults(fitOrig, fitList, TRUE)
##############################
