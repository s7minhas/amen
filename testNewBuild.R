rm(list=ls())
library(devtools) ; devtools::install_github('s7minhas/amen')
library(amen)
library(rbenchmark)

############################################################
############################################################
## New functions involving replicated data are labeled as xx_repL instead of xx_rep:
# ame_repL.R
# rbeta_ab_repL_fc.R
# rrho_mh_repL.R
# rs2_repL_fc.R
############################################################
############################################################

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
	burn=10,nscan=100,odens=1,plot=FALSE, print=FALSE)
############################################################
############################################################

############################################################
############################################################
## Do no harm: Make sure that results from modified functions are the same as original

##############################
# test to make sure they return the same results
runTests <- function(orig, modded){
	if(round(sum(orig$BETA-modded$BETA),10)!=0){ stop("BETA results don't match.") }
	if(round(sum(orig$VC-modded$VC),10)!=0){ stop("VC results don't match.") }
	if(round(sum(orig$APM-modded$APM),10)!=0){ stop("APM results don't match.") }
	if(round(sum(orig$BPM-modded$BPM),10)!=0){ stop("BPM results don't match.") }
	if(round(sum(orig$U-modded$U),10)!=0){ stop("U results don't match.") }
	if(round(sum(orig$V-modded$V),10)!=0){ stop("V results don't match.") }
	if(round(sum(orig$UVPM-modded$UVPM),10)!=0){ stop("UVPM results don't match.") }
	if(round(sum(orig$GOF-modded$GOF),10)!=0){ stop("GOF results don't match.") }
	N=dim(YX_bin_long$Y)[3]
	for(t in 1:N){ if(round(sum(orig$EZ[,,t]-modded$EZ[[t]]),10)!=0){ stop("EZ results don't match.") } }
	for(t in 1:N){ if(round(sum(orig$YPM[,,t]-modded$YPM[[t]],na.rm=TRUE),10)!=0){ stop("YPM results don't match.") } }
}
##############################

##############################
# load data
data(YX_bin_long) ; data(YX_bin_list) # same as other but replicates are stored in list

# binary ame_rep with symmetric=FALSE, R=0, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(YX_bin_long$Y,YX_bin_long$X,R=0,
	model='bin',
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL(YX_bin_list$Y,YX_bin_list$X,R=0,
	model='bin',
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=FALSE, R=2, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(YX_bin_long$Y,YX_bin_long$X,R=2,
	model='bin',
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL(YX_bin_list$Y,YX_bin_list$X,R=2,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=TRUE, R=2, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(YX_bin_long$Y,YX_bin_long$X,R=2,
	model='bin', symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL(YX_bin_list$Y,YX_bin_list$X,R=2,
	model='bin', symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.
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

# binary ame_rep with symmetric=FALSE, R=0, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(Y,Xdyad,Xnode,Xnode,R=0,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL(yL,xDyadL,xNodeL,xNodeL,R=0,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=FALSE, R=2, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(Y,Xdyad,Xnode,Xnode,R=2,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL(yL,xDyadL,xNodeL,xNodeL,R=2,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=FALSE, R=2, symmetric=TRUE, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(Y,Xdyad,Xnode,Xnode,R=2,
	model='bin', symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL(yL,xDyadL,xNodeL,xNodeL,R=2,
	model='bin', symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.
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
fitOrig<-ame_rep(Y,Xd,Xn,R=1,
	model="ord",symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL(yL,xDyadL,xNodeL,R=1,
	model="ord",symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.
##############################

############################################################
############################################################
