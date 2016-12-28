rm(list=ls())
library(amen)
library(rbenchmark)
library(ggplot2)
library(reshape2)
library(magrittr)

############################################################
# testing whether obs with NAs in Y end up affecting the likelihood
# and thus parameter estimates

data(YX_bin_long)

Y=YX_bin_long$Y ; Xdyad=YX_bin_long$X ; 
for(t in 1:dim(Y)[3]){diag(Y[,,t])=NA}
tmp = array(NA, dim=c(50,50,4,4),dimnames=list(NULL,NULL,c('dyadVar1','dyadVar2','dyadVar3','binomDyadVar',NULL)))
tmp[,,1:3,] = Xdyad
set.seed(6886) ; stuff = rbinom(50*50*4,1,.7)
# stuff = (stuff - mean(stuff))/sd(stuff)
tmp[,,4,] = array(stuff,dim=c(50,50,4))
Xdyad = tmp

set.seed(6886)
fit <- ame_rep(Y, Xdyad, R=2,
                   model='nrm',
                   burn=500,nscan=1000,odens=1,seed=6886,
                   plot=FALSE, print=FALSE)

# run with glm
yLong = melt(Y) ; yLong$id = paste(yLong$Var1, yLong$Var2, yLong$Var3, sep='_')
xLong = melt(Xdyad) %>% dcast(.,Var1+Var2+Var4~Var3,value.var='value')
xLong$id = paste(xLong$Var1, xLong$Var2, xLong$Var4, sep='_')
datLong = xLong ; datLong$dv = yLong$value[match(datLong$id, yLong$id)]
naive=lm(dv~dyadVar1 + dyadVar2 + dyadVar3 + binomDyadVar, data=datLong)

# add some NAs to Y array
Y_withMiss = Y
Y_withMiss[c(10,30,40),c(32,43,41),c(1,3)] = NA
Y_withMiss[c(21,18,14,1,2,48,47),c(5,13,50,9,7,6,21),c(2,4)] = NA
Y_withMiss[abs(c(21,18,14,1,2,48,47)-7),c(5,13,43,9,7,6,21)+7,c(1,4)] = NA
set.seed(6886) ; missRow = sample(1:50, 20) ; set.seed(7686) ; missCol = sample(1:50, 20) ; 
Y_withMiss[missRow,missCol,c(2,3)] = NA
set.seed(2132) ; missRow2 = sample(1:50, 20) ; set.seed(87676) ; missCol2 = sample(1:50, 20) ; 
Y_withMiss[missRow2,missCol2,c(1,2)] = NA
fit_withMiss <- ame_rep(Y_withMiss, Xdyad, R=2,
                   model='nrm',
                   burn=500,nscan=1000,odens=1,seed=6886,
                   plot=FALSE, print=FALSE)

# run with glm
yLong = melt(Y_withMiss) ; yLong$id = paste(yLong$Var1, yLong$Var2, yLong$Var3, sep='_')
xLong = melt(Xdyad) %>% dcast(.,Var1+Var2+Var4~Var3,value.var='value')
xLong$id = paste(xLong$Var1, xLong$Var2, xLong$Var4, sep='_')
datLong = xLong ; datLong$dv = yLong$value[match(datLong$id, yLong$id)]
naive_withMiss=lm(dv~dyadVar1 + dyadVar2 + dyadVar3 + binomDyadVar, data=datLong)

# change values in Xdyad
Xdyad_weird = Xdyad
Xdyad_weird[c(10,30,40),c(32,43,41),,c(1,3)] = NA
Xdyad_weird[c(21,18,14,1,2,48,47),c(5,13,50,9,7,6,21),,c(2,4)] = NA
Xdyad_weird[abs(c(21,18,14,1,2,48,47)-7),c(5,13,43,9,7,6,21)+7,,c(1,4)] = NA
Xdyad_weird[missRow,missCol,,c(2,3)] = NA
Xdyad_weird[missRow2,missCol2,,c(1,2)] = NA
fit_withMiss_weirdX <- ame_rep(Y_withMiss, Xdyad_weird, R=2,
                        model='nrm',
                        burn=500,nscan=1000,odens=1,seed=6886,
                        plot=FALSE, print=FALSE)

# change values in Xdyad
Xdyad_weird = Xdyad
Xdyad_weird[c(10,30,40),c(32,43,41),c(1,4),c(1,3)] = NA
Xdyad_weird[c(21,18,14,1,2,48,47),c(5,13,50,9,7,6,21),c(2,3),c(2,4)] = NA
Xdyad_weird[abs(c(21,18,14,1,2,48,47)-7),c(5,13,43,9,7,6,21)+7,1,c(1,4)] = NA
Xdyad_weird[missRow,missCol,c(2,3),c(2,3)] = NA
Xdyad_weird[missRow2,missCol2,4,c(1,2)] = NA
fit_withMiss_weirdX2 <- ame_rep(Y_withMiss, Xdyad_weird, R=2,
                               model='nrm',
                               burn=500,nscan=1000,odens=1,seed=6886,
                               plot=FALSE, print=FALSE)

# compare
sum(is.na(Y_withMiss))/length(Y_withMiss)
apply(fit$BETA, 2, mean)
coef(naive)
coef(naive_withMiss)
apply(fit_withMiss$BETA, 2, mean)
apply(fit_withMiss_weirdX$BETA, 2, mean)
apply(fit_withMiss_weirdX2$BETA, 2, mean)

summary(fit_withMiss)
summary(fit_withMiss_weirdX)
summary(fit_withMiss_weirdX2)

summary(naive_withMiss)

ugh = melt(fit_withMiss$BETA)
ugh$model='withMiss'
ugh2 = melt(fit_withMiss_weirdX$BETA)
ugh2$model = 'withMiss_weirdX'
ugh = rbind(ugh,ugh2)
ggplot(ugh, aes(x=Var1, y=value, color=model)) + geom_line() + facet_wrap(model~Var2,scales='free')

############################################################

############################################################
############################################################
## Running with changing actor composition

data(YX_bin_long)

Y=YX_bin_long$Y ; Xdyad=YX_bin_long$X ; 

# restructure Y
Y <- lapply(1:dim(YX_bin_long$Y)[3], function(t){YX_bin_long$Y[,,t]})
Xdyad <- lapply(1:dim(YX_bin_long$X)[4], function(t){YX_bin_long$X[,,,t]})

# add labels
set.seed(6886) ; actors <- as.character( sample(300:700,size=50,replace=FALSE) )
Y <- lapply(Y, function(y){ dimnames(y)[[1]] <- actors ; dimnames(y)[[2]] <- actors ; return(y) })
Xdyad <- lapply(Xdyad, function(x){ dimnames(x)[[1]] <- actors ; dimnames(x)[[2]] <- actors ; return(x) })

set.seed(6886) ; toRem = sample(1:length(actors), 5)
toRem = actors[toRem]
Y[[1]] = Y[[1]][-which(actors %in% toRem),-which(actors %in% toRem)]
Y[[2]] = Y[[2]][-which(actors %in% toRem),-which(actors %in% toRem)]
Xdyad[[1]] = Xdyad[[1]][-which(actors %in% toRem),-which(actors %in% toRem),]
Xdyad[[2]] = Xdyad[[2]][-which(actors %in% toRem),-which(actors %in% toRem),]

actorSet = actors
nActors = length(actors)
N = length(Y)
pdLabs = NULL

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

# run mod
fitList<-ame_repL2(Y,Xdyad,R=2,
	model='bin',
	burn=10,nscan=100,odens=1,plot=FALSE, print=FALSE)

# load data
data(YX_bin_long)
Y2=YX_bin_long$Y ; Xdyad2=YX_bin_long$X 

# add labels
set.seed(6886) ; actors2 <- as.character( sample(300:700,size=50,replace=FALSE) )
rownames(Y2) <- colnames(Y2) <- rownames(Xdyad2) <- colnames(Xdyad2) <- actors2
actors2 <- rownames(Y2)
set.seed(6886) ; toRem2 = sample(1:length(actors2), 5)
toRem2 = actors2[toRem2]

# add NAs
Y2[toRem2, ,1:2] <- NA
Y2[,toRem2 ,1:2] <- NA
Xdyad2[toRem2, , , 1:2] <- NA
Xdyad2[, toRem2, , 1:2] <- NA

fitOrig <- ame_rep(Y2, Xdyad2, R=2,
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
	# for(t in 1:N){ if(round(sum(orig$EZ[,,t]-modded$EZ[[t]]),10)!=0){ stop("EZ results don't match.") } }
	# for(t in 1:N){ if(round(sum(orig$YPM[,,t]-modded$YPM[[t]],na.rm=TRUE),10)!=0){ stop("YPM results don't match.") } }
	for(t in 1:N){ if(round(sum(orig$EZ[,,t]-modded$EZ[,,t]),10)!=0){ stop("EZ results don't match.") } }
	for(t in 1:N){ if(round(sum(orig$YPM[,,t]-modded$YPM[,,t],na.rm=TRUE),10)!=0){ stop("YPM results don't match.") } }
}
##############################

##############################
# load data
data(YX_bin_long) ; data(YX_bin_list) # same as other but replicates are stored in list

# binary ame_rep with symmetric=FALSE, R=0, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(YX_bin_long$Y,YX_bin_long$X,R=0,
	model='bin',
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL2(YX_bin_list$Y,YX_bin_list$X,R=0,
	model='bin',
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=FALSE, R=2, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(YX_bin_long$Y,YX_bin_long$X,R=2,
	model='bin',
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL2(YX_bin_list$Y,YX_bin_list$X,R=2,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=TRUE, R=2, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(YX_bin_long$Y,YX_bin_long$X,R=2,
	model='bin', symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL2(YX_bin_list$Y,YX_bin_list$X,R=2,
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
fitList<-ame_repL2(yL,xDyadL,xNodeL,xNodeL,R=0,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=FALSE, R=2, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(Y,Xdyad,Xnode,Xnode,R=2,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL2(yL,xDyadL,xNodeL,xNodeL,R=2,
	model='bin', 
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.

# binary ame_rep with symmetric=FALSE, R=2, symmetric=TRUE, rvar=TRUE, cvar=TRUE
fitOrig<-ame_rep(Y,Xdyad,Xnode,Xnode,R=2,
	model='bin', symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
fitList<-ame_repL2(yL,xDyadL,xNodeL,xNodeL,R=2,
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
fitList<-ame_repL2(yL,xDyadL,xNodeL,R=1,
	model="ord",symmetric=TRUE,
	burn=5,nscan=5,odens=1,plot=FALSE, print=FALSE)
runTests(fitOrig,fitList) # if nothing returned that's good.
##############################

############################################################
############################################################
