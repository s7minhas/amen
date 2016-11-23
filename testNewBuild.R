rm(list=ls())
library(devtools) ; devtools::install_github('s7minhas/amen')
library(amen)

# New functions involving replicated data are labeled as xx_repL instead of xx_rep

## Do no harm
# First make sure that results from modified functions are the same as original
data(YX_bin_long)
data(YX_bin_list) # same as above but replicates are stored in list

# original fn
fitOrig<-ame_rep(YX_bin_long$Y,YX_bin_long$X,burn=5,nscan=5,odens=1,model="bin",plot=FALSE)
fitList<-ame_repL(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,model="bin",plot=FALSE)

