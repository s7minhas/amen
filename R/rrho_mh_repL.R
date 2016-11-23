#' Metropolis update for dyadic correlation with independent replicate data
#' 
#' Metropolis update for dyadic correlation with independent replicate data. 
#' 
#' 
#' @usage rrho_mh_repL(E.T, rho, s2 = 1)
#' @param E.T List of square residual relational matrices. The list objects
#' contain different replicates. Each list object
#' is a square residual relational matrix. 
#' @param rho current value of rho
#' @param s2 current value of s2
#' @return a new value of rho
#' @author Peter Hoff, Yanjun He
#' @export rrho_mh_repL
rrho_mh_repL <-
  function(E.T,rho,s2=1)
  {
    N<-length(E.T)
    EM <- do.call('rbind', lapply(E.T, function(E){
      em<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] )/sqrt(s2)
      return(em) }) )
    
    emcp<-sum(EM[,1]*EM[,2])
    emss<-sum(EM^2)
    
    m<- nrow(EM)
    sr<- 2*(1-cor(EM)[1,2]^2)/sqrt(m)
    
    rho1<-rho+sr*qnorm( runif(1,pnorm( (-1-rho)/sr),pnorm( (1-rho)/sr)))
    
    lhr<-(-.5*(m*log(1-rho1^2)+(emss-2*rho1*emcp)/(1-rho1^2)))-
      (-.5*(m*log(1-rho^2 )+(emss-2*rho*emcp )/(1-rho^2 )))   +
      ( (-.5*log(1-rho1^2)) - (-.5*log(1-rho^2)) )
    
    if(log(runif(1))<lhr) { rho<-rho1 }
    rho
  }
