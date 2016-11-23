#' Gibbs update for dyadic variance with independent replicate relational data
#' 
#' Gibbs update for dyadic variance with independent replicate relational data
#' 
#' 
#' @usage rs2_repL_fc(E.T, rho)
#' @param E.T List of square residual relational matrices. The list
#' is for different replicates. Each object in the list
#' is a square residual relational matrix
#' @param rho current value of rho
#' @return a new value of s2
#' @author Peter Hoff, Yanjun He
#' @export rs2_repL_fc
rs2_repL_fc <-
  function(E.T,rho)
  {
    N<-length(E.T)
    H<-mhalf( solve(matrix(c(1,rho,rho,1),2,2)) )
    EM<-do.call('rbind', lapply(E.T, function(E){
      em <- cbind(E[upper.tri(E)],t(E)[upper.tri(E)] ) %*%H
      return(em) }) )
    1/rgamma(1, (length(EM)+1)/2 , (sum(EM^2)+1)/2 )
  }
