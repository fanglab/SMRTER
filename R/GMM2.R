#' GMM2 
#'
#' GMM2 estimates the Gaussian mixture models
#' 
#' @param x a vector of numeric values 
#' @param mu a vector of values indicating the means of mixture models
#' @param fixed a vector of values indicating the index of 
#'
#' @return a list of values comprising the following components
#' \itemize{
#'  \item {G} an integer value indicating the number of Gaussian mixture models
#'  \item {lamda} a vector of numeric values indicating the propotion of each mixture model
#'  \item {mu} a vector of numeric values indicating the mean of each mixture model
#'  \item {sd} a vector of numeric values indicating the standard deviation of each mixture model
#' }
#' 
#' @export
#'
#' @examples
#' 
#' x1 <- rnorm(300,mean=1,sd=0.5)
#' x2 <- rnorm(700,mean=4,sd=1)
#' x <- c(x1,x2)
#' plot(density(x))
#' 
#' # fix both two means
#' GMM2(x,mu=c(1,4),fixed=c(1,2) )
#' # fix one mean, and estimate the other 
#' GMM2(x,mu=c(1,10),fixed=c(1) )
#' # do not fix means, and estimate both
#' GMM2(x,mu=c(0,10),fixed=c() )
#' 
#' 
GMM2 <- function(x , mu , fixed=c(1:length(mu)) )
{
  G = length(mu) # G is the number of mixture model
  sdv = rep(1,G)
  lamda = rep(1,G)/G
  
  for(iter in 1:100)
  {
    ###  update P=p(x|q,mu,sdv)
    P_gaussian = lapply( 1:length(lamda) , function(i) { dnorm( x , mean=mu[i], sd=sdv[i] )+1e-300 }  )   
    
    ### update Z                                                                                                                                                                                              
    Z =  sapply( 1:length(lamda), function(i)  {  P_gaussian[[i]] * lamda[i] } )  
    rowsum = rowSums(Z)
    Z =  apply(Z,2,function(x){ x/rowsum })	 # p(q|x,theta)
    
    ### update mu
    for(i in setdiff( c(1:G) , fixed ))
    {	mu[i]  =  sum(x * Z[,i],na.rm=T) / sum(Z[,i],na.rm=T)   }
    
    ### update C and Theta
    sdv =  sapply( 1:length(mu), function(i) {  
      var = sum( (x-mu[i])^2 * Z[,i] ,na.rm=T ) / sum(Z[,i],na.rm=T ) ;  
      sqrt(var) 
      })  
    
    ### update lamda
    lamda = colSums(Z,na.rm=T)/nrow(Z) # p(j|theta); w
    
    #cat(lamda[1],lamda[2],mu[1],mu[2],sdv[1],sdv[2],'\n')
  }
  
  list(G=G , lamda=lamda , mu=mu , sd=sdv)
  
}

