#' callFDR
#'
#' callFDR calls False Discovery Rate by comparing the foreground values with the background values
#' 
#' @param foreground a vector of numeric values, either IPD ratio or Score from the native DNA
#' @param background a vector of numeric values, either IPD ratio or Score from the WGA
#' @param thresList a vector of numeric values indicating the threshold list for FDR calling
#' @param adjust a boolean value indicating whether the estimated FDRs are adjusted to be decreased with the increased thresholds
#' @param verbose a boolean value indicating whether output on screen
#'
#' @return a data frame values comprsing two columns
#' \itemize{
#'  \item {thres} a vector of numeric values indicating the threshold list 
#'  \item {fdr} a vector of numeric values indicating the estimated FDRs correpsonding to the thresholds
#' }
#' 
#' @export
#'
#' @examples
#' 
#' 
callFDR = function( foreground , background , thresList=NULL , adjust=T , verbose=F )
{
  lenF = length(foreground)
  lenB = length(background)
  if(is.null(thresList)) { thresList = seq(  min(foreground) , max(foreground) , length.out=300 ) }
  
  fdr = sapply( thresList , function(thres) 
    { 
      rf = ( sum(foreground>=thres)+1 ) / (lenF+1) 
      rb = ( sum(background>=thres)+1 ) / (lenB+1) 
      f=rb/rf 
      if(verbose) cat(thres,f,'\n') 
      f 
    } )
  
  if(adjust) 
  { 
    for(i in 2:length(fdr)) 
    { 
      if(fdr[i-1]>1) fdr[i-1]=1
      fdr[i] = ifelse( fdr[i-1]<fdr[i] , fdr[i-1] , fdr[i] )  
    } 
  }
  
  data.frame( thres=thresList, fdr=fdr )
  
}

