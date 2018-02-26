#' splitSmrt
#'
#' splitSmrt splits the smrt data into chunks based on both refName and strand
#' 
#' @param smrt a data frame of smrt data comprising at least refName, strand, tpl
#' @param sorted a boolean value indicating whether the output is sorted based on tpl
#'
#' @return a list comprising of all refNames on strand of either 0 or 1. The name of each item in the list is "refName.strand". 
#' 
#' @export
#'
#' @examples
#' 
#' 
splitSmrt = function( smrt , sorted=TRUE )
{
  res = split( smrt , as.character(smrt$refName) )
  do.call( c, lapply( res , function(r) {
    t0 = split( r , as.character(r$strand) )
    if( sorted ) { 
      t0 = lapply( t0 , function(t1) t1[order(t1$tpl),] )
    }
    t0
  } ) )
}
