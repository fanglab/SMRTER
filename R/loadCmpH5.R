#' loadCmpH5
#'
#' loadCmpH5 loads the cmp.h5 format data
#' 
#' @param path a character value indicating the file path of cmp.h5 file
#'
#' @return a list comprising the following components
#' \itemize{
#'  \item {cmpH5} a h5 data saving the IPDs and sequence for each subread
#'  \item {index} a data frame saving the basic information of each subread, e.g., chromosome, start/end location in genome, and start/end location in the read 
#'  \item {accuracy} an vector of read accuracy of one subread
#'  \item {readlength} a vector of integer values indicating the subread length
#' }

#' @export
#'
#' @examples
#' 
loadCmpH5 = function( path )
{
	system.time( cmpH5 <- PacBioCmpH5(path) )
	system.time( index <- alnIndex(cmpH5)  )
	index$alignedStrand = 1 -  index$alignedStrand 
	# use 1 minus, to covert the read space to template space; 
	# you can confirm by comparing with the csv file for the coverage of a specific position of either strand
	system.time( accuracy <- getAccuracy(cmpH5,denominator=getTemplateSpan) )
	system.time( readlength <- getReadLength(cmpH5) )
	
	list(  cmpH5=cmpH5 , index=index , accuracy=accuracy , readlength=readlength )
}


#library(foreach)
#library(doMC)
#library(  h5r , lib.loc='/hpc/users/zhus02/setup/my_own_R_package/R_3.3.1/')
#library( pbh5 , lib.loc='/hpc/users/zhus02/setup/my_own_R_package/R_3.3.1/')

