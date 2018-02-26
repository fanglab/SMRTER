#' readSmrtCsv
#'
#' readSmrtCsv reads the csv format smrt data (modifications.csv) 
#' 
#' @param file a character value indicating the file name of the csv format smrt data 
#'
#' @return a data frame of smrt data comprising the following columns
#' \itemize{
#'  \item {refName} a character value indicating the chromosomes
#'  \item {strand} a character value indicating the strand taking on either 0 (reference) or 1 (opposite)
#'  \item {tpl} an integer value indicating the genomic location 
#'  \item {base}	a character value indiating the cognate base at this position in the reference
#'  \item {score} an integer value indicating the Phred-transformed pvalue of modification detection, i.e. 10*(-log10(pvalue))
#'  \item {coverage} an integer value indicating the count of valid IPDs for the in-silico control mode or the mean of case and control coverage for the case-control mode
#'  \item {ipdRatio} a numeric value indicating the ipd Ratio
#'  \item {modelPrediction/controlMean} a numeric value indicating the normalized mean IPD for the in-silico control mode/case-control mode
#' }
#' 
#' @export
#'
#' @examples
#' 
#' @seealso \code{\link{readSmrtGff}}; 
#'
readSmrtCsv = function( file )
{
  RData = gsub( "csv$" , 'RData' , file )
  
  if( file.exists(RData) )
  {
    loadFile = load( RData )
    cat('loaded ',loadFile,'\n')
  } else {			
    smrt_portal <- read.table( file , sep=',' , header=T )
    save( smrt_portal , file=RData )
  }

  smrt_portal
  
}
