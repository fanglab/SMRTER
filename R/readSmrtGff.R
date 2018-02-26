#' readSmrtGff
#'
#' readSmrtGff reads the gff format smrt data (modifications.gff) 
#' 
#' @param file a character value indicating the file name of the gff format smrt data 
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
#'  \item {tag} a character value indicating the estimated modification types, including "modified_base","m6A","m5C" and "m4C"
#'  \item {identificationQv} an integer value indicating the score matching the modification signitures of "m6A","m5C" and "m4C"
#'  \item {context}	a character value indiating the bases in [-20,+20] of this position in the reference
#' }
#' 
#' @export
#'
#' @examples
#' 
#' @seealso \code{\link{readSmrtCsv}}; 
#'

readSmrtGff = function( file )
{
  
  RData = gsub( "gff$" , "RData" , file )
    
  if( file.exists(RData) )
  {
    loadFile = load( RData )
    cat('loaded ',loadFile,'\n')
  } else {	
    
    gff <- read.table( file , sep='\t' )
    
    split9 = strsplit( as.character(gff[,9]) , ';' )
    
    getFeature = function(feature) sapply( split9 , function(x) { ifelse( any( grepl(feature,x) ), x[grep(feature,x)] , NA ) } ) 
    
    features = c('coverage','context','IPDRatio','frac','fracLow','fracUp','identificationQv')
    info = sapply( features ,  function(x) { cat('getting feature: ', x,'\n'); getFeature(x)  } )
    info = sapply( 1:length(features) , function(i) { gsub( paste(features[i],'=',sep='') , '' , as.character(info[,i]) ) } )
    info = data.frame(info)
    colnames(info) = features
    
    strand = rep( 0 , nrow(gff) )
    strand[ as.character(gff[,7])=='-' ] = 1
    
    ic = function(x) as.integer(as.character(x))
    nc = function(x) as.numeric(as.character(x))
    
    base = substr( as.character(info$context) , 21 , 21 )
    
    gff = data.frame( refName = gff[,1] , tpl = ic(gff[,4]) , strand =ic(strand)  , score =ic(gff[,6]), ipdRatio = nc(info$IPDRatio), 
                      coverage = ic(info$coverage) , tag = gff[,3] , identificationQv = ic(info$identificationQv) , base = base , context = info$context )
    
    save( gff , file=RData )
  }
  
  gff
}

