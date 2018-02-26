#' hemiMethylation
#'
#' hemiMethylation checks strand-specific methylation from the Single-molecule level
#' 
#' @param SMipdRatio a list with the name of refName and strand "refName.strand". Each item is a data frame comprising the following columns:
#' \itemize{
#'  \item {tpl} an integer value indicating the genomic location 
#'  \item {ipdRatio} a numeric value indicating the IPD ratio at this position of one specific molecule
#' }
#' @param draw a boolean value indicating whehter draw the figure or not
#'
#' @return a data frame comprising the following columns:
#' \itemize{
#'  \item {pos_ipdRatio} a vector of numeric values indicating the IPD ratio on the positive strand
#'  \item {neg_ipdRatio} a vector of numeric values indicating the IPD ratio on the negative strand
#' }
#' 
#' @export
#'
#' @examples
#' 
#' @seealso \code{\link{smrterSyntheticControl}}; \code{\link{smrterTraditionalControl}}; 
#' 
#' 
hemiMethylation = function(SMipdRatio , draw=FALSE , point=c('plot','smoothScatter') )
{
  
  matchDoubleStrand = function( pos , neg )
  {
    pos_mol = rownames(pos)
    neg_mol = rownames(neg)
    pos_tpl = pos[,1]
    neg_tpl = neg[,1]-1
    
    pos_tag = paste( pos_mol , pos_tpl , sep='_' )
    neg_tag = paste( neg_mol , neg_tpl , sep='_' )
    
    pos_ipdRatio = pos[order(pos_tag) , 2 ]
    neg_ipdRatio = neg[order(neg_tag) , 2 ]
    pos_tag = sort(pos_tag)
    neg_tag = sort(neg_tag)
    
    pos_ipdRatio = pos_ipdRatio[ pos_tag %in% neg_tag ]
    neg_ipdRatio = neg_ipdRatio[ neg_tag %in% pos_tag ]
    
    cbind( pos_ipdRatio , neg_ipdRatio )
  }
  
  
  checkDS_fromSMipdRatio = function( SMipdRatio )
  {
    tag = names(SMipdRatio)
    split_tag = strsplit(tag,'[.]') 
    chr = sapply( split_tag , function(x){ x[length(x)-1] } )
    str = sapply( split_tag , function(x){ x[length(x)  ] } )
    
    matchedDS = lapply( unique(chr) , function(chri) {
      if( sum(chr==chri)==2 )
      {
        cat( 'mathcing ==> ',tag[chr==chri] , '\n' )
        pos = SMipdRatio[[ which( chr==chri & str==0 ) ]]
        neg = SMipdRatio[[ which( chr==chri & str==1 ) ]]
        matchDoubleStrand( pos, neg )
      } 
    } )
    
    do.call(rbind,matchedDS)
  }
  
  
	ipdRatio_2str = checkDS_fromSMipdRatio( SMipdRatio )
	
	if(draw)
	{
	# par(mfrow=c(2,2))
		
	#	if( 'plot' %in% point )
	#	plot(          log2(ipdRatio_2str[,1]) , log2(ipdRatio_2str[,2]) , xlab='Single molecule pos strand log2 ipdRatio' , 
	#					ylab='Single molecule neg strand log2 ipdRatio' ,main='Checking HemiMethylation')
						
	#	if( 'smoothScatter' %in% point )
		smoothScatter( log2(ipdRatio_2str[,1]) , log2(ipdRatio_2str[,2]) , xlab='Single molecule pos strand log2 ipdRatio' , 
						ylab='Single molecule neg strand log2 ipdRatio' ,main='Checking HemiMethylation')
	}
	
	ipdRatio_2str
}











