#' getConsensus 
#' 
#' getConsensus performs the DNA modification consensus analysis along the biological elements of interest
#'
#' @param smrt a data frame of smrt data, taking at least columns: refName, strand, tpl, ipdRatio/score
#' @param annot a data frame of annotation giving the genomic locations of the biological elements of interest
#' @param window an integer value indicating the size of flanking region (bp)
#' @param chr a character value indicating the name of the column representing the "chromosome" in the annotation
#' @param strand a character value indicating the name of the column representing the "strand" in the annotation
#' @param start a character value indicating the name of the column representing the "start" in the annotation
#' @param criterion a character value taking on either "ipdRatio" or "score", by default, "ipdRatio"
#' @param smooth an integer value indicating the size of sliding window for smoothing 
#' @param cores an integer value indicating the number of computer cores for parallel computing
#'
#' @return a vector of numeric value (2*window+1 long), representing the consensus of the criterion (by default, ipdRatio) around the start sites of the biological elements
#' @export
#'
#' @examples
#' 
#' @seealso \code{\link{readSmrtGff}}; \code{\link{readSmrtCsv}};
#' 
getConsensus =	function( smrt , annot , window=1000 , chr='chr' , strand='strand' , start='start', criterion="ipdRatio" , smooth=0 , cores=1 )
{	
  
  library(foreach)
  library(doMC)
  registerDoMC(cores)
  
  getVal = function( smrtI , startI , strandI )
  {
    f = function( starti , strandi )
    {	
      range1 = which( abs( smrtI$tpl - starti ) <= window )
      val = rep(0, 2*window+1 )
      if( strandi == '+')
      {	
        range2 = smrtI$tpl[range1] - starti + window + 1 
      } else {
        range2 = starti - smrtI$tpl[range1] + window + 1 
      }
      val[range2] = smrtI[ range1 , criterion ]
      val
    }
    t( mapply( f , startI , strandI ) )	
  }
  
  uni_chr = intersect( unique(as.character(smrt$refName)) , unique(as.character(annot[,chr])) )
  
  all_val = foreach( chrI = uni_chr ) %dopar% 
  {
    cat(chrI,'\n')
    smrtI    = subset(smrt , refName == chrI )
    annotI = annot[ as.character(annot[,chr])==chrI , ]
    strandI = as.character( annotI[,strand] )
    startI = annotI[ , start ]
    if( length(startI)>0 ) 
      getVal( smrtI , startI , strandI )
  }
  
  names(all_val) = uni_chr
  all_val = all_val[order(names(all_val))]
  consensus_matrix = do.call(rbind,all_val)
  
  annot = annot[ as.character(annot[,chr]) %in% uni_chr , ]
  split_annot = split( annot , as.character(annot[,chr])  )
  split_annot = split_annot[order(names(split_annot))]
  annot = do.call(rbind,split_annot)
  
  res = list( consensus_matrix = consensus_matrix , annot = annot  )
  
  consensus = apply( consensus_matrix , 2 , function(x) ifelse(sum(x>0)==0,0,mean(x[x>0])) )
  
  if(smooth>0)
  {
    tmp = consensus
    l = round( (smooth-1)/2 )
    for(i in 1:length(consensus))
    {
      s = max( 1 , i-l  )
      e = min( length(consensus) , i+l )
      tmp[i] = mean( consensus[s:e] ) 
    }
    consensus = tmp
  }
  
  consensus
  
}

