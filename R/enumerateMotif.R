#' enumerateMotif
#'
#' enumerateMotif enumerates all k-mer motifs
#' 
#' @param kmer an integer value indicating the motif length
#' @param meth_char a character value indicating the methylated position in the motif
#'
#' @return a data frame comprising the following components
#' \itemize{
#'  \item {motif} a character value indicating the methylation motifs
#'  \item {pos} an integer indicating the methylation position in the motif
#' }
#' 
#' @export
#'
#' @examples
#' 
enumerateMotif = function( kmer , meth_char )
{	
  alphabet = c('G','A','T','C')
  motifs = alphabet
  for(i in 1:(kmer-1) )
    motifs = c( sapply( motifs, function(x) { paste(x,alphabet,sep='') } ) )
  
  motifs = motifs[ grep( meth_char , motifs) ]
  pos = sapply( motifs , function(str) { grep(meth_char,strsplit(str,'')[[1]]) } )
  pos_num = sapply(pos,length)
  data.frame( motif = rep(motifs,pos_num) , pos = do.call(c,pos) , stringsAsFactors = FALSE )
  
}

