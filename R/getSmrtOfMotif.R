#' getSmrtOfMotif
#'
#' getSmrtOfMotif extracts the rows in a smrt data for the methylation position in a given methylation motif
#' 
#' @param smrt a data frame of smrt data comprising at least refName, strand and tpl
#' @param motif a character value indicating the methylation motif
#' @param methy_pos an integer value indicating the methylated position in the methylation motif
#' @param genome a DNAStringSet of the genome reference
#'
#' @return a data frame of smrt data comprising of the methylation positions
#' @export
#'
#' @examples
#' 
getSmrtOfMotif = function( smrt , motif ,  methy_pos , genome )
{
  cat('searching for ',motif,' at position ' ,methy_pos,' ...\n')
  
  # if( is.null(genome) ) genome = readDNAStringSet(ref)	
  
  motifF = motif
  motifB = as.character(reverseComplement(DNAString(motif)))
  positionF = vmatchPattern(  motifF , genome , fixed=FALSE)
  positionB = vmatchPattern(  motifB , genome , fixed=FALSE)
  
  cat('splitting data ...\n')
  split_smrt = splitSmrt(smrt,sorted=F)
  
  smrt_chr = list()
  chr_list = unique(gsub( ".[01]$" , "" , names(split_smrt) ))
  for( chr in chr_list )
  {
    cat( "  ->" , chr , "\n" )
    posF_chr = start(positionF[[chr]]) + methy_pos - 1
    posB_chr = start(positionB[[chr]]) + nchar(motif) - methy_pos
    smrt_chr_0 = subset( split_smrt[[ paste0(chr,'.0') ]] , tpl %in% posF_chr )
    smrt_chr_1 = subset( split_smrt[[ paste0(chr,'.1') ]] , tpl %in% posB_chr )
    smrt_chr[[chr]] = rbind(smrt_chr_0 , smrt_chr_1)
  }
  
  cat('combining results ...\n')
  do.call(rbind,smrt_chr)
  
}


getSmrtOfMotif2 = function( smrt , motif ,  methy_pos , genome=NULL , ref=NULL , non_motif_sites = F   )
{
   
  if( is.null(genome) ) genome = readDNAStringSet(ref)	
  
  baseOfInterest = substr(motif,methy_pos,methy_pos)
  cat('extracting sites at ',baseOfInterest,' ...\n')
  smrt = subset(smrt,base==baseOfInterest )
  
  motifF = motif
  motifB = as.character(reverseComplement(DNAString(motif)))
  positionF = vmatchPattern(  motifF , genome , fixed=FALSE)
  positionB = vmatchPattern(  motifB , genome , fixed=FALSE)
  
  cat('splitting data ...\n')
  split_smrt = splitSmrt(smrt,sorted=F)
  
  cat('searching for ',motif,' at position ' ,methy_pos,' ...\n')
  smrt_chr = list()
  non_smrt_chr = list()
  for( chr in names(split_smrt) )
  {
    cat( "  ->" , chr , "\n" )
    posF_chr = start(positionF[[chr]]) + methy_pos - 1
    posB_chr = start(positionB[[chr]]) + nchar(motif) - methy_pos
    smrt_chr_0 = subset( split_smrt[[chr]][['0']] , tpl %in% posF_chr )
    smrt_chr_1 = subset( split_smrt[[chr]][['1']] , tpl %in% posB_chr )
    smrt_chr[[chr]] = rbind(smrt_chr_0 , smrt_chr_1)
    
    if( non_motif_sites )
    {
      non_smrt_chr_0 = subset( split_smrt[[chr]][['0']] , ! tpl %in% posF_chr )
      non_smrt_chr_1 = subset( split_smrt[[chr]][['1']] , ! tpl %in% posB_chr )
      non_smrt_chr[[chr]] = rbind( non_smrt_chr_0 , non_smrt_chr_1)
    }
  }
  
  cat('combining results ...\n')
  smrt_motif = do.call(rbind,smrt_chr)
  smrt_non_motif = do.call(rbind,non_smrt_chr)
  
  list(smrt_motif=smrt_motif,smrt_non_motif=smrt_non_motif)
  
}



