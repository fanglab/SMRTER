#' overlapWithBed 
#' 
#' overlapWithBed checks the intersection of DNA modifications between SMRTseq and other detection methods
#'
#' @param smrt a data frame of SMRT smrt data, indicating the putative DNA modifications detected by SMRT
#' @param bed a data frame comprising at least chr, start, end, indicating the regions of DNA modifications detected by other techniques, e.g., DIPseq
#' @param genome_size an integer value indicating the size of genome
#' @param coverage_thres a vector of integer values, indicating the list of thresholds for read coverage
#' @param score_thres a vector of integer values, indicating the list of thresholds for modification detection scores
#' @param meth_base a character value, indicating the methylated base, by default, "A"
#' @param tags  a vector of character values, indicating the types of modifications. If NULL, overlapWithBed() checks all types of modifications for that base indicated by "meth_base"
#' @param cores an integer value indicating the number of computer cores for parallel computing
#'
#' @return a data frame comprising the following columns
#' \itemize{
#'  \item {total_smrt} an integer value indicating the number of putative modifications by SMRTseq
#'  \item {overlap_smrt} an integer value indicating the number of SMRTseq-detected modifications which are also supported by other techniques
#'  \item {observed_ratio} a numeric value indicating the observed fraction of SMRTseq-detected modifications supported by other techniques 
#'  \item {expected_ratio} a numeric value indicating the expected fraction of SMRTseq-detected modifications supported by other techniques 
#'  \item {fold_enrich} the odds ratio beteween observed_ratio and expected_ratio
#'  \item {pval} the p value for the significance of overlap between SMRTseq and other techniques compared with the random chance
#' }
#' 
#' 
#' @export
#'
#' @examples
#' 
#' @seealso \code{\link{readSmrtGff}}; \code{\link{readSmrtCsv}};
#' 
overlapWithBed = function( smrt , bed , chr='chr' , start='start' , end='end' , genome_size , coverage_thres=c(10,20,21,30) , score_thres=c(20,30,31,40) , meth_base='A' , tags = NULL , cores=10 )
{
  
  library(foreach)
  library(doMC)
  registerDoMC(cores)
  
  overlap = function( smrt_list , bed_list , coverageT , scoreT , expected_ratio )
  {
    
    smrt_of_interest = lapply( smrt_list , function(x) 
      subset( x , score>=scoreT[1] & score<=scoreT[2] & coverage>=coverageT[1] & coverage<=coverageT[2] ) )
    
    smrt_overlap = foreach( chr = names(bed_list) ) %dopar%
    {
      smrt_chr = smrt_of_interest[[chr]]
      bed_chr = bed_list[[chr]]
      do.call( c , lapply( smrt_chr$tpl , function(x) any( x>=bed_chr[,start] & x<=bed_chr[,end] ) ))
    }
    smrt_overlap = do.call( c , smrt_overlap )
    
    if(!is.null(smrt_overlap))
    {
      total_smrt = sum( sapply(smrt_of_interest,nrow) )
      overlap_smrt = sum( smrt_overlap )
      observed_ratio = mean( smrt_overlap )
      pval = binom.test( overlap_smrt , total_smrt , expected_ratio )$p.val
      fold_enrich = (overlap_smrt/total_smrt) / expected_ratio
      data.frame( total_smrt , overlap_smrt , observed_ratio , expected_ratio , fold_enrich , pval )
    } else{
      data.frame( total_smrt=0 , overlap_smrt=0 , observed_ratio=0 , expected_ratio , fold_enrich=0 , pval=NA )
    }
    
  }

  if( is.null(tags) )
  {
    smrt = subset( smrt , base==meth_base )
  } else {
    smrt = subset( smrt , base==meth_base & tag%in%tags )
  }
  smrt_list = split( smrt , as.character(smrt$refName) )
  bed_list = split( bed , as.character(bed[,chr]) )
  
  expected_ratio = sum(bed[,end]-bed[,start])/genome_size
  
  coverage_thres = matrix(coverage_thres,ncol=2,byrow=TRUE)
  score_thres = matrix(score_thres,ncol=2,byrow=TRUE)
  
  res = list()
  num = 0
  for( i in 1:nrow(coverage_thres) )
  {
    for( j in 1:nrow(score_thres) )
    {
      ci = coverage_thres[i,]
      sj = score_thres[j,]
      num=num+1
      thres = paste0('coverage in [', ci[1] ,",", ci[2] ,'] and score in [', sj[1] ,",", sj[2] ,']')
      cat(thres,'\n')
      mar_res = overlap( smrt_list , bed_list , ci , sj , expected_ratio )
      res[[num]] = data.frame( thres , mar_res )
    }
  }
  
  do.call(rbind,res)
  
}
