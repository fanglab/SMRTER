#' smrterSyntheticControl
#'
#' smrterSyntheticControl calls methylation fraction at single molecule level with in silico control
#' 
#' @param path a character value indicating the file path of .cmp.h5 file
#' @param sitesInfo_interest a data frame of smrt data indicating the sites of interest for the single molecule analysis. It should cover the basic information for those sites, including at least "refName", "strand", "tpl", and "modelPrediction"
#' @param split_chunk_num an integer value indicating the number of split chunks for parallel computing
#' @param cores an integer value indicating the number of computer cores for parallel computing
#' @param type a character values indicating the types of analysis. It takes on either "callMod" or "SMipdRatio", where "callMod" represents calling methylation fraction, and "SMipdRatio" represents collecting IPD ratio of single molecule for the following analysis
#' @param coverage_per_molecule an integer value indicating the threshold of number of subreads covering one single site from one molecule 
#' @param number_of_molecule an integer value indicating the threshold of number of molecules covering one single site 
#' @param coverage_thres an integer value indicating the threshold of number of subreads covering one single site from all molecules
#' @param score_thres an integer value indicating the threshold of modification score to filter sites of interest
#'
#' @return the type of "callMod" returns a list with the name of refName and strand "refName.strand". Each item is a data frame comprising the following columns:
#' \itemize{
#'  \item {tpl} an integer value indicating the genomic location 
#'  \item {tMean} a numeric value indicating the capped mean of normalized IPDs observed at this position (AGG)
#'  \item {tErr} a numeric value indicating the capped standard error of normalized IPDs observed at this position: standard deviation/sqrt(coverage)
#'  \item {coverage} an integer value indicating the count of valid IPDs 
#'  \item {molNum} an integer value indicating the count of valid molecules
#'  \item {modelPrediction} a numeric value indicating the normalized mean IPD for the in-silico control mode
#'  \item {ipdRatio} a numeric value indicating tMean/modelPrediction
#'  \item {score} an integer value indicating the Phred-transformed pvalue of modification detection, i.e. 10*(-log10(pvalue))
#'  \item {ipdRatio2} a numeric value indicating the mean ipd Ratio of modified molecules at a single site. 
#'  Two Gaussian mixture models are learned for unmodified and modified distributions, 
#'  where the mean for the unmodified distribution is fixed as 1, 
#'  and the mean for the modified is learned (i.e., ipdRatio2)
#'  \item {frac} a numeric value indicating the fraction of modified molecules at a single site estimated from the single molecule level
#' }
#' the type of "SMipdRatio" returns a list with the name of refName and strand "refName.strand". Each item is a data frame comprising the following columns:
#' \itemize{
#'  \item {tpl} an integer value indicating the genomic location 
#'  \item {ipdRatio} a numeric value indicating the IPD ratio at this position of one specific molecule
#' }
#' 
#' @export
#'
#' @examples
#' 
smrterSyntheticControl = function( path , sitesInfo_interest , split_chunk_num=2000 , cores=20 , type = 'callMod' ,
                                   coverage_per_molecule=5 , number_of_molecule=10 , coverage_thres=coverage_per_molecule*number_of_molecule , score_thres=20)
{
  # sitesInfo_interest is a smrtportal csv file with refname, strand, tpl, modelPrediction
  
  registerDoMC( cores )
  cat( "registered cores:" , getDoParWorkers(), "\n" )
  
  system.time( native <- loadCmpH5( path ) )
  
  num = 0
  result = list()
  uni_refName = intersect( unique(as.character( native$index$fullRefName )) , unique(as.character(sitesInfo_interest$refName)) )
  
  for( refI in uni_refName )
  {
    for( strandI in c(0,1) )
    {
      cat( refI ,':  strand: ',strandI, '.....\n' )
      cat( '    loading ipd .....\n' )
      
      system.time( chunk_info       <-   fetchChunks( native , refI , strandI ) )
      system.time( globalPercentile <-   getGlobalPercentile(chunk_info$normIpd) )
      system.time( chunk_insilico   <-   subset( sitesInfo_interest, refName==refI & strand==strandI ) )
      system.time( chunk_info       <-   subset( chunk_info , tpl %in% chunk_insilico$tpl & !is.na(normIpd) ) )
      
      if( nrow(chunk_info)>0 )
      {
        system.time( split_chunk_info <-   split( chunk_info , cut( chunk_info$tpl , split_chunk_num , labels=F ) ) )
        
        cat( '    calling modifications from chunks:',length(split_chunk_info),' .....\n' )
        num = num+1
        system.time( result[[num]] <- foreach( i = 1:length(split_chunk_info) , .inorder=FALSE ) %dopar% {
          cat(i,'\n')
          info      = split_chunk_info[[i]]
          region    = range( info$tpl ) 
          insilico  = subset( chunk_insilico, tpl>=region[1] & tpl<=region[2] ) 
          info      = info[ info$tpl %in% insilico$tpl & !is.na(info$normIpd) , ]
          
          if( type == 'callMod'    )   res <- MultiPositionSyntheticControl( info , insilico , globalPercentile ,coverage_per_molecule , number_of_molecule , coverage_thres , score_thres )
          if( type == 'SMipdRatio' )   res <- IpdRatioSyntheticControl_SM( info , insilico , globalPercentile ,coverage_per_molecule , number_of_molecule , coverage_thres )
          if( type == 'AGGipdRatio' )  res <- IpdRatioSyntheticControl_Aggregate( info , insilico , globalPercentile ,coverage_thres )
          
          res
          
        } )
        names(result)[num] = paste( refI,strandI,sep='.'  )
      } # nrow(chunk_info)
      
    } # strandI
  } # refI
  
  lapply( result , function(x) do.call(rbind,x) )
  
}


