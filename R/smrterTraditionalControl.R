#' smrterTraditionalControl
#'
#' smrterTraditionalControl calls methylation fraction at single molecule level with WGA control
#' 
#' @param case_path a character value indicating the file path of cmp.h5 file for native DNA
#' @param cont_path a character value indicating the file path of cmp.h5 file for WGA
#' @param sitesInfo_interest a data frame of smrt data indicating the sites of interest for the single molecule analysis. It should cover the basic information for those sites, including at least "refName", "strand", "tpl", and "modelPrediction"
#' @param split_chunk_num an integer value indicating the number of split chunks for parallel computing
#' @param cores an integer value indicating the number of computer cores for parallel computing
#' @param type a character values indicating the types of analysis. It takes on either "callMod" or "SMipdRatio", where "callMod" represents calling methylation fraction, and "SMipdRatio" represents collecting IPD ratio of single molecule for the following analysis
#' @param coverage_per_molecule an integer value indicating the threshold of number of subreads covering one single site from one molecule in native DNA
#' @param number_of_molecule an integer value indicating the threshold of number of molecules covering one single site in native DNA
#' @param case_coverage an integer value indicating the threshold of number of subreads covering one single site from all molecules in native DNA
#' @param cont_coverage an integer value indicating the threshold of number of subreads covering one single site from all molecules in WGA
#' @param score_thres an integer value indicating the threshold of modification score to filter sites of interest
#'
#' @return the type of "callMod" returns a list with the name of refName and strand "refName.strand". Each item is a data frame comprising the following columns:
#' \itemize{
#'  \item {tpl} an integer value indicating the genomic location 
#'  \item {caseMean} a numeric value indicating mean of normalized case IPDs observed at this position (AGG)
#'  \item {caseStd} a numeric value indicating the standard deviation of case IPDs observed at this position
#'  \item {controlMean} a numeric value indicating the mean of normalized control IPDs observed at this position (AGG)
#'  \item {controlStd} a numeric value indicating the standard deviation of control IPDs observed at this position
#'  \item {coverage} an integer value indicating the mean of case and control coverage
#'  \item {caseCoverage} an integer value indicating count of valid case IPDs at this position
#'  \item {controlCoverage} an integer value indicating count of valid control IPDs at this position
#'  \item {molNum} an integer value indicating the count of valid case molecules
#'  \item {ipdRatio} a numeric value indicating caseMean / controlMean
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
#' 
#' @export
#'
#' @examples
#' 
#' 
smrterTraditionalControl = function( case_path , cont_path , sitesInfo_interest=NULL , split_chunk_num=2000 , cores=20, type = 'callMod' ,
                                     coverage_per_molecule=5 , number_of_molecule=10 , case_coverage=coverage_per_molecule*number_of_molecule , cont_coverage=20 , score_thres=20)
{
  
  registerDoMC( cores )
  cat( "registered cores:" , getDoParWorkers(), "\n" )
  
  #split_chunk_num = 2000
  
  system.time( case <- loadCmpH5( case_path ) )
  system.time( control <- loadCmpH5( cont_path ) )
  
  uni_refName = intersect( unique(as.character( case$index$fullRefName )) , 
                           unique(as.character( control$index$fullRefName )) ) 
  
  num = 0
  result = list()
  
  for( refI in uni_refName )
  {	
    for( strandI in c(0,1) )
    {
      cat( refI ,':  strand: ',strandI, '.....\n' )
      cat( '  case loading ipd .....\n' )
      system.time( case_chunk_info <- fetchChunks( case , refI , strandI ) )
      system.time( case_globalPercentile <- getGlobalPercentile(case_chunk_info$normIpd) )
      
      cat( '  control loading ipd .....\n' )
      system.time( control_chunk_info <- fetchChunks( control , refI , strandI ) )
      system.time( control_globalPercentile <- getGlobalPercentile(control_chunk_info$normIpd) )
      
      if( !is.null(sitesInfo_interest) )
      {
        chunk_region = subset( sitesInfo_interest , refName==refI & strand==strandI )
        case_chunk_info = subset( case_chunk_info , tpl %in% chunk_region$tpl )
        control_chunk_info = subset( control_chunk_info , tpl %in% chunk_region$tpl )
      }
      
      if( nrow(case_chunk_info)>0 )
      {
        system.time( case_split_chunk_info <- split( case_chunk_info , cut( case_chunk_info$tpl , split_chunk_num , labels=F ) ) )
        
        cat( '    calling modifications from chunks:',length(case_split_chunk_info),' .....\n' )
        num = num+1
        system.time( result[[num]] <- foreach( i = 1:length(case_split_chunk_info) , .inorder=FALSE ) %dopar% {
          
          cat(i,'\n')
          case_info = case_split_chunk_info[[i]]
          region = range( case_info$tpl ) 
          control_info = subset( control_chunk_info, tpl>=region[1] & tpl<=region[2] ) 
          case_info = case_info[ case_info$tpl %in% control_info$tpl & !is.na(case_info$normIpd) , ]
          control_info = control_info[ control_info$tpl %in% case_info$tpl & !is.na(control_info$normIpd) , ]
          
          if( type == 'callMod'    )   res <- MultiPositionTraditionalControl(case_info , control_info , case_globalPercentile , control_globalPercentile , coverage_per_molecule , number_of_molecule , case_coverage , cont_coverage , score_thres)
          if( type == 'SMipdRatio' )   res <- IpdRatioTraditionalControl_SM(case_info , control_info , case_globalPercentile , control_globalPercentile , coverage_per_molecule , number_of_molecule , case_coverage , cont_coverage )
          if( type == 'AGGipdRatio' )  res <- IpdRatioTraditionalControl_Aggregate(case_info , control_info , case_globalPercentile , control_globalPercentile , case_coverage , cont_coverage )
          
          res
          
        } )
        names(result)[num] = paste( refI,strandI,sep='.'  )
      } # nrow(case_chunk_info)
      
    } # strandI
  } # refI
  
  lapply( result , function(x) do.call(rbind,x) )
  
}

