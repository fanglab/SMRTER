#' drawMotifHeatmap
#'
#' drawMotifHeatmap plots the heatmap of motif enrichment scores
#' 
#' @param motif_fdr a data frame comprising the follow components
#' \itemize{
#'  \item {motif} a character value indicating the methylation motifs
#'  \item {pos} an integer indicating the methylation position in the motif
#'  \item {score} a numeric value indicating the motif enrichment score for the corresponding motif
#' }
#' @param outputName a character value indicating the output pdf file name
#'
#' @return
#' @export
#'
#' @examples
#' 
drawMotifHeatmap = function( motif_fdr , outputName )
{
  
  library(ggplot2)
  library(RColorBrewer)
  
  motif_fdr2 = subset(motif_fdr , pos%in%c(2,3))
  
  centered_motifs = sapply( 1:nrow(motif_fdr2)  , function(i) {
    if(motif_fdr2[i,2]==2) m = paste0( "centered at [" , substr(motif_fdr2[i,1],2,2) , "]" , substr(motif_fdr2[i,1],3,3) )
    if(motif_fdr2[i,2]==3) m = paste0( "centered at "  , substr(motif_fdr2[i,1],2,2) , "[" , substr(motif_fdr2[i,1],3,3) , "]" )
    m
  } )
  
  first_base = substr(motif_fdr2[,1],1,1)
  fourth_base = substr(motif_fdr2[,1],4,4) 
  score_mod = -log10(motif_fdr2[,3])
  motif_results = data.frame( centered_motifs, first_base, fourth_base, score_mod=-log10(motif_fdr2[,3]) )
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  #myPalette <- colorRampPalette(c(  "green" , 'yellow' , "red" ))
  
  pdf( outputName )
  print( ggplot(motif_results) +
    geom_tile(aes(y=first_base, x=fourth_base, fill=score_mod)) +
    facet_wrap(~centered_motifs) +
    scale_fill_gradientn(name="Mod Score", colours=myPalette(100)) +
    labs(title="Methylation score for the 4-mer motifs") +
    labs(x="The 4th base in the 4-mer motif", y="The 1st base in the 4-mer motif") +
    theme_bw() )
  dev.off()
  
}






