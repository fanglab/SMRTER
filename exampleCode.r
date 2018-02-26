
	
library( Biostrings )
library( foreach )
library( doMC )
library( h5r , lib.loc='/hpc/users/zhus02/schzrnas/sjzhu/tmp/install_tmp')
library( pbh5 , lib.loc='/hpc/users/zhus02/schzrnas/sjzhu/tmp/install_tmp')
library( SMRTER , lib.loc='/hpc/users/zhus02/schzrnas/sjzhu/tmp/install_tmp')

####################################################################################################
######################  FDR calling                          ######################################
####################################################################################################
setwd("/hpc/users/zhus02/Perspective/data/Example_Algae_scaffold18/")

nat <- readSmrtCsv( "nat_modifications.csv" )
wga <- readSmrtCsv( "wga_modifications.csv"  )
natA <- subset( nat , base=='A' )
wgaA <- subset( wga , base=='A' ) 

scaffold_18 = readDNAStringSet( "Algae_scaffold_18.fasta" )
natVATB = getSmrtOfMotif( smrt=natA , motif="VATB" , methy_pos=2 , genome=scaffold_18 )

fdrAllA = callFDR( natA$ipdRatio , wgaA$ipdRatio , thresList=seq(1,5,0.1) )
fdrVATB = callFDR( natVATB$ipdRatio , wgaA$ipdRatio , thresList=seq(1,5,0.1)  )

####################################################################################################
######################   density plot for both native and WGA  ####################################
####################################################################################################
setwd("/sc/orga/projects/fangg03a/sjzhu/Projects/SMRT/Perspective/data/Example_Algae_scaffold18/test")
pdf("DensityPlot.pdf")
par(mfrow=c(2,2))

plot(density(natA$ipdRatio),col='red',lwd=1.5,xlim=c(0,8),main="Native As vs WGA As",ylab="Frequency",xlab="IPD ratio")
lines(density(wgaA$ipdRatio),col='black',lwd=1.5)
legend("topright",legend=c("Native As","WGA As"),col=c('red','black'),lty=1,lwd=1.5)

plot(density(natA$ipdRatio),col='red',lwd=1.5,xlim=c(0,8),ylim=c(0,0.01),main="Native As vs WGA As (Zoom in)",ylab="Frequency",xlab="IPD ratio")
lines(density(wgaA$ipdRatio),col='black',lwd=1.5)
abline(v=4.4,col="blue",lwd=1.5,lty=2)

plot(density(natVATB$ipdRatio),col='red',lwd=1.5,xlim=c(0,8),main="Native VATBs vs WGA As",ylab="Frequency",xlab="IPD ratio")
lines(density(wgaA$ipdRatio),col='black',lwd=1.5)
legend("topright",legend=c("Native VATBs","WGA As"),col=c('red','black'),lty=1,lwd=1.5)

plot(density(natVATB$ipdRatio),col='red',lwd=1.5,xlim=c(0,8),ylim=c(0,0.05),main="Native VATBs vs WGA As (Zoom in)",ylab="Frequency",xlab="IPD ratio")
lines(density(wgaA$ipdRatio),col='black',lwd=1.5)
abline(v=4.0,col="green",lwd=1.5,lty=2)

dev.off()



####################################################################################################
######################   call motif from aggregate level      ######################################
####################################################################################################
motifFdr = callFDR_eachMotif( smrt=natA , genome=scaffold_18 , kmer=4 , meth_char = 'A' , cores=10 , thres=4 ) #, outputFolder )
setwd("/sc/orga/projects/fangg03a/sjzhu/Projects/SMRT/Perspective/data/Example_Algae_scaffold18/test")
drawMotifHeatmap( motifFdr , outputName="motif_4mer.pdf" )

setwd("/sc/orga/projects/fangg03a/sjzhu/Projects/SMRT/Perspective/data/Example_Algae_scaffold18/test")
motifFdr = callMotif(  nat=natA , wga=wgaA , genome=scaffold_18 , kmer=4 , meth_char='A' , cores=10 , criterion="ipdRatio" , thres=4 ) 


####################################################################################################
######################   call fraction from single molecule level   ################################
####################################################################################################
fracInsilico <- smrterSyntheticControl( path="nat_aligned_reads.cmp.h5" , sitesInfo_interest=natVATB , split_chunk_num=2000 , cores=20 , type = 'callMod' , 
							coverage_per_molecule=3 , number_of_molecule=5 , coverage_thres=15 , score_thres=0 )
resInsilico <- do.call( rbind , fracInsilico )
quantile( subset(resInsilico,ipdRatio>4.5)$frac , seq(0,1,0.1) , na.rm=T )	


####################################################################################################
######################   call fraction from single molecule level   ################################
####################################################################################################
fracTraditional <- smrterTraditionalControl( case_path="nat_aligned_reads.cmp.h5" , cont_path="wga_aligned_reads.cmp.h5" , 
                sitesInfo_interest=natVATB , split_chunk_num=2000 , cores=20 , type = 'callMod' , 
                coverage_per_molecule=3 , number_of_molecule=5 , case_coverage=15 , cont_coverage=10 , score_thres=0 )
resTraditional <- do.call( rbind , fracTraditional )
quantile( subset(resTraditional,ipdRatio>4.5)$frac , seq(0,1,0.1) , na.rm=T )	


####################################################################################################
######################   check hemi-methylation                     ################################
####################################################################################################
SMipdRatio <- smrterSyntheticControl( path="nat_aligned_reads.cmp.h5" , sitesInfo_interest= subset(natVATB,ipdRatio>4.5) , split_chunk_num=2000 , cores=20 , type = 'SMipdRatio' , 
							coverage_per_molecule=3 , number_of_molecule=5 , coverage_thres=15 , score_thres=0 )
pdf("hemiMethylation.pdf")
ipdRatio.2str <- hemiMethylation(SMipdRatio,draw=TRUE)
dev.off()

####################################################################################################
#########   calling fraction at SM by fixing two GMMs with means of 1 and 4       ##################
####################################################################################################
fracSm <- tapply( SMipdRatio[[1]]$ipdRatio , SMipdRatio[[1]]$tpl , function(ipdRsm) 
{
  if(length(ipdRsm)>=5) GMM2( ipdRsm , c(1,4) , fixed=c(1,2) )$lamda[2]
})
fracSm = fracSm[ sapply(fracSm,length)>0 ]
fracFixed1n4 <- data.frame( tpl=as.integer(names(fracSm)) , frac=do.call(c,fracSm) )


####################################################################################################
######################   call motif from single molecule level   ###################################
####################################################################################################
SMipdRatio2 <- smrterSyntheticControl( path="nat_aligned_reads.cmp.h5" , sitesInfo_interest= natA , split_chunk_num=2000 , cores=20 , type = 'SMipdRatio' , 
                                      coverage_per_molecule=5 , number_of_molecule=1 , coverage_thres=3 , score_thres=0 )

motifFdrSM = callFDR_eachMotif( split_smrt=SMipdRatio2 , genome=scaffold_18 , kmer=4 , meth_char = 'A' , cores=10 , thres=4 ) #, outputFolder )


####################################################################################################
################################   call consensus   ################################################
####################################################################################################
setwd("/sc/orga/projects/fangg03a/sjzhu/Projects/SMRT/Perspective/data/Example_Algae_chr1")
chromosome_1 <- readDNAStringSet("Algae_chromosome_1.fasta") 
gene <- read.table( 'gene_annot_chromosome_1.txt' , sep='\t') 
dip <- read.table('m6dA_DIPseq_chromosome_1.txt',sep='\t') 
natChr1VATB  <- readSmrtCsv( 'nat_modifications_chromosome_1.csv' ) 
methVATB <- subset( natChr1VATB , ipdRatio>4.5 ) 
consensus = getConsensus( smrt=natChr1VATB , window=2000 , annot=gene , chr='V1' , strand='V7' , start='V4' , criterion="ipdRatio" , smooth=10 , cores=10 )
setwd("/sc/orga/projects/fangg03a/sjzhu/Projects/SMRT/Perspective/data/Example_Algae_chr1/test")
pdf('consensus.pdf',h=5,w=10)
plot(c(-2000:2000),consensus,type='l',main='Consensus analysis of V[A]TB at gene TSS',xlab='Distance to TSS (bp)',ylab='Normalized IPD ratio of V[A]TB')
dev.off()

####################################################################################################
################################   overlap between SMRTseq and DIP   ###############################
####################################################################################################
overlap = overlapWithBed( smrt=methVATB , bed=dip , chr='V1' , start='V2' , end='V3' , 
                genome_size=nchar(chromosome_1) , coverage_thres=c(0,10,11,20,21,30) , score_thres=c(0,20,21,40,41,100) , 
                meth_base='A' , tags = NULL , cores=10 )
subset(overlap,total_smrt>10 & !is.na(pval))  	




