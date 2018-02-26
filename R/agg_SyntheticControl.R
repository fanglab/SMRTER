
IpdRatioSyntheticControl_Aggregate = function(info , insilico , globalPercentile , coverage_thres=0)
{	
	tmp = lapply( sort(unique(info$tpl)) , function(site) {
		#cat(site,'\n')
		dat = info[ info$tpl == site , ]
		site_ipd = dat$normIpd
		
		
		if( 'modelPrediction' %in% colnames(insilico) ) 
		{
			modelPrediction = insilico[ insilico$tpl == site , 'modelPrediction' ]
		} else {
			modelPrediction = insilico[ insilico$tpl == site , 'controlMean' ]
		}
		
		if( length(site_ipd)>=coverage_thres ) 
		{
			res = PositionIpdRatioSyntheticControl_Aggregate( site_ipd , globalPercentile , modelPrediction , moleculeID=NULL  )
			data.frame( tpl=site , res)
		}  
	} )
	
	do.call(rbind,tmp)
}


PositionIpdRatioSyntheticControl_Aggregate = function( site_ipd , globalPercentile , modelPrediction  ,moleculeID=NULL )
{
	
	site_ipd = cappingSyntheticControl( site_ipd,globalPercentile,modelPrediction  )
	
	tMean = mean(site_ipd,na.rm=TRUE)
	tErr = sd(site_ipd,na.rm=T) / sqrt(length(site_ipd))
	ipdRatio = tMean / modelPrediction
	if(is.na(ipdRatio)) { ipdRatio = 1.0 }
	coverage = length(site_ipd)
	pvalue = computeObservationTstatistic(tMean, tErr , modelPrediction ,coverage)
	score = round(-10.0 * log10(pvalue))
	
	data.frame( tMean=tMean , tErr=tErr , coverage=coverage , modelPrediction=modelPrediction  , ipdRatio=ipdRatio , score=score )
	
}







IpdRatioSyntheticControl_Aggregate2 = function(info , insilico , globalPercentile , coverage_thres=0)
{	
	tmp = lapply( sort(unique(info$tpl)) , function(site) {
		#cat(site,'\n')
		dat = info[ info$tpl == site , ]
		site_ipd = dat$normIpd
		
		if( 'modelPrediction' %in% colnames(insilico) ) 
		{
			modelPrediction = insilico[ insilico$tpl == site , 'modelPrediction' ]
		} else {
			modelPrediction = insilico[ insilico$tpl == site , 'controlMean' ]
		}
		
		if( length(site_ipd)>=coverage_thres ) 
		{
			site_ipd = cappingSyntheticControl( site_ipd,globalPercentile,modelPrediction  )
			res = mean(site_ipd)/modelPrediction
			c( site , res )
		} 
	} )
	
	do.call(rbind,tmp)
}
