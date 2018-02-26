
MultiPositionSyntheticControl = function(info,insilico , globalPercentile , coverage_per_molecule , number_of_molecule , coverage_thres , score_thres)
{	
	tmp = lapply( sort(unique(info$tpl)) , function(site) {
		#cat(site,'\n')
		dat = info[ info$tpl == site , ]
		site_ipd = dat$normIpd
		moleculeID = dat$moleculeID
		
		if( 'modelPrediction' %in% colnames(insilico) ) 
		{
			modelPrediction = insilico[ insilico$tpl == site , 'modelPrediction' ]
		} else {
			modelPrediction = insilico[ insilico$tpl == site , 'controlMean' ]
		}
		
		if( length(site_ipd)>=coverage_thres ) 
		{
			res = computePositionSyntheticControl( site_ipd , globalPercentile , modelPrediction , moleculeID , 
					coverage_per_molecule , number_of_molecule , coverage_thres , score_thres)
			data.frame( tpl=site , res)
		} 
	} )
	
	do.call(rbind,tmp)
}


computePositionSyntheticControl = function( site_ipd , globalPercentile , modelPrediction  ,moleculeID=NULL ,
									coverage_per_molecule , number_of_molecule , coverage_thres , score_thres)
{
	
	site_ipd = cappingSyntheticControl( site_ipd,globalPercentile,modelPrediction  )
	
	tMean = mean(site_ipd,na.rm=TRUE)
	tErr = sd(site_ipd,na.rm=T) / sqrt(length(site_ipd))
	ipdRatio = tMean / modelPrediction
	if(is.na(ipdRatio)) { ipdRatio = 1.0 }
	coverage = length(site_ipd)
	pvalue = computeObservationTstatistic(tMean, tErr , modelPrediction ,coverage)
	score = round(-10.0 * log10(pvalue))
	
	frac = NA
	ipdRatio2 = NA
	molNum = NA
	
	if( !is.null(moleculeID) & score >= score_thres & length(site_ipd)>=coverage_thres )
	{  
		smol_ipd = tapply(site_ipd,moleculeID,mean)
		smol_len = tapply(moleculeID,moleculeID,length)
		smol_ipd = smol_ipd[ smol_len>=coverage_per_molecule & !is.na(smol_ipd) & !is.infinite(smol_ipd)  ]
		#smol_ipd = smol_ipd[ !is.na(smol_ipd) & !is.infinite(smol_ipd)  ]
		molNum = length(smol_ipd)
		if( molNum >= number_of_molecule )
		{
			gmm = GMM( smol_ipd ,modelPrediction)
			frac = gmm$lamda[2]
			ipdRatio2 = gmm$mu[2]/modelPrediction
		} 
	} 
	
	if( is.null(moleculeID) & score >= score_thres & length(site_ipd)>=coverage_thres ) 
	{
		frac = GMM( site_ipd ,modelPrediction)$lamda[2]
	}
	
	data.frame( tMean=tMean , tErr=tErr , coverage=coverage , molNum=molNum , modelPrediction=modelPrediction  , ipdRatio=ipdRatio , score=score , ipdRatio2=ipdRatio2, frac=frac )
	
}

