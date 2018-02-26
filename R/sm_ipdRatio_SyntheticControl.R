
IpdRatioSyntheticControl_SM = function(info , insilico , globalPercentile , 
										coverage_per_molecule , number_of_molecule , coverage_thres)
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
			res = IpdRatioPositionSyntheticControl_SM( site_ipd , globalPercentile , modelPrediction , moleculeID , coverage_per_molecule )
			rep_site = rep(site,length(res))
			data.frame( tpl=rep_site , ipdRatio=res )
		} 
	} )
	
	do.call(rbind,tmp)
}

IpdRatioPositionSyntheticControl_SM = function(site_ipd , globalPercentile , modelPrediction  ,moleculeID , coverage_per_molecule )
{
	site_ipd = cappingSyntheticControl( site_ipd,globalPercentile,modelPrediction  )
	
	smol_ipd = tapply(site_ipd,moleculeID,mean)
	smol_ipdRatio = smol_ipd/modelPrediction
	smol_len = tapply(moleculeID,moleculeID,length)
	smol_ipdRatio = smol_ipdRatio[ smol_len>=coverage_per_molecule & !is.na(smol_ipdRatio) & !is.infinite(smol_ipdRatio)  ]
	
	smol_ipdRatio
}

###############################################################################################################################



