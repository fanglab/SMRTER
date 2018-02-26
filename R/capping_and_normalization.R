
getGlobalPercentile = function(normIpd)
{
	ifelse( length(normIpd)<10 , 5.0 , quantile( normIpd , 0.99 ,na.rm=TRUE ) )
}
	
cappingSyntheticControl  = function( site_ipd,globalPercentile,modelPrediction  )
{
	percentile = min(0.9, (1-1/(length(site_ipd)-1)))
	localPercentile = quantile(site_ipd, percentile,na.rm=TRUE)
	capValue = max(globalPercentile, 4.0 * modelPrediction, localPercentile) # Capping 3
	#capValue = max(globalPercentile, localPercentile) # Capping 3.1; we can use 4.0 * modelPrediction later on, 
	site_ipd[site_ipd>capValue] = capValue
	site_ipd
}

cappingTraditionalControl = function( case_site_ipd , control_site_ipd , case_globalPercentile , control_globalPercentile )
{
	percentile = min(0.9, (1-1/(length(case_site_ipd)-1)))
	localPercentile = quantile(case_site_ipd, percentile,na.rm=TRUE)
	capValue = max(case_globalPercentile , 4.0 * median(case_site_ipd) , localPercentile) # Capping 3
	case_site_ipd[case_site_ipd>capValue] = capValue
	
	percentile = min(0.9, (1-1/(length(control_site_ipd)-1)))
	localPercentile = quantile(control_site_ipd, percentile,na.rm=TRUE)
	capValue = max(control_globalPercentile , 4.0 * median(control_site_ipd) , localPercentile) # Capping 3
	control_site_ipd[control_site_ipd>capValue] = capValue
	
	list( case = case_site_ipd , control = control_site_ipd  )
	
}

subreadNormalizationFactor = function(subread_ipd)
{
	capValue = min( 10 , quantile(subread_ipd,0.99,na.rm=TRUE) ); 
	subread_ipd[subread_ipd>capValue] = capValue; 
	if(length(subread_ipd)<2) { normalization = 0.1 } else {	normalization = mean(subread_ipd,na.rm=T) }
	normalization
}
