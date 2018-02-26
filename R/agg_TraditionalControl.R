
IpdRatioTraditionalControl_Aggregate = function(case_info  ,control_info , case_globalPercentile , control_globalPercentile , case_coverage , cont_coverage )
{
	uni_tpl = intersect( unique( case_info$tpl ) , unique( control_info$tpl ) )
	tmp = lapply( sort(uni_tpl) , function(site) {
		#cat(site,'\n')
		case_site_info = case_info[ case_info$tpl == site , ]
		case_site_ipd = case_site_info$normIpd
		
		control_site_info = control_info[ control_info$tpl == site , ]
		control_site_ipd = control_site_info$normIpd
		
		if( length(case_site_ipd)>=case_coverage & length(control_site_ipd)>=cont_coverage ) 
		{
			res = PositionIpdRatioTraditionalControl_Aggregate( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile )
			data.frame( tpl=site , res)
		} 
	} )
	
	do.call(rbind,tmp)
}

PositionIpdRatioTraditionalControl_Aggregate= function( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile )
{
	
	capping_site_ipd = cappingTraditionalControl( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile)
	case_site_ipd = capping_site_ipd[['case']]
	control_site_ipd = capping_site_ipd[['control']]
	
	coverage 		= 	round( (length(case_site_ipd)+length(control_site_ipd))/2 )
	caseCoverage    = 	length(case_site_ipd)
	controlCoverage = 	length(control_site_ipd)
	caseMean		= 	mean(case_site_ipd)
	caseStd 		= 	sd(case_site_ipd)
	controlMean 	=	mean(control_site_ipd)
	controlStd 		= 	sd(control_site_ipd)
	
	trim = floor( 0.001*length(control_site_ipd) ) : ceiling( (1-0.03)*length(control_site_ipd) )
	ctrlMean = mean( sort(control_site_ipd)[ trim ] )
	ipdRatio = ifelse( abs(ctrlMean) > 1e-3 , mean(case_site_ipd)/ ctrlMean , 1 )
	
	testResults = testProcedure(case_site_ipd, control_site_ipd)
    testStatistic = testResults[['testStatistic']]
    pvalue = testResults[['pvalue']]
	score = round(-10.0 * log10(pvalue))
	
	data.frame( caseMean=caseMean , caseStd=caseStd , controlMean=controlMean , controlStd=controlStd , 
				coverage=coverage , caseCoverage=caseCoverage , controlCoverage=controlCoverage , 
				ipdRatio=ipdRatio , score=score )
	
}
