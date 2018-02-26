
MultiPositionTraditionalControl = function(case_info,control_info,case_globalPercentile,control_globalPercentile,
											coverage_per_molecule , number_of_molecule , case_coverage , cont_coverage , score_thres )
{
	uni_tpl = intersect( unique( case_info$tpl ) , unique( control_info$tpl ) )
	tmp = lapply( sort(uni_tpl) , function(site) {
		#cat(site,'\n')
		case_site_info = case_info[ case_info$tpl == site , ]
		case_site_ipd = case_site_info$normIpd
		case_moleculeID = case_site_info$moleculeID
		
		control_site_info = control_info[ control_info$tpl == site , ]
		control_site_ipd = control_site_info$normIpd
		
		if( length(case_site_ipd)>=case_coverage & length(control_site_ipd)>=cont_coverage ) 
		{
			res = computePositionTraditionalControl( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile, case_moleculeID,
													coverage_per_molecule , number_of_molecule , case_coverage , cont_coverage , score_thres )
			data.frame( tpl=site , res)
		} 
	} )
	
	do.call(rbind,tmp)
}

computePositionTraditionalControl = function( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile, case_moleculeID=NULL, 
												coverage_per_molecule , number_of_molecule , case_coverage , cont_coverage , score_thres )
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
	frac = NA
	ipdRatio2 = NA
	molNum = NA
	
	if( !is.null(case_moleculeID) & score >= score_thres & length(case_site_ipd)>=case_coverage )
	{  
		smol_ipd = tapply(case_site_ipd,case_moleculeID,mean)
		smol_len = tapply(case_moleculeID,case_moleculeID,length)
		smol_ipd = smol_ipd[ smol_len>=coverage_per_molecule & !is.na(smol_ipd) & !is.infinite(smol_ipd)  ]
		#smol_ipd = smol_ipd[ !is.na(smol_ipd) & !is.infinite(smol_ipd)  ]
		molNum = length(smol_ipd)
		if( molNum>=number_of_molecule )
		{
			gmm = GMM( smol_ipd ,controlMean)
			frac = gmm$lamda[2]
			ipdRatio2 = gmm$mu[2]/controlMean
		} 
	} 
	
	if( is.null(case_moleculeID) & score >= score_thres & length(case_site_ipd)>=case_coverage & length(control_site_ipd)>=cont_coverage) 
	{
		frac = GMM( case_site_ipd ,controlMean)$lamda[2]
	}
	
	data.frame( caseMean=caseMean , caseStd=caseStd , controlMean=controlMean , controlStd=controlStd , 
				coverage=coverage , caseCoverage=caseCoverage , controlCoverage=controlCoverage , 
				molNum=molNum , ipdRatio=ipdRatio , score=score , ipdRatio2=ipdRatio2, frac=frac )
	
}
