
IpdRatioTraditionalControl_SM = function(case_info,control_info,case_globalPercentile,control_globalPercentile, 
														coverage_per_molecule , number_of_molecule , case_coverage , cont_coverage)
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
			res = IpdRatioPositionTraditionalControl_SM( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile , case_moleculeID, coverage_per_molecule )
			rep_site = rep(site,length(res))
			cbind( rep_site , res )
		} 
	} )
	
	do.call(rbind,tmp)
}



IpdRatioPositionTraditionalControl_SM = function( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile , case_moleculeID, coverage_per_molecule )
{
	capping_site_ipd = cappingTraditionalControl( case_site_ipd, control_site_ipd , case_globalPercentile, control_globalPercentile)
	case_site_ipd = capping_site_ipd[['case']]
	control_site_ipd = capping_site_ipd[['control']]
	
	case_smol_ipd = tapply(case_site_ipd,case_moleculeID,mean)
	trim = floor( 0.001*length(control_site_ipd) ) : ceiling( (1-0.03)*length(control_site_ipd) )
	ctrlMean = mean( sort(control_site_ipd)[ trim ] )
	
	case_smol_ipdRatio = case_smol_ipd/ctrlMean
	smol_len = tapply(case_moleculeID,case_moleculeID,length)
	case_smol_ipdRatio = case_smol_ipdRatio[ smol_len>=coverage_per_molecule & !is.na(case_smol_ipdRatio) & !is.infinite(case_smol_ipdRatio)  ]
	if( abs(ctrlMean) < 1e-3 ) { case_smol_ipdRatio = rep(1,length(case_smol_ipdRatio)) }
	
	case_smol_ipdRatio
}






