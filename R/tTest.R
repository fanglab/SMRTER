
computeObservationTstatistic = function(tMean, tErr ,  modelPrediction ,coverage)   ##  synthetic control
{
	em = 0.01 + 0.03 * modelPrediction + 0.06 * modelPrediction ^ (1.7)
	t = -(tMean - modelPrediction) / sqrt(em^2 + tErr^2)
    df = max(1, coverage-1)
	pvalue = pt(t,df)
	pvalue
}

testProcedure = function(x , y , exclude=0.95)  ##  traditional control
{
	cappedSlog = function(v)
	{
		q = quantile(v, exclude)
		v = v[!is.na(v)]
		v[v>q] = q
		v[v<=0] = 1/(75+1)
		log(v)
	}
	
	erfc = function(x) { 2*pnorm(-sqrt(2)*x) }
	
	x1 = cappedSlog(x)
    x2 = cappedSlog(y)
    sx1 = var(x1) / length(x1)
    sx2 = var(x2) / length(x2)
	totalSE = sqrt(sx1 + sx2)
    if (totalSE == 0)
	{	
		stat = 0	
	} else {
		stat = ( mean(x1) - mean(x2)) / totalSE
	}

	pval = 0.5 * erfc( stat / sqrt(2) )
	
	list( testStatistic=stat, pvalue=pval )
	
}



identificationQv = function( smrt_region , pos , positiveControl , insilicoControl   )
{

	logpdf = function(x, df)
	{
		r = df
		lPx = lgamma((r + 1) / 2) - lgamma(r / 2)
        lPx = lPx - 0.5 * log(r * pi) - (r + 1) / 2 * log(1 + (x ** 2) / r)
        lPx
	}
	
	singleScore = function( tMean, tErr ,  modelPrediction ,coverage )
	{
		#prior = self.modPriors[context[self.pre]]
		
		em = 0.01 + 0.03 * modelPrediction + 0.06 * modelPrediction ^ (1.7)
		t = -(tMean - modelPrediction) / sqrt(em^2 + tErr^2)
		df = max(1, coverage-1)
		#log( dt(t,df) )
		logpdf(t, df)
	}
	
	scoreRegion = function( smrt_region , start , end , controlModel )
	{
		sc = 0
		for( loci in start:(end+1) )
		{		
			index = which( smrt_region$tpl == loci )
			if( length(index)==1 )
			{
				likelihood = singleScore( smrt_region$tMean[index] , smrt_region$tErr[index] ,  controlModel[index] ,smrt_region$coverage[index] )
			} else {
				likelihood = 0
			}
			sc = sc + likelihood
		}
		sc
	}

	post = 4 
	pre  = 10
	
	modScore   =  scoreRegion( smrt_region , pos - post, pos + pre, positiveControl ) + log(0.05) # log(0.05) is prior for H J K
	noModScore =  scoreRegion( smrt_region , pos - post, pos + pre, insilicoControl )
	llr = modScore - noModScore
	qModScore = 10 * llr * log10(exp(1)) + 10 * log10( 1+exp(-llr) )
	qModScore
	
}



