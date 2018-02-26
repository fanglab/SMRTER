
GMM = function( ipd , modelPrediction=NULL )
{
	if( is.null(modelPrediction) )
	{
		mu = c(0,3)
	} else {
		mu = c(modelPrediction,modelPrediction+2)
	}
	
	sdv = c(1,1)
	lamda = c(0.5,0.5)                                                                                                                                                                                       
	for(iter in 1:100)
	{
		###  update P
		P_gaussian = lapply( 1:length(lamda) , function(i) { dnorm( ipd , mean=mu[i], sd=sdv[i] )+1e-300 }  )  # mixture gaussian: p(x|q,mu,sdv) 
		### update Z                                                                                                                                                                                              
		Z =  sapply( 1:length(lamda), function(i)  {  P_gaussian[[i]] * lamda[i] } )  
		rowsum = rowSums(Z)
		Z =  apply(Z,2,function(x){ x/rowsum })	 # p(q|x,theta)
		colnames(Z) = c('theta1','theta2')	
		### update C and Theta
		if( !is.null(modelPrediction) )
		{
			mu[2]  =  sum(ipd * Z[,2]) / sum(Z[,2])
		} else {
			mu  =  sapply( 1:length(mu), function(i) {  sum(ipd * Z[,i],na.rm=T) / sum(Z[,i],na.rm=T)  }  ) 
		}
		sdv =  sapply( 1:length(mu), function(i) {  var = sum( (ipd-mu[i])^2 * Z[,i] ,na.rm=T ) / sum(Z[,i],na.rm=T ) ;  sqrt(var) }    )  
		### update lamda
		lamda = colSums(Z,na.rm=T)/nrow(Z) # p(j|theta); w
		#cat(lamda[1],lamda[2],mu[1],mu[2],sdv[1],sdv[2],'\n')
	}
	
	if(mu[2]<mu[1])
	{
		lamda = c(1,0)
		mu = rep(mean(ipd),2)
		sdv = rep(sd(ipd),2)
	}
	list(lamda=lamda,mu=mu,sd=sdv)
}


EMM = function( ipd , modelPrediction=NULL )
{	

	expPdf = function( data , mu ) { exp( -(data/mu) ) / mu }
		
	estimateSingleFraction = function(mu1, data, mu0, L)
	{
		a0 = expPdf(data, mu0)
		a1 = expPdf(data, mu1)
		res = optimize(mixModelFn , c(0.01, 0.99) , tol=1e-02 , a0=a0 , a1=a1 )$minimum
		#res = fminbound(mixModelFn, 0.01, 0.99, args=(a0, a1), xtol=1e-02)
		#if (sum(a1/a0) <= L) res = 0.0
		#if (sum(a0/a1) <= L) res = 1.0
		res	
	}

	mixModelFn = function( p, a0, a1)
	{
		tmp = (1 - p) * a0 + p * a1
		sum( -log(tmp[tmp>0]) )
	}
		  
	optimalMixProportion = function(data, mu0, L)
	{   
		#mu1 = fminbound(estimateSingleFraction, mu0, 10.0 * mu0, args=(data, mu0, L, False), xtol=1e-01) # shijia: fminbound estimates the best mu1 between mu0 and 10.0 * mu0
		mu1 = optimize(estimateSingleFraction, c( mu0, 10* mu0), tol=1e-01 , data=data , mu0=mu0, L=L )$minimum
		estimateSingleFraction(mu1, data, mu0, L)
	}
	
	L = length(ipd)
	optimalMixProportion(ipd, modelPrediction, L)
	
}   
   
   
   
   
