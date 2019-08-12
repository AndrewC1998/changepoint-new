CROPS <- function(data, penalty="CROPS", pen.value, method="PELT", test.stat="Normal", class=TRUE, param.est=TRUE, minseglen, shape, nquantiles, func){
  if(method != "PELT"){stop('CROPS is a valid penalty choice only if method="PELT", please change your method or your penalty.')}
  if (test.stat == "empirical_distribution"){
    nonparametric.ed.sumstat = function(data,K=nquantiles){ # This now takes into account the integral transformation
      ##USE K points in integral
      n <- length(data)
      if(K>n) K=n
      Q <- matrix(0,K,n+1)
      x=sort(data)
      yK= -1 + (2*(1:K)/K-1/K)
      c=-log(2*n-1)
      pK=(1+exp(c*yK))^-1
      for (i in 1:K){
        j=as.integer((n-1)*pK[i] + 1)
        Q[i,-1] <- cumsum(data<x[j])+0.5*cumsum(data==x[j])
      }
      return(Q)
    }
    sumstat <- nonparametric.ed.sumstat(data, K = nquantiles)
  }else{
    nquantiles <- NA
    mu <- mean(data)
    sumstat <- cbind(c(0,cumsum(coredata(data))),c(0,cumsum(coredata(data)^2)),cumsum(c(0,(coredata(data)-mu)^2)))
  }

  switch(test.stat,
    "empirical_distribution" = {stat="ed"},
    "Normal" = {stat = "norm"},
    "Exponential" = {stat = "exp"},
    "Gamma" = {stat = "gamma"},
    "Poisson" = {stat = "poisson"},
    {stop("Only empirical_distribution, Normal, Exponential, Gamma and Poisson are valid test statistics")}
  )
  costfunc = paste0(func, ".", stat)

  out = range_of_penalties(sumstat, cost=costfunc, min_pen=pen.value[1], max_pen=pen.value[2], minseglen=minseglen, nquantiles=nquantiles)

  if(func=="nonparametric"){
    cpttype="nonparametric"
  }else if(func=="var"){
    cpttype="variance"
  }else if(func=="meanvar"){
    cpttype="mean and variance"
  }else{
    cpttype="mean"
  }

  if(class==TRUE){
      ans = class_input(data=data,cpttype=cpttype, method=method, test.stat=test.stat, penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.est, out=out,shape=shape)
      if(func=="var"){
        param.est(ans)=c(param.est(ans),mean=mu)
      }
    return(ans)
  }else{return(out)}
}
