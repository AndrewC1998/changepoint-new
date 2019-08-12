NPPELT = function(sumstat,pen=0, cost_func = "nonparametric.ed", minseglen = 1, nquantiles = 100){
  # function that uses the PELT method to calculate changes in mean where the segments in the data are assumed to be Normal
  n=dim(sumstat)[2]-1
  if(n<2){stop('Data must have at least 2 observations to fit a changepoint model.')}

  storage.mode(sumstat) = 'double'
  error=0
  
  if(cost_func=="nonparametric.ed"){
      MBIC = 0
  }else{
      MBIC = 1
  }

  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)

  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'

  answer=list()
  answer[[6]]=1
  on.exit(.C("FreeNPPELT",answer[[6]]))

  storage.mode(lastchangelike) = 'double'
  storage.mode(lastchangecpts) = 'integer'
  storage.mode(numchangecpts) = 'integer'

  answer=.C('NPPELT',cost_func = cost_func, sumstat = sumstat, n = as.integer(n), pen = as.double(pen), cptsout = cptsout, error = as.integer(error), minseglen = as.integer(minseglen), nquantiles = as.integer(nquantiles), lastchangelike = lastchangelike, lastchangecpts = lastchangecpts, numchangecpts = numchangecpts, MBIC = as.integer(MBIC))
  if(answer$error>0){
    stop("C code error:",answer$error,call.=F)
  }
  return(list(lastchangecpts = answer$lastchangecpts, cpts = sort(answer$cptsout[answer$cptsout>0]), lastchangelike = answer$lastchangelike, ncpts = answer$numchangecpts))

}
