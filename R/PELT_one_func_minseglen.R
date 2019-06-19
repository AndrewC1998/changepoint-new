PELT = function(sumstat,pen=0, cost_func = "norm.mean", shape = 1, minseglen = 1){
  # function that uses the PELT method to calculate changes in mean where the segments in the data are assumed to be Normal
  n = length(sumstat[,1]) - 1
  m = length(sumstat[1,]) - 1
  MBIC = 0
  tol = 0 
  if(n<2){stop('Data must have at least 2 observations to fit a changepoint model.')}
  
  storage.mode(sumstat) = 'double'
  error=0
  
  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)
  
  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'
  
  answer=list()
  answer[[6]]=1
  on.exit(.C("FreePELT",answer[[6]]))
  
  storage.mode(lastchangelike) = 'double'
  storage.mode(lastchangecpts) = 'integer'
  storage.mode(numchangecpts) = 'integer'
  
  min = 0
  optimal = 0
  max = 0
  
  # answer=.C('PELT',cost_func, y3, y2,y,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape))
  answer=.C('PELT', cost_func, sumstat, as.integer(n), as.integer(m), as.integer(min), as.integer(optimal), as.integer(max), as.double(pen), cptsout, as.integer(error), as.double(shape), as.integer(minseglen), as.double(tol), lastchangelike, lastchangecpts, numchangecpts, as.integer(MBIC))
  
  
  if(answer[[6]]>0){
    stop("C code error:",answer[[6]],call.=F)
  }
  return(list(lastchangecpts=answer[[14]],cpts=sort(answer[[8]][answer[[8]]>0]), lastchangelike=answer[[13]], ncpts=answer[[15]]))
  
}
