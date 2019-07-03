PELT = function(sumstat, pen=0, cost_func = "norm.mean", shape = 1, minseglen = 1){
  # function that uses the PELT method to calculate changes in mean where the segments in the data are assumed to be Normal
  n = length(sumstat[,1])-1
  m = length(sumstat[1,])
  tol = 0
  if(cost_func == "norm.mean" || cost_func == "var.norm" || cost_func == "meanvar.norm" || cost_func == "meanvar.exp" || cost_func == "meanvar.gamma" || cost_func == "meanvar.poisson"){
    MBIC = 0
  }else{
    MBIC = 1
  }
  if(n<2){stop('Data must have at least 2 observations to fit a changepoint model.')}

  storage.mode(sumstat) = 'double'
  error=0

  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)

  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'

  answer=list()
  answer[[7]]=1
  on.exit(.C("FreePELT",answer[[7]]))

  storage.mode(lastchangelike) = 'double'
  storage.mode(lastchangecpts) = 'integer'
  storage.mode(numchangecpts) = 'integer'

  min = 0
  optimal = 0
  max = 0

  answer=.C('PELT', cost_func, sumstat, as.integer(n), as.integer(m), as.double(pen), cptsout, as.integer(error), as.double(shape), as.integer(min), as.integer(optimal), as.integer(max), as.integer(minseglen), as.double(tol), lastchangelike, lastchangecpts, numchangecpts, as.integer(MBIC))


  if(answer[[7]]>0){
    stop("C code error:",answer[[7]],call.=F)
  }

  return(list(lastchangecpts=answer[[15]],cpts=sort(answer[[6]][answer[[6]]>0]), lastchangelike=answer[[14]], ncpts=answer[[16]]))

}
