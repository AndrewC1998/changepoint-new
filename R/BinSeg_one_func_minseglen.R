BINSEG = function(sumstat, pen = 0, cost_func = "norm.mean", shape = 1, minseglen = 2,  Q=5){

  n = length(sumstat[,1])
  m = length(sumstat[1,])
  tol = 0
  MBIC = 0
  if(n<2){stop('Data must have at least 2 observations to fit a changepoint model.')}
  if(Q>((n/2)+1)){stop(paste('Q is larger than the maximum number of segments',(n/2)+1))}
  if(Q<=0){stop(paste('Q is the maximum number of changepoints so should be greater than 0'))}
  shape = 1

  storage.mode(sumstat) = 'double'

  cptsout = rep(0,Q) # sets up null vector for changepoint answer
  likeout = rep(0,Q) # sets up null vector for likelihood of changepoints in cptsout
  storage.mode(cptsout)='integer'
  storage.mode(likeout)='double'
  op_cps = 0
  min = 0
  optimal = 0
  max = 0
  error = 0

  answer=.C('binseg', cost_func, sumstat, as.integer(n), as.integer(m), as.integer(min), as.integer(optimal), as.integer(max), as.double(pen), as.integer(Q), cptsout, as.integer(error), as.integer(minseglen), likeout, as.integer(op_cps), as.double(shape), as.double(tol), as.integer(MBIC))
  if(answer[[9]]==Q){warning('The number of changepoints identified is Q, it is advised to increase Q to make sure changepoints have not been missed.')}
  if(answer[[9]]==0){cpts=n}
  else{cpts=c(sort(answer[[6]][1:answer[[9]]]),n)}
  return(list(cps=rbind(answer[[6]],2*answer[[8]]),cpts=cpts,op.cpts=answer[[9]],pen=pen))
  ##answer[6] is cptsout, answer[8] is likeout ("beta value")
  }
