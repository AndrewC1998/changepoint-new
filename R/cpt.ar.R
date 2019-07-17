cpt.ar <- function(data, penalty = "MBIC", pen.value = 0, min.order = 1, max.order=5, method = "PELT", dist = "Normal", class = TRUE, param.estimates = TRUE, minseglen = 1, shape = 0, tol = 1e-07){

  ##Check arguments are valid
  if(any(!complete.cases(data))){stop("data has missing values, this function cannot handle missing values")}
  if(any(!is.numeric(data))){stop("data must be a numeric vector")}
  if(!is.character(penalty) || length(penalty)>1)
  stop("Argument 'penelty' is invalid.")
  #value of 'penalty' & 'pen.value' checked within changepoint::penalty_decision
  if(!is.character(method) || length(method)>1)
  stop("Argument 'method' is invalid.")
  if(method!="AMOC" && method != "PELT") ##RESTRICTION IN USE
  stop("Invalid method, must be AMOC or PELT.")
  if(!is.character(dist) || length(dist)>1)
  stop("Argument 'dist' is invalid.")
  if(dist != "Normal"){  ##RESTRICTION IN USE
      warning(paste0("dist = ",dist," is not supported. Converted to dist='Normal'"))
      dist <- "Normal"
  }
  if(!is.logical(class) || length(class)>1)
  stop("Argument 'class' is invalid.")
  if(!is.logical(param.estimates) || length(param.estimates)>1)
  stop("Argument 'param.estimates' is invalid.")
  if(!is.numeric(minseglen) || length(minseglen)>1)
  stop("Argument 'minseglen' is invalid.")
  if(minseglen <= 0 || minseglen%%1 != 0) ##Further checks applied later
  stop("Argument 'minseglen' must be positive integer.")
  if(!is.numeric(tol) || length(tol)!=1)
  stop("Argument 'tol' is invalid.")
  if(tol<0) stop("Argument 'tol' must be positive.") ##Argument shape is assessed by the command where it is to be used.

  #Generate the minimum summary statistic to send to C
  if(min.order==0){
      warning("The choice of AR should not realistically be used in the case. AR(1) should be set as the minimum")
      sumstat = cbind(data,rep(1,length(data)))
  }else{
      sumstat = design(data,min.order)
  }

  if(minseglen >= min.order){
      minseglen = length(sumstat[1,])
  }

  cost_func = "ar.norm"
  n = length(sumstat[,1])
  MBIC = 0
  if(penalty=="MBIC"){MBIC=1}

  answer=list()
  answer[[7]] = 1
  on.exit(.C("FreePELT",answer[[7]]))

  answer <- .C('PELT', cost_func=cost_func, sumstat=as.double(sumstat), n=as.integer(n), m=as.integer(length(sumstat[1,])), pen=as.double(pen.value), cptsout=vector("integer",n), error=as.integer(0), shape=as.double(shape), minorder=as.integer(min.order), optimalorder = as.integer(min.order), maxorder = as.integer(max.order), minseglen=as.integer(minseglen),
  tol=as.double(tol), lastchangelike=vector("double", n+1), lastchangecpts=vector("integer", n+1), numchangecpts=vector("integer", n+1), MBIC=as.integer(MBIC))
  if(max.order < answer$optimalorder){
    answer$optimalorder = max.order
  }
  if(answer$err!=0){
      stop("C code error:",answer$err,call.=F)
  }

  if(class==TRUE){
    #Convert to cpt.reg object
    ans <- new("cpt.reg")
    data.set(ans) <- design(data, answer$optimalorder)
    cpttype(ans) <- paste0("Autoregressive structure of order ", answer$optimalorder)
    method(ans) <- method
    distribution(ans) <- dist
    pen.type(ans) <- penalty
    pen.value(ans) <- answer$pen
    cpts(ans) <- list(lastchangecpts=answer$lastchangecpts[1:answer$n], cpts=sort(answer$cptsout[answer$cptsout>0]), lastchangelike=answer$lastchangelike[1:answer$n], ncpts=answer$numchangecpts[1:answer$n])
    if(method=="PELT") ncpts.max(ans) <- Inf
    if(param.estimates) ans = param(ans)
     return(ans)
  }else{
    return(list(order = answer$optimalorder, lastchangecpts=answer$lastchangecpts[1:answer$n], cpts=sort(answer$cptsout[answer$cptsout>0]), lastchangelike=answer$lastchangelike[1:answer$n], ncpts=answer$numchangecpts[1:answer$n]))
  }
}
