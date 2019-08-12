cpt.plot=function(data, type = "mean", penalty="MBIC", pen.value=0 ,method = "AMOC", Q=5, test.stat="Normal", class=TRUE, param.estimates=TRUE, shape=1, tol = 1e-07, know.mean = FALSE, mu = NA, minseglen=2){
    # Q = maximum number of changepoints
    if(type == "mean"){
        n = length(data)
        if(n > Q){
            for(i in (2*Q):n){
                Sys.sleep(0.1)
                plot(x=1:i,y=data[1:i],xlim=c(0,n),ylim=range(data), xlab = i, ylab = "", type ="l")  #type = "l" can be used for line plot
                ansvar=cpt.mean(data[1:i],penalty=penalty,pen.value=pen.value,method=method,Q=Q,test.stat=test.stat,class=class,param.estimates=param.estimates,minseglen=minseglen)
                change.values <- cpts(ansvar)
                if(length(ansvar@cpts)>0){
                    for(i in 1:Q){
                        abline(v=change.values, col = "blue")
                    }
                }
            }
        }else{
            stop("Q must be less than n")
        }
    }else if(type == "var"){
      n = length(data)
      if(n > Q){
          for(i in (2*Q):n){
              Sys.sleep(0.1)
              plot(x=1:i,y=data[1:i],xlim=c(0,n),ylim=range(data), xlab = i, ylab = "", type ="l")  #type = "l" can be used for line plot
              ansvar=cpt.var(data[1:i],penalty=penalty,pen.value=pen.value,know.mean=know.mean,mu=mu,method=method,Q=Q,test.stat=test.stat,class=class,param.estimates=param.estimates,minseglen=minseglen)
              change.values <- cpts(ansvar)
              if(length(ansvar@cpts)>0){
                  for(i in 1:Q){
                      abline(v=change.values, col = "blue")
                  }
              }
          }
      }else{
          stop("Q must be less than n")
      }
    }else if(type == "meanvar"){
      n = length(data)
      if(n > Q){
          for(i in (2*Q):n){
              Sys.sleep(0.1)
              plot(x=1:i,y=data[1:i],xlim=c(0,n),ylim=range(data), xlab = i, ylab = "", type ="l")  #type = "l" can be used for line plot
              ansvar=cpt.meanvar(data[1:i],penalty=penalty,pen.value=pen.value,method=method,Q=Q,test.stat=test.stat,class=class,param.estimates=param.estimates,shape=shape,minseglen=minseglen)
              change.values <- cpts(ansvar)
              if(length(ansvar@cpts)>0){
                  for(i in 1:Q){
                      abline(v=change.values, col = "blue")
                  }
              }
          }
      }else{
          stop("Q must be less than n")
      }
    }else if(type == "reg"){
      if(class(data) != "matrix"){
        stop("Must manually create the design matrix.")
      }
      n = length(data[,1])
      if(n > Q){
          for(i in Q:n){
              Sys.sleep(0.1)
              plot(x=1:i,y=data[1:i,1],xlim=c(0,n),ylim=range(data[,1]), xlab = i, ylab = "", type ="l")  #type = "l" can be used for line plot
              ansvar=cpt.reg(data[1:i,],penalty=penalty,pen.value=pen.value,method=method,dist=test.stat,class=class,param.estimates=param.estimates,minseglen=minseglen,shape=shape,tol=tol)
              change.values <- cpts(ansvar)
              if(length(ansvar@cpts)>0){
                  for(i in 1:Q){
                      abline(v=change.values, col = "blue")
                  }
              }
          }
      }else{
          stop("Q must be less than n")
      }
    }else{
        stop("Type must be mean, var, meanvar or reg")
    }
}
