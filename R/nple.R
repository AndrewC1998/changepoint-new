#Detect changes in variance using a localised estimate of the variance from LSW model
nple <- function(data, penalty = "Manual", pen.value = 6 * log(length(data) - crop), method = "PELT", class = TRUE, minseglen = 30, nquantiles = 10, crop = 16, family = "DaubExPhase", filter.number = 10, binwidth = 1){

    n = length(data)

    data.lacv <- my.lacf(data, family = family, filter.number = filter.number, binwidth = binwidth, lag.max = 0)[1:(n - crop)]

    cps <- cpt.np(data = data.lacv, penalty = penalty, pen.value = pen.value, method = method, class = class, minseglen = minseglen, nquantiles = nquantiles)

    return(cps)
  }


my.lacf <- function(x, filter.number = 10, family = "DaubLeAsymm", smooth.dev = var, AutoReflect = TRUE, lag.max = NULL, WPsmooth.type = "RM", binwidth, tol = 0.1, maxits = 5, ABBverbose = 0, verbose = FALSE, ...){
      # function taken from locits 1.7.3
      # RK edit to non-2^J length
      TT <- length(x)
      Jp <- ceiling(logb(TT, 2))
      add <- 2^Jp - TT
      # NvSK says Nh is # of non-zero elements in the filter.
      # The support of the discrete wavelets is something like this:
      # Lj<-(2^j -1)*(Nh-1)+1
      filter <- wavethresh::filter.select(filter.number, family)
      Nh <- length(filter$H != 0)
      Lj=(2^(1:Jp-1))*(Nh-1)+1

      xa <- c(rep(0, times = floor(add/2)), x,rep(0,times=ceiling(add/2)))
      lxa <- length(xa) #should be 2^(J+1) if TT not equal to  2^J


      dsname = deparse(substitute(x))
      if (WPsmooth.type == "RM") {
        if (missing(binwidth) || binwidth == 0)
          binwidth <- locits::AutoBestBW(x = xa, filter.number = filter.number, family = family, smooth.dev = smooth.dev, AutoReflect = AutoReflect, tol = tol, maxits = maxits, plot.it = FALSE, verbose = ABBverbose)
        if (verbose == TRUE)
          cat("Linear Smoothing. Bandwidth is: ", binwidth, "\n")
      }
      S <- locits::ewspec3(x = xa, filter.number = filter.number, family = family, smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = WPsmooth.type, binwidth = binwidth, ...)$S
      Smat <- matrix(S$D, nrow = length(xa), ncol = Jp)
      Psi <- wavethresh::PsiJmat(-Jp+(add!=0), filter.number = filter.number, family = family)
      nc <- ncol(Psi)
      L <- (nc - 1)/2
      dimnames(Psi) <- list(NULL, c(-L:0, 1:L))
      if (is.null(lag.max))
        lag.max <- floor(10 * (log10(length(x))))
      if (L + 1 + lag.max > ncol(Psi)) {
        warning(paste("lag.max too high. Have reset it to ",
                      ncol(Psi) - L - 1, ". Higher lags are zero"))
        lag.max <- ncol(Psi) - L - 1
      }

      if(add>0){
        Smat=Smat[,-Jp] # remove extra level
        Smat=Smat[-(1:floor(add/2)),] # remove start
        Smat=Smat[-(length(x)+1:ceiling(add/2)),] # remove end
      }
      the.lacf <- Smat %*% Psi[, (L + 1):(L + 1 + lag.max)]
      return(the.lacf)
  }
