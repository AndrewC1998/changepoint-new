context("plot(diagnostic = TRUE) tests")

# Generate cpt.range object
testdata <- changepoint::ftse100$V2
obj.cpt.range <- cpt.var(testdata, method = "PELT",
                         penalty = "CROPS", pen.value = c(5, 500))

# For code coverage
plot(obj.cpt.range, diagnostic = TRUE)
plot(obj.cpt.range, diagnostic = TRUE, type = "h")

set.seed(1) # Note: new data sets must be added at the end.
singmeandata <- c(rnorm(100,0,1),rnorm(100,10,1))
mulmeandata <- c(rnorm(100,0,1),rnorm(100,10,1),rnorm(100,20,1),rnorm(100,50,1))
nochangedata <- c(rnorm(200,0,1))
singvardata <- c(rnorm(100,10,1),rnorm(100,10,5))
mulvardata <- c(rnorm(100,20,10),rnorm(100,20,15),rnorm(100,20,20),rnorm(100,20,25))
singmeanvardata <- c(rnorm(50,0,1),rnorm(50,3,10))
mulmeanvardata <- c(rnorm(50,0,1),rnorm(50,5,3),rnorm(50,10,1),rnorm(50,3,10))
mulmeanvarexpdata <- c(rexp(50,1), rexp(50,3), rexp(50,5), rexp(50,7)) #rate values correct
mulmeanvarpoisdata <- c(rpois(50,1), rpois(50,2), rpois(50,3), rpois(50,5)) #lambda values correct?
constantdata <- rep(1, 200)
shortdata <- c(2)
negativedata <- jitter(rep(-100, 200) )
characterdata <- rep("ert", 200)
#NAdata - creates 10 random NA within singmeandata
NAdata <- singmeandata
rn <- sample(1:length(singmeandata), 10, replace=F)
for(i in rn){
  NAdata[i] <- NA
}
NAdata[1] <- NA
data <- list(singmeandata,mulmeandata, nochangedata, singvardata, mulvardata, mulmeanvardata, mulmeanvarexpdata, mulmeanvarpoisdata, constantdata)

type <- c("mean", "var", "meanvar", "reg", "other")

for(i in 1:length(data)){
  for( j in 1:length(type)){
    if(type[j] == "mean" || type[j] == "var" || type[j] == "meanvar"){
      expect_error(cpt.plot(data[[i]], type = type[j], Q = length(data[[i]]) + 1), "Q must be less than n")
    }else if(type[j] == "reg"){
      expect_error(cpt.plot(data[[i]], type = type[j], Q = length(data[[i]]) + 1), "Must manually create the design matrix.")
      expect_error(cpt.plot(design(data[[i]],i), type = type[j], Q = length(data[[i]]) + 1), "Q must be less than n")
    }else{
      expect_error(cpt.plot(data[[i]], type = type[j], Q = length(data[[i]]) + 1), "Type must be mean, var, meanvar or reg")
    }
  }
}

# Tests for plots
#These functions are still somewhat experimental and even though the generated plots match the true plots
#there still seems to be some issues being caused. https://github.com/r-lib/vdiffr/issues
#we will readd this feature once the original issue is resolved.

#vdiffr::expect_doppelganger("Diagnostic plot (default)", plot(obj.cpt.range, diagnostic = TRUE))
#vdiffr::expect_doppelganger("Diagnostic plot (histogram)", plot(obj.cpt.range, diagnostic = TRUE, type = "h"))
