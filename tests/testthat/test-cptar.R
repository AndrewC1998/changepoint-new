context("cpt.ar tests")

# testing functions, aim to get 100% test coverage on exported code
# testing for cpt.reg function

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
data <- list(singmeandata,mulmeandata, nochangedata, singvardata, mulvardata, mulmeanvardata, mulmeanvarexpdata, mulmeanvarpoisdata, negativedata)

expect_error(cpt.ar(NAdata), "data has missing values, this function cannot handle missing values")

expect_error(cpt.ar(characterdata), "data must be a numeric vector")

for(i in 1:length(data)){
  expect_error(cpt.ar(data[[i]], penalty = 1), "Argument 'penalty' is invalid.")

  expect_error(cpt.ar(data[[i]], method = 1), "Argument 'method' is invalid.")

  expect_error(cpt.ar(data[[i]], method = "other method"), "Invalid method, must be PELT.")

  expect_error(cpt.ar(data[[i]], dist = 1), "Argument 'dist' is invalid.")

  #expect_warning(cpt.ar(data[[i]], dist = "Exponential", class = FALSE), "dist = Exponential is not supported. Converted to dist='Normal'")

  expect_error(cpt.ar(data[[i]], class = 1), "Argument 'class' is invalid.")

  expect_error(cpt.ar(data[[i]], param.estimates = 1), "Argument 'param.estimates' is invalid.")

  expect_error(cpt.ar(data[[i]], minseglen = "character"), "Argument 'minseglen' is invalid.")

  expect_error(cpt.ar(data[[i]], minseglen = -2), "Argument 'minseglen' must be positive integer.")

  expect_error(cpt.ar(data[[i]], tol = "character"), "Argument 'tol' is invalid.")

  expect_error(cpt.ar(data[[i]], tol = -2), "Argument 'tol' must be positive.")

  expect_error(cpt.ar(data[[i]], min.order = 0), "The choice of AR should not realistically be used in the case. AR with order 1 should be set as the minimum")

  expect_error(cpt.ar(data[[i]], max.order = length(data[[i]])), "The order trying to be fit is unrealistic and will lead to errors. If you wish to try an explicitly large model manually use cpt.reg")
}
