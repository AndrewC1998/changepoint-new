context("cpt.np tests")

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
data <- list(singmeandata,mulmeandata, nochangedata, singvardata, mulvardata, mulmeanvardata, mulmeanvarexpdata, mulmeanvarpoisdata, constantdata, negativedata)

method <- c("AMOC", "BinSeg", "SegNeigh")
method2 <- c("BinSeg", "BinSeg", "SegNeigh")

for(i in 1:length(data)){
  for(j in 1:length(method)){
    suppressWarnings(expect_error(cpt.np(data[[i]], test.stat = "CUSUM"), "Invalid Method, must be AMOC, SegNeigh or BinSeg"))

    suppressWarnings(expect_error(cpt.np(data[[i]], test.stat = "CUSUM", method = method[j]), "MBIC penalty is not valid for nonparametric test statistics."))

    suppressWarnings(expect_error(cpt.np(data[[i]], test.stat = "CUSUM", method = method[j], penalty = "Asymptotic"), "Asymptotic penalty values must be > 0 and <= 1"))

    suppressWarnings(expect_error(cpt.np(data[[i]], test.stat = "CUSUM", method = method[j], penalty = "Asymptotic", pen.value = 0.01), "Asymptotic penalties have not been implemented yet for CUSUM"))

    expect_warning(cpt.np(data[[i]], test.stat = "CUSUM", method = method[j], penalty = "Manual"), "Traditional penalty values are not appropriate for the CUSUM test statistic")

    suppressWarnings(expect_error(cpt.np(data[[i]], test.stat = "CSS"), "CSS does not satisfy the assumptions of PELT, use SegNeigh or BinSeg instead."))

    suppressWarnings(expect_error(cpt.np(data[[i]], test.stat = "CSS", method = method2[j]), "MBIC penalty is not valid for nonparametric test statistics."))

    expect_warning(cpt.np(data[[i]], test.stat = "CSS", method = method2[j], penalty = "Asymptotic", pen.value = 0.01), "Asymptotic penalty value is not accurate for multiple changes, it should be treated the same as a manual penalty choice.")
  }
}
