\name{cpt.ar}

\alias{cpt.ar}

\title{

Identifying and selecting changes Regression

}

\description{

Calculates the optimal positioning and (potentially) number of changepoints for data using the user specified method. The model effectively checks several choices of AR models and fits the most appropriate choice.

}

\usage{
cpt.ar(data, penalty = "MBIC", pen.value = 0, min.order = 1, max.order=5,
method = "PELT", dist = "Normal", class = TRUE, param.estimates = TRUE,
minseglen = 1, shape = 0, tol = 1e-07)
}

\arguments{

  \item{data}{

	A vector or matrix of the time series data.

}

  \item{penalty}{

	Choice of "None", "SIC", "BIC", "MBIC", "AIC", "Hannan-Quinn", "Asymptotic", "Manual" and "CROPS" penalties.  If Manual is specified, the manual penalty is contained in the pen.value parameter. If Asymptotic is specified, the theoretical type I error is contained in the pen.value parameter.  If CROPS is specified, the penalty range is contained in the pen.value parameter; note this is a vector of length 2 which contains the minimum and maximum penalty value.  Note CROPS can only be used if the method is "PELT". The predefined penalties listed DO count the changepoint as a parameter, postfix a 0 e.g."SIC0" to NOT count the changepoint as a parameter.

}

  \item{pen.value}{

	The theoretical type I error e.g.0.05 when using the Asymptotic penalty.  A vector of length 2 (min,max) if using the CROPS penalty.  The value of the penalty when using the Manual penalty option - this can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.  

}

  \item{min.order}{

  A user defined minimum order p for the AR(p) models. The user may stipulate any value between 1 and max.order. Note that a choice of AR(0) would be a case of change in variance and so is not allowed.

}

  \item{max.order}{

  A user defined maximum order p for the AR(p) models. This is allowed so a user can stop the choice of any model beyond a certain order.

}
   \item{method}{

  Choice of "AMOC", "PELT" or "BinSeg".

}

  \item{dist}{

	This denotes the distribution that the data follows.

}

 \item{class}{

	Logical.  If TRUE then an object of class \code{cpt} is returned.

}

  \item{param.estimates}{

	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.

}

  \item{minseglen}{

  Positive integer giving the minimum segment length (no. of observations between changes), default is the minimum allowed by theory.

}

  \item{shape}{

  Value of the assumed known shape parameter required when test.stat="Gamma".

}

  \item{tol}{

  The tolerance at which model decisions are made.

}
}

\details{

	This function is used to find changes in regression for data. The changes are found using the PELT method supplied which can be single changepoint (AMOC) or multiple changepoints using exact (PELT).  A changepoint is denoted as the last observation of the segment / regime.

}

\value{

  If \code{class=TRUE} then an object of S4 class "cpt.reg" is returned.  The slot \code{cpts} contains the changepoints that are returned.  For \code{class=FALSE} the structure is as follows.



A list is returned with elements: 

  \item{order}{The functions decision on the optimal AR order for the data.}

	\item{lastchangecpts}{Vector of identified last changepoint locations.}

  \item{cpts}{Locations of the changepoints.}

  \item{lastchangelike}{Vector of likelihoods for the data.}

  \item{ncpts}{The number of changepoints that have been found at each given data point stored as a vector.}

}

\references{

Change in Normal mean: Hinkley, D. V. (1970) Inference About the Change-Point in a Sequence of Random Variables, \emph{Biometrika} \bold{57}, 1--17



PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598



CROPS: Haynes K, Eckley IA, Fearnhead P (2014) Efficient penalty search for multiple changepoint problems (in submission), arXiv:1412.3617



Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512



MBIC: Zhang, N. R. and Siegmund, D. O. (2007) A Modified Bayes Information Criterion with Applications to the Analysis of Comparative Genomic Hybridization Data. \emph{Biometrics} \bold{63}, 22-32.

}

\author{

Andrew Connell and Rebecca Killick

}





\seealso{

\code{\link{cpt.var}},\code{\link{cpt.meanvar}},\code{\link{plot-methods}},\code{\linkS4class{cpt}}

}

\examples{
set.seed(1)
  x <- 1:200
  data <- c()
  for(i in x){
    data <- c(data,rnorm(1,sample(1:i, 1, replace=TRUE),0.1*sample(1:i, 1, replace=TRUE)))
}

  #out <- cpt.ar(data)
  #cpts(out)      ##changepoints
  #param.est(out) ##parameter estimates (rows: beta estimates per segment)
  #plot(out)      ##plot of fit
  

}



\keyword{methods}

\keyword{univar}

\keyword{models}

\keyword{ts}
