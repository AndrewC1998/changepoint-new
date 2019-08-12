\name{cpt.nple}
\alias{cpt.nple}
\title{
Detects changes in local stationary wavelet autocorelation.
}
\description{
Plots the optimal positioning of changepoints for data using the user specified method. 
}
\usage{
cpt.nple(data, penalty = "Manual", pen.value = 6 * log(length(data) - crop), 
         method = "PELT", class = TRUE, minseglen = 30, nquantiles = 10, 
        crop = 16, family = "DaubExPhase", filter.number = 10, binwidth = 1)
}
\arguments{
  \item{data}{
	The data that the user wants to plot. 
  }
  \item{penalty}{
  Choice of "None", "SIC", "BIC", "MBIC", AIC", "Hannan-Quinn", "Manual" and "CROPS" penalties.  If Manual is specified, the manual penalty is contained in the pen.value parameter. If CROPS is specified, the penalty range is contained in the pen.value parameter; note this is a vector of length 2 which contains the minimum and maximum penalty value.  Note CROPS can only be used if the method is "PELT". The predefined penalties listed DO count the changepoint as a parameter, postfix a 0 e.g."SIC0" to NOT count the changepoint as a parameter.
  }
  \item{pen.value}{
  The value of the penalty when using the Manual penalty option.  A vector of length 2 (min,max) if using the CROPS penalty.
  }
  \item{method}{
  Currently the only method is "PELT".
  }
  \item{class}{
  Logical. If TRUE then an object of class cpt is returned.
  }
  \item{minseglen}{
  Positive integer giving the minimum segment length (number of observations between changes), default is the minimum allowed by theory.
  }
  \item{nquantiles}{
  The number of quantiles to calculate.
  }
  \item{crop}{
  The amount of data pruned.
  }
  \item{family}{
   Option passed through to wavethresh functions. The option allows selection from which of a range of families the wavelet will come from. Please read wavethresh documentation for more details.
  }
  \item{filter.number}{
   Option passed through to wavethresh functions. Allows user to choose a filter. Please read wavethresh documentation for more details.
  }
  \item{binwidth}{
   If the periodogram smoothing is "RM" then this is the number of consecutive observations used in the running mean smooth. Currently the only periodogram smoothing supported is "RM".
  }
}
\details{
 A nonparametric approach to detecting changes in variance within a time series which we demonstrate is resilient to de- partures from the assumption of Normality or presence of outliers. Our method is founded on a local estimate of the variance provided by the Locally Stationary Wavelet (LSW) framework. Within this setting, the structure of this local estimate of the variance will be piecewise constant if a time series has piecewise constant variance. Consequently, changes in
the variance of a time series can be detected in a non-parametric setting.
}
\value{
 If \code{class=TRUE} then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are returned.  For \code{class=FALSE} the structure is as follows.

If data is a vector (single dataset) then a vector/list is returned depending on the value of method.  If data is a matrix (multiple datasets) then a list is returned where each element in the list is either a vector or list depending on the value of method.

If method is PELT then a vector is returned containing the changepoint locations for the penalty supplied. If the penalty is CROPS then a list is returned with the elements:

\item{cpt.out}{A data frame containing the value of the penalty value where the number of segmentations chages, the number of segmentations and the value of the cost at that penalty value.}
\item{changepoints}{The optimal changepoints for the different penalty values startings with the lowest penalty value.}
}
\references{
Change in Normal mean: Hinkley, D. V. (1970) Inference About the Change-Point in a Sequence of Random Variables, \emph{Biometrika} \bold{57}, 1--17

PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598

MBIC: Zhang, N. R. and Siegmund, D. O. (2007) A Modified Bayes Information Criterion with Applications to the Analysis of Comparative Genomic Hybridization Data. \emph{Biometrics} \bold{63}, 22-32.
}
\author{
Andrew Connell, Rebecca Killick
}


\seealso{
\code{\link{cpt.np}},\code{\linkS4class{cpt}}
}
\examples{
 x <- c(rnorm(100,50,1), rnorm(100,50,3))
 cpt.nple(x)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
