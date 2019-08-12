\name{NPPELT}
\alias{NPPELT}
\title{
NPPELT (Non-Parametric Pruned Exact Linear Time) - Only intended for developer use.
}
\description{
Implements the Non-parametric PELT method for identifying changepoints in a given set of summary statistics for a specified cost function and penalty.

This function is called by \code{cpt.np}and \code{cpt.nple} when \code{method="PELT"}.  This is not intended for use by regular users of the package.  It is exported for developers to call directly for speed increases or to fit alternative cost functions.

WARNING: No checks on arguments are performed!
}
\usage{
NPPELT(sumstat,pen=0, cost_func = "nonparametric.ed", minseglen = 1, nquantiles = 100)
}
\arguments{
  \item{sumstat}{
	A matrix containing the summary statistics of data within which you wish to find a changepoint.  Currently assumes 3 columns and uses the number of rows as the length of the data +1 (initial value of 0).
}
  \item{pen}{
  Default choice is 0, this should be evaluated elsewhere and a numerical value entered.  This should be positive - this isn't checked but results are meaningless if it isn't.
}
  \item{cost_func}{
  The friendly name of the cost function to be called in C.  If using your own cost function, this must be the name of the C function to use.
}
  \item{minseglen}{
  Positive integer giving the minimum segment length (no. of observations between changes), default is 1.  No checks are performed on the input value so it could be larger than feasible to have changes in the data.
  }
  \item{nquantiles}{
  The number of quantiles to calculate when test.stat = "empirical_distribution".
  }
}
\details{
	This function is used as a wrapper function to implement the non-parametric PELT algorithm in C. It simply creates the necessary worker vectors, ensures all inputs are the correct type, and passes everything to the C function.
	
	This function is exported for developer use only.  It does not perform any checks on inputs (other than type coersion) and is simply a wrapper function for the C code.
}
\value{
	A list is returned with elements:
	\item{lastchangecpts}{Vector of length n containing the last changepoint prior to each timepoint.}
	\item{cpts}{Ordered list of optimal number of changepoints ending with n.}
	\item{lastchangelike}{Vector of lenght n containing the likelihood of the optimal segmentation up to each timepoint.}
	\item{ncpts}{Number of changes identified.}
}
\references{
PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598

CROPS: Haynes K, Eckley IA, Fearnhead P (2014) Efficient penalty search for multiple changepoint problems, \emph{Springer} \bold{27(05)}, 1293--1305
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.np}},\code{\link{cpt.nple}},\code{\link{plot-methods}},\code{\linkS4class{cpt}}
}
\examples{
#This function should only be used by developers, see its use in cpt.np.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
