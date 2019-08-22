\name{BICvalues}
\alias{BICvalues}
\title{
Generic Function - BICvalues
}
\description{
Generic function
}
\usage{
BICvalues(object)
}
\arguments{
  \item{object}{
	Depending on the class of \code{object} depends on the values used (and if one exists)
}
}
\details{
	Generic Function
}
\value{
	Depends on the class of \code{object}, see individual methods
}
\author{
Andrew Connell
}


\seealso{
\code{\link{BICvalues-methods}}
}
\examples{
x=new("cpt.ar") # new cpt object
BICvalues(x) # retrieves the BICvalues slot from x
}

\keyword{methods}
\keyword{cpt}
\keyword{internal}