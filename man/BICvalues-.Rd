\name{BICvalues<-}
\alias{BICvalues<-}
\title{
Generic Function - BICvalues<-
}
\description{
Generic function
}
\usage{
BICvalues(object)<-value
}
\arguments{
  \item{object}{
	Depending on the class of \code{object} depends on the method used (and if one exists)
}
  \item{value}{
	Replacement value
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
\code{\link{BICvalues<--methods}}
}
\examples{
x=new("cpt.ar") # new cpt object
BICvalues(x)<-c(0) # replaces the current BICvalues slot of x with empty vector.
}

\keyword{methods}
\keyword{cpt}
\keyword{internal} 