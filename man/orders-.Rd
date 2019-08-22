\name{orders<-}
\alias{orders<-}
\title{
Generic Function - orders<-
}
\description{
Generic function
}
\usage{
orders(object)<-value
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
\code{\link{orders<--methods}}
}
\examples{
x=new("cpt.ar") # new cpt object
orders(x)<-array(0) # replaces the current orders slot of x with empty vector.
}

\keyword{methods}
\keyword{cpt}
\keyword{internal} 