\name{orders}
\alias{orders}
\title{
Generic Function - orders
}
\description{
Generic function
}
\usage{
orders(object)
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
\code{\link{orders-methods}}
}
\examples{
x=new("cpt.ar") # new cpt object
orders(x) # retrieves the orders slot from x
}

\keyword{methods}
\keyword{cpt}
\keyword{internal}