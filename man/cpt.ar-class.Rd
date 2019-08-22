\name{cpt.ar-class}
\Rdversion{1.1}
\docType{class}
\alias{cpt.ar-class}
\alias{cpt.ar-method,cpts}
\alias{cpt.ar-method,cpttype}
\alias{cpt.ar-method,data.set}
\alias{cpt.ar-method,data.set.ts}
\alias{cpt.ar-method,test.stat}
\alias{cpt.ar-method,ncpts.max}
\alias{cpt.ar-method,method}
\alias{cpt.ar-method,minseglen}
\alias{cpt.ar-method,param.est}
\alias{cpt.ar-method,pen.type}
\alias{cpt.ar-method,pen.value}
\alias{cpt.ar-method,cpts<-}
\alias{cpt.ar-method,cpttype<-}
\alias{cpt.ar-method,data.set<-}
\alias{cpt.ar-method,test.stat<-}
\alias{cpt.ar-method,ncpts.max<-}
\alias{cpt.ar-method,method<-}
\alias{cpt.ar-method,minseglen<-}
\alias{cpt.ar-method,param.est<-}
\alias{cpt.ar-method,pen.type<-}
\alias{cpt.ar-method,pen.value<-}
\alias{cpt.ar-method,print}
\alias{cpt.ar-method,summary}
\alias{cpt.ar-method,plot}
\alias{cpt.ar-method,param}
\alias{cpt.ar-method,logLik}

\title{Class "cpt.ar"}
\description{
A class for changepoint objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cpt.ar", ...)}.
  \describe{
    \item{\code{new("cpt.ar", ...)}:}{creates a new object with class cpt.ar }
  }
}
\section{Slots}{
  \describe{
  \item{\code{BICvalues}:}{Object of class \code{"numeric"}, a coerced time series of the original data. }
  \item{\code{orders}:}{Object of class \code{"array"}, a coerced time series of the original data. }
    \item{\code{data.set}:}{Object of class \code{"matrix"}, a coerced time series of the original data. }
    \item{\code{cpttype}:}{Object of class \code{"character"}, the type of changepoint that was identified. }
    \item{\code{method}:}{Object of class \code{"character"}, the method that was used to search for changepoints. }
    \item{\code{test.stat}:}{Object of class \code{"character"}, the test statistic for the analysis of the data. }
    \item{\code{pen.type}:}{Object of class \code{"character"}, the penalty type specified in the analysis.}
    \item{\code{pen.value}:}{Object of class \code{"numeric"}, the value of the penalty used in the analysis.}
    \item{\code{minseglen}:}{Object of class \code{"numeric"}, the minimum segment length (no. of observations between changepoints) used in the analysis.}
    \item{\code{cpts}:}{Object of class \code{"numeric"}, vector of changepoints identified. }
    \item{\code{ncpts.max}:}{Object of class \code{"numeric"}, maximum number of changepoint that can be identified. }
    \item{\code{param.est}:}{Object of class \code{"list"}, list where each element is a vector of parameter estimates, if requested. }
    \item{\code{date}:}{Object of class \code{"character"}, date and time the changepoint analysis was run. }
        \item{\code{version}:}{Object of class \code{"character"}, version number of the package used when the analysis was run.}

  }
}
\section{Methods}{
  \describe{
    \item{BICvalues}{\code{signature(object = "cpt.ar")}: retrieves BIC values slot }
    \item{orders}{\code{signature(object = "cpt.ar")}: retrieves segment orders slot }
    \item{cpts}{\code{signature(object = "cpt.ar")}: retrieves cpts slot }
    \item{cpttype}{\code{signature(object = "cpt.ar")}: retrieves cpttype slot }
    \item{data.set}{\code{signature(object = "cpt.ar")}: retrieves vector version of data.set slot }
    \item{data.set.ts}{\code{signature(object = "cpt.ar")}: retrieves time series version of data.set slot }
    \item{test.stat}{\code{signature(object = "cpt.ar")}: retrieves test.stat slot }
    \item{ncpts.max}{\code{signature(object = "cpt.ar")}: retrieves ncpts.max slot }
    \item{method}{\code{signature(object = "cpt.ar")}: retrieves method slot }
    \item{minseglen}{\code{signature(object = "cpt.ar")}: retrieves minseglen slot }
    \item{param.est}{\code{signature(object = "cpt.ar")}: retrieves param.est slot }
    \item{pen.type}{\code{signature(object = "cpt.ar")}: retrieves pen.type slot }
    \item{pen.value}{\code{signature(object = "cpt.ar")}: retrieves pen.value slot }
    \item{BICvalues<-}{\code{signature(object = "cpt.ar")}: replaces BICvalues slot }
    \item{orders<-}{\code{signature(object = "cpt.ar")}: replaces orders slot }
    \item{cpts<-}{\code{signature(object = "cpt.ar")}: replaces cpts slot }
    \item{cpttype<-}{\code{signature(object = "cpt.ar")}: replaces cpttype slot }
    \item{data.set<-}{\code{signature(object = "cpt.ar")}: replaces data.set slot }
    \item{test.stat<-}{\code{signature(object = "cpt.ar")}: replaces test.stat slot }
    \item{ncpts.max<-}{\code{signature(object = "cpt.ar")}: replaces ncpts.max slot }
    \item{method<-}{\code{signature(object = "cpt.ar")}: replaces method slot }
    \item{minseglen<-}{\code{signature(object = "cpt.ar")}: replaces minseglen slot }
    \item{param.est<-}{\code{signature(object = "cpt.ar")}: replaces param.est slot }
    \item{pen.type<-}{\code{signature(object = "cpt.ar")}: replaces pen.type slot }
    \item{pen.value<-}{\code{signature(object = "cpt.ar")}: replaces pen.value slot }
    \item{print}{\code{signature(object = "cpt.ar")}: prints details of the cpt object including summary}
    \item{summary}{\code{signature(object = "cpt.ar")}: prints a summary of the cpt object }
    \item{plot}{\code{signature(object = "cpt.ar")}: plots the cpt object with changepoints highlighted}
    \item{param}{\code{signature(object = "cpt.ar")}: calculates the parameter estimates for the cpt object}
    \item{logLik}{\code{signature(object = "cpt.ar")}: returns the overall log-likelihood of the cpt object}
  }
}

\author{
Andrew Connell
}

\seealso{
\code{\link{data.set-methods}},\code{\link{cpts-methods}},\code{\link{cpt.ar}}
}
\examples{
showClass("cpt.ar") # shows the structure of the cpt class

x=new("cpt.ar") # creates a new object with the cpt class defaults
cpts(x) # retrieves the cpts slot from x
cpts(x)<-c(10,50,100) # replaces the cpts slot from x with c(10,50,100)

# Example of a change in variance at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,0,10))
#ans=cpt.ar(x) 
#print(ans) # prints details of the analysis including a summary
#summary(ans)
}
\keyword{classes}
\keyword{cpt}
\keyword{internal} 