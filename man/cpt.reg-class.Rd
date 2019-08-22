\name{cpt.reg-class}
\Rdversion{1.1}
\docType{class}
\alias{cpt.reg-class}
\alias{cpt.reg-method,cpts}
\alias{cpt.reg-method,cpttype}
\alias{cpt.reg-method,data.set}
\alias{cpt.reg-method,data.set.ts}
\alias{cpt.reg-method,test.stat}
\alias{cpt.reg-method,ncpts.max}
\alias{cpt.reg-method,method}
\alias{cpt.reg-method,minseglen}
\alias{cpt.reg-method,param.est}
\alias{cpt.reg-method,pen.type}
\alias{cpt.reg-method,pen.value}
\alias{cpt.reg-method,cpts<-}
\alias{cpt.reg-method,cpttype<-}
\alias{cpt.reg-method,data.set<-}
\alias{cpt.reg-method,test.stat<-}
\alias{cpt.reg-method,ncpts.max<-}
\alias{cpt.reg-method,method<-}
\alias{cpt.reg-method,minseglen<-}
\alias{cpt.reg-method,param.est<-}
\alias{cpt.reg-method,pen.type<-}
\alias{cpt.reg-method,pen.value<-}
\alias{cpt.reg-method,print}
\alias{cpt.reg-method,summary}
\alias{cpt.reg-method,plot}
\alias{cpt.reg-method,param}
\alias{cpt.reg-method,logLik}

\title{Class "cpt.reg"}
\description{
A class for changepoint objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cpt.reg", ...)}.
  \describe{
    \item{\code{new("cpt.reg", ...)}:}{creates a new object with class cpt.reg }
  }
}
\section{Slots}{
  \describe{
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
    \item{cpts}{\code{signature(object = "cpt.reg")}: retrieves cpts slot }
    \item{cpttype}{\code{signature(object = "cpt.reg")}: retrieves cpttype slot }
    \item{data.set}{\code{signature(object = "cpt.reg")}: retrieves vector version of data.set slot }
    \item{data.set.ts}{\code{signature(object = "cpt.reg")}: retrieves time series version of data.set slot }
    \item{test.stat}{\code{signature(object = "cpt.reg")}: retrieves test.stat slot }
    \item{ncpts.max}{\code{signature(object = "cpt.reg")}: retrieves ncpts.max slot }
    \item{method}{\code{signature(object = "cpt.reg")}: retrieves method slot }
    \item{minseglen}{\code{signature(object = "cpt.reg")}: retrieves minseglen slot }
    \item{param.est}{\code{signature(object = "cpt.reg")}: retrieves param.est slot }
    \item{pen.type}{\code{signature(object = "cpt.reg")}: retrieves pen.type slot }
    \item{pen.value}{\code{signature(object = "cpt.reg")}: retrieves pen.value slot }
    \item{cpts<-}{\code{signature(object = "cpt.reg")}: replaces cpts slot }
    \item{cpttype<-}{\code{signature(object = "cpt.reg")}: replaces cpttype slot }
    \item{data.set<-}{\code{signature(object = "cpt.reg")}: replaces data.set slot }
    \item{test.stat<-}{\code{signature(object = "cpt.reg")}: replaces test.stat slot }
    \item{ncpts.max<-}{\code{signature(object = "cpt.reg")}: replaces ncpts.max slot }
    \item{method<-}{\code{signature(object = "cpt.reg")}: replaces method slot }
    \item{minseglen<-}{\code{signature(object = "cpt.reg")}: replaces minseglen slot }
    \item{param.est<-}{\code{signature(object = "cpt.reg")}: replaces param.est slot }
    \item{pen.type<-}{\code{signature(object = "cpt.reg")}: replaces pen.type slot }
    \item{pen.value<-}{\code{signature(object = "cpt.reg")}: replaces pen.value slot }
    \item{print}{\code{signature(object = "cpt.reg")}: prints details of the cpt object including summary}
    \item{summary}{\code{signature(object = "cpt.reg")}: prints a summary of the cpt object }
    \item{plot}{\code{signature(object = "cpt.reg")}: plots the cpt object with changepoints highlighted}
    \item{param}{\code{signature(object = "cpt.reg")}: calculates the parameter estimates for the cpt object}
    \item{logLik}{\code{signature(object = "cpt.reg")}: returns the overall log-likelihood of the cpt object}
  }
}

\author{
Rebecca Killick
}

\seealso{
\code{\link{cpts-methods}},\code{\link{cpt.reg}}
}
\examples{
showClass("cpt.reg") # shows the structure of the cpt class

x=new("cpt.reg") # creates a new object with the cpt class defaults
cpts(x) # retrieves the cpts slot from x
cpts(x)<-c(10,50,100) # replaces the cpts slot from x with c(10,50,100)

# Example of a change in variance at 100 in simulated normal data
set.seed(1)
x=design(c(rnorm(100,0,1),rnorm(100,0,10)),1)
ans=cpt.reg(x) 
print(ans) # prints details of the analysis including a summary
summary(ans)
plot(ans) # plots the data with change (vertical line) at 100
logLik(ans) # raw likelihood of the data with changepoints, second value is likelihood + penalty
}
\keyword{classes}
\keyword{cpt}
\keyword{internal} 