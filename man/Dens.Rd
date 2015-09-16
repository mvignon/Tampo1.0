\name{Dens}
\alias{Dens}
\title{
Relative density probability
}
\description{
Description of any distribution based on relative density probability. The function can be used alone for descriptive purpose but was written for input in the \code{Dens2Hab} function. 
}
\usage{
Dens(x, IC = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector containing concentrations to be described} 
  \item{IC}{Confidence interval to be used}
}

\value{
A list including the following elements:
\item{Fun}{Extrapolated density function.}
\item{NFun}{Extrapolated relative [0-1] density function.}
\item{Min}{Min value in \code{x}.}
\item{Max}{Min value in \code{x}.}
\item{Mode}{Mode for the highest density value (standardized to 1).}
\item{LP}{Left side relative density for the provided \code{IC}.}
\item{RP}{Rigth side relative density for the provided \code{IC}.}
\item{IC}{Remind the confidence interval used.}
}

\note{
"Relative" refers to the fact that the highest density probability is arbitrarily set to 1.
}
\seealso{
\code{\link{Dens2Hab}}
}
\examples{
# Random data
x<-c(rnorm(500),rnorm(200,mean=2,sd=2))
# Density description
D<-Dens(x)
# Left side density threshold
D$LP
# Use extrapolation function
plot(D$Fun(seq(from=0,to=5,by=0.1)),type="l")
}

\keyword{Density}
\keyword{univariate}
