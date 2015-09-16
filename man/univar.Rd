\name{univar}
\alias{univar}
\docType{data}
\title{
Univariate data
}
\description{
Elemental sequence (Sr:Ca) obtained using LA-ICPMS from one otolith of \emph{Salmo trutta} in the Kerguelen Islands.
}
\usage{data(univar)}
\format{
  A data frame with 688 observations on the following 2 variables.
  \describe{
    \item{\code{Distance}}{Distance to the nucleus (micrometers).}
    \item{\code{Sr:Ca}}{Sr:Ca ratio}
  }
}

\examples{
data(univar)
plot(univar[1:688,])
}
\keyword{datasets}
\keyword{univariate}
