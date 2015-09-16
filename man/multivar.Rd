\name{multivar}
\alias{multivar}
\docType{data}
\title{
Multivariate data
}
\description{
Elemental/isotopic sequences obtained using LA-ICPMS from one otolith of \emph{Salmo trutta} in the Kerguelen Islands.
}
\usage{data(multivar)}
\format{
  A data frame with 435 observations on the following 11 variables.
  \describe{
     \item{\code{Distance}}{Distance to the nucleus (micrometers).}
    \item{\code{Mg24}}{Mg24:Ca ratio}
    \item{\code{Mg25}}{Mg25:Ca ratio}
    \item{\code{Mn55}}{Mn55:Ca ratio}
    \item{\code{Zn66}}{Zn66:Ca ratio}
    \item{\code{Zn68}}{Zn68:Ca ratio}
    \item{\code{Rb85}}{Rb85:Ca ratio}
    \item{\code{Sr86}}{SR86:Ca ratio}
    \item{\code{Sr88}}{Sr88:Ca ratio}
    \item{\code{Ba135}}{Ba135:Ca ratio}
    \item{\code{Ba138}}{Ba138:Ca ratio}
  }
}

\examples{
data(multivar)
matplot(multivar$Distance,scale(multivar[2:11]),type="l",lty="solid")
pairs(multivar[2:11])

}
\keyword{datasets}
\keyword{multivariate}
