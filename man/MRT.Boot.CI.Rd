\name{MRT.Boot.CI}
\alias{MRT.Boot.CI}
\title{
 Assessing the robustness of an univariate partition process by bootstrapping a given bagging fraction
}
\description{
Original data are bootstraped to ensure that the models produced robust detection of time series discontinuities. The function is identical to \code{MRT.Boot} except that  histogram of the bootstrapped values is not plotted in the lower panel of the graph. Instead, confident intervals appear directly on the data.
}
\usage{
MRT.Boot.CI(mrt.x, mrt.y, data, Size, Bagging = 0.5, N = 100, IC = 0.8)
}
\arguments{
  \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column name for the elemental concentration to be used. Data in this colum must be numeric.}
  \item{data}{ Dataframe to be analysed.}
  \item{Size}{
Size of the requested tree to be obtained (number og groups).
}
  \item{Bagging}{
The in-bag fraction [0-1]. Default is 0.5.
}
  \item{N}{
The number of boostraps. Default is 100.
}
  \item{IC}{
Confidence interval to be used for quantiles. Default is 0.8.
}
}

\value{
A list including the following elements:
\item{Mean}{Mean threshold values for each split, ordered sequentially.}
\item{Boot}{A matrix of dimension [\code{N},\code{Size}] containing threshold values for each split (columns) and each bootstrap (rows).}
\item{Qant}{A matrix containing qantiles for all splits.}
\item{OutSize}{The number of boostraps that did not produced tree of size \code{Size}.}
}

\note{
Quantiles are defined as [(1 - \code{IC})/2 ; 1-(1-\code{IC})/2)]. Default quantiles are [0.1 ; 0.9] using \code{IC}=0.8.
1-\code{OutSize} gives the true number of bootstraps used for means and qantiles calculation. If \code{OutSize}>0, \code{Boot} is of dimension [1-\code{OutSize},\code{Size}].
}

\seealso{
\code{\link{MRT.Boot}}
}
\examples{
data(univar)
model<-MRT.Boot.CI("Distance","Sr:Ca",univar,Size=4,Bagging=0.5,N=200)
}

\keyword{univariate}
\keyword{bootstrap}
