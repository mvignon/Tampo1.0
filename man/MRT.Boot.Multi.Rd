\name{MRT.Boot.Multi}
\alias{MRT.Boot.Multi}
\title{
Assessing the robustness of a multivariate partition process by bootstrapping a given bagging fraction
}
\description{
Multivariate expansion of \code{MRT.Boot}. Original data are bootstraped to ensure that the models produced robust detection of time series discontinuities.
}
\usage{
MRT.Boot.Multi(mrt.x, mrt.y, data, Size, Bagging = 0.5, N = 100, IC = 0.8, Scale = FALSE)
}
\arguments{
   \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column names for the elemental concentrations to be used. Data in these columns must be numeric. Columns are aggregated using \code{c()}.}
  \item{data}{ Dataframe to be analysed.}
  \item{Size}{Size of the requested tree to be obtained (number og groups).}
  \item{Bagging}{
The in-bag fraction [0-1]. Default is 0.5.
}
  \item{N}{The number of boostraps. Default is 100.}
  \item{IC}{Confidence interval to be used for quantiles. Default is 0.8.}
\item{Scale}{Do concentrations need to be scaled?}
}
\details{
The exact threshold value selected at each node is intrinsically case-sensitive
(instability to small perturbations of the data). This is of primarily importance because elemental ratios extracted from the same otolith can slightly differs according to the transect position or the machine laser resolution. To examine properties of the groups formed by chronological clustering and ensure that the models produced robust detection of time series discontinuities, the measurements were further resampled for model building using a given ("in-bag") fraction of the original data (caracterized by the \code{Bagging} argument). Characteristics of a large ensemble of MRT models provide a more accurate characterization of time series discontinuities. From the ensemble, the function return the number of time a given node appears in a particular zone. This allows quantifying the robustness of structural changes in the data. The associated values (mean and quantiles) are recorded, and an histogram of the bootstrapped values is plotted in the lower panel of the graph.
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
\code{\link{MRT.Size.Multi}}
}
\examples{
data(multivar)
model<-MRT.Boot.Multi("Distance", c("Mg24","Zn66"),multivar,Size=4,Scale=TRUE,N=100)


}

\keyword{multivariate}
\keyword{bootstrap}
