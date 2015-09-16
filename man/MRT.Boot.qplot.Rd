\name{MRT.Boot.qplot}
\alias{MRT.Boot.qplot}
\title{
Assessing the partition congruence between elements.
}
\description{
The function identify and quantify which elements alone are intrinsically clustered as in the multivariate analysis, and which one exhibit the same clustering patterns. This is the core basis for \code{\link{plotMRT}}, \code{\link{Mapp}}, \code{\link{Dendro2}}, \code{\link{Dendro3}} and \code{\link{MMCompare}}.

}
\usage{
MRT.Boot.qplot(mrt.x, mrt.y, data, Size, Bagging = 0.5, N = 100, IC = 0.8, Scale = FALSE)
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
The result from the multivariate analysis appears as a clustering compromise that maximizes the homogeneity within each group for all the elements. The function have been created to specifically assess the congruence between elements. To do this, the bootstrap process is applied independently on 1) the overall dataset containing all the elements and 2) on each element separately.
}
\value{
Values are intended to be analyzed/plotted using either \code{\link{plotMRT}}, \code{\link{Mapp}}, \code{\link{Dendro2}}, \code{\link{Dendro3}} or \code{\link{MMCompare}}.

A list including the following elements:
\item{BootMulti}{A matrix containing threshold values for each split (columns) and each bootstrap (rows) using multivariate data.}
\item{BootMono}{A list of matrix containing threshold values for each split (columns) and each bootstrap (rows) using univariate data (one matrix per element analyzed).}
\item{Outsize}{The number of boostraps that did not produced tree of size \code{Size} using multivariate data.}
\item{N}{Remind the requested number of bootstrap}
\item{Data.Hist}{A dataframe containing all thresholds values (from both univariates and multivariate analysis) reorganized by values, Elements and Splits (respectively in columns named "Dat", "Element" and "Split".}
}


\seealso{
\code{\link{plotMRT}}, \code{\link{Mapp}}, \code{\link{Dendro2}}, \code{\link{Dendro3}}, \code{\link{MMCompare}}.
}

\examples{
# Simple comparison of two isotopes using default arguments
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Sr86","Sr88"),
multivar[which(multivar$Distance>-100),], Size=3,Scale=TRUE,N=50)
plotMRT(model)
}


\keyword{multivariate}
\keyword{univariate}
\keyword{comparison}
\keyword{bootstrap}
