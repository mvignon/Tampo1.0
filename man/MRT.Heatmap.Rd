\name{MRT.Heatmap}
\alias{MRT.Heatmap}

\title{
Compute comparison indices based on unions and intersections of the splitting regions among elements and/or splits
}
\description{
Quantitative analysis of results from \code{MRT.Boot.qplot}. Internal function used by \code{\link{Mapp}}, \code{\link{Dendro2}}, \code{\link{Dendro3}} and \code{\link{MMCompare}}. Not to be used alone !
}
\usage{
MRT.Heatmap(model, IC = 0.9, Index = "S")
}
\arguments{
  \item{model}{Results from \code{MRT.Boot.qplot}.}
  \item{IC}{Confidence interval to be used for quantiles. Default is 0.9.}
  \item{Index}{The index to be used ("S", "SS" or "J"). Default is "S".}
}
\details{
Three main indexes are available for computation:
\enumerate{
\item{\emph{S}}{ : \emph{Sij} is the simple ratio between the intersection of the two regions \emph{i}-\emph{j} and the specific area of \emph{i}.
Please note that the index is asymmetrical (i.e. \emph{Sij} is different from \emph{Sji}).}
\item{\emph{SS}}{ : \emph{SSij} is a weighted version of \emph{Sij} that take into consideration the width of the \emph{i} region using (\emph{Tot}-\emph{i})/\emph{Tot} as the weigth coefficient, \emph{Tot} being the maximum width allowed in the data set (for a non localized splitting region occurring through the entire range of values, from the nucleus to the otolith edge). This avoids some bias when using a high number of splits that results in widely unlocalized splitting regions.}
\item{\emph{J}}{ : \emph{Jij} corresponds to the classical Jaccard index (ration between the intersection and the union of two regions \emph{i}-\emph{j} ). This index is symmetrical (i.e. \emph{Jij} equal \emph{Jji}).}
}
}

\value{
A list including the following elements:
\item{Index}{A square matrix containing index for univariate analyses. Values are sorted by split and then by element. The matrix is symmetrical only when using Index="J".}
\item{heat}{A dataframe containing all values reorganized in rows. To be used by \code{Mapp}.}
\item{method}{Remind the Index used for the analysis}
\item{Clust}{A dataframe listing splits and elements separately. To be used by \code{Mapp}, \code{Dendro2} and \code{Dendro3}}
\item{Multi}{A list of object of class \code{interval} corresponding to splitting regions from multivariate bootstrap and according to \code{IC}. Intervals are ordered sequentially by split.}
}


\seealso{
\code{\link{MRT.Boot.qplot}}, code{\link{Mapp}}, \code{\link{Dendro2}}, \code{\link{Dendro3}}, \code{\link{MMCompare}}.}

\examples{
# Simple comparison of two isotopes using default arguments
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Sr86","Sr88"),
multivar[which(multivar$Distance>-100),], Size=3,Scale=TRUE,N=50)
results<-MRT.Heatmap(model)
names(results)
}

\keyword{multivariate}
\keyword{univariate}
\keyword{comparison}
