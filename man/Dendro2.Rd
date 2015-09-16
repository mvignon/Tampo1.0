\name{Dendro2}
\alias{Dendro2}
\title{
Dendrogram for comparing element and split overlaps
}
\description{
The function summaries the distance between splits/elements in one dimension using a complete linkage based on the similarity matrix as depicted by the \code{Map} function and computed by \code{MRT.Heatmap}. 
}
\usage{
Dendro2(model, Index = "S", IC = 0.9, Mean = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{model}{Results from \code{MRT.Boot.qplot}.}
\item{Index}{The index to be used ("S", "SS" or "J"). Default is "S".}
\item{IC}{The confident interval to be used. Default is 0.9.}
\item{Mean}{Should the mean value be used for asymmetrical matrix.}
}
\details{
For practical reasons, inputting a non-symmetric dissimilarity matrix ("S" and "SS" only) in most algorithm would probably generate irrelevant dendrogram. As a consequence, by default, only the lower triangle of the overlapping matrix is used, the rest is ignored. This has no specific effect when using the symmetrical Jaccard distance but caution must be paid when using asymmetrical \emph{Sij} and \emph{SSij} indexes (if a measure is asymmetric, then it no longer fits the definition of a distance). Alternatively, turning "TRUE" the \code{Mean} argument will create a dendrogram based on the mean distance value from the asymmetrical overlapping matrix.
Please also note that the associated cophenetic coefficient correlation will automatically appear in the console.
}
\value{
Multiple objects of class \code{gg} created from the \pkg{ggplot2} package.
}


\seealso{
\code{\link{Dendro3}}, \code{\link{MRT.Heatmap}}, \code{\link{MRT.Boot.qplot}}
}
\examples{
# Simple comparison of two isotopes using default arguments
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Sr86","Sr88"),
multivar[which(multivar$Distance>-100),], Size=4,Scale=TRUE,N=30)
Dendro2(model)
}

\keyword{multivariate}
\keyword{univariate}
\keyword{comparison}
