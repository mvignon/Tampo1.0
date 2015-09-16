\name{Dendro3}
\alias{Dendro3}
\title{
Dendrogram for comparing elemental overlaps
}
\description{
The function summaries the distance between elements (over all splits) in one dimension using a complete linkage based on the similarity matrix computed by \code{MRT.Heatmap}. 
}
\usage{
Dendro3(model, Index = "S", IC = 0.9, Mean = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{model}{Results from \code{MRT.Boot.qplot}.}
\item{Index}{The index to be used ("S", "SS" or "J"). Default is "S".}
\item{IC}{The confident interval to be used. Default is 0.9.}
\item{Mean}{Should the mean value be used for asymmetrical matrix.}
}

\value{
An object of class \code{gg} created from the \pkg{ggplot2} package.
}


\seealso{
\code{\link{Dendro2}}, \code{\link{MRT.Heatmap}}, \code{\link{MRT.Boot.qplot}}
}
\examples{
# Simple comparison using default arguments
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Mg24","Mn55","Sr88"),
multivar[which(multivar$Distance>-100),], Size=4,Scale=TRUE,N=25)
Dendro3(model)

# Comparing more elements
# X<-c(colnames(multivar[,-1]))
# model2<-MRT.Boot.qplot("Distance",X,
# multivar[which(multivar$Distance>-100),], Size=4,Scale=TRUE,N=100)
# Dendro3(model2,Index="J")
}

\keyword{multivariate}
\keyword{univariate}
\keyword{comparison}
