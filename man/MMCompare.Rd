\name{MMCompare}
\alias{MMCompare}
\title{
Congruency between mono- and multi-elemental signals
}
\description{
The function quantify which elements alone (univariate analysis) are intrinsically clustered as in the multivariate analysis. To do so, one of the three indexes (see \code{MRT.Heatmap} for details) is calculated for each element between its own splitting regions and those from the multivariate dataset. Comparison is performed regardless of split numbers, over all splits.
}
\usage{
MMCompare(model, IC = 0.9, Index = "S")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{model}{Results from \code{MRT.Boot.qplot}.}
\item{IC}{The confident interval to be used. Default is 0.9.}
\item{Index}{The index to be used ("S", "SS" or "J"). Default is "S".}
}

\value{
\item{Over1}{A matrix containing the index values \emph{Indexmono-multi} for each element.}
\item{Over1}{An optional matrix containing the index values \emph{Indexmulti-mono} for each element. Only for asymmetrical indexex (i.e. \emph{S} and \emph{SS}).}
}

\note{
By default, the "S"" index is used. Because of its asymmetrical nature, two complementary values are reported (i.e. \emph{Smono-multi} and \emph{Smulti-mono}).
}


\seealso{
\code{\link{Dendro3}}, \code{\link{MRT.Heatmap}}, \code{\link{MRT.Boot.qplot}}
}

\examples{
# Simple comparison of two isotopes
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Mg24","Mg25"),
multivar[which(multivar$Distance>-100),], Size=3,Scale=TRUE,N=50)
MMCompare(model,Index="J")

# More exhaustive comparison omong elements
# X<-c(colnames(multivar[,-1]))
# model2<-MRT.Boot.qplot("Distance",X,
# multivar[which(multivar$Distance>-100),], Size=4,Scale=TRUE,N=100)
# MMCompare(model2,Index="J")
# MMCompare(model2)

}

\keyword{multivariate}
\keyword{univariate}
\keyword{comparison}
