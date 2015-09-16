\name{Mapp}
\alias{Mapp}
\title{
Complete investigation of all paiwise comparisons among elements and/or splits overlaps
}
\description{
Quantitative analysis of results from \code{MRT.Boot.qplot} using either "S", "SS" or "J" index as computed by \code{MRT.Heatmap}. Results are visualized using a colour gradient.
}
\usage{
Mapp(model, Index = "S", what = "Total", IC = 0.9, Values = T, write = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Results from \code{MRT.Boot.qplot}.}
 \item{Index}{The index to be used ("S", "SS" or "J"). Default is "S".}
  \item{what}{What to compare in panels. Either "Total", "Split" or "Element".}
  \item{IC}{The confident interval to be used. Default is 0.9.}
  \item{Values}{Plot index values on the map.}
  \item{write}{To write results.}
}
  
\value{
\enumerate{
\item{}An object of class \code{gg} created from the \pkg{ggplot2} package.
\item{}An optional tab-delimited file (.txt) named "Index" within the current directory.
}
}


\seealso{
\code{\link{MRT.Heatmap}}, \code{\link{MRT.Boot.qplot}}
}

\examples{
# Simple comparison of two isotopes using default arguments
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Sr86","Sr88"),
multivar[which(multivar$Distance>-100),], Size=3,Scale=TRUE,N=50)
Mapp(model)

# Changing index
# Mapp(model, Index="J")
# Comparing splits (one panel per element)
# Mapp(model, what="Split")
# Comparing elements (one panel per split)
# Mapp(model, what="Element")

# Map for more splits and all elements
# X<-c(colnames(multivar[,-1]))
# model2<-MRT.Boot.qplot("Distance",X,
# multivar[which(multivar$Distance>-100),], Size=10,Scale=TRUE,N=100)
# Mapp(model2)
}

% R documentation directory.
\keyword{multivariate}
\keyword{univariate}
\keyword{comparison}
