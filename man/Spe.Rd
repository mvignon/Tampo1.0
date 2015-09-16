\name{Spe}
\alias{Spe}
\title{
Specificity of splitting regions
}
\description{
The function plot the splitting range (\code{IC} only) for each element (univariate analysis) through the different splits, as well as for the multivariate analysis. This helps distinguishing non localized splitting region occurring through the entire range of values, from the nucleus to the otolith edge. This avoids some bias when using a high number of splits that results in widely unlocalized splitting regions.
}
\usage{
Spe(model, IC = 0.9)
}
\arguments{
\item{model}{Results from \code{MRT.Boot.qplot}.}
\item{IC}{The confident interval to be used. Default is 0.9.}
}
\details{
\enumerate{
\item{First panel}{ Splitting ranges for multivariate analysis ("All") and for each element from the \code{model}.}
\item{Second panel}{ Splitting ranges (both uni- and multivariate) through the different splits.}
}
}

\value{
Multiple objects of class \code{gg} created from the \pkg{ggplot2} package.
}

\note{
Ranges are in same units than \code{mrt.x} from the \code{MRT.Boot.qplot} function.
}

\seealso{
\code{\link{MRT.Heatmap}}, \code{\link{MRT.Boot.qplot}}
}
\examples{

# Range comparison of two isotopes through 3 splits
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Sr86","Sr88"),
multivar[which(multivar$Distance>-100),], Size=4,Scale=TRUE,N=50)
Spe(model)

# Range comparison of three elements through 6 splits
# model2<-MRT.Boot.qplot("Distance",c("Mg24","Mn55","Zn68"),
# multivar[which(multivar$Distance>-100),], Size=7,Scale=TRUE,N=100)
# Spe(model2)
}

\keyword{multivariate}
\keyword{univariate}
\keyword{specificity}