\name{plotMRT}
\alias{plotMRT}
\title{
Plot results from \code{MRT.Boot.qplot} using ggplot2.
}
\description{
A graphical description of bootstrap results from \code{MRT.Boot.qplot}. Allow assessing the congruency between splits and elements.
}
\usage{
plotMRT(aa, type = "density", by = "Split", col = 1, palette = "Set1", scale = "free")
}
\arguments{
  \item{aa}{Results from \code{MRT.Boot.qplot}.}
  \item{type}{Type of graphic to be plotted. "density"" for smooth density plot or "count" for raw data.}
  \item{by}{Plots can be organized by "Split" or "Element". "Total" also available when using type="count" and display results irrespective of split numbers.}
  \item{col}{Numeric. Number of columns to display facets.}
  \item{palette}{Palette colors for splits or elements. "Blues" "Set2" "Set3" "Dark2" "Pastel1" "Pastel2" "Paired" and "Accent". Default is "Set1". See ggplot2 manual for more palettes.}
  \item{scale}{Are scales shared across all facets ("fixed"), or do they vary across rows ("free_x"), columns ("free_y"), or both rows and columns ("free").}
}

\value{An object of class \code{gg} created from the \pkg{ggplot2} package.}


\seealso{
\code{\link{MRT.Boot.qplot}}
}
\examples{
# Simple comparison of two isotopes using default arguments
data(multivar)
model<-MRT.Boot.qplot("Distance",c("Sr86","Sr88"),
multivar[which(multivar$Distance>-100),], Size=3,Scale=TRUE,N=50)
plotMRT(model)

# Changing graphics
# - Count data instead of density plot
# plotMRT(model,type="count")
# - Counts irrespective of split number
# plotMRT(model,type="count",by="Total")
# - Density plots organized by elements
# plotMRT(model,type="density",by="Element")
# - Identical using 2 columns
# plotMRT(model,type="density",by="Element",col=2)

# Axis modification
# - Free x axis and fixed y axis
# plotMRT(model,scale="free_x")
# - Free y axis and fixed x axis
# plotMRT(model,scale="free_y")
# - Fixed axis
# plotMRT(model,scale="fixed")
}

\keyword{multivariate}
\keyword{univariate}
\keyword{comparison}
