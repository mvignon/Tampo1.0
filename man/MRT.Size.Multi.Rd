\name{MRT.Size.Multi}
\alias{MRT.Size.Multi}
\title{
Multivariate regression tree with constrain on the group number
}
\description{
Partitioning (with a time contiguity constraint) the data using an arbitrary tree size number. The function uses \code{mvpart} as a basic function to persue the multivariate tree analyses (from the \pkg{mvpart} package). Multivariate expension of \code{MRT.Size}.
}
\usage{
MRT.Size.Multi(mrt.x, mrt.y, data, Size, Scale = FALSE, Graph = FALSE)
}
\arguments{
    \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column names for the elemental concentrations to be used. Data in these columns must be numeric. Columns are aggregated using \code{c()}.}
  \item{data}{ Dataframe to be analysed.}
  \item{Size}{Size of the requested tree to be obtained (number og groups).}
  \item{Scale}{Do concentrations need to be scaled?}
  \item{Graph}{
Optional supplementary graphs depicting apparent and relative (cross-validated) R-squares according to tree sizes.}
}
\details{
\bold{Scaling} \cr
Elements can appear with different scales so that one of them can appear extremely squashed (the plot is automatically scaled to contain all elements. Because each split minimizes the total sums of squares of the distances of measurements from the centroids of the groups to which they belong, scaling the data (using the basic \code{Scale} argument function) would result in changing the times series discontinuities. If data are not scaled, it is likely that the discontinuities obtained by regression tree correlate changes in the variable having the largest scale variation, irrespective of other variables. While it is possible to perform both analyses, scaling the data appears as the more reasonable choice since in a discriminant context users are interested in relative variations in different element, not in their absolute values. In doing so, scaled values plotted in the graphic cannot be interpreted directly. However, the \code{Mean} and \code{SD} extracted from the
function are automatically corrected in order to be readable at the appropriate scale.

\bold{Pruning} \cr
Omitting the \code{Size} argument will result in pruning the lower branches on the basis of v-fold cross-validation error (the selected tree size minimizes the cross-validated error based on the '1se' rule, see \code{\link{mvpart}}).
}
\value{
A list including the following elements:
\item{Change}{Values of internal splits (as appearing in the consol).}
\item{Mean}{Means for each group (groups are organized chronologically, from 1 to \code{Size}).}
\item{SD}{Standard deviation for each group (still organized chronologically).}
\item{Model}{An object of class \code{rpart} created by the \code{mvpart} function from the \pkg{mvpart} package.}
\item{Belong}{The group belonging for each initial observations (same length as \code{mrt.x}).}
\item{N}{The sequential number of observations in each group.}
\item{Predict}{The predicted values (stepwise constant model).}
\item{Size}{The real tree size created by the \code{mvpart} function.}

\item{outpout}{Text for \code{summary.MRT}}
\item{X}{Original \code{mrt.x} values}
\item{Y}{Original \code{mrt.y} values}
}



\seealso{
\code{\link{MRT.Size}}
}
\examples{
data(multivar)
# For unscaled elements
model<-MRT.Size.Multi ("Distance", c("Mg24","Zn66"),multivar,Size=4,Scale=FALSE)

# For scaled elements
model2<-MRT.Size.Multi ("Distance", c("Mg24","Zn66"),multivar,Size=4,Scale=TRUE)

# The raw data correspond to a transect passing through the primordium and starting
# 200 micrometer before. Because data are supposed to be roughly symmetrical apart 
# from the primordium, we decided to start the analysis at only -100 to have enough
# representation of the maternal effect.
model3<-MRT.Size.Multi ("Distance",c("Mg24","Mn55","Sr88"),
multivar[which(multivar$Distance>-100),],Size=4,Scale=TRUE)
plot(model3$Model)
text(model3$Model)

# Identical with unspecified tree size
model4<-MRT.Size.Multi ("Distance",c("Mg24","Mn55","Sr88"),
multivar[which(multivar$Distance>-100),],Scale=TRUE)
}

\keyword{multivariate}