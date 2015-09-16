\name{MRT.All}
\alias{MRT.All}
\title{
Univariate regression tree with sequential constrain on the group number
}
\description{
The function is a manual alternative of \code{MRT.AllSize} that do not create
PDF file but instead show the sequential results.
}
\usage{
MRT.All(mrt.x, mrt.y, data)
}
\arguments{
\item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column name for the elemental concentration to be used. Data in this colum must be numeric.}
  \item{data}{ Dataframe to be analysed.}
}
\details{
 The first graph starts with two splits (three groups) and the user can repeatedly add other splits by pressing the <enter> key. There is practically no upper limit but, because the number of observation in terminal nodes tends to decrease in a monotonic manner in relation to the tree size, it is best to keep a reasonable number of splits. When the partition process is ended (by pressing the <N> key), the R-square and relative error are shown.
}
\value{
Display trees of increasing complexity.
\item{Cptable}{The table of optimal prunings based on a complexity parameter. See \code{\link{rpart.object}}.}
\item{FullTree}{An object of class \code{rpart} created by the \code{mvpart} function from the \pkg{mvpart} package and creating a tree of size \code{Size}.
}
}



\seealso{
\code{\link{MRT.Size}}
}
\examples{
data(univar)
#MRT.All("Distance", "Sr:Ca", univar, Size=50)
#  Do you want to split again? Press <enter> to continue or N to stop:
#  Splitting continue until <N> is typed.
}
\keyword{univariate}
