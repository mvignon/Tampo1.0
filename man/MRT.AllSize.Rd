\name{MRT.AllSize}
\alias{MRT.AllSize}
\title{Univariate regression tree using a whole range of constrain on the group number}
\description{
Identical to \code{MRTR.Size} except that data are partitioned using a full range of sizes and results are sequentially printed in a PDF file.}
\usage{
MRT.AllSize(mrt.x, mrt.y, data, Size = 10, GraphName = "")
}
\arguments{
   \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column name for the elemental concentration to be used. Data in this colum must be numeric.}
  \item{data}{ Dataframe to be analysed.}
  \item{Size}{Size of the requested tree to be obtained (number og groups).}
  \item{GraphName}{Name of the PDF file to be created}
}

\details{
Instead of fitting separately models of various sizes, the \code{MRT.AllSize} function may be useful to automatically explore the whole range of sizes from 2 to \code{Size}. In an exploratory context, it is therefore recommended to use substantial range. One specificity of this function relies on the fact that models (both tree graphics and categorized data) are not displayed in R but are instead 'printed' in a PDF file. Remember that, by default, the PDF file will be created in the working directory and will be named "MRT.AllSizeGraphics". The created file consists in tree and categorized graphics alternatively. The last slide of the PDF file typically exhibits the R-square and Relative Error versus the number of splits.
}
\value{
Print a PDF file within the current working directory.
\item{Cptable}{The table of optimal prunings based on a complexity parameter. See \code{\link{rpart.object}}.}
\item{FullTree}{An object of class \code{rpart} created by the \code{mvpart} function from the \pkg{mvpart} package and creating a tree of size \code{Size}.}
}


\seealso{
\code{\link{MRT.Size}}
}
\examples{
data(univar)
#MRT.AllSize("Distance", "Sr:Ca", univar, Size=20)
#to close the PDF file (if not automatically)
#dev.off() 
}

\keyword{univariate}
