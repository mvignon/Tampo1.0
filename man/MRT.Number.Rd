\name{MRT.Number}
\alias{MRT.Number}
\title{Univariate regression tree with constrain on the number of measures in groups.
}
\description{
The function try to generate a tree (whose size depend on the \code{Size} argument) by fixing a minimum number of consecutive measures within each cluster generated through the partition process.
}
\usage{
MRT.Number(mrt.x, mrt.y, data, Size, Number = 100, Graph = FALSE, 
           Habitat1 = FALSE, Habitat2 = FALSE)
}
\arguments{
 \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column name for the elemental concentration to be used. Data in this colum must be numeric.}
  \item{data}{ Dataframe to be analysed.}
  \item{Size}{Size of the requested tree to be obtained (number og groups). May be limited accounting for the \code{Number} argument.}
  \item{Number}{
Minimum number of observations requested in a group to be created.
}
  \item{Graph}{
Optional supplementary graphs depicting apparent and relative (cross-validated) R-squares according to tree sizes.}
\item{Habitat1}{Plot habitat use according to the given threshold.}
\item{Habitat2}{Plot habitat use according to the two given thresholds.}
}
\details{
This process removes short clusters corresponding to biologically irrelevant/uninterpretable information. Because the number of consecutive  measures in terminal nodes tends to decrease in a monotonic manner in relation to the tree size, the biggest tree that fully met the conditions is returned. If the provided minimum number of consecutive measures within each cluster is not limitative, then, the function creates a tree of size \code{Size}.
}
\value{
A list including the following elements:
\item{Change}{Values of internal splits (as appearing in the consol).}
\item{Model}{An object of class \code{rpart} created by the \code{mvpart} function from the \pkg{mvpart} package. The corresponding tree accounts for both \code{Number} and \code{Size}}
\item{Mean}{Means for each group, organized chronologically.}
\item{SD}{Standard deviation for each group, organized chronologically.}
\item{Belong}{The group belonging for each initial observations (same length as \code{mrt.x}).}
\item{N}{The sequential number of observations in each group.}
\item{Predict}{The predicted values (stepwise constant model).}
\item{Size}{The real tree size created by the \code{mvpart} function.}   
  }


\note{
Omitting the \code{Size} argument will result in pruning the lower branches on the basis of v-fold cross-validation error (the selected tree size minimizes the cross-validated error based on the '1se' rule, see \code{\link{mvpart}}).
}


\seealso{
\code{\link{MRT.Size}}, \code{\link{MRT.2C}}, \code{\link{MRT.Mean}}
}
\examples{
data(univar)
# Unconstrained model with small groups
model <- MRT.Size("Distance", "Sr:Ca", univar,Size=20)
model$Size # 20 groups
model$N # A lot of groups have very few measures

# Constrained model with groups having at least 40 measures
model2 <- MRT.Number("Distance", "Sr:Ca", univar,Number=40,Size=20)
model2$Size # Number of groups filling the "number" argument 
model2$N # The "number" argument was respected.

# 'Number' constrain with automatic optimal pruning
model3 <- MRT.Number("Distance", "Sr:Ca", univar,Number=40)

# Other 'Number' constrains with automatic optimal pruning
model4 <- MRT.Number("Distance", "Sr:Ca", univar,Number=100)
model5 <- MRT.Number("Distance", "Sr:Ca", univar,Number=150)

# Habitat(s) can also be viewed
MRT.Number("Distance", "Sr:Ca", univar,Number=40,Habitat1=3)
MRT.Number("Distance", "Sr:Ca", univar,Number=40,Habitat2=c(2.5,4))
}

\keyword{univariate}
\keyword{constrain}