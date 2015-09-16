\name{MRT.2C}
\alias{MRT.2C}
\title{
Univariate regression tree with double constrain on the difference between the means of two consecutive clusters and the minimum number of consecutive measures within each cluster.
}
\description{
The function integrates the two constraints from \code{\link{MRT.Number}} and \code{\link{MRT.Mean}} simultaneously. It allows automatically estimating the pruning characteristics that generate the number of groups fitting these criteria.
}
\usage{
MRT.2C(mrt.x, mrt.y, data, Cond = 1, Number = 100, Size, Graph = FALSE,
       Habitat1 = FALSE, Habitat2 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column name for the elemental concentration to be used. Data in this colum must be numeric.}
  \item{data}{ Dataframe to be analysed.}
  \item{Cond}{The minimal difference requested between the mean elemental concentrations of two consecutive clusters. Same units that \code{mrt.y}
}
 \item{Number}{
Minimum number of observations requested in a group to be created.
}
 \item{Size}{
Size of the requested tree to be obtained (number og groups). May be limited accounting for the \code{Cond} and \code{Number} arguments. Default is 10.
}
 \item{Graph}{
Optional supplementary graphs depicting apparent and relative (cross-validated) R-squares according to tree sizes.
}
\item{Habitat1}{Plot habitat use according to the given threshold.}
\item{Habitat2}{Plot habitat use according to the two given thresholds.}
}

\value{
A list including the following elements:
\item{Change}{Values of internal splits (as appearing in the consol).}
\item{Model}{An object of class \code{rpart} created by the \code{mvpart} function from the \pkg{mvpart} package.  The corresponding tree accounts simultaneously for \code{Cond}, \code{Number} and \code{Size}}
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
\code{\link{MRT.Number}}, \code{\link{MRT.Mean}}
}
\examples{
data(univar)
# Unconstrained model
model <- MRT.Size("Distance", "Sr:Ca", univar,Size=20)
model$Size # 20 groups
abs(diff(model$Mean)) # A lot of consecutive groups have limited mean differences
model$N # A lot of groups have very few measures

# model constrained twice
model2 <- MRT.2C("Distance", "Sr:Ca", univar,Cond=0.5,Number=38,Size=20)
abs(diff(model2$Mean)) # The "Cond" argument was respected.
model2$N # The "number" argument was also respected.

# Double constrain with automatic optimal pruning
model3 <- MRT.2C("Distance", "Sr:Ca", univar,Cond=0.1,Number=38,Size=20)
model3$Size # the model results in 12 groups
model4 <- MRT.2C("Distance", "Sr:Ca", univar,Cond=0.1,Number=38)
model4$Size 
# the model now results in only 5 groups on the basis of the v-fold CV error

# Habitat(s) can also be viewed
MRT.2C("Distance", "Sr:Ca", univar,Cond=0.1,Number=38,Size=20,Habitat1=3.5)
MRT.2C("Distance", "Sr:Ca", univar,Cond=0.1,Number=38,Size=20,Habitat2=c(2.5,4))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univariate}
\keyword{constrain}

