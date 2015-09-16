\name{MRT.Mean}
\alias{MRT.Mean}
\title{Univariate regression tree with constrain on the difference between the means of two consecutive clusters.
}
\description{
The function prunes a fully grown tree by fixing the minimal difference requested between the mean elemental concentrations of two consecutive clusters. This process removes clusters that do not exhibit enough consecutive variation to be interpretable. If the provided cluster differences is not limitative, then, the function creates a tree of size \code{Size}.
}
\usage{
MRT.Mean(mrt.x, mrt.y, data, Cond = 1, Size, Graph = FALSE, Habitat1 = FALSE,
         Habitat2 = FALSE)
}
\arguments{
 \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column name for the elemental concentration to be used. Data in this colum must be numeric.}
  \item{data}{ Dataframe to be analysed.}
  \item{Cond}{The minimal difference requested between the mean elemental concentrations of two consecutive clusters. Same units that \code{mrt.y}
}
  \item{Size}{
Size of the requested tree to be obtained (number og groups). May be limited accounting for the \code{Cond} argument.
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
\item{Model}{An object of class \code{rpart} created by the \code{mvpart} function from the \pkg{mvpart} package.  The corresponding tree accounts for both \code{Cond} and \code{Size}}
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
\code{\link{MRT.Size}}, \code{\link{MRT.2C}}, \code{\link{MRT.Number}}
}
\examples{
data(univar)
# Unconstrained model with limited mean differences
model <- MRT.Size("Distance", "Sr:Ca", univar,Size=20)
model$Size # 20 groups
abs(diff(model$Mean)) # A lot of consecutive groups have limited mean differences

# Constrained model with means of consecutive groups having at 1 unit of difference
model2 <- MRT.Mean("Distance", "Sr:Ca", univar,Cond=1,Size=20)
model2$Size # Number of groups filling the "Cond" argument 
abs(diff(model2$Mean)) # The "number" argument was respected.
# Other constrain
model3 <- MRT.Mean("Distance", "Sr:Ca", univar,Cond=0.5,Size=20)

# 'Cond' constrain with automatic optimal pruning
model3$Size # the model results in 8 groups
model4 <- MRT.Mean("Distance", "Sr:Ca", univar,Cond=0.5)
model4$Size
# model4 results in 6 groups and therefore discard 2 possibles groups
#   on the basis of the v-fold cross-validation error

# Habitat(s) can also be viewed
MRT.Mean("Distance", "Sr:Ca", univar,Cond=0.5,Habitat1=3)
MRT.Mean("Distance", "Sr:Ca", univar,Cond=0.5,Habitat2=c(2.5,4))
MRT.Mean("Distance", "Sr:Ca", univar,Cond=1.5,Habitat2=c(2.5,4))
}

\keyword{univariate}
\keyword{constrain}
