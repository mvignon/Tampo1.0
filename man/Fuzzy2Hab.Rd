\name{Fuzzy2Hab}
\alias{Fuzzy2Hab}
\title{
Fuzzy extraction of habitat use from density probability
}
\description{
The function uses relative density probability (likelihood) from two habitats to perform a fuzzy extraction of habitat use. Density probabilities must previously be extracted from the \code{Dens} function. 
}
\usage{
Fuzzy2Hab(mrt.x, mrt.y, data, Dens1, Dens2, loess = FALSE, span = 0.05,
          Graph = FALSE, W = TRUE)
}
\arguments{
  \item{mrt.x}{Column name for the chronological values to be used, typically the distance to the nucleus. Data in this colum must be numeric.
}
  \item{mrt.y}{Column name for the elemental concentrations to be used. Data in this column must be numeric.
}
  \item{data}{Dataframe to be analysed.}
  \item{Dens1}{Results of \code{Dens} for habitat 1.}
  \item{Dens2}{Results of \code{Dens} for habitat 2.}
  \item{loess}{Display fuzzy clustering using Loess instead of raw data.}
  \item{span}{\emph{Alpha} parameter controling the degree of smoothing of the loess. See \code{\link[stats]{loess}}.
}
  \item{Graph}{To (additionally) plot the probability distributions of the two habitats (similar to \code{Dens} graphics).}
  \item{W}{To open a new windows. Must be turned off when \code{Dens2Hab} is used to print in a PDF file.}
}
\details{
In a fuzzy clustering approach, each data point belongs to all habitats, but only to a certain degree. Each data point is assigned a probability of belonging to each habitat (based on the empirical distribution extracted by the \code{Dens} function). Fuzzy clustering provides a way of capturing quantitative uncertainty and even visualizing it. While the function mainly rely on the graphical result, attribution is expressed as the percentage of \code{mrt.y} values falling exclusively inside the confidence interval of one habitat or falling outside any of the two intervals. These values are considered respectively as \emph{Assigned} and \emph{Unassigned}. Alternatively, a value falling inside both confidence intervals may be equally likely to belong to both clusters and clustering is considered as \emph{Equivocal}. This only appear when confidence intervals overlap.
}
\value{
A list including the following elements:
\item{Equivocal}{Percentage of \code{mrt.y} values equivocally attributed to both habitats.}
\item{Unassigned}{Percentage of \code{mrt.y} values remaining unassigned to any habitat.}
\item{Assigned}{Percentage of \code{mrt.y} values specifically assigned to one habitat.}
%% RAJOUTER DATA.FRAME AVEC PROBA RESPECTIVE !!!!!
}


\seealso{
\code{\link{Dens}}
\code{\link{Dens2Hab}}
}
\examples{
data(univar)

# hypothetical distribution for habitat 1
x1<-rnorm(1000,mean=4.5,sd=1)
D1<-Dens(x1)
# hypothetical distribution for habitat 2
x2<-rnorm(1000,mean=2,sd=1)
D2<-Dens(x2)

# Habitat extraction
Fuzzy2Hab("Distance","Sr:Ca",univar,D1,D2,Graph=FALSE,loess=FALSE)
Fuzzy2Hab("Distance","Sr:Ca",univar,D1,D2,Graph=FALSE,loess=TRUE)
Fuzzy2Hab("Distance","Sr:Ca",univar,D1,D2,Graph=FALSE,loess=TRUE,span=0.4)

}
\keyword{Density}
\keyword{univariate}
\keyword{habitat}
