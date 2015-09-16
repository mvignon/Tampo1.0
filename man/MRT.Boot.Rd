\name{MRT.Boot}
\alias{MRT.Boot}
\title{
 Assessing the robustness of an univariate partition process by bootstrapping a given bagging fraction
}
\description{
Original data are bootstraped to ensure that the models produced robust detection of time series discontinuities.
}
\usage{
MRT.Boot(mrt.x, mrt.y, data, Size, Bagging = 0.5, N = 100, IC = 0.8)
}
\arguments{
 \item{mrt.x}{
Column name for the chronological values to be used as a time contiguity constraint, typically the distance to the nucleus. Data in this colum must be numeric.}
  \item{mrt.y}{
Column name for the elemental concentration to be used. Data in this colum must be numeric.}
  \item{data}{ Dataframe to be analysed.}
  \item{Size}{
Size of the requested tree to be obtained (number og groups).
}
  \item{Bagging}{
The in-bag fraction [0-1]. Default is 0.5.
}
  \item{N}{
The number of boostraps. Default is 100.
}
  \item{IC}{
Confidence interval to be used for quantiles. Default is 0.8.
}
}
\details{
The exact threshold value selected at each node is intrinsically case-sensitive
(instability to small perturbations of the data). This is of primarily importance because elemental ratios extracted from the same otolith can slightly differs according to the transect position or the machine laser resolution. To examine properties of the groups formed by chronological clustering and ensure that the models produced robust detection of time series discontinuities, the measurements were further resampled for model building using a given ("in-bag") fraction of the original data (caracterized by the \code{Bagging} argument). Characteristics of a large ensemble of MRT models provide a more accurate characterization of time series discontinuities. From the ensemble, the function return the number of time a given node appears in a particular zone. This allows quantifying the robustness of structural changes in the data. The associated values (mean and quantiles) are recorded, and an histogram of the bootstrapped values is plotted in the lower panel of the graph.
}
\value{
A list including the following elements:
\item{Mean}{Mean threshold values for each split, ordered sequentially.}
\item{Boot}{A matrix of dimension [\code{N},\code{Size}] containing threshold values for each split (columns) and each bootstrap (rows).}
\item{Qant}{A matrix containing qantiles for all splits.}
\item{OutSize}{The number of boostraps that did not produced tree of size \code{Size}.}

}

\note{
Quantiles are defined as [(1 - \code{IC})/2 ; 1-(1-\code{IC})/2)]. Default quantiles are [0.1 ; 0.9] using \code{IC}=0.8.
1-\code{OutSize} gives the true number of bootstraps used for means and qantiles calculation. If \code{OutSize}>0, \code{Boot} is of dimension [1-\code{OutSize},\code{Size}].
}



\seealso{
\code{\link{MRT.Boot.CI}}
}
\examples{
# Real data example
data(univar)
model<-MRT.Boot("Distance","Sr:Ca",univar,Size=4,Bagging=0.5,N=100)
mean(model$Boot[,1])==model$Mean[1]

# Simulated data example
# Generating a random time series exhibiting two apparent discontinuities,
# the first being abrupt and the second being obviously more smoothly.
# Rand<-c(rep(-0.5,length=250)+rnorm(250,0,0.25), rep(1,length= 500)+
#     rnorm(500,0,0.25), sin(seq(from=1.5,to=5,length=500))+
#     rnorm(500,0,0.25),rep(-1,length=200)+rnorm(200,0,0.25))
# xx<-1:length(Rand)
# Dat<-data.frame(X=xx,Y=Rand*100);plot(Dat$X,Dat$Y)
# MRT.Boot("X","Y", Dat, Size=3, Bagging=0.5, N=200)

}

\keyword{univariate}
\keyword{bootstrap}
