Dens <-
function(
  x,          # Data to investigate.
  IC=0.95     # IC to be used.
)
{
  # P<-c((1-IC)/2,1-(1-IC)/2)
  # Calculate density probability.
  d<-density(x)
  Dmax<-max(d$y)
  dd<-1/Dmax*d$y # To standardize the density function between [0-1]
  # Plot histogramme.
  hist(x,freq=FALSE)
  lines(d)
  # Extrapolate density function.
  approxfun(d$x,d$y,yleft=0,yright=0)->f1
  approxfun(d$x,dd,yleft=0,yright=0)->f2
  # Cumulative density probability.
  CumD<-numeric(length=length(d$x))
  CumD[1]<-d$y[1]
  for(i in 2:length(d$x)) {CumD[i]<-CumD[i-1]+d$y[i]}
  CumDD<-1/max(CumD)*CumD
  Min<-d$x[max(which(CumDD<=(1-IC)/2))]
  Max<-d$x[min(which(CumDD>=1-(1-IC)/2))]
  Mode<-d$x[which(dd==max(dd))]
  abline(v=c(Min,Max))
  # Advanced Plot.
  maxyplot<-max(dd)*1.05
  minxplot<-min(d$x)-0.1*(range(d$x)[2]-range(d$x)[1])
  maxxplot<-max(d$x)+0.1*(range(d$x)[2]-range(d$x)[1])
  x<-seq(minxplot,maxxplot,0.1)
  plot(d$x,dd,ylim=c(0, maxyplot),xlim=c(minxplot,maxxplot),type="n",ylab="relative density probability",xlab="Elemental concentration")
  title(main=paste("Density probability at ",IC*100,"%")) 
  dim<-par("usr")
  polygon(d$x, dd, col="#FF404075", border=NA)
  xxmin<-d$x[which(d$x<=Min)]
  xxmax<-d$x[which(d$x>=Max)]
  polygon(c(min(xxmin), xxmin, max(xxmin)), c(0, dd[which(d$x<=Min)], 0), col="#FF4040", border=NA)
  polygon(c(min(xxmax), xxmax, max(xxmax)), c(0, dd[which(d$x>=Max)], 0), col="#FF4040", border=NA)
  abline(v=c(Min,Max,d$x[which(d$y==max(d$y))]),col="white")
  lines(d$x, dd)
  text(d$x[which(d$y==max(d$y))],max(dd),pos=3, paste(IC*100,"% IC"), col="black")
  text(c(Min,Max),c(dd[which(d$x==Min)],dd[which(d$x==Max)]),c(paste(round(Min,digits = 2)),paste(round(Max,digits = 2))),pos=3,)
  segments(dim [1], dd[which(d$x==Min)], Min, dd[which(d$x==Min)], col="black")
  segments(dim [2], dd[which(d$x==Max)], Max, dd[which(d$x==Max)], col="black")
  text(dim[1]+(Min-dim [1])/3,dd[which(d$x==Min)],pos=3, paste(round(1/Dmax*f1(Min),digits = 2)), col="black")
  text(dim[2]-(dim [2]-Max)/3,dd[which(d$x==Max)],pos=3, paste(round(1/Dmax*f1(Max),digits = 2)), col="black")

  ans<-list(Fun=f1,NFun=f2,Min=Min,Max=Max,Mode=Mode,LP=1/Dmax*f1(Min),RP=1/Dmax*f1(Max),IC=IC)
}
