MRT.All <-
function(
  mrt.x,                  # Column name for the distance to the nucleus.
  mrt.y,                  # Column name for the elemental concentration.
  data                    # Name of the dataframe to analyse (must be charged).
)
{
  # Preliminary (package & arguments) 
  if(missing(data)){
    stop("Data is missing")
  }

  
  # Data acquisition & Modelling
  y.data<-na.omit(data[,mrt.y])
  x.data<-data[which(data[,mrt.y]!="NA"),mrt.x]
  if (length(x.data)!=length(y.data)){
    cat(paste("ERROR: 'x' and 'y' differ in length","\n"))
  }
  Data<-as.data.frame(cbind(x.data,y.data))
  names(Data)<-c("X","Y")
 

Size<-2
while(Size!=0){
    # Model
    j<-Size
	if(Size==2) {kii<-0}
	model<-mvpart(Y~X,data=Data,size=j,cp=0.0000001)
    #model<-mvpart(Y~X,data=Data,size=j,cp=0.0001,xv="1se")
    P<-predict(model)
    # Renumber clusters sequentially
    gr<-model$where
    aa<-1;gr2<-rep(1,length(gr))
    for(k in 2:length(gr)) {
      if (gr[k]!=gr[k-1]) aa <-aa+1
      gr2[k] <- aa
    }
    levels(as.factor(gr2))
    ki<-length(levels(as.factor(gr2)))
    # Liste des Standard deviation + Meanpar groupe (apres groupes sequentiels)
    SD<-rep(1,length(levels(as.factor(gr2))))
    for(k in 1:length(SD)) {
      SD[k]<-sd(Data$Y[which(gr2==k)])
    }
    ME<-rep(1,length(levels(as.factor(gr2))))
    for(k in 1:length(ME)) {
      ME[k]<-mean(Data$Y[which(gr2==k)])
    }
    # Plot avec Mean+-SD/groupes
    if (ki>kii) {
		plot(Data$X,Data$Y,main=paste("Sample:",mrt.y,"\n",paste("Number of splits=",ki-1)),col="white",xlab="Distance from primordium",ylab=expression(paste("Sr:Ca (x",10^-3,")")))
		grid()
		x<-par("usr"); rect(x[1],x[3],x[2],x[4],col="#6495ED25")
		suite<-rep(1,k)
		for(k in 1:length(suite)) {
		suite[k]<-length(Data$Y[which(gr2==k)])
		}
		qq<-rep(SD,suite)
		for(k in 1:ki) {
		abline(v=min(Data$X[which(gr2==k)]), col="red", lty="dotted")
		points(min(Data$X[which(gr2==k)]),x[3], col="red",pch=17,cex=2)
		rect(min(Data$X[which(gr2==k)]),min(P[which(gr2==k)]-qq[which(gr2==k)]),max(Data$X[which(gr2==k)]),max(P[which(gr2==k)]+qq[which(gr2==k)]),col="#6495ED50")
		palette(rainbow(ki))
		points (Data$X[which(gr2==k)],Data$Y[which(gr2==k)],col="red",pch=21,bg=k+1)
		}
		lines(Data$X,P)
		kii<-ki
	
		yorn <- readline(prompt = "Do you want to split again? Press <enter> to continue or N to stop: ")
		if(yorn == 'N' || yorn =='n') {
			cat('analysis stopped with',ki,"groups","\n")
			Size<-0
		} else {
			Size<-Size+1
			cat('regression tree updated',"\n")
		}
	} else {Size<-Size+1}
	}
  rsqrpart<-function (x) 
  {
    if (!inherits(x, "rpart")) 
      stop("Not legitimate rpart")
    p.rpart <- x$cptable
    xstd <- p.rpart[, 5]
    xerror <- p.rpart[, 4]
    rel.error <- p.rpart[, 3]
    nsplit <- p.rpart[, 2]
    method <- x$method
    if (!method == "anova") 
      cat("May not be applicable for this method\n")
    plot(nsplit, 1 - rel.error, xlab = "Number of Splits", ylab = "R-square",ylim = c(0, 1), type = "l")
    grid()
    A<-par("usr"); rect(A[1],A[3],A[2],A[4],col="#6495ED25")
    points(nsplit, 1 - rel.error,pch=21,col="black",bg="white")
    par(new = TRUE)
    plot(nsplit, 1 - xerror, xlab = "", ylab = "",ylim = c(0, 1), type = "l",lty=2)
    points(nsplit, 1 - xerror,pch=21,col="black",bg="grey")
    legend(2/3*(A[2]-A[1]),1/4*(A[4]-A[3]) , c("Apparent", "X Relative"), lty = 1:2, col="black",pch=21,cex=1,box.col="black",pt.bg=c("white","grey"),bg="white",bty="0")
    ylim <- c(min(xerror - xstd) - 0.1, max(xerror + xstd) + 
      0.1)
    plot(nsplit, xerror, xlab = "Number of Splits", ylab = "X Relative Error", ylim = ylim, type = "l")
    grid()
    A<-par("usr"); rect(A[1],A[3],A[2],A[4],col="#6495ED25")
    points(nsplit, xerror,pch=21,col="black",bg="grey")
    segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
    invisible()
  }
   dev.new();par(mfrow=c(1,2));rsqrpart(model)
  
  # Function returns
  ans<-list(Cptable=model$cptable,FullTree=model)
}
