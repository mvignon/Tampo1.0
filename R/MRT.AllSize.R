MRT.AllSize <-
function(
  mrt.x,                  # Column name for the distance to the nucleus.
  mrt.y,                  # Column name for the elemental concentration.
  data,                   # Name of the dataframe to analyse (must be charged).
  Size=10,                # Size of the requested tree to be obtained.
  GraphName=""            # Name of the PDF file to be created.
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
  names(Data)<-c("Distance","Y")
  
  if(GraphName==""){
	NAMEInit<-NAME<-"MRT.AllSizeGraphics"
	nu<-1
	while(file.exists(paste(NAME,".pdf",sep=""))){
		NAME<-paste(NAMEInit,"(",nu,")",sep="")
		nu<-nu+1
	}
	pdf(paste(NAME,".pdf",sep=""))
  } else {  pdf(paste(GraphName,".pdf")) }
  
  
  for(j in 2:Size){ # Mettre 3 si l'on veut supprimer seulement 2 groupes avec 1 split
    # Model
    model<-mvpart(Y~Distance,data=Data,size=j,cp=0.0001)
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
    plot(Data$Distance,Data$Y,main=paste("Sample:",mrt.y,"\n",paste("Number of splits=",ki-1)),col="white",xlab="Distance from primordium",ylab=expression(paste("Sr:Ca (x",10^-3,")")))
    grid()
    x<-par("usr"); rect(x[1],x[3],x[2],x[4],col="#6495ED25")
    suite<-rep(1,k)
    for(k in 1:length(suite)) {
      suite[k]<-length(Data$Y[which(gr2==k)])
    }
    qq<-rep(SD,suite)
    for(k in 1:ki) {
      abline(v=min(Data$Distance[which(gr2==k)]), col="red", lty="dotted")
      points(min(Data$Distance[which(gr2==k)]),x[3], col="red",pch=17,cex=2)
      rect(min(Data$Distance[which(gr2==k)]),min(P[which(gr2==k)]-qq[which(gr2==k)]),max(Data$Distance[which(gr2==k)]),max(P[which(gr2==k)]+qq[which(gr2==k)]),col="#6495ED50")
      palette(rainbow(ki))
      points (Data$Distance[which(gr2==k)],Data$Y[which(gr2==k)],col="red",pch=21,bg=k+1)
    }
    lines(Data$Distance,P)
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
  rsqrpart(model)
  dev.off()
  
  # Results typing
  cat("-----------------------------------------------------","\n")
  cat("PDF file successfully created in the following directory: ","\n")
  cat(getwd(),"\n")

  # Function returns
  ans<-list(Cptable=model$cptable,FullTree=model)
}
