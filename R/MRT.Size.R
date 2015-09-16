if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"rsq.rpart"))
}

MRT.Size <-
function(
  mrt.x,                  # Column name for the distance to the nucleus.
  mrt.y,                  # Column name for the elemental concentration.
  data,                   # Name of the dataframe to analyse (must be charged).
  Size,                   # Size of the requested tree to be obtained.
  Graph=FALSE,             # Supplementary graphs (optional).
  Habitat1=FALSE,         # Plot habitat use according to the given threshold.
  Habitat2=FALSE          # Plot habitat use according to the two given thresholds.
  )
{
  
  # Preliminary (package & arguments) 
  default.pars=par(no.readonly = TRUE)
  if(missing(data)){
    stop("Data is missing")
  }

   if (is.numeric(Habitat1)&is.numeric(Habitat2)) {
       stop("Can't use simultaneously Habitat1 and Habitat2 arguments. Please consider using only one.")
     }
     if (is.numeric(Habitat1)=="FALSE" & is.logical(Habitat1)=="FALSE") {
       stop("Habitat1 must be numeric.")
     }
     if (is.numeric(Habitat2)=="FALSE" & is.logical(Habitat2)=="FALSE") {
       stop("Habitat2 must be numeric.")
     }
     if (is.numeric(Habitat2)) {
       if(length(Habitat2)!=2) {
         stop("Habitat2 must contain two thresholds values")
       }
     }
     if(length(Habitat2)==2) {
       if(Habitat2[1]>Habitat2[2]) {
         Habitat2[1]->a;Habitat2[1]<-Habitat2[2];Habitat2[2]<-a    
       }
     }
  
  # Data acquisition & Modelling
  y.data<-na.omit(data[,mrt.y])
  x.data<-data[which(data[,mrt.y]!="NA"),mrt.x]
  if (length(x.data)!=length(y.data)){
    cat(paste("ERROR: 'x' and 'y' differ in length","\n"))
  }
  Data<-as.data.frame(cbind(x.data,y.data))
  names(Data)<-c("X","Y")
  if(missing(Size)) {
    model<-mvpart(Y~X,data=Data,size=Size,xv="1se")
  } else {
    model<-mvpart(Y~X,data=Data,size=Size,cp=0.0001,xv="1se")
  }
  #   model<-mvpart(Y~X,data=Data,size=Size,minsplit=100)
  #   model<-mvpart(Y~X,data=Data,minsplit=200)
  P<-predict(model)
  
  # Cluster attributes
  gr<-model$where
  aa<-1;gr2<-rep(1,length(gr))      # Renumber clusters sequentially.
  for(i in 2:length(gr)) {
    if (gr[i]!=gr[i-1]) aa <-aa+1
    gr2[i] <- aa
  }
  levels(as.factor(gr2))
  k<-length(levels(as.factor(gr2)))
  SD<-rep(1,length(levels(as.factor(gr2))))   # SD for each group.
  for(i in 1:length(SD)) {
    SD[i]<-sd(Data$Y[which(gr2==i)])
  }
  ME<-rep(1,length(levels(as.factor(gr2))))   # Mean for each group.
  for(i in 1:length(ME)) {
    ME[i]<-mean(Data$Y[which(gr2==i)])
  }
  
    # Comparing data with the habitat threshold(s) provided
     if (is.numeric(Habitat1)) {
       if (xor(as.numeric(Habitat1)>max(Data$Y) , as.numeric(Habitat1)<min(Data$Y))) {
         cat("Note: Threshold value for habitat out of range. Fish always in the same habitat");cat("\n")
       } 
     }
     if (is.numeric(Habitat2)) {
       if (Habitat2[2]>max(Data$Y) | Habitat2[1]<min(Data$Y)){
         cat("Note: Threshold value(s) for habitat out of range. Fish always in the same habitat");cat("\n")
       } 
     }
	 
  # Plotting
  if (any(is.numeric(Habitat1),is.numeric(Habitat2))) {
        layout(matrix(c(1,1,1,2),4,1))
        par(oma=c(0,0,0,0))
        #par(oma=c(0,0,0,0),mar=c(1,1,1,1))
        } else {layout(1,1,1)} 
  plot(Data$X,Data$Y,main=paste("Sample:",mrt.y,"\n",paste("Number of splits=",k-1)),col="white",xlab="Distance from primordium",ylab=expression(paste("Sr:Ca (x",10^-3,")")))
  grid()
  x<-par("usr"); rect(x[1],x[3],x[2],x[4],col="#6495ED25")  # Backgrond color can be changed here using any RGB Hex. color.
  suite<-rep(1,k)   # Suite: Number of measures in each cluster.
  for(i in 1:length(suite)) {
    suite[i]<-length(Data$Y[which(gr2==i)])
  }
  qq<-rep(SD,suite) # qq: SD of the belonging cluster at each point
  Change<-rep(0,length=k-1)
  N<-rep(0,length=k)
  abline(v=range(Data$X), col="red")
  points(range(Data$X),rep(x[3],2), col="red",pch=24,cex=2,bg="white")
  for(i in 1:k) {
    if(i!=k){
      abline(v=min(Data$X[which(gr2==i)]), col="red", lty="dotted")
      if(i!=1) {points(min(Data$X[which(gr2==i)]),x[3], col="red",pch=17,cex=2)}
      Change[i]<-max(Data$X[which(gr2==i)])
      N[i]<-length(Data$X[which(gr2==i)])
      
    } else{
      abline(v=min(Data$X[which(gr2==i)]), col="red", lty="dotted")
      points(min(Data$X[which(gr2==i)]),x[3], col="red",pch=17,cex=2)
      N[i]<-length(Data$X[which(gr2==i)])
    }
    
    rect(min(Data$X[which(gr2==i)]),min(P[which(gr2==i)]-qq[which(gr2==i)]),max(Data$X[which(gr2==i)]),max(P[which(gr2==i)]+qq[which(gr2==i)]),col="#6495ED50")
    palette(rainbow(k))
    points (Data$X[which(gr2==i)],Data$Y[which(gr2==i)],col="red",pch=21,bg=i+1)
  }
  for(i in 1:k) {
    lines(Data$X[which(gr2==i)],rep(ME[i],length(Data$X[which(gr2==i)])))
  }
  lines(Data$X,P)
  
  
  ######################################
         # 3c) Plotting habitat use (Optional)
         if (Habitat1!="FALSE") {
           abline(h=Habitat1,lwd=4,col="White")
           abline(h=Habitat1,lty=6,lwd=2)
           par(mar=c(5,2,0,2))
           plot(Data$X,Data$X,axes=F,ann=F,ylim=c(0,2),col="white",main="Habitat use");axis(1);title(xlab="Distance from primordium")
           legend("top",c("Habitat 1","Habitat 2","Unassigned"),horiz=T,fill=c("red","blue","grey"))
           for(l in 1:k) {
             if(min(P[which(gr2==l)]-qq[which(gr2==l)])>Habitat1) {
               #rect(min(xx[which(gr2==l)]),0,max(xx[which(gr2==l)]),1,col="red")
               #             reds <- colorRampPalette(c('coral2', 'red3')) 
               reds <- colorRampPalette(c("#8B0000", "#FFAEB9", "#FF0000"))
               h<-seq(from=0,to=1,length=50)
               for(m in 1:50) {
                 polygon(c(min(Data$X[which(gr2==l)]),max(Data$X[which(gr2==l)])),c(h[m],h[m]),border=reds(50)[m],col=NA)
               }
               #             gradient.rect (min(xx[which(gr2==l)]),0,max(xx[which(gr2==l)]),1,col=reds(50))
             } else {
               if(max(P[which(gr2==l)]+qq[which(gr2==l)])<Habitat1){
                 #               blues <- colorRampPalette(c('dark blue', 'light blue')) 
                 #               blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 h<-seq(from=0,to=1,length=50)
                 for(m in 1:50) {
                   polygon(c(min(Data$X[which(gr2==l)]),max(Data$X[which(gr2==l)])),c(h[m],h[m]),border=blues(50)[m],col=NA)
                 }
                 #rect(min(xx[which(gr2==l)]),0,max(xx[which(gr2==l)]),1,col="blue")
                 #               gradient.rect (min(xx[which(gr2==l)]),0,max(xx[which(gr2==l)]),1,col=blues(50))
               } else {
                 #               greys <- colorRampPalette(c('slategray4', 'slategray1')) 
                 greys<-colorRampPalette(c("#3D3D3D", "#CCCCCC", "#6B6B6B"))
                 h<-seq(from=0.1,to=0.9,length=50)
                 for(m in 1:50) {
                   polygon(c(min(Data$X[which(gr2==l)]),max(Data$X[which(gr2==l)])),c(h[m],h[m]),border=greys(50)[m],col=NA)
                 }
                 #               gradient.rect (min(xx[which(gr2==l)]),0.1,max(xx[which(gr2==l)]),0.9,col=greys(50))
               }
             }
           }
         }
         if (is.numeric(Habitat2)) {
           abline(h=Habitat2,lwd=4,col="White")
           abline(h=Habitat2,lty=6,lwd=2)
           par(mar=c(5,2,0,2))
           plot(Data$X,Data$X,axes=F,ann=F,ylim=c(0,2),col="white",main="Habitat use");axis(1);title(xlab="Distance from primordium")
           legend("top",c("Habitat 1","Habitat 2","Unassigned"),horiz=T,fill=c("red","blue","grey"))
           for(l in 1:k) {
             if(min(P[which(gr2==l)]-qq[which(gr2==l)])>Habitat2[2]){
               #rect(min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col="red")
               #             reds <- colorRampPalette(c('coral2', 'red3')) 
               reds <- colorRampPalette(c("#8B0000", "#FFAEB9", "#FF0000"))
               h<-seq(from=0,to=1,length=50)
               for(m in 1:50) {
                 polygon(c(min(Data$X[which(gr2==l)]),max(Data$X[which(gr2==l)])),c(h[m],h[m]),border=reds(50)[m],col=NA)
               }
               #             gradient.rect (min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col=reds(50))
             } else {
               if(max(P[which(gr2==l)]+qq[which(gr2==l)])<Habitat2[1]){
                 #               blues <- colorRampPalette(c('dark blue', 'light blue')) 
                 #               blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 h<-seq(from=0,to=1,length=50)
                 for(m in 1:50) {
                   polygon(c(min(Data$X[which(gr2==l)]),max(Data$X[which(gr2==l)])),c(h[m],h[m]),border=blues(50)[m],col=NA)
                 }
                 #rect(min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col="blue")
                 #               gradient.rect (min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col=blues(50))
               } else {
                 #               greys <- colorRampPalette(c('slategray4', 'slategray1')) 
                 greys<-colorRampPalette(c("#3D3D3D", "#CCCCCC", "#6B6B6B"))
                 h<-seq(from=0.1,to=0.9,length=50)
                 for(m in 1:50) {
                   polygon(c(min(Data$X[which(gr2==l)]),max(Data$X[which(gr2==l)])),c(h[m],h[m]),border=greys(50)[m],col=NA)
                 }
                 #               gradient.rect (min(xx[which(gr2==k)]),0.1,max(xx[which(gr2==k)]),0.9,col=greys(50))
               }
             }
           }
         }
  
  
  
  if(Graph=="TRUE"){
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
      legend(2,0.2, c("Apparent", "X Relative"), lty = 1:2, col="black",pch=21,cex=1,box.col="black",pt.bg=c("white","grey"),bg="white",bty="0")
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
  }
  par(new=TRUE)
  on.exit(par(default.pars))
  # Results typing
  output<-""
  output<-paste(output,"---------------------------------","\n")
  output<-paste(output,"Internal partitions identified at: ","\n")
  for(i in 1:k-1){
    if(i!=k-1){
      output<-paste(output,"    ",Change[i],"\n")
    } else {
      output<-paste(output,"    ",Change[i]," (from primordium)")
    }
  }
  cat(output)
  cat("\n")
  # Function returns
  output2<-""
  output2<-paste(output2,mrt.y,"~",mrt.x,", Data =",deparse(substitute(model.data)),", using MRT.Size","\n")
  output2<-paste(output2,"Size=",k,"\n")
  output2<-paste(output2,output,"\n")
  ans<-list()
  class(ans) <- "MRT"
  ans<-list(Change=Change,Model=model,Mean=ME,SD=SD,Belong=gr2,N=N,Predict=P,Size=k,output=output2,X=Data$X,Y=Data$Y)
}
