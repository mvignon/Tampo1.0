MRT.Size.Multi <-
function(
     mrt.x,                  # Column name for the distance to the nucleus.
     mrt.y,                  # Column name for the elemental concentration.
     data,                   # Name of the dataframe to analyse (must be charged).
     Size,
     Scale=FALSE,            # To scale the elemental data.
     Graph=FALSE
   )
   {
     # Preliminary (package & arguments) 
     if(missing(data)){
       stop("Data is missing")
     }

     default.pars=par(no.readonly = TRUE)
     # Data acquisition & Modelling
     L<-length(mrt.y)+1
     # A FAIRE : verifier longueur des L vecteurs et prendre le plus petit? avec notes d'info qui s'affiche
     #####
     y.data<-as.data.frame(na.omit(data[,mrt.y]))
     x.data<-data[which(y.data[,1]!="NA"),mrt.x]
     if (length(x.data)!=dim(y.data)[1]){
       cat(paste("ERROR: 'x' and 'y' differ in length","\n"))
     }
     Data<-as.data.frame(cbind(x.data,y.data))
     names(Data)[1]<-"X"
     if(Scale==TRUE) {
       NonScaleData<-Data
       Data[,-1]<-scale(Data[,-1],scale = TRUE,center=TRUE)
     }
     model<-mvpart(data.matrix(Data[,2:L])~X,data=Data,size=Size,cp=0.0001)
     P<-as.data.frame(predict(model));names(P)<-mrt.y
     
     # Cluster attributes
     gr<-model$where
     aa<-1;gr2<-rep(1,length(gr))      # Renumber clusters sequentially.
     for(i in 2:length(gr)) {
       if (gr[i]!=gr[i-1]) aa <-aa+1
       gr2[i] <- aa
     }
     k<-length(levels(as.factor(gr2)))
     SD<-matrix(1,ncol=L-1,nrow=k)               # SD for each group (lines) and element (column).
     ME<-matrix(1,ncol=L-1,nrow=k)               # Mean for each group (lines) and element (column).
     for(i in 1:dim(SD)[1]) {
       for(j in 1:dim(SD)[2]) {
         SD[i,j]<-sd(Data[which(gr2==i),j+1])
         ME[i,j]<-mean(Data[which(gr2==i),j+1])
       }
     }
     colnames(ME)<-colnames(SD)<-mrt.y
     rownames(ME)<-rownames(SD)<-rownames(ME,do.NULL=FALSE,prefix="Group ")
     
     # Plotting
     if(Scale=="FALSE") {
	 	matplot(Data$X,Data[,2:L],main=paste("Number of splits = ",k-1),col="white",xlab="Distance from primordium",ylim=Ymax<-c(min(Data[,-1]),max(Data[,-1])),ylab="Original Units",type="l",lty=1)
     } else matplot(Data$X,Data[,2:L],main=paste("Number of splits = ",k-1),col="white",xlab="Distance from primordium",ylim=Ymax<-c(min(Data[,-1]),max(Data[,-1])),ylab="Scaled Units",type="l",lty=1)
	 grid()
     x<-par("usr"); rect(x[1],x[3],x[2],x[4],col="#6495ED25")  # Backgrond color can be changed here using any RGB Hex. color.
     suite<-rep(1,k)   # Suite: Number of measures in each cluster.
     for(i in 1:length(suite)) {
       suite[i]<-dim(Data[which(gr2==i),])[1]
     }
     qq<-matrix(1,ncol=dim(y.data)[2],nrow=length(x.data)) # qq: SD of the belonging cluster at each point
     for(i in 1:dim(y.data)[2]) {
       qq[,i]<-rep(SD[,i],suite)
     }
     Change<-rep(0,length=k-1)
     N<-rep(0,length=k)
     matlines(Data$X,Data[,2:L],ylim=Ymax<-c(min(Data[,-1]),max(Data[,-1])),type="l",lty=1,col="darkgrey")
     abline(v=range(Data$X), col="red")
     points(range(Data$X),rep(x[3],2), col="red",pch=24,cex=2,bg="white")
     for(i in 1:k) {
       if(i==1){
         Change[i]<-max(Data$X[which(gr2==i)])
         N[i]<-length(Data$X[which(gr2==i)])
       } else {
         abline(v=min(Data$X[which(gr2==i)]), col="red", lty="dotted")
         points(min(Data$X[which(gr2==i)]),x[3], col="red",pch=17,cex=2)
         if(i!=k){  Change[i]<-max(Data$X[which(gr2==i)])}
         N[i]<-length(Data$X[which(gr2==i)])
       }
       palette(rainbow(k))
       for(j in 1:dim(y.data)[2]) {
         points (Data$X[which(gr2==i)],Data[which(gr2==i),j+1],col="red",pch=21,bg=i+1)
         rect(min(Data$X[which(gr2==i)]),min(P[which(gr2==i),j]-qq[which(gr2==i),j]),max(Data$X[which(gr2==i)]),max(P[which(gr2==i),j]+qq[which(gr2==i),j]),col="#6495ED50")
         lines(Data$X[which(gr2==i)],rep(ME[i,j],length(Data$X[which(gr2==i)])),col=j,lwd=2,lty=j)
       }
     }
     #legend("top",inset=-0.05,mrt.y,horiz=T,fill=1:dim(y.data)[2])
     #legend("top",inset=-0.05,bg="white",mrt.y,horiz=T, col=1:dim(y.data)[2], lty=1:dim(y.data)[2])
	 legend("top",inset=-0.005,bg="white",mrt.y,horiz=T, col=1:dim(y.data)[2], lty=1:dim(y.data)[2])
     
     if(Graph=="TRUE"){
       dev.new();par(mfrow=c(1,2));rsq.rpart(model)
     }
     cat("-----------------------------","\n")
     cat("Internal partitions identified at: ","\n")
     for(i in 1:k-1){
       if(i!=k-1){
         cat("    ");cat(Change[i],"\n")
       } else {
         cat("    ");cat(Change[i]);cat(" (from primordium)")
       }
     }
	 cat("\n")
     par(default.pars)
     #   return(Change)
     if(Scale==TRUE) {
       k<-length(levels(as.factor(gr2)))
       SD<-matrix(1,ncol=L-1,nrow=k)         
       ME<-matrix(1,ncol=L-1,nrow=k)
       for(i in 1:dim(SD)[1]) {
         for(j in 1:dim(SD)[2]) {
           SD[i,j]<-sd(NonScaleData[which(gr2==i),j+1])
           ME[i,j]<-mean(NonScaleData[which(gr2==i),j+1])
         }
       }
       colnames(ME)<-colnames(SD)<-mrt.y
       rownames(ME)<-rownames(SD)<-rownames(ME,do.NULL=FALSE,prefix="Group ")
     }
     ans<-list(Change=Change,Mean=ME,SD=SD,Model=model,Belong=gr2,N=N,Predict=P,Size=k)
   }
