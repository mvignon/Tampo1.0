MRT.Boot.Multi <-
function(
     mrt.x,                  # Column name for the distance to the nucleus.
     mrt.y,                  # Column name for the elemental concentration.
     data,                   # Name of the dataframe to analyse (must be charged).
     Size,                   # Size of the requested tree.
     Bagging=0.5,            # The in-bag fraction [0-1].
     N=100,                  # The number of boostrap (depending on your data, start using moderate number)
     IC=0.8,
     Scale=FALSE
   )
   {    #x=data,i=colonne,Size=tree size,replace=T/F ???,Bagging [0-1]
     
     ######################################
     # 1) Preliminary (package & arguments) 
     ######################################
     
     if(require("boot")=="FALSE"){
       print("trying to install boot")
       install.packages("boot")
       if(require("boot")){
         print("boot installed and loaded")
       } else {
         stop("could not install boot. Please try manually")
       }
     }
	 
	 default.pars=par(no.readonly = TRUE)
     t1 <- unclass(Sys.time())
	 
     
     ######################################
     # 2)  Data acquisition & Modelling
     ######################################
     L<-length(mrt.y)+1
     y.data<-as.data.frame(na.omit(data[,mrt.y]))
     x.data<-data[which(y.data[,1]!="NA"),mrt.x]
     if (length(x.data)!=dim(y.data)[1]){
       cat(paste("ERROR: 'x' and 'y' differ in length","\n"))
     }
     Data<-as.data.frame(cbind(x.data,y.data))
     names(Data)[1]<-"X"
     if(Scale==TRUE) {
       Data[,-1]<-scale(Data[,-1],scale = TRUE,center=TRUE)
     }
     model0<-mvpart(data.matrix(Data[,2:L])~X,data=Data,size=Size,cp=0.0001)
     P<-as.data.frame(predict(model0));names(P)<-mrt.y
     
     # To renumber clusters sequentially.
     gr<-model0$where
     aa<-1;gr2<-rep(1,length(gr))
     for(k in 2:length(gr)) {
       if (gr[k]!=gr[k-1]) aa <-aa+1
       gr2[k] <- aa
     }
     ki<-length(levels(as.factor(gr2)))
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
     
     ######################################
     # 3)  Bootstrapping
     ######################################
     # Create the bagging fraction.
     Max<-length(Data$X)*Bagging
     a<-Data
     Results<-matrix(ncol=ki,nrow=N)       # Size ?? remplacer par ki-N ?
	 progress<-txtProgressBar(min = 1, max = N, style = 3)
     for(i in 1:N){
       b <- sample(1:nrow(a), replace=FALSE)
       cc<-a[b[1:Max],]
       dd<-cc[order(cc[,1]),]
       modelOpt<-mvpart(data.matrix(dd[,2:L])~X,data=dd,size=ki,cp=0.0001)
       POpt<-predict(modelOpt,dd)
       grb<-modelOpt$where
       aa<-1;gr2b<-rep(1,length(grb))
       for(k in 2:length(grb)) {
         if (grb[k]!=grb[k-1]) aa <-aa+1
         gr2b[k] <- aa
       }
       Condition<-0
       v<-vector(length=ki)
       if(length(levels(as.factor(gr2b)))==ki){ # Verify that the boostrapped tree is of requested size
         for(k in 1:ki) {
           v[k]<-min(dd$X[which(gr2b==k)]) # v: min valu for each group
         }
       } else {
         Condition<-Condition+1
       }
       tt<-as.vector(t(v[1:ki]))
       Results[i,]<-tt
	   setTxtProgressBar(progress, i)
     }
     close(progress)
	 Resultss<-Results[which(Results[,2]!=0),-1] #resultats des boostraps.
     
     ######################################
     # 4)  Quantiles calculus
     ######################################
     quant<-function(x){
       quantile(x,c(0.1,0.9))
     }
     z<-apply(Resultss,2,quant)
     colnames(z)<-colnames(z,do.NULL=FALSE,prefix="Split ")
     m<-apply(Resultss,2,mean)
     m<-t(as.matrix(m))
     colnames(m)<-colnames(m,do.NULL=FALSE,prefix="Split ")
     rownames(m)<-""
     #colnames(z)<-c("1","2","3","4","5")
     #colnames(Results[i,])<-cat(paste("SpN??",i))
     
     ######################################
     # 5)  Plot avec Mean+-SD/groupes
     ######################################
     # Raw data and background analysis.
     layout(matrix(c(1,1,2,2),4,1))
     matplot(Data$X,Data[,2:L],main=paste("Sample",names(Data[,2]),"\n",paste("Number of splits=",k-1)),col="white",xlab="Distance from primordium",ylim=Ymax<-c(min(Data[,-1]),max(Data[,-1])),ylab=expression(paste("Sr:Ca (x",10^-3,")")),type="l",lty=1)
     grid()
     x<-par("usr"); rect(x[1],x[3],x[2],x[4],col="#6495ED25")  # Backgrond color can be changed here using any RGB Hex. color.
     suite<-rep(1,ki)   # Suite: Number of measures in each cluster.
     for(i in 1:length(suite)) {
       suite[i]<-dim(Data[which(gr2==i),])[1]
     }
     qq<-matrix(1,ncol=dim(y.data)[2],nrow=length(x.data)) # qq: SD of the belonging cluster at each point
     for(i in 1:dim(y.data)[2]) {
       qq[,i]<-rep(na.omit(SD[,i]),suite)
     }
     Change<-rep(0,length=k-1)
     matlines(Data$X,Data[,2:L],ylim=Ymax<-c(min(Data[,-1]),max(Data[,-1])),type="l",lty=1,col="darkgrey")
     abline(v=range(Data$X), col="red")
     points(range(Data$X),rep(x[3],2), col="red",pch=24,cex=2,bg="white")
     for(i in 1:k) {
       if(i==1){
         Change[i]<-max(Data$X[which(gr2==i)])
       } else {
         abline(v=min(Data$X[which(gr2==i)]), col="red", lty="dotted")
         points(min(Data$X[which(gr2==i)]),x[3], col="red",pch=17,cex=2)
         if(i!=k){  Change[i]<-max(Data$X[which(gr2==i)])}
       }
       palette(rainbow(k))
       for(j in 1:dim(y.data)[2]) {
         points (Data$X[which(gr2==i)],Data[which(gr2==i),j+1],col="red",pch=21,bg=i+1)
         rect(min(Data$X[which(gr2==i)]),min(P[which(gr2==i),j]-qq[which(gr2==i),j]),max(Data$X[which(gr2==i)]),max(P[which(gr2==i),j]+qq[which(gr2==i),j]),col="#6495ED50")
         lines(Data$X[which(gr2==i)],rep(ME[i,j],length(Data$X[which(gr2==i)])))
       }
     }
     
     hist(Resultss,breaks=seq(from=min(Data$X),to=max(Data$X),by=5),xlim=c(min(Data$X),max(Data$X)),main=paste("Bootstrap",N),col="white",xlab="")
     grid();p<-par("usr"); rect(p[1],p[3],p[2],p[4],col="#6495ED25")
     abline(v=m, col="red", lty="dotted")
     par(new=T)
     hist(Resultss,breaks=seq(from=min(Data$X),to=max(Data$X),by=5),xlim=c(min(Data$X),max(Data$X)),main=paste("Bootstrap",N),border="black",col="darkblue",xlab="Distance from primordium")
     abline(v=range(Data$X), col="red")
     par(new=F)
     par(default.pars)
     
     t2 <- unclass(Sys.time())
     elapsed.time.minutes <- round((t2-t1),3)  #calculate the total elapsed time
     cat("-----------------------------","\n")
     cat("(elapsed time - ",round(elapsed.time.minutes,1),"seconds)","\n")
     cat("-----------------------------","\n")
     cat("Mean results:","\n")
     print(m)
     cat("\n")
     cat("Quantile results:","\n")
     print(z)
     layout(matrix(1,1))
     par(default.pars)
     ans<-list(Mean=m,Boot=Resultss,Qant=z,OutSize=N-dim(Resultss)[1])
   }
