MRT.Boot.CI <-
function(
     mrt.x,                  # Column name for the distance to the nucleus.
     mrt.y,                  # Column name for the elemental concentration.
     data,                   # Name of the dataframe to analyse (must be charged).
     Size,                   # Size of the requested tree.
     Bagging=0.5,            # The in-bag fraction [0-1].
     N=100,                  # The number of boostrap (depending on your data, start using moderate number)
     IC=0.8                  # Confidence interval to be used for quantiles.
)
   {  	#x=data,i=colonne,Size=tree size,replace=T/F ???,Bagging [0-1]
     
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
     y.data<-na.omit(data[,mrt.y])
     x.data<-data[which(data[,mrt.y]!="NA"),mrt.x]
     xx<-as.data.frame(x.data)
     if (length(x.data)!=length(y.data)){
       cat(paste("ERROR: 'x' and 'y' differ in length","\n"))
     }
     Data<-as.data.frame(cbind(x.data,y.data))
     names(Data)<-c("X","Y")
     model0<-mvpart(Y~X,data=Data,size=Size)
     P<-predict(model0)
     
     # To renumber clusters sequentially.
     gr<-model0$where
     aa<-1;gr2<-rep(1,length(gr))
     for(k in 2:length(gr)) {
       if (gr[k]!=gr[k-1]) aa <-aa+1
       gr2[k] <- aa
     }
     ki<-length(levels(as.factor(gr2)))
     
     ######################################
     # 3)  Bootstrapping
     ######################################
     # Create the bagging fraction.
     Max<-length(Data$X)*Bagging
     a<-data.frame(A=Data$X,B=as.matrix(Data$Y))
     Results<-matrix(ncol=ki,nrow=N)       # Size ?? remplacer par ki-N ?
	 progress<-txtProgressBar(min = 1, max = N, style = 3)
     for(i in 1:N){
       b <- sample(1:nrow(a), replace=FALSE)
       cc<-a[b[1:Max],]
       dd<-cc[order(cc$A),]
       modelOpt<-mvpart(B~A,data=dd,size=ki)
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
           v[k]<-min(dd$A[which(gr2b==k)]) # v: min valu for each group
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
     P<-c((1-IC)/2,1-(1-IC)/2)
     quant<-function(x){
       quantile(x,P)
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
     par(mfrow=c(1,1))
     plot(Data$X,Data$Y,col="white",ylab=expression(paste("Sr:Ca (x",10^-3,")")),xlab="Distance from primordium")
     grid();p<-par("usr"); rect(p[1],p[3],p[2],p[4],col="#6495ED25")
     title(paste("Sample ",mrt.y,"Tree Size = ",length(levels(as.factor(gr2)))))
     p<-par("usr"); rect(p[1],p[3],p[2],p[4],col="#6495ED25")
     suite<-rep(1,k)
     # Confidence interval. #arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
     abline(v=min(Data$X[which(gr2==1)]), col="red")
     palette(rainbow(ki))
     for(k in 1:(ki-1)) {
       if(k!=(ki-1)) {
         if(k!=1) { abline(v=min(Data$X[which(gr2==k)]), col="red", lty="dotted")}
         arrows(z[1,k],p[3]+0.1*(p[4]-p[3]),z[2,k], p[3]+0.1*(p[4]-p[3]), angle=90,code=3,col="black",lwd=1,length = 0.2)
         points (Data$X[which(gr2==k)],Data$Y[which(gr2==k)],pch=21,bg=k+1,col="black")
       } else {
         abline(v=min(Data$X[which(gr2==k)]), col="red", lty="dotted")
         arrows(z[1,k],p[3]+0.1*(p[4]-p[3]),z[2,k], p[3]+0.1*(p[4]-p[3]), angle=90,code=3,col="black",lwd=1,length = 0.2)
         abline(v=max(Data$X[which(gr2==k)]), col="red", lty="dotted")
         points (Data$X[which(gr2==k)],Data$Y[which(gr2==k)],pch=21,bg=k+1,col="black")
         abline(v=max(Data$X[which(gr2==k+1)]), col="red")
         points (Data$X[which(gr2==k+1)],Data$Y[which(gr2==k+1)],pch=21,bg=k+2,col="black")
       }
       rug(Resultss,ticksize=-0.05,col="darkgray")
     }
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
	 cat("\n")
     par(default.pars)
     #ans<-list(Mean=m,Boot=Resultss,Qant=z,OutSize=N-dim(Resultss)[1])
     ans<-list(Mean=m,Boot=Resultss,Qant=z,OutSize=N-dim(Resultss)[1])
   }
