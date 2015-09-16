if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"gam"))
}

Gam2Hab <-
function(
  mrt.x,                  # Column name for the distance to the nucleus.
  mrt.y,                  # Column name for the elemental concentration.
  data,                   # Name of the dataframe to analyse (must be charged).
  IC=0.95,                # Confidence interval for the fitted model [0-1].
  Habitat1=FALSE,         # Threshold value to be used for discriminating between 2 habitats.
  Habitat2=FALSE,         # Threshold values to be used for discriminating between 2 habitats.
  strict="FALSE",         # Indicate whether IC is ignored to define habitat changes. If "TRUE", no equivocal zone appears and habitat changes correspond to intersection of threshold and fitted gam.
  point="TRUE",           # Plot original data. 
  grid="FALSE",           # Plot a background grid.
  ylab="Elemental ratio" # Legend to be used for the Y axis.
  )
{
  # Preliminary (package & arguments) 
  if (IC>1|IC<0){
    stop("IC must be between 0-1")
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
 if (is.numeric(Habitat1)) {
       if(length(Habitat1)>1) {
         stop("Habitat1 must contain only one threshold value")
       }
     }
     if(length(Habitat2)==2) {
       if(Habitat2[1]>Habitat2[2]) {
         Habitat2[1]->a;Habitat2[1]<-Habitat2[2];Habitat2[2]<-a    
       }
     }
  if(require("mgcv")=="FALSE"){
    print("trying to install mgcv")
    install.packages("mgcv")
    if(require("mgcv")){
      print("mgcv installed and loaded")
    } else {
      stop("could not install mgcv. Please try manually")
    }
  }
  # Data acquisition & Modelling  
  y.data<-na.omit(data[,mrt.y])
  x.data<-data[which(data[,mrt.y]!="NA"),mrt.x]
  if (length(x.data)!=length(y.data)){
    stop(paste("ERROR: 'x' and 'y' differ in length","\n"))
  }
    Data<-as.data.frame(cbind(x.data,y.data))
    names(Data)<-c("X","Y")
    REG=gam(Y~s(X,k=30),data=Data,family=gaussian(link="identity"))
    Predict<-predict(REG,se=TRUE)
    fit<-Predict$fit
    se<-Predict$se.fit
    lcl<-fit-(qnorm(IC+((1-IC)/2)))*se
    ucl<-fit+(qnorm(IC+((1-IC)/2)))*se
    # Basic Plottting
if (any(is.numeric(Habitat1),is.numeric(Habitat2))) {
        layout(matrix(c(1,1,1,2),4,1))
        par(oma=c(0,0,0,0))
        #par(oma=c(0,0,0,0),mar=c(1,1,1,1))
    } else {layout(1,1,1)} 
if(point=="TRUE"){
    plot(0,type="n",xlab="Distance from primordium",ylab="Elemental ratio",xlim=c(min(x.data),max(x.data)),ylim=c(min(y.data),max(y.data)))
    }else{
plot(0,type="n",xlab="Distance from primordium",ylab="Elemental ratio",xlim=c(min(x.data),max(x.data)),ylim=c(min(lcl),max(ucl)))
}
x<-par("usr")
    if (grid=="TRUE"){
      grid()  
    }
    polygon(c(x.data[order(x.data)],x.data[order(x.data,decreasing=T)]),c(ucl[order(x.data)],lcl[order(x.data,decreasing=T)]),col="#6495ED25",border=NA)
    lines(x.data,fit,lwd=1)
if(is.numeric(Habitat1)){
abline(h=Habitat1,lty=2,col="red")
} else{
abline(h=Habitat2,lty=2,col="red")
}
    lines(x.data,ucl,lwd=0.3,lty=3)
    lines(x.data,lcl,lwd=0.3,lty=3)
    if (point=="TRUE"){
      points(Data$X,Data$Y)  
    }
    # Habitat
    # Change
    Change<-rep(0,length=length(Data$X));b<-length(Data$X)-1
Hab<-rep(0,length=b)
if(is.numeric(Habitat1)){
if(strict=="TRUE") {
  for(i in 1:b){
if (fit[i]>Habitat1){Hab[i]<-1
} else {if (fit[i]<Habitat1){Hab[i]<-2}
}
if (fit[i]>Habitat1 & fit[i+1]<Habitat1){
  Change[i+1]<-1
}
if (fit[i]<Habitat1 & fit[i+1]>Habitat1){
  Change[i+1]<-1
}
  }
} else {
  Slidingthreshold<-as.data.frame(cbind(Habitat1-(fit-lcl),Habitat1+(ucl-fit)));names(Slidingthreshold)<-c("Low","Up")
  for(i in 1:b){
if (fit[i]>Slidingthreshold$Up[i]){Hab[i]<-1
} else {if(fit[i]<Slidingthreshold$Low[i]){Hab[i]<-2}
}
if (fit[i]<Slidingthreshold$Low[i] & fit[i+1]>Slidingthreshold$Low[i+1]){
  Change[i+1]<-1
}
if (fit[i]<Slidingthreshold$Up[i] & fit[i+1]>Slidingthreshold$Up[i+1]){
  Change[i+1]<-1
}
if (fit[i]>Slidingthreshold$Up[i] & fit[i+1]<Slidingthreshold$Up[i+1]){
  Change[i+1]<-1
}
if (fit[i]>Slidingthreshold$Low[i] & fit[i+1]<Slidingthreshold$Low[i+1]){
  Change[i+1]<-1
}
  }
}
rec<-rep(0,length(Change))
init<-1
for (i in 1:length(Change)){
if(Change[i]==1){
rec[i]<-init+1
init<-init+1
} else {rec[i]<-init}}
}
if(is.numeric(Habitat2)){
if(strict=="TRUE") {
  for(i in 1:b){
if (fit[i]>Habitat2[2]){Hab[i]<-1
} else {if (fit[i]<Habitat2[1]){Hab[i]<-2}
}
  }
} else {
  Slidingthreshold<-as.data.frame(cbind(Habitat2[1]-(fit-lcl),Habitat2[2]+(ucl-fit)));names(Slidingthreshold)<-c("Low","Up")
  for(i in 1:b){
if (fit[i]>Slidingthreshold$Up[i]){Hab[i]<-1
} else {if(fit[i]<Slidingthreshold$Low[i]){Hab[i]<-2}
}
  }
}
for(i in 1:(b-1)){
if(Hab[i+1]!=Hab[i]){
Change[i+1]<-1
}
}



rec<-rep(0,length(Change))
init<-1
for (i in 1:length(Change)){
if(Change[i]==1){
rec[i]<-init+1
init<-init+1
} else {rec[i]<-init}}
}


############
    # Plotting #
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
 
if (any((Habitat1!="FALSE"|Habitat2!="FALSE"))) {
if(is.numeric(Habitat1)){
abline(h=Habitat1,lwd=4,col="White")
abline(h=Habitat1,lty=6,lwd=2)
   } else{
   abline(h=Habitat2,lwd=4,col="White")
abline(h=Habitat2,lty=6,lwd=2)
   }
           par(mar=c(5,2,0,2))
           plot(x.data,y.data,axes=F,ann=F,ylim=c(0,2),col="white",main="Habitat use");axis(1);title(xlab="Distance from primordium")
           legend("top",c("Habitat 1","Habitat 2","Unassigned"),horiz=T,fill=c("red","blue","grey"))
           for(l in 1:nlevels(as.factor(rec))) {
             if(any(Hab[which(rec[-length(rec)]==l)]==1)) {
               reds <- colorRampPalette(c("#8B0000", "#FFAEB9", "#FF0000"))
               h<-seq(from=0,to=1,length=50)
               for(m in 1:50) {
                 polygon(c(min(Data$X[which(rec==l)]),max(Data$X[which(rec==l)])),c(h[m],h[m]),border=reds(50)[m],col=NA)
               }
             } else {
               if(any(Hab[which(rec[-length(rec)]==l)]==2)){
                 blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 h<-seq(from=0,to=1,length=50)
                 for(m in 1:50) {
                   polygon(c(min(Data$X[which(rec==l)]),max(Data$X[which(rec==l)])),c(h[m],h[m]),border=blues(50)[m],col=NA)
                 }
               } else {
                 greys<-colorRampPalette(c("#3D3D3D", "#CCCCCC", "#6B6B6B"))
                 h<-seq(from=0.1,to=0.9,length=50)
                 for(m in 1:50) {
                   polygon(c(min(Data$X[which(rec==l)]),max(Data$X[which(rec==l)])),c(h[m],h[m]),border=greys(50)[m],col=NA)
                 }
               }
             }
           }
         }


  an<-as.numeric(x.data[which(Change!=0)])
  cat("Changes identified at: ","\n")
  for(i in 1:length(an)){
    if(i!=length(an)){
    cat("    ");cat(an[i],"\n")
    } else {
      cat("    ");cat(an[i]);cat(" ( from primordium)","\n")
      }
    }
  ans<-list(Change=as.numeric(x.data[which(Change!=0)]),Habitat=cbind(x.data[-length(x.data)],Hab))
}
