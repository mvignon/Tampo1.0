if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"snip.rpart",
	"prune"))
}
 
MRT.Mean<-function(
     mrt.x,                  # Column name for the distance to the nucleus.
     mrt.y,                  # Column name for the elemental concentration.
     data,                   # Name of the dataframe to analyse (must be charged).
     Cond=1,                 # Min. difference requested between two consecutive groups.
     Size,                   # Maximum tree size (number of groups) to be tested.
     Graph=FALSE,            # Supplementary graphs (optional).
     Habitat1=FALSE,         # Plot habitat use according to the given threshold.
     Habitat2=FALSE          # Plot habitat use according to the two given thresholds.
   )
   {
     
     ######################################
     # 1) Preliminary (package & arguments) 
     ######################################
     default.pars=par(no.readonly = TRUE)
     t1 <- unclass(Sys.time())
     
     if(missing(data)){
       stop("Data is missing")
     }
     if (Cond <= 0){
       stop("'Cond' must be > 0")
     }
     if(missing(Size)){
	} else {
        if (Size <= 0){
       stop("'Size' must be > 0")
	}
	}
     if (Cond> max(na.omit(data[,2]))-min(na.omit(data[,2]))){
       stop("WARNING: Condition out of range !!! Please, try another reduced value")
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

     
     #windows() ## Peut etre ?? mettre ailleurs (plus tard) ????
     
     ######################################
     # 2)  Data acquisition & Modelling
     ######################################
     
     ######################################
     # 2a) Basic model with maximum tree size (according to the 'Size' argument (default=10).
     # Data acquisition & Modelling
     y.data<-na.omit(data[,mrt.y])
     x.data<-data[which(data[,mrt.y]!="NA"),mrt.x]
     if (length(x.data)!=length(y.data)){
       cat(paste("ERROR: 'x' and 'y' differ in length","\n"))
     }
     Data<-as.data.frame(cbind(x.data,y.data))
     names(Data)<-c("X","Y")

     if(missing(Size)) {miss<-1
	} else {miss<-0}
     model<-mvpart(Y~X,data=Data,size=Size,cp=0.00001)
	Size<-max(model$cptable[,2])+1      # Verify the tree size obtained (should be equal to the 'Size' argument).
     
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
     ######################################
     # 2b) Computing the whole range of models from 2 to 'Size' by pruning the fully grown basic tree.
     ME<-matrix(NA,nrow=Size,ncol=Size)  # Matrix for group means (ncol=nrow= max. number of groups to be tested). Note that the first column will remain empty (only one group).
     SD<-matrix(NA,nrow=Size,ncol=Size)  # Matrix for group SD (identical dimension).
     MEE<-rep(NA,length=Size)				    # Matrix to verify the number of groups (nb. splits +1) within each pruned tree (should normally range from 2 to 'Size').
     for(j in 2:Size){
       if(j==2){
         Prunemodel<- snip.rpart(model, toss=c(2,3))   # Only keep the first node (i.e. first two groups).
       } else {
         CP<-model$frame[which(round(model$frame[,"complexity"],5)==round(model$cptable[min(which(model$cptable[,2]>=j-1)),1],5)),6]   # Needed because 'model$cptable' and 'model$frame' are not rounded using the same decimal number.
         Prunemodel<- prune(model, cp=CP[1])   # Prune the basic model at the requested tree size.
       }
       P<-predict(Prunemodel)          # Predicted values.
       gr<-Prunemodel$where                    # To renumber clusters sequentially.
       nb<-1;gr2<-rep(1,length(gr))
       for(k in 2:length(gr)) {
         if (gr[k]!=gr[k-1]) nb <-nb+1
         gr2[k] <- nb
       }
       ki<-max(as.numeric(gr2))                # Verify the number of groups (should be equal to 'Size').
       for(k in 1:ki) {                        # Fill the ME, SD & MEE matrices for each pruned model.
         ME[k,j]<-mean(Data$Y[which(gr2==k)])
         SD[k,j]<-sd(Data$Y[which(gr2==k)])
         MEE[j]<-length(na.omit(ME[,j]))
       }
     }
     
     ######################################
     # 2c) Testing conditions (min. mean difference as expressed in the 'Cond' argument) between pairs of consecutive groups.
     Limit<-matrix(NA,nrow=Size,ncol=Size)     # Condition matrix to be filled.
     for(k in 2:Size){
       M<-MEE-1
       for(l in 1:M[k]){
         Limit[l,k]<-ifelse(abs(ME[l+1,k]-ME[l,k])>=Cond,"TRUE","FALSE")
       }
     }
     c<-rep(NA,length=Size)                    # Find which tree size always met the condition.
     for(k in 2:Size){
       c[k]<-all(na.omit(Limit[,k])=="TRUE")
     }
     if(all(na.omit(c)=="TRUE")){
       TreeSize<-Size;Condition<-1} else{
         if(all(na.omit(c)=="FALSE")){
           TreeSize<-0;Condition<-2} else{
             i<-3;while ( length(which(na.omit(Limit[,i]=="FALSE")))<1 ) {
               TreeSize<-i;i<-i+1;Condition<-0
             }
           }
       }
     # Condition=1 when the condition is never restrictive (TreeSize=Size).
     # Condition=2 when the condition is always restrictive (no groups possible, TreeSize=0).
     # Condition=0 when the condition is sometime restrictive (within the 2:Size range, TreeSize=i).
     TreeSize;#(NbSplits<-TreeSize-1).
     
     ######################################
     # 2d) Computing optimized model (size='TreeSize') that met the conditions.
     if(Condition==2){
       graphics.off()
       stop("WARNING: No tree size met the condition! Please, consider decreasing the 'Cond' argument")} else{
	if(miss==1) {
       modelOpt<-mvpart(Y~X,data=Data,size=TreeSize,xv="1se")
      } else {
	modelOpt<-mvpart(Y~X,data=Data,size=TreeSize,cp=0.00001)
	}
         POpt<-predict(modelOpt)
         gr<-modelOpt$where
         aa<-1;gr2<-rep(1,length(gr))
         for(k in 2:length(gr)) {
           if (gr[k]!=gr[k-1]) aa <-aa+1
           gr2[k] <- aa
         }
         ki<-length(levels(as.factor(gr2)))            # Verify the number of groups (should be equal to 'TreeSize').
         SDOpt<-rep(1,length(levels(as.factor(gr2))))  # Matrix for group SD.
         for(k in 1:length(SDOpt)) {
           SDOpt[k]<-sd(Data$Y[which(gr2==k)])
         }
         MEOpt<-rep(1,length(levels(as.factor(gr2))))  # Matrix for group means.
         for(k in 1:length(MEOpt)) {
           MEOpt[k]<-mean(Data$Y[which(gr2==k)])
         }
         
         ######################################
         # 3)  Plotting the optimized model
         ######################################
         
         ######################################
         # 3a) Plotting background and title.
         #windows()
         if (any(is.numeric(Habitat1),is.numeric(Habitat2))) {
         layout(matrix(c(1,1,1,2),4,1))
         par(oma=c(0,0,0,0))
           #par(oma=c(0,0,0,0),mar=c(1,1,1,1))
         } else {layout(1,1,1)}
         plot(Data$X,Data$Y,col="white",ylab=expression(paste("Sr:Ca (x",10^-3,")")),xlab="Distance from primordium")
         grid();x<-par("usr"); rect(x[1],x[3],x[2],x[4],col="#6495ED25")
         if(Condition==1){
           title(paste("Sample",mrt.y,"\n",paste("Tree Size = ",length(levels(as.factor(gr2))))))
           mtext(paste("Note: The condition is never restrictive"),side=1,line=3,padj=0,adj=1,cex=0.65)
         } else{
           title(paste("Sample",mrt.y,"\n",paste("Optimal Tree Size = ",TreeSize)))
           mtext(paste("IMPORTANT: The condition was restrictive in creating this optimized tree"),side=1,line=3,padj=0,adj=1,cex=0.65)
         }
         #abline(v=(seq(0,100,25)), col="lightgray", lty="dotted") 
         #abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
         x<-par("usr"); rect(x[1],x[3],x[2],x[4],col="#6495ED25")
         
         ######################################
         # 3b) Plotting raw+predicted values & Mean+-SD for each group.
         suite<-rep(1,k)
         for(k in 1:length(suite)) {
           suite[k]<-length(Data$Y[which(gr2==k)])
         }
         qq<-rep(SDOpt,suite)
		 N<-rep(0,length=k)
         Change<-rep(0,length=ki-1)
         #		polygon(c(min(xx),xx,max(xx)),c(0,P+qq,0),col="grey")
         #		polygon(c(min(xx),xx,max(xx)),c(0,P-qq,0),col="white")
         for(k in 1:ki) {
           if(k!=ki){
             abline(v=min(Data$X[which(gr2==k)]), col="red", lty="dotted")
             points(min(Data$X[which(gr2==k)]),x[3], col="red",pch=17,cex=2)
             Change[k]<-max(Data$X[which(gr2==k)])
			 N[k]<-length(Data$X[which(gr2==k)])
           } else{
             abline(v=min(Data$X[which(gr2==k)]), col="red", lty="dotted")
             points(min(Data$X[which(gr2==k)]),x[3], col="red",pch=17,cex=2)
             abline(v=max(Data$X[which(gr2==k)]), col="red", lty="dotted")
             points(max(Data$X[which(gr2==k)]),x[3], col="red",pch=17,cex=2)
			 N[k]<-length(Data$X[which(gr2==k)])
           }
           rect(min(Data$X[which(gr2==k)]),min(POpt[which(gr2==k)]-qq[which(gr2==k)]),max(Data$X[which(gr2==k)]),max(POpt[which(gr2==k)]+qq[which(gr2==k)]),col="#6495ED50")
           palette(rainbow(ki))
           points (Data$X[which(gr2==k)],Data$Y[which(gr2==k)],col="red",pch=21,bg=k+1)
           #polygon(c(min(xx[which(gr2==i)]),xx[which(gr2==i)],max(xx[which(gr2==i)])),c(0,P[which(gr2==i)]+qq[which(gr2==i)],0),col="grey",border=NA)
           #polygon(c(min(xx[which(gr2==i)]),xx[which(gr2==i)],max(xx[which(gr2==i)])),c(0,P[which(gr2==i)]-qq[which(gr2==i)],0),col="white",border=NA)
           #lines(xx[which(gr2==k)],rep(ME[i],length(xx[which(gr2==k)])),col=k)
         }
         lines(Data$X,POpt)
         
         ######################################
         # 3c) Plotting habitat use (Optional)
         if (Habitat1!="FALSE") {
           abline(h=Habitat1,lwd=4,col="White")
           abline(h=Habitat1,lty=6,lwd=2)
           par(mar=c(5,2,0,2))
           plot(Data$X,Data$X,axes=F,ann=F,ylim=c(0,2),col="white",main="Habitat use");axis(1);title(xlab="Distance from primordium")
           legend("top",c("Habitat 1","Habitat 2","Unassigned"),horiz=T,fill=c("red","blue","grey"))
           for(k in 1:ki) {
             if(min(POpt[which(gr2==k)]-qq[which(gr2==k)])>Habitat1) {
               #rect(min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col="red")
               #             reds <- colorRampPalette(c('coral2', 'red3')) 
               reds <- colorRampPalette(c("#8B0000", "#FFAEB9", "#FF0000"))
               h<-seq(from=0,to=1,length=50)
               for(l in 1:50) {
                 polygon(c(min(Data$X[which(gr2==k)]),max(Data$X[which(gr2==k)])),c(h[l],h[l]),border=reds(50)[l],col=NA)
               }
               #             gradient.rect (min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col=reds(50))
             } else {
               if(max(POpt[which(gr2==k)]+qq[which(gr2==k)])<Habitat1){
                 #               blues <- colorRampPalette(c('dark blue', 'light blue')) 
                 #               blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 h<-seq(from=0,to=1,length=50)
                 for(l in 1:50) {
                   polygon(c(min(Data$X[which(gr2==k)]),max(Data$X[which(gr2==k)])),c(h[l],h[l]),border=blues(50)[l],col=NA)
                 }
                 #rect(min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col="blue")
                 #               gradient.rect (min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col=blues(50))
               } else {
                 #               greys <- colorRampPalette(c('slategray4', 'slategray1')) 
                 greys<-colorRampPalette(c("#3D3D3D", "#CCCCCC", "#6B6B6B"))
                 h<-seq(from=0.1,to=0.9,length=50)
                 for(l in 1:50) {
                   polygon(c(min(Data$X[which(gr2==k)]),max(Data$X[which(gr2==k)])),c(h[l],h[l]),border=greys(50)[l],col=NA)
                 }
                 #               gradient.rect (min(xx[which(gr2==k)]),0.1,max(xx[which(gr2==k)]),0.9,col=greys(50))
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
           for(k in 1:ki) {
             if(min(POpt[which(gr2==k)]-qq[which(gr2==k)])>Habitat2[2]){
               #rect(min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col="red")
               #             reds <- colorRampPalette(c('coral2', 'red3')) 
               reds <- colorRampPalette(c("#8B0000", "#FFAEB9", "#FF0000"))
               h<-seq(from=0,to=1,length=50)
               for(l in 1:50) {
                 polygon(c(min(Data$X[which(gr2==k)]),max(Data$X[which(gr2==k)])),c(h[l],h[l]),border=reds(50)[l],col=NA)
               }
               #             gradient.rect (min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col=reds(50))
             } else {
               if(max(POpt[which(gr2==k)]+qq[which(gr2==k)])<Habitat2[1]){
                 #               blues <- colorRampPalette(c('dark blue', 'light blue')) 
                 #               blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 blues <- colorRampPalette(c("#27408B", "#D1EEEE", "#4169E1"))
                 h<-seq(from=0,to=1,length=50)
                 for(l in 1:50) {
                   polygon(c(min(Data$X[which(gr2==k)]),max(Data$X[which(gr2==k)])),c(h[l],h[l]),border=blues(50)[l],col=NA)
                 }
                 #rect(min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col="blue")
                 #               gradient.rect (min(xx[which(gr2==k)]),0,max(xx[which(gr2==k)]),1,col=blues(50))
               } else {
                 #               greys <- colorRampPalette(c('slategray4', 'slategray1')) 
                 greys<-colorRampPalette(c("#3D3D3D", "#CCCCCC", "#6B6B6B"))
                 h<-seq(from=0.1,to=0.9,length=50)
                 for(l in 1:50) {
                   polygon(c(min(Data$X[which(gr2==k)]),max(Data$X[which(gr2==k)])),c(h[l],h[l]),border=greys(50)[l],col=NA)
                 }
                 #               gradient.rect (min(xx[which(gr2==k)]),0.1,max(xx[which(gr2==k)]),0.9,col=greys(50))
               }
             }
           }
         }
         
         ######################################
         # 3d) R-square and Relative Error plot, versus the number of splits (for the optimized tree only) 
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
             legend(1,0.2 , c("Apparent", "X Relative"), lty = 1:2, col="black",pch=21,cex=1,box.col="black",pt.bg=c("white","grey"),bg="white",bty="0")
             ylim <- c(min(xerror - xstd) - 0.1, max(xerror + xstd) + 
               0.1)
             plot(nsplit, xerror, xlab = "Number of Splits", ylab = "X Relative Error", ylim = ylim, type = "l")
             grid()
             A<-par("usr"); rect(A[1],A[3],A[2],A[4],col="#6495ED25")
             points(nsplit, xerror,pch=21,col="black",bg="grey")
             segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
             invisible()
           }
           dev.new();par(mfrow=c(1,2));rsqrpart(modelOpt) 
         }
		par(new=TRUE)
		on.exit(par(default.pars))
    }

     if(Condition!=2){
	 t2 <- unclass(Sys.time())
     elapsed.time.minutes <- round((t2-t1),3)  #calculate the total elapsed time
     cat("------------------------------","\n")
     cat("(elapsed time - ",round(elapsed.time.minutes,1),"seconds)","\n")
     cat("---------------------------------","\n")
     cat("Internal partitions identified at: ","\n")
     for(i in 1:ki-1){
       if(i!=k-1){
         cat("    ");cat(Change[i],"\n")
       } else {
         cat("    ");cat(Change[i]);cat(" (from primordium)")
       }
     }
	 }
	 cat("\n")
     ans<-list(Change=Change,Model=modelOpt,Mean=MEOpt,SD=SDOpt,Belong=gr2,N=N,Predict=POpt,Size=ki)
   }
   
  