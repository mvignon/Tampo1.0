if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"interval",
	"interval_union",
	"interval_intersection",
	"interval_is_empty",
	"interval_measure"))
}

MRT.Heatmap <-
function(
     model,                  # The "MRT" object as obtained by the "MRT.Boot.qplot" function 
     IC=0.90,                # The confident interval to be used.
	 Index="S"				 # The index to be used ("S","SS","J")
   )
   { 
	if(require("sets")=="FALSE"){
		print("trying to install sets")
		install.packages("sets")
		if(require("sets")){
		  print("sets installed and loaded")
		} else {
		  stop("could not install sets. Please try manually")
		}
	}
	t1 <- unclass(Sys.time())
	
	P<-c((1-IC)/2,1-(1-IC)/2)  # Percentiles as in IC
	quant<-function(x){		   # Quantile function
       quantile(x,P)
    }
	
	a<-model$BootMono			# Data
	L<-length(a)*dim(a[[1]])[2] # Side length of the heatmap to be created (map=L*L)
	R11<-matrix(0,ncol=2,nrow=L) # 2 colonnes (Quantile min-max) et L lignes organisées Split1(élément1,élé2,élé3,...)+Split2(élé1,élé2,...)+Split3(élé1,...
	b<-1
	for (j in 1:dim(a[[1]])[2]) {
		for (i in 1:length(a)) {
		R11[b,]<-quant(a[[i]][,j])
		b<-b+1
		}
	}
	
	R111<-matrix(0,ncol=2,nrow=L) # idem El1(Spl1,Spl2,...), El2(Spl1,...) utilisé pour Rmulti
	b<-1
	for (i in 1:length(a)) {
		for (j in 1:dim(a[[1]])[2]) {
		R111[b,]<-quant(a[[i]][,j])
		b<-b+1
		}
	}
	
	R11<-R11+abs(min(R11))	# Pour n'avoir que des valeurs positives
	colnames(R11)<-c("Qmin","Qmax")
	
	Rmulti<-list(NA)
	b<-1
	for (i in 1:length(a)) {
		for (j in 1:dim(a[[1]])[2]) {
			if(j==1){transit<-interval(R111[b,1],R111[b,2])
			b<-b+1
			}else{transit<-interval_union(transit,interval(R111[b,1],R111[b,2]))
			b<-b+1}
		}
		Rmulti[[i]]<-transit
	}
	
	progress<-txtProgressBar(min = 1, max = L, style = 3,char="#")
	
	################### Ordre des splits + éléments dans la matrice Overlap2
	Sp<-vector(length=dim(a[[1]])[2]) # Sequantial split
	for (i in 1:length(Sp)) {
	Sp[i]<-paste("Split",i)
	}
	Spp<-rep(Sp,each=length(a))
	El<-vector(length=length(a))	#Sequantial elements
	for (i in 1:length(El)) {
	El[i]<-names(a)[i]
	}
	Ell<-rep(El,times=length(Sp))
	
	RR<-matrix(0,ncol=2,nrow=L)	# store splits and elements Col1 for Splits and Col2 for Element
	RR<-cbind(Spp,Ell)
	colnames(RR)<-c("Splits","Elements")

	RRR<-matrix(0,ncol=1,nrow=L)
	for (i in 1:dim(RR)[1]) {
	RRR[i]<-paste(Spp[i],Ell[i])
	}
	
	###################### 
	Overlap<-matrix(NA,nrow=L,ncol=L)	# matrice L*L Overlap
	for (i in 1:L) {
		for (j in 1:L) {
			aa<-interval_intersection(interval(R11[i,1],R11[i,2]),interval(R11[j,1],R11[j,2]))
			if (interval_is_empty(aa)){
				Overlap[i,j]<-0
			} else {
			aa<-as.vector(aa)
			Overlap[i,j]<-range(aa)[2]-range(aa)[1]
			}
		}
		setTxtProgressBar(progress, i)
	}

	Range<-rep(NA,length=L)		# Range of values for each split and each element
	for (i in 1:L) {
		bb<-range(interval(R11[i,1],R11[i,2]))
		Range[i]<-range(bb)[2]-range(bb)[1]
	}
	
	if (Index=="S") {
		Overlap2<-matrix(NA,nrow=L,ncol=L)	# Indice asymétrique I(a,b)=overlap(a,b)/range(a)
		for (i in 1:L) {
				Overlap2[i,]<-Overlap[i,]/Range[i]
				}
	}
	if (Index=="J") {
		Overlap2<-matrix(NA,nrow=L,ncol=L)	# matrice L*L Overlap as defined by the jaccard index
		for (i in 1:L) {
			for (j in 1:L) {
				inter<-interval_intersection(interval(R11[i,1],R11[i,2]),interval(R11[j,1],R11[j,2]))
				union<-interval_union(interval(R11[i,1],R11[i,2]),interval(R11[j,1],R11[j,2]))
				if (interval_is_empty(inter)){
					Overlap2[i,j]<-0
				} else {
				Overlap2[i,j]<-interval_measure(inter)/interval_measure(union)
				#OverJac[i,j]<-interval_division(inter/union)
				}
			}
		}
	}

	if (Index=="SS") {
				Overlap2<-matrix(NA,nrow=L,ncol=L)	# Indice asymétrique I(a,b)=overlap(a,b)/range(a)
			min<-min(model$Data.Hist$Dat)
		mina<-min+abs(min)
		maxa<-max(model$Data.Hist$Dat)+abs(min)
		TOT<-maxa-mina
		for (i in 1:L) {
				Overlap2[i,]<-((Overlap[i,]/Range[i])*((TOT-Range[i])/TOT)) # Standardisation par la largeur de l'élément (évite fort index lorsque peut préçi eté tendu)
				}
	}
	colnames(Overlap2)<-rownames(Overlap2)<-paste(Spp,Ell)
	
	### Vector for Heatmap using ggplot2
	V<-vector()
	for (i in 1:dim(Overlap2)[1]) {
		V<-c(V,Overlap2[i,])
	}
	
	### Store heatmap data for using with ggplot2 ###
	heat<-data.frame(X=rep(RRR,each=dim(Overlap2)[1]),Y=rep(RRR,times=dim(Overlap2)[1]),V=V,Sp1=rep(Sp,each=length(V)/dim(a[[1]])[2]),Sp2=rep(Spp,times=dim(Overlap2)[1]),El1=rep(Ell,each=dim(Overlap2)[1]),El2=rep(Ell,times=dim(Overlap2)[1]),f=dim(Overlap2)[1]/nlevels(as.factor(Sp))+0.5,t=dim(Overlap2)[1]+0.5,b=dim(Overlap2)[1]/nlevels(as.factor(Sp)))
	### Store data for hierarchical clustering ###
	clust<-data.frame(Splits=Spp,Elements=Ell,SE=paste(Spp,Ell))
	
	t2 <- unclass(Sys.time())
    elapsed.time.minutes <- round((t2-t1),3)  #calculate the total elapsed time
	cat("-----------------------------","\n")
    cat("(elapsed time - ",round(elapsed.time.minutes,1),"seconds)","\n")
	ans<-list(Index=Overlap2,heat=heat,method=Index,Clust=clust,Multi=Rmulti)
}
