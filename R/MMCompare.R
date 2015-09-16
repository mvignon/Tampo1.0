if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"interval",
	"interval_union",
	"interval_intersection",
	"interval_is_empty",
	"interval_measure"))
}

MMCompare <-
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
	names(Rmulti)<-names(model$BootMono)
	
	c<-model$BootMulti
	R22<-matrix(0,ncol=2,nrow=dim(c)[2]) # 2 colonnes (Quantile min-max)
	for (i in 1:dim(c)[2]) {
		R22[i,]<-quant(c[,i])
	}

	R222<-c()
	for (i in 1:dim(R22)[1]) {
			if(i==1) {transit<-interval(R22[i,1],R22[i,2])
			}else{transit<-interval_union(transit,interval(R22[i,1],R22[i,2]))}
		R222<-transit
	}
	names(R222)<-"Multi"
	
	Mono<-Rmulti
	Multi<-R222
	Overlap1<-matrix(NA,nrow=length(Mono),ncol=1)	# Mono/Multi
	Overlap2<-matrix(NA,nrow=length(Mono),ncol=1)	# Multi/Mono (for asymmetrical matrix)

	if (Index=="S") { # Indice asymétrique I(a,b)=overlap(a,b)/range(a)
		for (i in 1:length(Mono)) {
			aa<-interval_intersection(Mono[i],Multi)
			if (interval_is_empty(aa)){
				Overlap1[i,1]<-0
				Overlap2[i,1]<-0
			} else {
			aa<-interval_measure(aa)
			Overlap1[i,1]<-aa/interval_measure(Mono[[i]])
			Overlap2[i,1]<-aa/interval_measure(Multi)
			}
		}
		rownames(Overlap1)<-rownames(Overlap2)<-names(model$BootMono)
		c<-data.frame(SMonoMulti=Overlap1,SMultiMono=Overlap2)
		barplot(as.matrix(t(c)),beside=T,legend = c(expression(paste(S[Mono-Multi])),expression(paste(S[Multi-Mono]))),ylab="S",main="Congruency between Mono and Multi signals")
		mtext(paste("based on",model$N,"boostraps, for",nlevels(model$Data.Hist$Split),"splits"))
	}
	
	if (Index=="J") { # Indice symétrique Overlap as defined by the jaccard index
		for (i in 1:length(Mono)) {
			aa<-interval_intersection(Mono[i],Multi)
			if (interval_is_empty(aa)){
				Overlap1[i,1]<-0
				Overlap2[i,1]<-0
			} else {
			aa<-interval_measure(aa)
			Overlap1[i,1]<-aa/interval_measure(interval_union(Mono[i],Multi))
			}
		}
		rownames(Overlap1)<-names(model$BootMono)
		c<-data.frame(J=Overlap1)
		barplot(as.matrix(t(c)),ylab="J",main="Congruency between Mono and Multi signals")
		mtext(paste("based on",model$N,"boostraps, for",nlevels(model$Data.Hist$Split),"splits"))
	}
	
	if (Index=="SS") { # Indice asymétrique I(a,b)=overlap(a,b)/range(a)
		min<-min(model$Data.Hist$Dat)
		mina<-min+abs(min)
		maxa<-max(model$Data.Hist$Dat)+abs(min)
		TOT<-maxa-mina
		for (i in 1:length(Mono)) {
			aa<-interval_intersection(Mono[i],Multi)
			if (interval_is_empty(aa)){
				Overlap1[i,1]<-0
				Overlap2[i,1]<-0
			} else {
			aa<-interval_measure(aa)
			Overlap1[i,1]<-(aa/interval_measure(Mono[[i]]))*((TOT-interval_measure(Mono[[i]]))/TOT)
			Overlap2[i,1]<-(aa/interval_measure(Multi))*((TOT-interval_measure(Mono[[i]]))/TOT)
			}
		}
		rownames(Overlap1)<-rownames(Overlap2)<-names(model$BootMono)
		c<-data.frame(SSMonoMulti=Overlap1,SSMultiMono=Overlap2)
		barplot(as.matrix(t(c)),beside=T,legend = c(expression(paste(SS[Mono-Multi])),expression(paste(SS[Multi-Mono]))),ylab="SS",main="Congruency between Mono and Multi signals")
		mtext(paste("based on",model$N,"boostraps, for",nlevels(model$Data.Hist$Split),"splits"))
	}
	c
	ans<-list(Over1=Overlap1,Over2=Overlap2)
}
