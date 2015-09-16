if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"ggplot",
	"aes",
	"geom_point",
	"facet_wrap",
	"scale_x_discrete",
	"geom_line",
	"geom_boxplot"))
}

Spe <-
function(
	model,		# The "MRT" object as obtained by the "MRT.Boot.qplot" function 
	IC=0.90    # The confident interval to be used.
	)
{	
	if(require("ggplot2")=="FALSE"){
			print("trying to install ggplot2")
			install.packages("ggplot2")
			if(require("ggplot2")){
			  print("ggplot2 installed and loaded")
			} else {
			  stop("could not install ggplot2. Please try manually")
			}
	}
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
	R11<-matrix(0,ncol=dim(a[[1]])[2],nrow=length(a)) # elements in lines and splits in columns
	for (i in 1:length(a)) {
		for (j in 1:dim(a[[1]])[2]) {
			b<-quant(a[[i]][,j])
			b<-b+abs(min(b))
			bb<-range(b)
			R11[i,j]<-range(bb)[2]-range(bb)[1]
		}
	}
	
	Sp<-vector(length=dim(a[[1]])[2]) # Sequantial split
	for (i in 1:length(Sp)) {
	Sp[i]<-paste("Split",i)
	}
	Sp2<-vector(length=dim(a[[1]])[2]) # Sequantial split
	for (i in 1:length(Sp)) {
	Sp2[i]<-i
	}
	Spp<-rep(Sp,each=length(a))
	Spp2<-rep(Sp2,each=length(a))
	El<-vector(length=length(a))	#Sequantial elements
	for (i in 1:length(El)) {
	El[i]<-names(a)[i]
	}
	Ell<-rep(El,times=length(Sp))
	
	b<-model$BootMulti	
	R22<-matrix(0,ncol=dim(b)[2],nrow=1) #  splits in columns
	for (i in 1:dim(b)[2]) {
			c<-quant(b[,i])
			c<-c+abs(min(c))
			bb<-range(c)
			R22[1,i]<-range(bb)[2]-range(bb)[1]
		}
	
	E<-rep("All",times=dim(b)[2])
	S<-vector(length=dim(b)[2])
	for (i in 1:dim(b)[2]) {
		S[i]<-paste("Split",i)
	}
	S2<-vector(length=dim(b)[2])
	for (i in 1:dim(b)[2]) {
		S2[i]<-i
	}
	# RR<-matrix(0,ncol=2,nrow=L)	# store splits and elements Col1 for Splits and Col2 for Element
	# RR<-cbind(Spp,Ell)
	# colnames(RR)<-c("Splits","Elements")

	# RRR<-matrix(0,ncol=1,nrow=L)
	# for (i in 1:dim(RR)[1]) {
	# RRR[i]<-paste(Spp[i],Ell[i])
	# }
	Range=NULL
	model<-data.frame(Range=c(c(R11),c(R22)),El=as.factor(c(Ell,E)),Sp=as.factor(c(Spp2,S2)))
	#dotchart(R11)
	# p=ggplot(qq, aes(Range,color=El))
	# p+geom_point()
	#p1<-ggplot(model) + geom_point(aes(x=Sp,y=Range)) + facet_wrap(~El) 
	p1<-ggplot(model,aes(x=Sp,y=Range)) + geom_point() + facet_wrap(~El) + scale_x_discrete(name="Splits")+ geom_line(aes(group = 1))#+ stat_smooth(aes(group = 1),method=loess)# 
	p2<-ggplot(model,aes(x=Sp,y=Range)) +geom_boxplot(aes(group=Sp))+ scale_x_discrete(name="Splits")
	multiplot(p1,p2, cols=1)
	ans<-list(data=model)
}
