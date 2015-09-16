if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"ggplot",
	"aes",
	"geom_tile",
	"geom_text",
	"scale_fill_gradientn",
	"geom_hline",
	"geom_vline",
	"labs",
	"facet_wrap",
	"scale_fill_gradient",
	"ggtitle"))
}

Mapp <-
function(
	model,		# The "MRT" object as obtained by the "MRT.Boot.qplot" function 
	Index="S",# The index to be used ("S","SS","J")
	what="Total",	# what to print. either "Total", "Split" or "Element".
	IC=0.90,                # The confident interval to be used.
	Values=T,	# Plot index values on the map
	write=F
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

	MRT<-MRT.Heatmap(model,Index=Index,IC=IC)
	model <-MRT$heat
	Overlap<-MRT$Index
	if(write==T) {
		write.table(Overlap,file="Overlapping Index.txt",sep="\t", col.names = factor(levels(model$X)), row.names = factor(levels(model$X)))
		cat("-----------------------------------------------------","\n")
		cat("text file successfully created in the following directory: ","\n")
		cat(getwd(),"\n")
	}
	X=Y=V=Sp1=Sp2=El1=El2=NULL
	if(what=="Total"){
		if(Values==T){
		p<-ggplot(model, aes(x=X, y=Y, fill=V))+geom_tile(colour = "white") +geom_text(aes(label = round(V, 1)))+ scale_fill_gradientn(colours = rev(heat.colors(10)),name = Index)+geom_hline(yintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+geom_vline(xintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+labs(x = "", y = "")
		} else p<-ggplot(model, aes(x=X, y=Y, fill=V))+geom_tile(colour = "white") + scale_fill_gradientn(colours = rev(heat.colors(10)),name = Index)+geom_hline(yintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+geom_vline(xintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+labs(x = "", y = "")
		#p<-ggplot(model, aes(x=X, y=Y, fill=V))+geom_tile(colour = "white")+ scale_fill_gradient(low = "red",high = "yellow",name = Index)+geom_hline(yintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+geom_vline(xintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+labs(x = "", y = "")
		#p<-ggplot(model, aes(x=X, y=Y, fill=V))+geom_tile(colour = "white")+ geom_text(aes(fill=V,label = round(V, 1)),size=model$V*5)+ scale_fill_gradient(low = "red",high = "yellow",name = Index)+geom_hline(yintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+geom_vline(xintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+labs(x = "", y = "")
		#p<-ggplot(model, aes(x=X, y=Y, fill=V))+geom_tile(colour = "white")+ geom_text(aes(label = round(V, 1)))+ scale_fill_gradient(low = "red",high = "yellow",name = Index)+geom_hline(yintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+geom_vline(xintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+labs(x = "", y = "")
		#p<-ggplot(model, aes(x=X, y=Y, fill=V))+geom_tile(colour = "white")+ geom_text(aes(label = round(V, 1)))+ scale_fill_gradient2(low = "white",mid="yellow",  high = "red",midpoint=0.5,name = Index)+geom_hline(yintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+geom_vline(xintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+labs(x = "", y = "")
		#p<-ggplot(model, aes(x=X, y=Y, fill=V))+geom_tile(colour = "white")+ geom_text(aes(label = round(V, 1)))+ scale_fill_gradientn(colours=c("cornflowerblue","darkmagenta","darkorange1"),name = Index)+geom_hline(yintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+geom_vline(xintercept=seq(from=model$f[1],to=model$t[1],by=model$b[1]),colour="black",size=1)+labs(x = "", y = "")
	
	#low = "white", high = "steelblue"
	#low = "blue",	  high = "yellow"
	#low="green", high="red"
	#scale_fill_gradient2(low = "white",mid="yellow",  high = "red",midpoint=0.5
	#scale_fill_gradientn(colours=c("navyblue","darkmagenta","darkorange1")
	

	}
	if(what=="Split"){
		p<-ggplot(model[which(model$El1==model$El2),], aes(x=Sp1, y=Sp2, fill=V))+ geom_tile()+geom_text(aes(label = round(V, 1)))+facet_wrap(~El1)+labs(x = "Splits", y = "Splits")+ scale_fill_gradient(low = "red",high = "yellow",name = Index)
	}
	if(what=="Element"){
		p<-ggplot(model[which(model$Sp1==model$Sp2),], aes(x=El1, y=El2, fill=V))+ geom_tile()+geom_text(aes(label = round(V, 1)))+facet_wrap(~Sp1)+labs(x = "Elements", y = "Elements")+ scale_fill_gradient(low = "red",high = "yellow",name = Index)
	}
	p+ggtitle("Overlapping Index")
	}
