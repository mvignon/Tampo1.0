if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"melt",
	"ggplot",
	"aes",
	"geom_histogram",
	"facet_wrap",
	"xlab",
	"ggtitle",
	"scale_fill_brewer",
	"geom_density"))
}

plotMRT <-
function(	aa,				# a model from the 'MRT.Boot.qplot' function
					type="density",	# "density for smooth density plot or "count" for raw data
					by="Split",		# Plot organized by "Split" or "Element". "Total" also available when using "count"
					col=1,			# Number of column to display
					palette="Set1",	# Colors for splits or elements. "Blues" "Set2" "Set3" "Dark2" "Pastel1" "Pastel2" "Paired" and "Accent"
					scale="free"	# Are scales shared across all facets ("fixed"), or do they vary across rows ("free_x"), columns ("free_y"), or both rows and columns ("free")
				)
{
  if(require("reshape")=="FALSE"){
    print("trying to install reshape")
    install.packages("reshape")
    if(require("reshape")){
      print("reshape installed and loaded")
    } else {
      stop("could not install reshape. Please try manually")
    }
  }
    if(require("reshape2")=="FALSE"){
    print("trying to install reshape2")
    install.packages("reshape2")
    if(require("reshape2")){
      print("reshape2 installed and loaded")
    } else {
      stop("could not install reshape2. Please try manually")
    }
  }
  if(require("ggplot2")=="FALSE"){
    print("trying to install ggplot2")
    install.packages("ggplot2")
    if(require("ggplot2")){
      print("ggplot2 installed and loaded")
    } else {
      stop("could not install ggplot2. Please try manually")
    }
  } 
  
  a<-aa$Data.Hist
  b<-melt(a,id=c("Element","Split"))
  b$Element<-as.factor(b$Element)
  b$Split<-as.factor(b$Split)
 value=Split=Element=NULL
  if(type=="count"){
   if (by=="Element"){
	ggplot(b,aes(x=value,fill=Split))+geom_histogram(alpha=0.6,binwidth=10)+facet_wrap(~Element,scales=scale,ncol=col)+xlab("Distance from primordium")+ ggtitle(paste("Counts for",aa$N,"Bootstraps"))+scale_fill_brewer(palette=palette)
  } else{
	if(by=="Split"){
		ggplot(b,aes(x=value,fill=Element))+geom_histogram(alpha=0.6,binwidth=10)+facet_wrap(~Split,scales=scale,ncol=col)+xlab("Distance from primordium")+ ggtitle(paste("Counts for",aa$N,"Bootstraps"))+scale_fill_brewer(palette=palette)
		} else{ #by should then be "Total"
			#qplot(Dat,data=aa$Data.Hist[which(aa$Data.Hist$Element!="All"),],geom='histogram',fill=Element,binwidth=10,main = paste("Bootstrap",aa$N), xlab = "Distance from primordium")
			ggplot(b[which(b$Element!="All"),],aes(x=value,fill=Element))+geom_histogram(alpha=0.6,binwidth=10)+xlab("Distance from primordium")+ ggtitle(paste("Counts for",aa$N,"Bootstraps"))+scale_fill_brewer(palette=palette)
		}
  
  }
  
  } else {
	if (by=="Element") { ggplot(b,aes(x=value,fill=Split))+geom_density(alpha=0.3)+facet_wrap(~Element,scales=scale,ncol=col)+xlab("Distance from primordium")+ ggtitle(paste("Density based on",aa$N,"Bootstraps"))+scale_fill_brewer(palette=palette)}
	else {ggplot(b,aes(x=value,fill=Element))+geom_density(alpha=0.3)+facet_wrap(~Split,scales=scale,ncol=col)+xlab("Distance from primordium")+ ggtitle(paste("Density based on",aa$N,"Bootstraps"))+scale_fill_brewer(palette=palette)}
	}
}
