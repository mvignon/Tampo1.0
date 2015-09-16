if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"dendro_data",
	"ggplot",
	"segment",
	"geom_segment",
	"aes",
	"scale_y_continuous",
	"scale_x_continuous",
	"label",
	"geom_point",
	"geom_text",
	"theme",
	"element_text",
	"grid.newpage",
	"pushViewport",
	"viewport",
	"grid.layout"))
}

Dendro2 <-
function(
	model,		# The "MRT" object as obtained by the "MRT.Boot.qplot" function 
	Index="S",  # The index to be used ("S","SS","J")
	IC=0.90,    # The confident interval to be used.
	Mean=F		# Compute and display mean from asymmetrical matrix ("S" and "SS" only)
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
	
	if(require("ggdendro")=="FALSE"){
			print("trying to install ggdendro")
			install.packages("ggdendro")
			if(require("ggdendro")){
			  print("ggdendro installed and loaded")
			} else {
			  stop("could not install ggdendro. Please try manually")
			}
	}

	MRT<-MRT.Heatmap(model,Index=Index,IC=IC)
	if(Index!="J" & Mean==T){fit<-hclust(as.dist(1-MRT$Index),method="complete")
	} else fit<-hclust(as.dist(1-((MRT$Index+t(MRT$Index))/2)),method="complete")
	coph <- cor(cophenetic(fit), as.dist(1-MRT$Index))
	b<-as.dendrogram(fit)
	model <- dendro_data(b, type="rectangle")


	### FOR SPLITS
	rain<-rainbow(length(fit$label))
	c<-rep(0,length=length(fit$labels))
	Splits<-rep(0,length=length(fit$labels))
	for(i in 1:nlevels(MRT$heat$Sp1)){
	l<-grepl(levels(MRT$heat$Sp1)[i],fit$labels[fit$order])
	c[l]<-rain[i]
	Splits[l]<-levels(MRT$heat$Sp1)[i]
	}

	### FOR Elements
	rain<-rainbow(length(fit$labels))
	c<-rep(0,length=length(fit$labels))
	Elements<-rep(0,length=length(fit$labels))
	for(i in 1:nlevels(MRT$heat$El1)){
	l<-grepl(levels(MRT$heat$El1)[i],fit$labels[fit$order])
	c[l]<-rain[i]
	Elements[l]<-levels(MRT$heat$El1)[i]
	}
	x=y=xend=yend=NULL
	p<-ggplot(segment(model)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))  + scale_y_continuous(name=paste("Index",Index),breaks=seq(0, 1, 0.1)) + scale_x_continuous(name="Splits/Elements",breaks=NULL)
	#+coord_flip()
	ddata_x<-dendro_data(b)
	model<-label(ddata_x)
	model<-cbind(model,Splits,Elements)
	
	p1<-p+ geom_point(data=model, size=5,aes(label=Splits, x=x, y=-0.1,col=Splits),show_guide =T)+ geom_text(data=model,size=4, aes(label=Elements, x=x, y=-0.2,col=Splits),show_guide =F)
	p2<-p+ geom_point(data=model, size=5,aes(label=Elements, x=x, y=-0.1,col=Elements),show_guide =T)+ geom_text(data=model,size=4, aes(label=Splits, x=x, y=-0.2,col=Elements),show_guide =F)

	p11<- p1+theme(axis.text.x  = element_text(angle=45))
	p22<- p2+theme(axis.text.x  = element_text(angle=45))
	multiplot(p11, p22, cols=2)
	cat("Cophenetic coefficient correlation: ","\n")
	cat(coph,"\n")
	
	}


multiplot <-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

