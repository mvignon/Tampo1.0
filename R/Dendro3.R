if(getRversion() >= "2.15.1") {utils::globalVariables(c(
	"interval_is_empty",
	"interval_intersection",
	"interval_measure",
	"interval_union",
	"dendro_data",
	"ggplot",
	"segment",
	"geom_segment",
	"aes",
	"scale_y_continuous",
	"scale_x_continuous",
	"label",
	"geom_text",
	"label",
	"grid.newpage",
	"pushViewport",
	"viewport",
	"grid.layout"))
}

Dendro3 <-
function(
	model,		# The "MRT" object as obtained by the "MRT.Boot.qplot" function 
	Index="S",  # The index to be used ("J")
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
	Multi<-MRT$Multi
	n<-nlevels(MRT$heat$El1)
	m<-matrix(ncol=n,nrow=n)
	progress<-txtProgressBar(min = 1, max = n, style = 3,char="#")
	for (i in 1:n){
		for (j in 1:n){
			if (interval_is_empty(interval_intersection(Multi[[i]],Multi[[j]]))){
				m[i,j]<-0
				} else {
				if(Index=="J") m[i,j]<-interval_measure(interval_intersection(Multi[[i]],Multi[[j]]))/interval_measure(interval_union(Multi[[i]],Multi[[j]]))
				if(Index=="S") m[i,j]<-interval_measure(interval_intersection(Multi[[i]],Multi[[j]]))/interval_measure(Multi[[i]])
				if(Index=="SS") {
					min<-min(model$Data.Hist$Dat)
					mina<-min+abs(min)
					maxa<-max(model$Data.Hist$Dat)+abs(min)
					TOT<-maxa-mina
					m[i,j]<-(interval_measure(interval_intersection(Multi[[i]],Multi[[j]]))/interval_measure(Multi[[i]]))*((TOT-interval_measure(Multi[[i]]))/TOT)
				}
			}
		}
		setTxtProgressBar(progress, i)
}
	close(progress)
	colnames(m)<-rownames(m)<-levels(MRT$heat$El1)
	
	if(Index!="J" & Mean==T){fit<-hclust(as.dist(1-((m+t(m))/2)),method="complete")
	} else fit<-hclust(as.dist(1-m),method="complete")
	coph <- cor(cophenetic(fit), as.dist(1-m))
	b<-as.dendrogram(fit)
	model <- dendro_data(b, type="rectangle")
	x=y=xend=yend=NULL
	p<-ggplot(segment(model)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))  + scale_y_continuous(name=paste("Index 1-",Index),breaks=seq(0, 1, 0.1)) + scale_x_continuous(name="Elements",breaks=NULL)
	#+coord_flip()
	ddata_x<-dendro_data(b)
	model<-label(ddata_x)
	# model<-cbind(model,Elements)
	p2<-p+ geom_text(data=model,size=4, aes(label=label, x=x, y=-0.05),show_guide =F)
	multiplot(p2, cols=1)
	#p2
	cat("Cophenetic coefficient correlation: ","\n")
	cat(coph,"\n")
	ans<-list(dist=m)
	}

	
	multiplot <-
function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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