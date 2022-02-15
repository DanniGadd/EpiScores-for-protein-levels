Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

### read in the supplementary table - worksheet saved as a CSV file from Danni's supp files ###
d = read.csv("U:/Datastore/IGMM/marioni-lab/Riccardo/Danni_paper_EpiScores_08March2021/Supplementary_Tables_040321_DG_fully_vs_wbc_vs_grim.csv")

### multiplot function ###
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

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

### end of multiplot function ###



library(ggplot2)

plot1 <- ggplot(d, aes(log(Hazard.Ratio), log(WBC.sensitivity.Hazard.Ratio)), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("Fully-adjusted log(Hazard Ratio)") +
      ylab("Full + WBC-adjusted log(Hazard Ratio)") +
      theme(plot.title = element_text(size = 11)) +
  geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=d[d$WBC.result!="attenuated",], shape=21, fill="white", size=2)+
  geom_point(data=d[d$WBC.result=="attenuated",], colour='red')
 
plot2 <- ggplot(d, aes(log(Hazard.Ratio), log(GrimAge.sensitivity.Hazard.Ratio)), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("Fully-adjusted log(Hazard Ratio)") +
      ylab("Full + WBC + GrimAge-adjusted log(Hazard Ratio)") +
      theme(plot.title = element_text(size = 11)) +
  geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=d[d$GrimAge.result!="attenuated",], shape=21, fill="white", size=2)+
  geom_point(data=d[d$GrimAge.result=="attenuated",], colour='red')

multiplot(plot1, plot2, cols=2)
