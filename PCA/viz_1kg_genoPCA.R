multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
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
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

ggbiplot = function(eigenvec,x,y,savepath){
  pcx = paste0('PC',as.character(x))
  pcy = paste0('PC',as.character(y))
  g = ggplot(eigenvec,aes_string(x=pcx,y=pcy,color='ancestry')) +
  geom_point(alpha = 0.75, shape = 16, size=2.5) +
  xlab(pcx) +
  ylab(pcy) +
  scale_color_manual(values=colscale,name='Ancestry') +
  theme(
    title = element_text(color='black',size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    axis.line.x = element_line(color='black', size = 0.4),
    axis.line.y = element_line(color='black', size = 0.4),
    axis.title.x = element_text(color='black', size = 14),
    axis.title.y = element_text(color='black', size = 14),
    axis.text = element_text(color = 'black', size = 11),
    plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
  ggsave(paste(savepath,'ancestryPCA_',pcx,'vs',pcy,'.png',sep=''), width=pwidth, plot = g)
} 

# Biplots of genotype PCA results
library(reshape2)

outpath='/sc/arion/projects/EPIASD/splicingQTL/PCA/'
setwd(outpath)
options(scipen=100, digits=3)

basepath='/sc/arion/projects/EPIASD/splicingQTL/PCA/studydata/' # This is the path where the PCA output lives
basefile='ASD-EPI_Plates1-2.1kg_phase3' # This is the base name of the .eigenval & .eigenvec files
meta1kgfile='/sc/arion/projects/EPIASD/splicingQTL/PCA/1kg_phase3_samplesuperpopinferreddata.txt'

# Plot params
pwidth = 13

# read in the eigenvectors, produced in PLINK
eigenvec = data.frame(read.table(paste0(basepath,basefile,'.eigenvec'), header=FALSE, skip=0, sep=" "))
rownames(eigenvec) = eigenvec[,2]
eigenvec = eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) = paste('PC', c(1:20), sep="")

# read in the 1kg metadata & read population from 3rd column (col name varies between metadata files I prepared)
ped = data.frame(read.table(meta1kgfile, header=TRUE, skip=0, sep="\t"))
colnames(ped)[3] = 'ancestry'
ped = ped[which(ped$IID %in% rownames(eigenvec)), ]
rownames(ped) = ped$IID

#ped = ped[match(rownames(eigenvec), ped$IID),]
#all(ped$IID == rownames(eigenvec)) == TRUE

summary(eigenvec)

# Read in metadata
# Add ancestry column to eigenvalues dataframe
eigenvec = cbind(eigenvec, ped[, 'ancestry'][match(rownames(eigenvec), rownames(ped))])
colnames(eigenvec)[21] = 'ancestry'

# Our samples do not match IDs in the 1kg ped file, so where Ancestry is NA, set Ancestry to 'This study'
levels(eigenvec[,'ancestry']) = c(levels(eigenvec[,'ancestry']),'This study')
eigenvec[is.na(eigenvec['ancestry']),'ancestry'] = 'This study'

# Prepare manual color scale
library(RColorBrewer)
colscale = brewer.pal(n = length(levels(eigenvec[,'ancestry'])), "Dark2")
colscale[which(levels(eigenvec[,'ancestry'])=='This study')] = 'black'

#Determine the proportion of variance of each component
#proportionvariances = ((apply(eigenvec, 1, sd)^2) / (sum(apply(eigenvec, 1, sd)^2)))*100

par(mfrow=c(1,2))
plot(eigenvec[,1], eigenvec[,2], xlab='PC1', ylab='PC2')
plot(eigenvec[,1], eigenvec[,3], xlab='PC1', ylab='PC3')

legend("topleft", bty="n", cex=1.5, title="", c("African","Hispanic","East Asian","Caucasian","South Asian"), fill=c("yellow","forestgreen","grey","royalblue","black"))

# ggplot it
library(ggplot2)

npcs = 4 # n PCs to visualize. All pairs 1:n will be plotted
PCs = seq(1,4,1)

lapply(PCs,function(x){
  lapply(PCs, function(y,x){
    if(x!=y){
      ggbiplot(eigenvec,x,y,outpath)
    }
  },x)
})
