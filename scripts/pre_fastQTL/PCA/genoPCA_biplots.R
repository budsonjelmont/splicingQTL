multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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

# Biplots of genotype PCA results
library(reshape2)

setwd('~/Documents/sqtl/PCA/')
options(scipen=100, digits=3)

#Read in the eigenvectors
eigen <- data.frame(read.table("Capstone4.sel.idsync.2allele.maf5.mind1.geno1.hwe1e-5.deduped.highLDexcl.indep1500_150_.2.eigenvec", header=TRUE, skip=0, sep=" ", stringsAsFactors = FALSE))
rownames(eigen) <- eigen[,2]
eigen <- eigen[,3:ncol(eigen)]

summary(eigen)

# Read in metadata
meta <- data.frame(read.table('../ID reformatting/meta_matchedIDs.csv', header=TRUE, comment.char='@', sep=",", stringsAsFactors = FALSE))
rownames(meta) = meta$PSIcountID
# Add ancestry column to eigenvalues dataframe
eigen = cbind(eigen, meta[, "ethnicity"][match(rownames(eigen), rownames(meta))])
colnames(eigen)[21] = 'ancestry'

#Determine the proportion of variance of each component
proportionvariances <- ((apply(eigen, 1, sd)^2) / (sum(apply(eigen, 1, sd)^2)))*100

#plot(eigen[,1], eigen[,2], xlab='PC1', ylab='PC2')

par(mfrow=c(1,2))
plot(eigen[,1], eigen[,2], xlab='PC1', ylab='PC2')
plot(eigen[,1], eigen[,3], xlab='PC1', ylab='PC3')


legend("topleft", bty="n", cex=1.5, title="", c("African","Hispanic","East Asian","Caucasian","South Asian"), fill=c("yellow","forestgreen","grey","royalblue","black"))

# ggplot it
library(ggplot2)

pc12 = ggplot(eigen,aes(x=PC1,y=PC2,color=ancestry)) +
  geom_point(alpha = 0.75, shape = 16, size=2.5) +
  xlab('PC1') +
  ylab('PC2') +
  scale_color_brewer(palette="Dark2",name="Ancestry",labels=c("AA","AS","CAUC","HISP","Other")) +
  theme(
    title = element_text(color='black',size=18, family='Arial'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    axis.line.x = element_line(color='black', size = 0.4),
    axis.line.y = element_line(color='black', size = 0.4),
    axis.title.x = element_text(color='black', size = 14, family='Arial'),
    axis.title.y = element_text(color='black', size = 14, family='Arial'),
    axis.text = element_text(color = 'black', size = 11, family='Arial'),
    legend.position='none',
    plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )

plot(pc12)
ggsave(paste('~/Documents/sqtl/PCA/ancestryPCA_PC1vs2.png',sep=''), plot = pc12)

pc23 = ggplot(eigen,aes(x=PC2,y=PC3,color=ancestry)) +
  geom_point(alpha = 0.75, shape = 16, size=2.5) +
  xlab('PC2') +
  ylab('PC3') +
  scale_color_brewer(palette="Dark2",name="Ancestry",labels=c("AA","AS","CAUC","HISP","Other")) +
  theme(
    title = element_text(color='black',size=18, family='Arial'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    axis.line.x = element_line(color='black', size = 0.4),
    axis.line.y = element_line(color='black', size = 0.4),
    axis.title.x = element_text(color='black', size = 14, family='Arial'),
    axis.title.y = element_text(color='black', size = 14, family='Arial'),
    axis.text = element_text(color = 'black', size = 11, family='Arial'),
    legend.position='none',
    plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )

plot(pc23)
ggsave(paste('~/Documents/sqtl/PCA/ancestryPCA_PC2vs3.png',sep=''), plot = pc23)

pc14 = ggplot(eigen,aes(x=PC1,y=PC4,color=ancestry)) +
  geom_point(alpha = 0.75, shape = 16, size=2.5) +
  xlab('PC1') +
  ylab('PC4') +
  scale_color_brewer(palette="Dark2",name="Ancestry",labels=c("AA","AS","CAUC","HISP","Other")) +
  theme(
    title = element_text(color='black',size=18, family='Arial'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    axis.line.x = element_line(color='black', size = 0.4),
    axis.line.y = element_line(color='black', size = 0.4),
    axis.title.x = element_text(color='black', size = 14, family='Arial'),
    axis.title.y = element_text(color='black', size = 14, family='Arial'),
    axis.text = element_text(color = 'black', size = 11, family='Arial'),
    legend.position='none',
    plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )

plot(pc14)
ggsave(paste('~/Documents/sqtl/PCA/ancestryPCA_PC1vs4.png',sep=''), plot = pc14)

pc13 = ggplot(eigen,aes(x=PC1,y=PC3,color=ancestry)) +
  geom_point(alpha = 0.75, shape = 16, size=2.5) +
  xlab('PC1') +
  ylab('PC3') +
  scale_color_brewer(palette="Dark2",name="Ancestry",labels=c("AA","AS","CAUC","HISP","Other")) +
  theme(
    title = element_text(color='black',size=18, family='Arial'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    axis.line.x = element_line(color='black', size = 0.4),
    axis.line.y = element_line(color='black', size = 0.4),
    axis.title.x = element_text(color='black', size = 14, family='Arial'),
    axis.title.y = element_text(color='black', size = 14, family='Arial'),
    axis.text = element_text(color = 'black', size = 11, family='Arial'),
    plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )

plot(pc13)
ggsave(paste('~/Documents/sqtl/PCA/ancestryPCA_PC1vs3.png',sep=''), plot = pc13)

p = multiplot(pc12, pc13, cols=2)

ggsave(paste('~/Documents/sqtl/PCA/ancestryPCA_multiplot.png',sep=''), plot = p)
