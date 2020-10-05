ggbiplot = function(eigenvec,x,y,savepath){
  pcx = paste0('PC',as.character(x))
  pcy = paste0('PC',as.character(y))
  g = ggplot(eigenvec,aes_string(x=pcx,y=pcy,color='ancestry',fill='ancestry')) +
  geom_point(alpha = 0.25, shape = 16, size=2.5) +
  xlab(pcx) +
  ylab(pcy) +
  scale_color_manual(values=colscale,name='Ancestry') +
  scale_fill_manual(values=colscale,name='Ancestry') +
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
    panel.border = element_rect(colour = 'black', fill=NA, size=2)
  )
  ggsave(paste0(savepath,'ancestryPCA_',pcx,'vs',pcy,'.png'), width=pwidth, plot = g)
  return(g)
} 

est_ancestry = function(eigenvec,g,x,y,conflevel,savepath){
  pcx = paste0('PC',as.character(x))
  pcy = paste0('PC',as.character(y))

  for (z in unique(eigenvec['ancestry'][eigenvec['ancestry']!='This study'])){
     dat = eigenvec[which(eigenvec['ancestry']==z),]
     # Add ellipse to plot to capture region containing x proportion of data points 
     g = g + stat_ellipse(data=dat,method='t',alpha=0.55,level=conflevel,na.rm=T)
     # Extract plot components
     build = ggplot_build(g)$data
     points = build[[1]]
     ell = build[[length(build)]]
     # Find which points are inside the ellipse, and add this to the data
     eigenvec[,paste0('est_',z)]=as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y))
   } 
  #) 
  ggsave(paste0(savepath,'ancestryPCAwithEllipse_',pcx,'vs',pcy,'.png'), width=pwidth, plot = g)
  return(eigenvec) 
}

# Make loadings heatmap (reqs library ComplexHeatmap)
plot_loadings = function(eigenvec,eigenval,outpath,basefile){
 pc_cols = grep('^PC[.]*', colnames(eigenvec), value = T, perl=T)
 product = sweep(eigenvec, 2, sqrt(eigenval), FUN='*')
 png(paste0(outpath,basefile,'.pca_loadings.png'))
 Heatmap(eigenvec[,pc_cols], name = 'PCA loadings 1:20',
        show_row_names = FALSE, show_column_names = TRUE, 
        row_dend_reorder = FALSE, column_dend_reorder = FALSE, 
        cluster_rows = FALSE, cluster_columns = FALSE, column_order = paste0('PC', c(1:20)),
        use_raster = TRUE, raster_device = 'png', raster_quality = 16)
 dev.off()
}

# Biplots of genotype PCA results
library(reshape2)

outpath='/sc/arion/projects/EPIASD/splicingQTL/PCA/'
setwd(outpath)
options(scipen=100, digits=3)

basepath='/sc/arion/scratch/belmoj01/splicingQTL/' # This is the path where the PCA output lives
basefile='Capstone4.sel.idsync.2allele.1kg_phase3' # This is the base name of the .eigenval & .eigenvec files
meta1kgfile='/sc/arion/projects/EPIASD/splicingQTL/PCA/1kg_phase3_samplesuperpopinferreddata.txt'
#estimate_ancestry=TRUE # Should an ancestry estimation be returned? This will also add ancestry ellipses to the plots written by ggbiplot()

conflevel=0.95 # Confidence level to use when drawing ellipse to estimate sample ancestry

# Plot params
pwidth = 13

# read in the eigenvectors, produced in PLINK
eigenvec = data.frame(read.table(paste0(basepath,basefile,'.eigenvec'), header=FALSE, skip=0, sep=' '))
eigenval = t(data.frame(read.table(paste0(basepath,basefile,'.eigenval'), header=FALSE, skip=0, sep=' ')))
rownames(eigenvec) = eigenvec[,2]
eigenvec = eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) = colnames(eigenval) = paste0('PC', c(1:20))

# read in the 1kg metadata & read population from 3rd column (col name varies between metadata files I prepared)
ped = data.frame(read.table(meta1kgfile, header=TRUE, skip=0, sep='\t'))
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
colscale = brewer.pal(n = length(levels(eigenvec[,'ancestry'])), 'Dark2')
colscale[which(levels(eigenvec[,'ancestry'])=='This study')] = 'black'

par(mfrow=c(1,2))
plot(eigenvec[,1], eigenvec[,2], xlab='PC1', ylab='PC2')
plot(eigenvec[,1], eigenvec[,3], xlab='PC1', ylab='PC3')

legend('topleft', bty='n', cex=1.5, title='', c('African','Hispanic','East Asian','Caucasian','South Asian'), fill=c('yellow','forestgreen','grey','royalblue','black'))

# ggplot it
library(ggplot2)
library(sp) # Used to extract plot parameters when reading which points are inside ellipse

npcs = 4 # n PCs to visualize. All pairs 1:n will be plotted
PCs = seq(1,4,1)

pc_pairs = lapply(PCs,function(x){
  lapply(PCs, function(y,x){
    if(x!=y){
      g=ggbiplot(eigenvec,x,y,outpath)
      if(x==2 & y==1){
        annot_eigenvec=est_ancestry(eigenvec,g,x,y,conflevel,outpath)
        write.table(annot_eigenvec,paste0(outpath,basefile,'.eigenvec_ancestry_est'),sep='\t',quote=F)
        return(annot_eigenvec)
      }
    }
  },x)
})

# Wrangle ancestry-annotated eigenvector to produce a single column listing all the the 1kg superpopulations 
# NOTE: solution below does not preserve multiple assigned ancestries if there is an overlap between 2 ellipses--instead only the first is retained in the results. e.g. 'CMC_PITT_178'
library(dplyr)

annot_eigenvec = pc_pairs[[2]][[1]] 
est_cols = grep('^est_[.]*', colnames(annot_eigenvec), value = T, perl=T)
anc = do.call(rbind,
    lapply(est_cols,function(x){
    annot_x = annot_eigenvec[which(annot_eigenvec[,x]),]
    tryCatch(
      {
        annot_x['ethnicity_est'] = gsub('est_','',x,perl=T)
        return(annot_x)
      },error = function(e) {
        return(annot_eigenvec[FALSE,]) 
      }
    )
  })
)

anc$ID = row.names(anc)
anc_est = anc %>% group_by(ID) %>% summarize(ancestry_est = paste0(sort(unique(ethnicity_est)),collapse=','))

annot_eigenvec$ID = row.names(annot_eigenvec)
ancestry_estimated = filter(annot_eigenvec,ancestry=='This study') %>% left_join(anc_est,by='ID')

write.table(ancestry_estimated,paste0(outpath,basefile,'.ancestry_est'),sep='\t',quote=F)
