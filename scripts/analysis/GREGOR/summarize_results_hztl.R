library(qvalue)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

resdir=args[1] # Dir containing the GREGOR results
resfile='StatisticSummaryFile.txt'
plottitle=args[2] # Title of plot
pngout=args[3] # Name of outfile

# Optional: path to a list of IDs for features to include in the plot (default is all)
#includesfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/eQTLs/TFsTested.txt'
#includesfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/sQTLs/RBPsTested.txt'
#includesfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Raj_NatComm_alzh2018/RBPsTested.txt'
#####

df = read.table(paste0(resdir,resfile),stringsAsFactors=FALSE,header=TRUE)

# Make label for each tested feature
df$feature = str_match(df$Bed_File,'.*-(.+).bed')[,2]

# If an includes file was passed, read the list of features to include & filter the data frame
if(exists("includesfile")){
  includes = read.table(includesfile,stringsAsFactors=FALSE,header=FALSE)
  df = df[which(df$feature %in% includes[,1]),]
}

# Calculate fold-change over expected # of SNPs
df$FC = df$InBed_Index_SNP/df$ExpectNum_of_InBed_SNP 

# FDR correction (pi0=1 for Benjamini-Hochberg FDR correction)
df$qval = qvalue(df$PValue, pi0=1)$qvalues

#df$label = round(df$FC,digits=2)
df$barlabel = '*'
df$barlabel[which(df$qval>=.05)] = ''

df$enrichment = 'QTL enrichment'

# Make label for each tested feature
df$feature = str_match(df$Bed_File,'.*-(.+).bed')[,2]

# Drop rows where FC is null
df = df[!is.na(df$FC),]

# Reorder by ascending FC as Walker 2018 did 
df = df[order(df$FC,decreasing=TRUE),] # This only orders the data frame; ggplot cares about the factor levels
df$feature = factor(df$feature, levels = df$feature[order(df$FC,decreasing=TRUE)])

# Write out results tables
write.table(df,paste0(resdir,resfile,'-withqvals'),sep='\t',row.names=F,col.names=T,quote=F)
write.table(df[which(df$qval<.05),],paste0(resdir,resfile,'-qval05'),sep='\t',row.names=F,col.names=T,quote=F)

# Plot results
library(ggplot2)

p = ggplot(df, aes(x=feature,y=FC,enrichment, fill=-log10(qval))) +
  geom_bar(stat='identity') +
#  coord_flip() + 
  ggtitle(plottitle) +
  geom_text(aes(label=barlabel), color = 'red', nudge_y=0.05, size = 5) +
  scale_y_continuous(limits = c(0,3.5)) + 
  scale_fill_gradient(low = 'mediumpurple4', high = 'mediumpurple1', limits=c(0,10)) +
#  scale_fill_gradient(low = 'darkred', high = 'red', limits=c(0,100)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
#    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color='black', size=10),
    axis.line.x = element_line(),
    axis.text.x = element_text(angle = 90, margin=margin(0,22,0,0,'pt'), color='black', size=10, hjust=1),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
#    legend.text.align = 0.5,
#    legend.title = element_text(hjust = -2.5, vjust = 1.6),
    legend.position = 'right',
    legend.direction = 'vertical'
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 4,
    title.position = 'top', title.hjust = 0.5))

plot(p)

aspect_ratio = 3 
height = 3
ggsave(paste0(resdir,pngout,'.png'), plot = p, dpi = 300, height=height, width=height*aspect_ratio)
ggsave(paste0(resdir,pngout,'.pdf'), plot = p, dpi = 300, height=height, width=height*aspect_ratio)
