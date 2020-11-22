library(ggplot2)

plot_dist_vs_y = function(df,x,xlab,xmin,xmax,y,ylab,ymin,ymax,color,colorlab,colormin,colormax,savepath){
  df[,x] = dat[,x]/1000
  g = ggplot(df,aes_string(x=x, y=y, color=color)) +
  geom_point(alpha = 0.33, shape = 16, size=2.5) +
  geom_vline(xintercept=0, linetype='dashed') +
  xlim(c(xmin,xmax)) + 
  xlab(xlab) +
  ylim(c(ymin,ymax)) + 
  ylab(ylab) +
  scale_colour_gradient(name = colorlab, low=colormin, high=colormax, guide='colourbar') +
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
    legend.title = element_text(color = 'black', size = 12),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
  ggsave(paste(savepath,'/distance_',x,'vs',y,'.png',sep=''), width=pwidth, plot = g)
}

#infile = '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl+ensg+nearest_TSS.txt'
infile = '/sc/hydra/projects/pintod02c/datasets-external/PEC/DER-10b_hg19_isoQTL.FPKM5.all.txt'
dat = read.table(infile, sep='\t', header=T)
pwidth = 14
#x = 'dist'
#x = 'nearest_TSS_dist'
x='SNP_distance_to_TSS'
#xlab = 'Distance to AS event (kb)'
xlab = 'Distance to nearest TSS (kb)'
xlow = -1000
xhigh = 1000 
#y = 'slope'
y = 'regression_slope'
ylab = 'Effect Size'
ylow = -2
yhigh = 2
#dat$qval_log10 = -log10(dat$qval)
dat$qval_log10 = -log10(dat$FDR)
colorcol = 'qval_log10'
colorlab = '-log10(FDR)' 
#colorlow = 'skyblue'
colorlow = 'violetred'
#colorhigh = 'slateblue4'
colorhigh = 'red'
plot_dist_vs_y(dat, x, xlab, xlow, xhigh, y, ylab, ylow, yhigh, colorcol, colorlab, colorlow, colorhigh, dirname(infile)) 
