library(ggplot2)

suppressMessages(library(argparser))
p = arg_parser('Plot the distribution of QTL effect sizes as a function of SNP-phenotype distance & colored by statistical significance based on values contained in the qtls.txt output')
p = add_argument(p, 'infile', help='qtls.txt results file containing a distance & statistical measure for each QTL-intron pair')
p = add_argument(p, 'xcol', help='Column name in qtls.txt containing a distance measure for each QTL (e.g. dist). Divided by 1000 before plotting.')
p = add_argument(p, '--xlab', default='dist', help='Label for x axis')
p = add_argument(p, '--xlo', type='numeric', default=-100, help='Min value for x axis')
p = add_argument(p, '--xhi', ,type='numeric', default=100, help='Max value for x axis')
p = add_argument(p, 'ycol', help='Column name in qtls.txt containing an effect size measure for each QTL (e.g. slope)')
p = add_argument(p, '--ylab', default='slope', help='Label for y axis')
p = add_argument(p, '--ylo', type='numeric', default=-2, help='Min value for y axis')
p = add_argument(p, '--yhi', type='numeric', default=2, help='Max value for y axis')
p = add_argument(p, 'colorcol', help='Column name in qtls.txt containing a continuous measure to color each point by (e.g. qval)')
p = add_argument(p, '--logit', flag=T, help='Take -log10 of [colorcol] values before plotting')
p = add_argument(p, '--colorlab', default='qval', help='Label to display on gradient legend')
p = add_argument(p, '--colorlo', default='skyblue', help='Start of color gradient')
p = add_argument(p, '--colorhi', default='slateblue4', help='End of color gradient')

args = parse_args(p)

infile = args$infile
distcol = args$distcol
distcol = args$distlab

# Plotting params
pwidth = 14

# Plotting function
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
  ggsave(paste(savepath,'/distance_',x,'vs',y,'.pdf',sep=''), width=pwidth, plot = g)
}

dat = read.table(infile, sep='\t', header=T)
x = args$xcol 
xlab = args$distlab
xlow = args$xlo
xhigh = args$xhi 
y = args$ycol 
ylab = args$ylab 
ylow = args$ylo 
yhigh = args$yhi
colorcol = args$colorcol 
colorlab = args$colorlab 
colorlow = args$colorlo 
colorhigh = args$colorhi 

if(args$logit){
  paste0('-log10 transforming column ',colorcol,'...')
  dat[colorcol] = -log10(dat[colorcol])
}

plot_dist_vs_y(dat, x, xlab, xlow, xhigh, y, ylab, ylow, yhigh, colorcol, colorlab, colorlow, colorhigh, dirname(infile)) 
