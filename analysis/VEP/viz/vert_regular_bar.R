vert_regular_bar = function(dat,fill,filltitle,ycat,yname,yaxlim,xcounts,xname,barlabel,errbarlow,errbarhigh,outpath){
  p=ggplot(data = dat, aes_string(fill = fill, x = xcounts, y = ycat), width = 0.6) +
    geom_bar(colour='black', stat = 'identity', position='dodge', width=0.7) +
    geom_text(aes_string(label=barlabel), color = 'red', nudge_y=0.075, size = 5) +
    geom_hline(yintercept=1, linetype = 'dashed') +
    geom_errorbar(aes_string(ymax = errbarhigh, ymin = errbarlow), position = position_dodge(width=0.5), width = 0.2) +
    #    scale_fill_gradient(low = muted('red'),mid = 'white',high = muted('blue'),midpoint = 0) +
    scale_fill_gradient(name = filltitle, low = 'darkslategrey',high = 'darkslategray1') +
    xlab(xname) +
    scale_y_continuous(name = yname, limits=c(0,yaxlim), breaks=seq(0,yaxlim,by=0.5), expand = c(0,.1)) +
    #	geom_text(data = label.df, label = '*') +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(color='black', size=0.6),
      axis.line.y = element_line(color='black', size=0.6),
      axis.text.x = element_text(colour='black', size=11, angle=90, hjust=1),
      axis.text.y = element_text(colour='black', size=11),
      axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=12),
      axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=12),
      axis.ticks = element_blank(),
      legend.position = 'right',
      legend.text = element_text(size = 12),
      plot.margin=unit(c(0.25,1,0.25,1), 'cm')
    )
  plot(p)
  aspect_ratio = 1.5
  height = 5
  ggsave(outpath,p, height = 7 , width = 7 * aspect_ratio)
}
