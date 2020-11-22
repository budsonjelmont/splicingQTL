#pie_chart = function(dat,fill,filltitle,ycat,yname,xcounts,xname,barlabel,outpath){
#  p = ggplot(data = dat, aes_string(fill = fill, x = xcounts, y = ycat), width = 1) +
#    geom_bar(colour='black', stat = 'identity', position='fill', width=1)
#  p = p + coord_polar(theta='x',start=0) +
#    xlab('') +
#    ylab('') +
#    theme(
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      panel.background = element_blank(),
#      axis.line.x = element_blank(),
#      axis.line.y = element_blank(),
#      axis.text.x = element_blank(),
#      axis.text.y = element_blank(),
#      axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=12),
#      axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=12),
#      axis.ticks = element_blank(),
#      legend.position = 'right',
#      legend.text = element_text(size = 12),
#      plot.margin=unit(c(0.25,1,0.25,1), 'cm')
#    )
#  plot(p)
#  aspect_ratio = 1.5
#  height = 5
#  ggsave(outpath,p, height = 7 , width = 7 * aspect_ratio)
#}

pie_chart = function(dat,palette,outpath){
  p = ggplot(data = dat, aes(fill = annotation, x = '', y = n_snps), width = 1) +
    geom_bar(colour='black', stat = 'identity', position='fill', width=1)
  p = p + coord_polar(theta='y',start=0) +
    xlab('') +
    ylab('') +
#    scale_fill_manual(values=colorRampPalette(brewer.pal(8, palette))(nrow(dat))) +
    scale_fill_manual(values=colorRampPalette(palette)(nrow(dat))) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'right',
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 17),
      plot.margin=unit(c(0.25,1,0.25,1), 'cm')
    )
  plot(p)
  aspect_ratio = 1.25
  height = 8
  ggsave(outpath,p, height = 5 , width = 8 * aspect_ratio)
}
