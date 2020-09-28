library(ggplot2)
library(RColorBrewer)

hztl_group_bar = function(dat,fill,ycat,yname,xcounts,xname,palette,outpath){
  p=ggplot(data = dat, aes_string(fill = fill, x = ycat, y = xcounts), width = 0.6) +
    geom_bar(colour='black', stat = 'identity', position="dodge", width=0.7) +
    coord_flip() +
    scale_fill_brewer('div', palette=palette) +
    xlab(yname) +
    scale_y_continuous(name = xname,expand = c(0,.1)) +
    #	geom_text(data = label.df, label = '*') +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(color='black', size=0.6),
      axis.line.y = element_line(color='black', size=0.6),
      axis.text.x = element_text(colour='black', size=11),
      axis.text.y = element_text(colour='black', size=11),
      axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=12),
      axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=12),
      axis.ticks = element_blank(),
      legend.position = c(.75,.875),
      legend.text = element_text(size = 12),
      legend.title=element_blank(),
      plot.margin=unit(c(0.25,1,0.25,1), 'cm')
    )
  plot(p)
  aspect_ratio = 1.5
  height = 5
  ggsave(outpath,p, height = 7 , width = 7 * aspect_ratio)
}
