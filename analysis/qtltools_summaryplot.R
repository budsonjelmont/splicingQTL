library(ggplot2)

df = data.frame(region=c('introns_only','exons_only','coding_exons_only','500bp_downstream_only','500bp_upstream_only','genes_only','3UTR','5UTR'),OR=c(0.9592288,0.8888918,0.5867461,1.05753,1.267125,0.96544,0.963563,1.1042),pval=c(0.009148,0.0008403,2.2e-16,0.3359,3.174e-05,0.02153,0.4692,0.1181),stringsAsFactors=FALSE)

df$label = round(df$pval,digits=2)
df$label[which(df$label==0)] = '<.001'

df$enrichment = 'sQTL enrichment'

p = ggplot(df, aes(region,enrichment, fill=OR)) + 
  ggtitle('sQTL enrichments by genomic region (GENCODE v19)') +
  geom_tile() +
  geom_text(aes(label=label), color = 'black', size = 6) + 
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 1, limits=c(0,2)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, margin=margin(0,22,0,0,'pt'), color='black', size=10, hjust=1),
    axis.text.y = element_text(color='black', size=10),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
#    legend.text.align = 0.5,
#    legend.title = element_text(hjust = -2.5, vjust = 1.6),
    legend.position = c(0.04, -0.25),
    legend.direction = 'horizontal'
  ) +
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1,
    title.position = 'top', title.hjust = 0.5))

plot(p)

ggsave('qtltools_res_summary.png', plot = p, dpi = 300)
