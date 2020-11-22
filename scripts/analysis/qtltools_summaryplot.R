library(ggplot2)

res = list('introns_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19intronsOnly_qtltools_enrich.txt',
  'exons_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19exonsOnly_qtltools_enrich.txt',
  'coding_exons_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19codingexonsOnly_qtltools_enrich.txt',
  '5000bp_downstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV195000downstreamOnly_qtltools_enrich.txt',
  '2000bp_downstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV192000downstreamOnly_qtltools_enrich.txt',
  '1000bp_downstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV191000downstreamOnly_qtltools_enrich.txt',
  '750bp_downstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19750downstreamOnly_qtltools_enrich.txt',
  '500bp_downstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19500downstreamOnly_qtltools_enrich.txt',
  '250bp_downstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19250downstreamOnly_qtltools_enrich.txt',
  '5000bp_upstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV195000upstreamOnly_qtltools_enrich.txt',
  '2000bp_upstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV192000upstreamOnly_qtltools_enrich.txt',
  '1000bp_upstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV191000upstreamOnly_qtltools_enrich.txt',
  '750bp_upstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19750upstreamOnly_qtltools_enrich.txt',
  '500bp_upstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19500upstreamOnly_qtltools_enrich.txt',
  '250bp_upstream_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19250upstreamOnly_qtltools_enrich.txt',
  'genes_only'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV19genesOnly_qtltools_enrich.txt',
  '3UTR'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV193UTRexonsOnly_qtltools_enrich.txt',
  '5UTR'='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/qtlTools_results/gencodeV195UTRexonsOnly_qtltools_enrich.txt')

stats = lapply(names(res),function(x){
  infile = res[[x]]
  print(infile)
  qtlres = read.table(infile, head=FALSE, stringsAsFactors=FALSE)
  ftestres = fisher.test(matrix(c(qtlres$V1, qtlres$V2, round(qtlres$V3), qtlres$V2), ncol=2))
  return(data.frame(region=x,pval=ftestres$p.value, OR=ftestres$estimate))
})

df = do.call(rbind,stats)

df$label = round(df$OR,digits=2)
df$label[which(df$pval>=.05)] = ''

df$enrichment = 'sQTL enrichment'

p = ggplot(df, aes(region,enrichment, fill=-log10(pval))) + 
  ggtitle('sQTL enrichments by genomic region (GENCODE v19)') +
  geom_tile() +
  geom_text(aes(label=label), color = 'black', size = 4) + 
  #scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 1, limits=c(0,20)) +
  scale_fill_gradient(low = 'white', high = 'red', limits=c(0,20)) +
  theme(
    plot.title = element_text(hjust = 0.5),
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
    legend.position = c(0.02, -0.25),
    legend.direction = 'vertical'
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 6,
    title.position = 'top', title.hjust = 0.5))

plot(p)

aspect_ratio = 1.5
height = 7
ggsave('qtltools_res_summary.png', plot = p, dpi = 300, height=height, width=height*aspect_ratio)
