library(argparse)
library(data.table)
library(RColorBrewer)
library(ggplot2)

hztl_facet_group_bar = function(dat,fill,ycat,yname,xcounts,xname,barlabel,groupby,facetby,palette,outpath){
  p=ggplot(data = dat, aes_string(fill = fill, x = ycat, y = xcounts, group=groupby), width = 0.6) +
    geom_bar(colour='black', stat = 'identity', position='dodge', width=0.7) +
    geom_text(aes_string(label=barlabel, group=groupby), color = 'black', position = position_dodge(width= 0.7), hjust = -0.2, size = 4) +
    coord_flip() +
    scale_fill_brewer('div', palette=palette) +
    xlab(yname) +
    scale_y_continuous(name = xname, expand = expansion(mult = c(0, 0.3))) +
    facet_wrap(as.formula(paste('~', facetby))) +
    #	geom_text(data = label.df, label = '*') +
    theme_bw() +
    theme(
      strip.text = element_text(size=16),
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
      legend.position = 'right',
      legend.text = element_text(size = 12),
      legend.title=element_blank(),
      plot.margin=unit(c(0.25,1,0.25,1), 'cm')
    )
  plot(p)
  aspect_ratio = 1.5
  height = 5
  ggsave(outpath,p, height = 7 , width = 7 * aspect_ratio)
}


parser = ArgumentParser(description='Plot map overlaps detected by find_QTL_overlaps.sh')
parser$add_argument('--overlapsfile', metavar='overlapsfile',help='Output from report_unique.py in find_QTL_overlaps.sh pipeline')
parser$add_argument('--outpathbase', metavar='outpathbase',help='Directory to write plots to')

args=parser$parse_args(c('--overlapsfile','~/Documents/sqtl/pfcIntersect/overlap_stats.tsv','--outpath','~/Documents/sqtl/pfcIntersect/')) 
#args = parser$parse_args()

overlapsfile=args$overlapsfile
outpathbase=args$outpathbase

res = fread(overlapsfile)

# Split annotation string
res[, c('tissue','known','const','dataset','datatype') := tstrsplit(res$comparison,'_',fixed=T)]

# Proportions -> percentages (don't need to this now--report_unique.py handles it)
#perccols = colnames(res)[grepl('perc',colnames(res))]
#res[,:= res[,..perccols]*100]
# Add label columns (percentage of overlap)
makePercentStr = function(x){
  return(paste0(sprintf(x, fmt = '%#.1f'),'%'))
}

res[,'perc_transcripts_str'] = makePercentStr(res$perc_transcripts)
res[,'perc_genes_str'] = makePercentStr(res$perc_transcripts)
res[,'perc_snps_str'] = makePercentStr(res$perc_transcripts)

##### Plots by data type #####
# GWAS plots
gwas_palette='Purples'
hztl_facet_group_bar(res[res$datatype=='GWAS'],'known','dataset','GWAS dataset','perc_transcripts','Transcripts (%)','perc_transcripts_str','known','tissue',gwas_palette,paste0(outpathbase,'GWAS_overlaps-transcripts.png'))
hztl_facet_group_bar(res[res$datatype=='GWAS'],'known','dataset','GWAS dataset','perc_genes','Genes (%)','perc_genes_str','known','tissue',gwas_palette,paste0(outpathbase,'GWAS_overlaps-genes.png'))
hztl_facet_group_bar(res[res$datatype=='GWAS'],'known','dataset','GWAS dataset','perc_snps','GWAS SNPs (%)','perc_snps_str','known','tissue',gwas_palette,paste0(outpathbase,'GWAS_overlaps-snps.png'))
# eQTL plots
eqtl_palette='Reds'
hztl_facet_group_bar(res[res$datatype=='eQTLs'],'known','dataset','eQTL dataset','perc_transcripts','Transcripts (%)','perc_transcripts_str','known','tissue',eqtl_palette,paste0(outpathbase,'eQTL_overlaps-transcripts.png'))
hztl_facet_group_bar(res[res$datatype=='eQTLs'],'known','dataset','eQTL dataset','perc_genes','Genes (%)','perc_genes_str','known','tissue',eqtl_palette,paste0(outpathbase,'eQTL_overlaps-genes.png'))
hztl_facet_group_bar(res[res$datatype=='eQTLs'],'known','dataset','eQTL dataset','perc_snps','eQTL SNPs (%)','perc_snps_str','known','tissue',eqtl_palette,paste0(outpathbase,'eQTL_overlaps-snps.png'))
# sQTL plots
sqtl_palette='Greens'
hztl_facet_group_bar(res[res$datatype=='sQTLs'],'known','dataset','sQTL dataset','perc_transcripts','Transcripts (%)','perc_transcripts_str','known','tissue',sqtl_palette,paste0(outpathbase,'sQTL_overlaps-transcripts.png'))
hztl_facet_group_bar(res[res$datatype=='sQTLs'],'known','dataset','sQTL dataset','perc_genes','Genes (%)','perc_genes_str','known','tissue',sqtl_palette,paste0(outpathbase,'sQTL_overlaps-genes.png'))
hztl_facet_group_bar(res[res$datatype=='sQTLs'],'known','dataset','sQTL dataset','perc_snps','sQTL SNPs (%)','perc_snps_str','known','tissue',sqtl_palette,paste0(outpathbase,'sQTL_overlaps-snps.png'))
# isoQTL plots
isoqtl_palette='Oranges'
hztl_facet_group_bar(res[res$datatype=='isoQTLs'],'known','dataset','isoQTL dataset','perc_transcripts','Transcripts (%)','perc_transcripts_str','known','tissue',isoqtl_palette,paste0(outpathbase,'isoQTL_overlaps-transcripts.png'))
hztl_facet_group_bar(res[res$datatype=='isoQTLs'],'known','dataset','isoQTL dataset','perc_genes','Genes (%)','perc_genes_str','known','tissue',isoqtl_palette,paste0(outpathbase,'isoQTL_overlaps-genes.png'))
hztl_facet_group_bar(res[res$datatype=='isoQTLs'],'known','dataset','isoQTL dataset','perc_snps','isoQTL SNPs (%)','perc_snps_str','known','tissue',isoqtl_palette,paste0(outpathbase,'isoQTL_overlaps-snps.png'))
# tQTL plots
tqtl_palette='YlOrBr'
hztl_facet_group_bar(res[res$datatype=='tQTLs'],'known','dataset','tQTL dataset','perc_transcripts','Transcripts (%)','perc_transcripts_str','known','tissue',tqtl_palette,paste0(outpathbase,'tQTL_overlaps-transcripts.png'))
hztl_facet_group_bar(res[res$datatype=='tQTLs'],'known','dataset','tQTL dataset','perc_genes','Genes (%)','perc_genes_str','known','tissue',tqtl_palette,paste0(outpathbase,'tQTL_overlaps-genes.png'))
hztl_facet_group_bar(res[res$datatype=='tQTLs'],'known','dataset','tQTL dataset','perc_snps','tQTL SNPs (%)','perc_snps_str','known','tissue',tqtl_palette,paste0(outpathbase,'tQTL_overlaps-snps.png'))


##### Plots source #####
# Do additional data wrangling first

# GWAS plots
gwas_palette='Purples'
hztl_facet_group_bar(res[res$dataset=='PEC2018'],'known','datatype','QTL SNP source','perc_transcripts','Transcripts (%)','perc_transcripts_str','known','tissue',gwas_palette,paste0(outpathbase,'PEC_overlaps-transcripts.png'))
hztl_facet_group_bar(res[res$dataset=='PEC2018'],'known','datatype','QTL SNP source','perc_snps','PEC SNPs (%)','perc_snps_str','known','tissue',gwas_palette,paste0(outpathbase,'PEC_overlaps-snps.png'))
# eQTL plots
eqtl_palette='Reds'
hztl_facet_group_bar(res[res$datatype=='eQTLs'],'known','dataset','eQTL dataset','perc_transcripts','Transcripts (%)','perc_transcripts_str','known','tissue',eqtl_palette,paste0(outpathbase,'eQTL_overlaps-transcripts.png'))
hztl_facet_group_bar(res[res$datatype=='eQTLs'],'known','dataset','eQTL dataset','perc_genes','Genes (%)','perc_genes_str','known','tissue',eqtl_palette,paste0(outpathbase,'eQTL_overlaps-genes.png'))
hztl_facet_group_bar(res[res$datatype=='eQTLs'],'known','dataset','eQTL dataset','perc_snps','eQTL SNPs (%)','perc_snps_str','known','tissue',eqtl_palette,paste0(outpathbase,'eQTL_overlaps-snps.png'))
# sQTL plots
