library(argparse)
library(data.table)
library(ggplot2)

parser = ArgumentParser(description='Read the HCP_titration_summary.tsv file created by make_stats_table.py & plot yields from the HCP titrations.')

parser$add_argument('--basepath', metavar='basepath', type='character', nargs='?', help='Dir containing the HCP_titration_summary.tsv file')

args = parser$parse_args()
#args = parser$parse_args('/sc/arion/projects/EPIASD/splicingQTL/output/fqtl_output_wasp_nominal_normal/minCovars+seqPC9/4genoPCs/')

basepath=args$basepath

res = fread(paste0(basepath,'/HCP_titration_summary.tsv'))

gg_linechart = function(df,xcol,xname,ycol,yname,title){
 max_hc = max(df[,..xcol])
 p = ggplot(df, aes_string(x = xcol, y = ycol)) +
    geom_point() + 
    geom_line(linetype='dotted',color='red') + 
    ggtitle(title) +
    scale_x_continuous(name = xname, limits=c(0,max_hc+5), breaks=seq(0,max_hc+5,by=5), expand = c(0,.75)) +
    scale_y_continuous(name = yname) +
#    scale_x_continuous(limits = ) + 
    theme_bw() + theme(
      plot.title = element_text(size=16, face='bold'),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size=14,color='black'),
      axis.title.x = element_text(size=16,color='black'),
      axis.text.y = element_text(size=14,color='black'),
      axis.title.y = element_text(size=16,color='black')
    )
 plot(p)
 ggsave(paste0(basepath,'/sQTL_titration-',xcol,'-vs-',ycol,'.png'),p)
 ggsave(paste0(basepath,'/sQTL_titration-',xcol,'-vs-',ycol,'.pdf'),p)
}

gg_linechart(res,'n_hcps','Hidden covariates','n_sgenes','Genes containing a QTL','Genes containing a QTL')
gg_linechart(res,'n_hcps','Hidden covariates','n_phenogenes','sGenes','sGenes')
gg_linechart(res,'n_hcps','Hidden covariates','n_phenogenes_jmb','sGenes','sGenes (bedtools)')
gg_linechart(res,'n_hcps','Hidden covariates','n_qtl_introns','sQTL-intron pairs','sQTL-intron pairs')
gg_linechart(res,'n_hcps','Hidden covariates','n_snps','sQTLs','sQTLs')
gg_linechart(res,'n_hcps','Hidden covariates','n_introns','Introns','Introns')
