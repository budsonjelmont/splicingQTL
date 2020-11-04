library(data.table)
library(ggplot2)

basepath='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/4genoPCs_nogenoInHCP/'
#basepath='/Users/judsonbelmont/Documents/sqtl/post_fastQTL_wrapup/'
res = fread(paste0(basepath,'HCP_titration_summary.tsv'))

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
      title = element_text(size=16, face='bold'),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size=13,color='black'),
      axis.title.x = element_text(size=15,color='black',face='bold'),
      axis.text.y = element_text(size=13,color='black'),
      axis.title.y = element_text(size=15,color='black',face='bold')
    )
 plot(p)
 ggsave(paste0(basepath,'sQTL_titration-',xcol,'-vs-',ycol,'.png'),p)
}

gg_linechart(res,'n_hcps','Hidden covariates','n_sgenes','sGenes','sGenes')
gg_linechart(res,'n_hcps','Hidden covariates','n_qtls','sQTLs','sQTLs')
gg_linechart(res,'n_hcps','Hidden covariates','n_snps','SNPs','SNPs')
gg_linechart(res,'n_hcps','Hidden covariates','n_introns','Introns','Introns')
