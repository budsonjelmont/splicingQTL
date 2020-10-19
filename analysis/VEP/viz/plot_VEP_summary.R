setwd('~/Documents/sqtl/analysis/VEP/')

library(data.table)   
library(qvalue)
library(RColorBrewer)
source('hztl_grouped_bar.R')
source('vert_regular_bar.R')
source('abovebelow0_bar.R')

dat = fread('~/Documents/sqtl/analysis/VEP/QTLSNPs.VEPregionAnnotation.summary',header=T)

# Clean up labels by removing '_variant'
dat$annotation = gsub('_variant','',dat$annotation)

# Remove totals & make separate data frame to calculate percentages
dat = dat[-which(dat$annotation=='total'),]
percdat = dat
percentages = apply(percdat[,-1],2,function(x){(x/sum(x,na.rm = T))*100})
percdat[,colnames(percentages)] = data.table(percentages)

# Melt data frames
meltdat = melt.data.table(dat, id.vars='annotation',variable.name='dataset',value.name='n_snps')
meltpercdat = melt.data.table(percdat, id.vars='annotation',variable.name='dataset',value.name='n_snps')

# Define palettes here
pecisoqtlpalette='Blues'
oursqtlpalette='Reds'

pecisoqtl_lohicolor=c('darkblue','white','firebrick4')
oursqtl_lohicolor=c('darkred','white','aquamarine')

# sQTL datasets compared
hztl_group_bar(meltdat[which(meltdat$dataset %in% c('our_sSNPs','Walker_sSNPs','Raj_sSNPs')),],'dataset','annotation','','n_snps','# SNPs','RdYlBu','sQTL_datasets_compared.png')
hztl_group_bar(meltdat[which(meltdat$dataset %in% c('our_sSNPs','Walker_sSNPs','Raj_sSNPs')),],'dataset','annotation','','n_snps','# SNPs','RdYlBu','sQTL_datasets_compared.pdf')
hztl_group_bar(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs','Walker_sSNPs','Raj_sSNPs')),],'dataset','annotation','','n_snps','% SNPs','RdYlBu','sQTL_datasets_compared-percent.png')
hztl_group_bar(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs','Walker_sSNPs','Raj_sSNPs')),],'dataset','annotation','','n_snps','% SNPs','RdYlBu','sQTL_datasets_compared-percent.pdf')
# QTL datasets compared
hztl_group_bar(meltdat[-which(meltdat[,'dataset']=='our_nonQTLSNPs'),],'dataset','annotation','','n_snps','# SNPs','RdYlBu','QTL_datasets_compared.png')
hztl_group_bar(meltdat[-which(meltdat[,'dataset']=='our_nonQTLSNPs'),],'dataset','annotation','','n_snps','# SNPs','RdYlBu','QTL_datasets_compared.pdf')
hztl_group_bar(meltpercdat[-which(meltpercdat[,'dataset']=='our_nonQTLSNPs'),],'dataset','annotation','','n_snps','% SNPs','RdYlBu','QTL_datasets_compared-percent.png')
hztl_group_bar(meltpercdat[-which(meltpercdat[,'dataset']=='our_nonQTLSNPs'),],'dataset','annotation','','n_snps','% SNPs','RdYlBu','QTL_datasets_compared-percent.pdf')
# our sSNPs vs nonQTL SNPs compared
hztl_group_bar(meltdat[which(meltdat$dataset %in% c('our_sSNPs','our_nonQTLSNPs')),],'dataset','annotation','','n_snps','# SNPs',oursqtlpalette,'ourQTLvsNonSNPs_compared.pdf')
hztl_group_bar(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs','our_nonQTLSNPs')),],'dataset','annotation','','n_snps','% SNPs',oursqtlpalette,'ourQTLvsNonSNPs_compared-percent.pdf')
# PEC isoSNPs vs non-isoSNPs compared
hztl_group_bar(meltdat[which(meltdat$dataset %in% c('PEC_isoSNPs','PEC_nonisoSNPs')),],'dataset','annotation','','n_snps','# SNPs',pecisoqtlpalette,'PECisoQTLvsNonSNPs_compared.pdf')
hztl_group_bar(meltpercdat[which(meltpercdat$dataset %in% c('PEC_isoSNPs','PEC_nonisoSNPs')),],'dataset','annotation','','n_snps','% SNPs',pecisoqtlpalette,'PECisoQTLvsNonSNPs_compared-percent.pdf')
# Pie chart of our sSNP categories
#pie_chart(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs')),] %>% filter(!is.na(n_snps)),'annotation','Category','annotation','Category','n_snps','% SNPs','annotation','ourQTLvsNonSNPs_compared-percent_pie.pdf')
pie_chart(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs')),] %>% filter(!is.na(n_snps)),oursqtl_lohicolor,'ourQTLSNPs_compared-percent_pie.pdf')
pie_chart(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs')),] %>% filter(!is.na(n_snps)),oursqtl_lohicolor,'ourQTLSNPs_compared-percent_pie.png')
# Pie chart of PEC isoSNP categories
#pie_chart(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs')),] %>% filter(!is.na(n_snps)),'annotation','Category','annotation','Category','n_snps','% SNPs','annotation','ourQTLvsNonSNPs_compared-percent_pie.pdf')
pie_chart(meltpercdat[which(meltpercdat$dataset %in% c('PEC_isoSNPs')),] %>% filter(!is.na(n_snps)),pecisoqtl_lohicolor,'PECisoQTLSNPs_compared-percent_pie.pdf')
pie_chart(meltpercdat[which(meltpercdat$dataset %in% c('PEC_isoSNPs')),] %>% filter(!is.na(n_snps)),pecisoqtl_lohicolor,'PECisoQTLSNPs_compared-percent_pie.png')

########################################################################################################################################################################
########################################################################################################################################################################

# Test proportions of QTL vs. non-QTL SNPs in each category - functionized
testEnrichment = function(dat,cat1,cat2,comparisontag,lohicolor){
  getme = c('annotation',cat1,cat2)
  comp = dat[,..getme]
  setkey(comp,'annotation')
  
  # Remove all rows containing only NAs
  comp = comp[which(!apply(is.na(comp[,-1]),1,all)),]
  # Convert remaining NAs to 0
  comp[which(is.na(comp[,..cat1])),cat1]=0
  comp[which(is.na(comp[,..cat2])),cat2]=0
  
  # Iterate through rows and test for enrichment
  fisher.res = lapply(comp$'annotation',function(x)
  {
    print(x)
    isQTLinAnnot=comp[x][[cat1]]
    isNotQTLinAnnot=comp[x][[cat2]]
    isQTLNotinAnnot=sum(comp[!x][[cat1]])
    isNotQTLNotinAnnot=sum(comp[!x][[cat2]])
    cont.table=data.frame(matrix(c(isQTLinAnnot,isQTLNotinAnnot,isNotQTLinAnnot,isNotQTLNotinAnnot),byrow=T,nrow=2))
    print(cont.table)
    test=fisher.test(cont.table, or=1, alternative='two.sided', conf.int=TRUE, conf.level=0.95, workspace=1e09)
    return(data.frame(annotation=x,pval=test$p.value,OR=test$estimate,confintlow=test$conf.int[1],confinthigh=test$conf.int[2]))
  }
  )
  # Make results table
  fisher.res.df = do.call(rbind,fisher.res)
  
  # FDR correction (Benjamini-Hochberg)
  #fisher.res.df$qval = qvalue(fisher.res.df$pval,pi0=1)$qvalues
  
  # FDR correction (Bonferroni)
  fisher.res.df$qval = p.adjust(fisher.res.df$pval,method='bonferroni')
  
  # Take log2 of OR
  fisher.res.df$ORlog2 = log2(fisher.res.df$OR)
  
  # Take -log10 of pvals
  fisher.res.df$log10qval = -log10(fisher.res.df$qval) 
  
  # Change factor order of annotation column in decreasing order of OR
  fisher.res.df$annotation = factor(fisher.res.df$annotation, levels = fisher.res.df$annotation[order(fisher.res.df$OR,decreasing=TRUE)])
  
  # Add label for significant tests
  fisher.res.df$barlabel = '*'
  fisher.res.df$barlabel[which(fisher.res.df$qval>=.05)] = ''
  
  # Remove categories with OR=Inf
  fisher.res.df = fisher.res.df[which(is.finite(fisher.res.df$OR)),]
  
  # Make plots
  abovebelow0_bar(fisher.res.df,'qval','adj.pvalue','ORlog2','log2(OR)','annotation','Variant category','barlabel',paste0(comparisontag,'_categories_tested.png'))
  abovebelow0_bar(fisher.res.df,'qval','adj.pvalue','ORlog2','log2(OR)','annotation','Variant category','barlabel',paste0(comparisontag,'_categories_tested.pdf'))
  
  # Non-log2 transformed
  vert_regular_bar(fisher.res.df,'log10qval',lohicolor[1],lohicolor[2],'right','-log10(pval)','OR','OR',5,'annotation','','barlabel','confintlow','confinthigh',paste0(comparisontag,'categories_tested-bonferroni.png'))
  vert_regular_bar(fisher.res.df,'log10qval',lohicolor[1],lohicolor[2],'right','-log10(pval)','OR','OR',5,'annotation','','barlabel','confintlow','confinthigh',paste0(comparisontag,'_categories_tested-bonferroni.pdf'))
  # Plot version w/o coloring by pval
  fisher.res.df['color']=0
  vert_regular_bar(fisher.res.df,'color',lohicolor[1],lohicolor[2],'none','color','OR','OR',12,'annotation','','barlabel','confintlow','confinthigh',paste0(comparisontag,'_categories_tested-bonferroni_nocolor.png'))
  vert_regular_bar(fisher.res.df,'color',lohicolor[1],lohicolor[2],'none','color','OR','OR',12,'annotation','','barlabel','confintlow','confinthigh',paste0(comparisontag,'_categories_tested-bonferroni_nocolor.pdf'))
  return(fisher.res.df)
}
# Run enrichment tests
oursnps_test = testEnrichment(dat,'our_sSNPs','our_nonQTLSNPs','PEC_isoQTLSNPs',oursqtl_lohicolor)
pecisosnps_test = testEnrichment(dat,'PEC_isoSNPs','PEC_nonisoSNPs','PEC_isoQTLSNPs',pecisoqtl_lohicolor)

# Tag the enrichment test results with the data source & combine them
oursnps_test[,'dataset'] = 'OurStudy_sSNPs'
pecisosnps_test[,'dataset'] = 'PEC_isoSNPs'

# To find colors:
# Call vert_regular_bar as above with desired low/high colors and save plot to p
# Call ggplot_build(p)$data[[4]]['fill'] to see color assigned
combined_palette = c('OurStudy_sSNPs'='#D18978','PEC_isoSNPs'='#977DC6')

enr.df=rbind(oursnps_test,pecisosnps_test)

# Set factor levels manually to make categories appear in 5'->3' order
enr.df$annotation = ordered(enr.df$annotation, levels = c('regulatory_region','TF_binding_site','intergenic_variant','intergenic','upstream_gene','5_prime_UTR','start_lost','synonymous','inframe_insertion','missense','splice_acceptor','splice_donor','splice_region','intron','stop_lost','stop_gained','stop_retained','incomplete_terminal_codon','3_prime_UTR','downstream_gene','non_coding_transcript_exon','mature_miRNA'))

# Plot params
yaxlim=12

p=ggplot(data = enr.df, aes_string(fill = 'dataset', x = 'annotation', y = 'OR'), width = 0.6) +
  geom_bar(colour='black', stat = 'identity', position='dodge', width=0.7) +
  geom_hline(yintercept=1, linetype = 'dashed') +
  geom_errorbar(aes_string(ymax = 'confinthigh', ymin = 'confintlow'), position = position_dodge(width=0.5), width = 0.2) +
  geom_text(aes_string(label='barlabel'), color = 'red', position = position_dodge(width=0.5), size = 8) +
  scale_fill_manual(values=combined_palette) +
  xlab('Annotation') +
  scale_y_continuous(name = 'OR', limits=c(0,yaxlim), breaks=seq(0,yaxlim,by=1), expand = c(0,.1)) +
  #	geom_text(data = label.df, label = '*') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color='black', size=0.6),
    axis.line.y = element_line(color='black', size=0.6),
    axis.text.x = element_text(colour='black', size=16, angle=90, hjust=1),
    axis.text.y = element_text(colour='black', size=16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=17),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.position = c(.85,.9),
    legend.text = element_text(size=15),
    plot.margin=unit(c(0.25,1,0.25,1), 'cm')
  )
plot(p)
aspect_ratio = 1.5
height = 5
ggsave('VEP_cat_enrichments_combined.png',p, height = 7 , width = 7 * aspect_ratio)
ggsave('VEP_cat_enrichments_combined.pdf',p, height = 7 , width = 7 * aspect_ratio)
