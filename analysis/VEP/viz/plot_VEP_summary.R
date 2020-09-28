setwd('~/Documents/sqtl/analysis/VEP/')

library(data.table)   
library(qvalue)
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

# sQTL datasets compared
hztl_group_bar(meltdat[which(meltdat$dataset %in% c('our_sSNPs','Walker_sSNPs','Raj_sSNPs')),],'dataset','annotation','','n_snps','# SNPs','RdYlBu','sQTL_datasets_compared.png')
hztl_group_bar(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs','Walker_sSNPs','Raj_sSNPs')),],'dataset','annotation','','n_snps','% SNPs','RdYlBu','sQTL_datasets_compared-percent.png')
# QTL datasets compared
hztl_group_bar(meltdat[-which(meltdat[,'dataset']=='our_nonQTLSNPs'),],'dataset','annotation','','n_snps','# SNPs','RdYlBu','QTL_datasets_compared.png')
hztl_group_bar(meltpercdat[-which(meltpercdat[,'dataset']=='our_nonQTLSNPs'),],'dataset','annotation','','n_snps','% SNPs','RdYlBu','QTL_datasets_compared-percent.pdf')
# our sSNPs vs nonQTL SNPs compared
hztl_group_bar(meltdat[which(meltdat$dataset %in% c('our_sSNPs','our_nonQTLSNPs')),],'dataset','annotation','','n_snps','# SNPs','RdYlBu','ourQTLvsNonSNPs_compared.pdf')
hztl_group_bar(meltpercdat[which(meltpercdat$dataset %in% c('our_sSNPs','our_nonQTLSNPs')),],'dataset','annotation','','n_snps','% SNPs','RdYlBu','ourQTLvsNonSNPs_compared-percent.png')

# Test proportions of QTL vs. non-QTL SNPs in each category
ourSNPs = dat[,c('annotation','our_nonQTLSNPs','our_sSNPs')]
setkey(ourSNPs,'annotation')

# Remove all rows containing only NAs
ourSNPs = ourSNPs[-which(apply(is.na(ourSNPs[,-1]),1,all)),]
# Convert remaining NAs to 0
ourSNPs[is.na(ourSNPs$our_sSNPs),'our_sSNPs']=0
ourSNPs[is.na(ourSNPs$our_nonQTLSNPs),'our_nonQTLSNPs']=0

# Iterate through rows and test for enrichment
fisher.res = lapply(ourSNPs$'annotation',function(x)
  {
    print(x)
    isQTLinAnnot=ourSNPs[x]$'our_sSNPs'
    isNotQTLinAnnot=ourSNPs[x]$'our_nonQTLSNPs'
    isQTLNotinAnnot=sum(ourSNPs[!x,'our_sSNPs'])
    isNotQTLNotinAnnot=sum(ourSNPs[!x,'our_nonQTLSNPs'])
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
fisher.res.df = fisher.res.df[-which(!is.finite(fisher.res.df$OR)),]

# Make plots
abovebelow0_bar(fisher.res.df,'qval','adj.pvalue','ORlog2','log2(OR)','annotation','Variant category','barlabel','ourSNPs_categories_tested.png')
abovebelow0_bar(fisher.res.df,'qval','adj.pvalue','ORlog2','log2(OR)','annotation','Variant category','barlabel','ourSNPs_categories_tested.pdf')

# Non-log2 transformed
vert_regular_bar(fisher.res.df,'log10qval','-log10(pval)','OR','OR',5,'annotation','','barlabel','confintlow','confinthigh','ourSNPs_categories_tested-bonferroni.png')
vert_regular_bar(fisher.res.df,'log10qval','-log10(pval)','OR','OR',5,'annotation','','barlabel','confintlow','confinthigh','ourSNPs_categories_tested-bonferroni.pdf')
