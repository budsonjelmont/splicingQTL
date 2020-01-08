library(qvalue)

resdir = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/'
infile = 'chrAll_combined'

d = read.table(paste(resdir,infile,sep=''), hea=F, stringsAsFactors=F)

colnames(d) = c('pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval', 'ppval', 'bpval')
png(paste(resdir,'nominal pvalue vs Beta-approximated pvalue.png',sep=''))
plot(d$ppval, d$bpval, xlab='Direct method', ylab='Beta approximation', main='Check plot')
abline(0, 1, col='red')
dev.off()

d$qval = qvalue(d$bpval)$qvalues

write.table(d[which(d$qval <= 0.05),], paste(resdir,'significant_sqtl.csv',sep=''), quote=F, row.names=F, col.names=T)

# count of unique SNPs
print(paste('# of unique SNPs: ',as.character(length(unique(d[which(d$qval < .05),'sid']))),sep=''))

# count of unique splicing events
print(paste('# of unique splicing events: ',as.character(length(unique(d[which(d$qval < .05),'pid']))),sep=''))
