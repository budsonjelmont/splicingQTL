library(qvalue)

args = commandArgs(trailingOnly=TRUE)

resdir = args[1]
#resdir = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/'
infile = 'chrAll_combined'

d = read.table(paste(resdir,infile,sep=''), header=F, stringsAsFactors=F)

colnames(d) = c('pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval','slope', 'ppval', 'bpval')
png(paste(resdir,'nominal pvalue vs Beta-approximated pvalue.png',sep=''))
plot(d$ppval, d$bpval, xlab='Direct method', ylab='Beta approximation', main='Check plot')
abline(0, 1, col='red')
dev.off()

#d$qval = qvalue(d$bpval)$qvalues #Storey
d$qval = p.adjust(d$bpval, method='fdr') #Benjamini-Hochberg 

d[,c('chr','start','end','clust')]= do.call(rbind,strsplit(d$pid,':'))

write.table(d, paste(resdir,'chrAll_combined+qval',sep=''), quote=F, row.names=F, col.names=T, sep='\t')
write.table(d[which(d$qval <= 0.05),], paste(resdir,'significant_sqtl.txt',sep=''), quote=F, row.names=F, col.names=T, sep='\t')

# count of significant sQTLs
print(paste('# significant sQTLs: ',as.character(dim(d[which(d$qval < .05),])[1])),sep='')

# count of unique clusters
print(paste('# of unique clusters: ',as.character(length(unique(d[which(d$qval < .05),'clust']))),sep=''))

# count of unique SNPs
print(paste('# of unique SNPs: ',as.character(length(unique(d[which(d$qval < .05),'sid']))),sep=''))

# count of unique splicing events
print(paste('# of unique splicing events: ',as.character(length(unique(d[which(d$qval < .05),'pid']))),sep=''))
