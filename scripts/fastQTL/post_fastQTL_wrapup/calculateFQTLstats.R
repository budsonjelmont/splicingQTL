library(qvalue)

args = commandArgs(trailingOnly=TRUE)

resdir = args[1] # e.g. '/sc/arion/projects/EPIASD/splicingQTL/output/fqtl_output_wasp_nominal_normal/minCovars+seqPC9/20genoPCs/55HCPs/'
#isnominal = args[2] # e.g. TRUE if fastQTL was run without the --permute flag (not used currently, count columns to determine if results are nominal or permuted)
fdrmethod = args[2] # {'st','bh'}
infile = 'chrAll_combined'

d = read.table(paste(resdir,'/',infile,sep=''), header=F, stringsAsFactors=F)

# Count columns to determine what kind of results these are
if(ncol(d)==5){ # Results are from nominal pass
  colnames(d) = c('pid', 'sid', 'dist', 'bpval','slope')
} else if (ncol(d)==11){ # Results are from permutation pass
  colnames(d) = c('pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval','slope', 'ppval', 'bpval')
  png(paste(resdir,'/nominal pvalue vs Beta-approximated pvalue.png',sep=''))
  plot(d$ppval, d$bpval, xlab='Direct method', ylab='Beta approximation', main='Check plot')
  abline(0, 1, col='red')
  dev.off()
} else {
  print('ERROR: unexpected number of columns in fastQTL output (should be 5 or 11)')
  quit(status=0)
}

# FDR correct using whichever method  was specified
if(fdrmethod=='st'){ #Storey
  d$qval = qvalue(d$bpval)$qvalues
} else if(fdrmethod=='bh'){ #Benjamini-Hochberg 
  d$qval = p.adjust(d$bpval, method='fdr') 
} else {
  print('ERROR: arg2 must specify Storey (\'st\') or Benjamini-Hochberg (\'bh\') correction')
  quit(status=0)
}

d[,c('chr','start','end','clust')]= do.call(rbind,strsplit(d$pid,':'))

write.table(d, paste(resdir,'/','chrAll_combined+qval',sep=''), quote=F, row.names=F, col.names=T, sep='\t')
write.table(d[which(d$qval <= 0.05),], paste(resdir,'/','significant_sqtl.txt',sep=''), quote=F, row.names=F, col.names=T, sep='\t')

# count of significant sQTLs
print(paste('# significant sQTLs: ',as.character(dim(d[which(d$qval < .05),])[1])),sep='')

# count of unique clusters
print(paste('# of unique clusters: ',as.character(length(unique(d[which(d$qval < .05),'clust']))),sep=''))

# count of unique SNPs
print(paste('# of unique SNPs: ',as.character(length(unique(d[which(d$qval < .05),'sid']))),sep=''))

# count of unique splicing events
print(paste('# of unique splicing events: ',as.character(length(unique(d[which(d$qval < .05),'pid']))),sep=''))
