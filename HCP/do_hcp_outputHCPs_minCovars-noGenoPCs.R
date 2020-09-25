args = commandArgs(trailingOnly=TRUE)

.libPaths( c( .libPaths(), "/hpc/users/belmoj01/.Rlib") )
library(Rhcpp)

# Read input files
#countsfile = '/sc/hydra/scratch/belmoj01/splicingQTL/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.05_perind.counts.idsync.deduped.gz.qqnorm_AllCombined'
countsfile = args[1] #'/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.deduped.gz.qqnorm_allCombined'
#countsfile = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_input/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.gz.qqnorm.AllCombined'
#covarsfile = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/covar/leafcutter-input_covar.20genoPCs.idsync.10splicingPCs.deduped.minCovars.txt'
covarsfile = args[2] # '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/covar_wasp/leafcutter-input_covar_WASP.20genoPCs.idsync.deduped.minCovars+seqPC9+20genoPCs.txt'
numhcps = as.numeric(args[3]) #5

counts = read.table(countsfile, header=TRUE, comment.char='@', check.names=FALSE, stringsAsFactors=FALSE)
# If concatenated file contains multiple header rows due to concatenation, drop all them
counts = counts[!counts$ID=='ID',]
# Drop first row after setting row names
#counts = counts[-1,]
rownames(counts) = counts$ID
counts = counts[ , !(colnames(counts) %in% colnames(counts)[1:4])]

covars_og = read.table(covarsfile, sep='\t', header=TRUE, comment.char='',stringsAsFactors=FALSE, check.names=FALSE)
covars = covars_og
rownames(covars) = covars[,1]
covars = covars[,-1]

# Transpose data frames
counts = as.data.frame(t(counts))
covars = as.data.frame(t(covars))

# Convert non-categorical columns from factors -> numbers, and dummy code categorical columns
catcols = c('study','sex','tissue')

for (i in colnames(covars)){
  if(!i %in% catcols){
   covars[,i] = as.numeric(as.character(covars[,i]))
  } else {
   f = as.formula(paste('~',i,sep=''))
   mm = model.matrix(f, data = covars)
   mm = mm[,2:ncol(mm),drop=F]
   covars[,colnames(mm)] = mm
  }
}
# Drop non-coded categorical columns
covars = covars[,!(colnames(covars) %in% catcols)]

# HCP params
excluded = c(paste('SV',1:4,sep=''),paste('splicingPC',1:20,sep=''), paste('genotypePC',1:20,sep=''))

Z = as.matrix(covars[,!(colnames(covars) %in% excluded)])
#Y = as.matrix(counts)
Y = apply(as.matrix(counts),2,as.numeric)
k = numhcps  # Walker uses 20 HCPs
lambda1 = 1
lambda2 = 1
lambda3 = 1

iter = 100

# Read expected result
res = hcp(Z, Y, k, lambda1, lambda2, lambda3, iter, fast=FALSE, stand=TRUE, log=FALSE)

summary(res)
hcps = data.frame(t(res$W))

# Add x HCPs to covariate data frame
hcps['id'] = paste('HCP',1:nrow(hcps),sep='')
write.table(hcps,paste(as.character(numhcps),'HCPs.txt',sep=''),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
