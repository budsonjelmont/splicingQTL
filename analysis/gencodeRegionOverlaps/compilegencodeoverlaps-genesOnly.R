# Read all BEDTOOLs overlaps for each Gencode region and annotate which overlaps each SNP is found in. De-redundinate as necessary, then summarize results
args = commandArgs(trailingOnly=TRUE)

snpfile = args[1] #'/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.uniqueSNPs.bed'
overlapdir = args[2] #'/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/gencodeOverlaps/'
# OR /sc/orga/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoqtls/gencodeOverlaps/

snps = read.table(snpfile, sep='\t',col.names=c('chr','start','end','snp_id','clust_id','strand'))

files = list.files(path = overlapdir, full.names = TRUE, recursive = FALSE)

# Get suffix of each file
getSuffix = function(x){
  split = strsplit(x,'_')
  suffix = strsplit(split[[1]][length(split[[1]])],'\\.')[[1]][1]
  return(suffix)
}

categories = unlist(lapply(files, getSuffix))

# Initialize data frame that will hold category counts 
catdf = data.frame(cbind(files,categories)) 

# Add columns to SNPs data frame summarizing genomic locations
getOverlap = function(cats,snpdf){
  bedfile = cats['files']
  bed = read.table(bedfile, sep='\t',col.names=c('chr','start','end','snp_id','clust_id','strand'))
  return(snpdf[,'snp_id'] %in% bed[,'snp_id'])
}

res = apply(catdf,1,getOverlap,snps)
colnames(res) = catdf$categories
rownames(res) = 1:nrow(res)

# Merge
snps = merge(snps, as.data.frame(res), by='row.names', all=TRUE)

# Change column names
categories = c('inGene','outsideGene')

colnames(snps)[8] = categories[1]

# Add outside gene column
snps['outsideGene']=!snps['inGene']

# Make summary table
regiondf = data.frame(region=categories)

regiondf['nsnps']=unlist(lapply(categories,function(x){return(sum(snps[,x]))}))
regiondf['nsnps%']=(regiondf['nsnps']/nrow(snps))*100

rownames(regiondf)=regiondf[,'region']

write.table(regiondf, paste0(overlapdir,'../gencode_inoutgene_summary.csv'), sep=',', col.names=T, row.names=F, quote=F)
 
# Make table of region counts
# Make summary plot
library(ggplot2)
