# Read all BEDTOOLs overlaps for each Gencode region and annotate which overlaps each SNP is found in. De-redundinate as necessary, then summarize results
args = commandArgs(trailingOnly=TRUE)

snpfile = args[1] #'/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.uniqueSNPs.bed'
overlapdir = args[2] #'/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/gencodeOverlaps_genesOnly-sGene/'
# OR /sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoqtls/gencodeOverlaps_genesOnly-sGene/
ensginresfile = args[3] #TRUE if argument 4 is fastQTL results summary file to get sGene ID from 
sigsqtlfile = args[4] #'/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl+ensg.csv'

# ENSG -> ENST map to assign genes to Gencode BED file
ensemblfile = '/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/gencode_ENSG_to_ENST.txt'
ensembl = read.table(ensemblfile, sep='\t', header=T)

snps = read.table(snpfile, sep='\t', col.names=c('chr','start','end','snp_id','clust_id','strand'))

# Join significant sQTLs summary file to get the associated sGene ID
# NOTE: after merge, nrow(SNPs) == # sQTLs
if(ensginresfile){
  sigsqtl = read.table(sigsqtlfile, sep='\t',header=T)
  snps = merge(snps, sigsqtl, by.x='snp_id', by.y='sid', all.x=TRUE)
}

files = list.files(path = overlapdir, full.names = TRUE, recursive = FALSE)

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
  bed = read.table(bedfile, sep='\t',col.names=c('chr','start','end','snp_id','clust_id','strand','chr_gc','start_gc','end_gc','id_gc','score_gc','strand_gc','cdsstart_gc','cdsend_gc','rgb_gc','bcount_gc','bstart_gc','blen_gc','len_overlap'))
  bed$id_gc = do.call(rbind,strsplit(as.character(bed$id_gc),'\\.'))[,1] 
  bed = merge(bed, ensembl, by.x='id_gc', by.y='Transcript.stable.ID', all.x=TRUE)
  if(ensginresfile){
    return(as.character(snpdf[,'ensg'])==as.character(bed[match(snpdf[,'snp_id'],bed[,'snp_id']),'Gene.stable.ID'])) 
  } else {
    return(snpdf[,'snp_id'] %in% bed[,'snp_id'])
  }
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
