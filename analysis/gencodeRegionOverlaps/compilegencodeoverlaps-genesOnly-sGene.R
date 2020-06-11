# Read all BEDTOOLs overlaps for each Gencode region and annotate which overlaps each SNP is found in. De-redundinate as necessary, then summarize results
args = commandArgs(trailingOnly=TRUE)

sigsqtlfile = args[1] # '/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/WalkerCell_significantSQTLS+ensg.txt'
# '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl+ensg.csv'
overlapdir = args[2] # '/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/gencodeOverlaps_genesOnly-sGene/'
# '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/gencodeOverlaps_genesOnly-sGene/'
snpidcolname = args[3] # 'variant_id' for Walker data, 'sid' for our sQTLs

# ENSG -> ENST map to assign genes to Gencode BED file
ensemblfile = '/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/gencodeV19/gencode_ENSG_to_ENST.txt'
ensembl = read.table(ensemblfile, sep='\t', header=T)

snps = read.table(sigsqtlfile, sep='\t',header=T)

# Join significant sQTLs summary file to get the associated sGene ID
# NOTE: after merge, nrow(SNPs) == # sQTLs
if (!'snp_id' %in% colnames(snps)){
  colnames(snps)[colnames(snps)==snpidcolname] = 'snp_id'
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
rownames(catdf)=catdf$categories

# Add columns to SNPs data frame summarizing genomic locations
getOverlap = function(cats,snpdf){
  bedfile = cats['files']
  category = as.character(cats['categories'][[1]])
  bed = read.table(bedfile, sep='\t',col.names=c('chr','start','end','snp_id','clust_id','strand','chr_gc','start_gc','end_gc','id_gc','score_gc','strand_gc','cdsstart_gc','cdsend_gc','rgb_gc','bcount_gc','bstart_gc','blen_gc','len_overlap'))
  bed$id_gc = do.call(rbind,strsplit(as.character(bed$id_gc),'\\.'))[,1]
  # Merge BEDTOOLS output with ENST->ENSG map to get gene IDs
  bed = merge(bed, ensembl, by.x='id_gc', by.y='Transcript.stable.ID', all.x=TRUE)
  m=merge(snpdf, bed[!duplicated(bed[,c('snp_id','Gene.stable.ID')]),c('snp_id','Gene.stable.ID','len_overlap')], by.x=c('snp_id','ensg'), by.y=c('snp_id','Gene.stable.ID'), all.x=T, suffixes = c('', '_bed'))
  m[category] = !is.na(m[,'len_overlap'])
  return(m[,c('snp_id','tscrid',category)])
  #return(!is.na(m[,'len_overlap'])) 
}

res = apply(catdf,1,getOverlap,snps)
# Add new columns to SNP data frame
snps=do.call(
 cbind,lapply(names(res),
  function(category){  
    merge(snps,res[[category]], by.x=c('snp_id','tscrid'), by.y=c('snp_id','tscrid'))
  }
 )
)
snps = cbind(snps,newcols)

# Change column names
categories = c('inGene','outsideGene')

colnames(snps)[colnames(snps) %in% catdf$categories] = categories[1]

# Add outside gene column
snps['outsideGene']=!snps['inGene']

# Now make summary columns of sSNPs falling within their sGene only
snps[is.na(snps['inGene']),'inGene']=FALSE
snps[is.na(snps['outsideGene']),'outsideGene']=TRUE
fnamesplit = unlist(strsplit(sigsqtlfile,'\\.'))
write.table(snps, paste0(fnamesplit[1],'+sSNPloc.',fnamesplit[2]), sep='\t', col.names=T, row.names=F, quote=F)

# Make summary table of sSNP locations
regiondf = data.frame(region=categories)
regiondf['nsnps']=unlist(lapply(categories,function(x){return(sum(snps[,x]))}))
regiondf['nsnps%']=(regiondf['nsnps']/nrow(snps))*100
write.table(regiondf, paste0(overlapdir,'../gencode_inoutgene_summary_insidesGeneONLY.csv'), sep=',', col.names=T, row.names=F, quote=F)
 
# Make table of region counts
# Make summary plot
library(ggplot2)
