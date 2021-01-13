library(dplyr)
library(stringr)
library(argparser)
library('doParallel') # I do this to make %dopar% work
source('/hpc/packages/minerva-common/leafcutter/1.0/leafcutter/leafcutter/R/utils.R')

p = arg_parser('Read list of phenotypes (chr:start:end:cluster IDs) and make map of clusters->gene_ids using functions from the leafcutter utils script')
#p = ArgumentParser('Read list of phenotypes (chr:start:end:cluster IDs) and make map of clusters->gene_ids using functions from the leafcutter utils script')
#p$add_argument('--phenos', metavar='phenos', type='character', help='List of cluster IDs') # e.g. '/sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts_clustids'
p = add_argument(p, '--phenos', help='List of cluster IDs to annotate', type='character')
#p$add_argument('--gtf', metavar='gtf', type='character', help='Path to GENCODE GTF')
p = add_argument(p, '--gtf', help='Path to GENCODE GTF', type='character')
#p$add_argument('--out', metavar='out', type='character', help='Name of output file to write')
p = add_argument(p, '--out', help='Name of output file to write', type='character')

#args = p$parse_args()
args = parse_args(p)
phenos = args$phenos
gtf = args$gtf
out = args$out

clu = read.table(phenos,header=T)
gc.gtf = read.table(gtf,header=F,comment.char='#',sep='\t',col.names=c('chr','source','type','start','end','score','strand','frame','attr'))

# Keep only exon features
gc.gtf = gc.gtf[gc.gtf$type=='exon',]

# Extract gene ID & add it as new column
gene_id_pat = 'gene_id ([^;]*);'
gc.gtf$gene_name = str_match(gc.gtf$attr, gene_id_pat)[,2]

# Parse clusters to data frame
clu.df = get_intron_meta(as.character(clu[,1]))

# Parse .gtf to data frame
#chr start end strand gene_name
#chr1 11869 12227 + DDX11L1
#chr1 12613 12721 + DDX11L1
#chr1 13221 14409 + DDX11L1
#chr1 11872 12227 + DDX11L1

clu.df.gc = map_clusters_to_genes(clu.df,gc.gtf)

# Write out the annotated output
write.table(clu.df.gc, out, sep='\t', row.names=F, col.names=T, quote=F) 
