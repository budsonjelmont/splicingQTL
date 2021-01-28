#! /usr/bin/env Rscript

# From FastQTL nominal pass
# Call eGene and eQTL/isoTx and isoQTL/sQTL and sQTL(intron)/sSNP

suppressMessages(library(argparser))
p = arg_parser('Calls eGene/eQTL or isoTx/isoQTL or sQTLsSNP from FastQTL nominal output')
p = add_argument(p, '--resdir', help='Dir containing the chrAll_combined file')
args = parse_args(p)
resdir = args$resdir

suppressMessages(library(data.table))
results = fread(paste0(resdir,'/chrAll_combined'), data.table=F)
colnames(results) = c('pid','sid','dist','npval','slope')
results$fdr = p.adjust(results$npval,method='fdr')
significant = results[results$fdr<.05,]

num_pheno = length(unique(significant$pid))
result_df = data.frame(num_pheno = c(num_pheno))
write.table(significant, paste0(resdir,'/chrAll_combined.FDR05'), quote=F, sep='\t', col.names=T, row.names=F)
write.table(result_df, paste0(resdir,'/chrAll_combined.FDR05_count'), quote=F, sep='\t', col.names=F, row.names=F)
