#! /usr/bin/env Rscript

# Adapted from https://urldefense.proofpoint.com/v2/url?u=https-3A__github.com_francois-2Da_fastqtl_blob_master_python_run-5FFastQTL-5Fthreaded.py&d=DwIGaQ&c=shNJtf5dKgNcPZ6Yh64b-A&r=OC3vB9j9iFQIMjLzUNsmpHPMq_MZjb0HpAie3TakjH8&m=o0sL6tRyRgWPli9u-FjVA8MaM8xA5sMpQWiVcmXMVXI&s=0-PVjvD9XtSSEvRnxC0BS1Ynf0PiIuL6hUUW8RyN6XU&e= 
# and https://urldefense.proofpoint.com/v2/url?u=https-3A__github.com_hobrien_GENEX-2DFB2_blob_master_R_calulateNominalPvalueThresholds.R&d=DwIGaQ&c=shNJtf5dKgNcPZ6Yh64b-A&r=OC3vB9j9iFQIMjLzUNsmpHPMq_MZjb0HpAie3TakjH8&m=o0sL6tRyRgWPli9u-FjVA8MaM8xA5sMpQWiVcmXMVXI&s=cdzP5HZOVDDCrKiqnS8cEB774T63XerCu_KbJU8VqGI&e= 

# Call egene/isotx/sintron(sqtl) form permutation FastQTL pass;
# and calculate nominal p-val threshold for qtl

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(dplyr))

p = arg_parser('Annotates FastQTL permutation output and runs qvalue')
p = add_argument(p, '--resdir', help='Dir containing the chrAll_combined file')
p = add_argument(p, '--level', help='egene, isotx, sqtl')
# p = add_argument(p, '--lambda', type='numeric', help='', default=NULL)
args = parse_args(p)
resdir = args$resdir
level = args$level

results = fread(paste0(resdir,'/chrAll_combined'), data.table=F)
colnames(results) = c('pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval', 'slope', 'ppval', 'bpval')
cat('-----> FastQTL results file dimension is: ', dim(results)[1], ' * ', dim(results)[2], '\n', sep='')
# 1. remove phenotypes w/o cis-variants
results = results[complete.cases(results),]
cat('-----> Complete dimension is: ', dim(results)[1],' * ', dim(results)[2], '\n', sep='')
# 2. remove duplicates
# keep phenotype with lowest nominal p-value (only applies when top_snp not tested in one dup)
# if nominal p identical (because top SNP tested in both), keep lowest corrected p-value (because window includes more SNPs)
results = arrange(results, npval, desc(bpval)) %>% group_by(pid) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()
cat('-----> Dup-removed dimension is: ', dim(results)[1], ' * ', dim(results)[2], '\n', sep='')

# 3. calculate q-values
# ‘signif’ rounds the values in its first argument to the specified number of significant digits
Q = qvalue(results[,'bpval'])
results$qval = signif(Q$qvalues, 6)
cat('-----> Proportion of significant phenotypes (1-pi0): ' , round((1 - Q$pi0), 2), '\n', sep='')
cat('-----> Number of ', level, ' @ FDR 0.05: ', sum(results[, 'qval']<0.05), '\n', sep='')

# 4. Determine globala min(p) significance threshold and calculate nominal p-val threshold for each phenotype
ub = sort(results[results$qval > .05, 'bpval'])[1]  # smallest p-value above FDR
lb = -sort(-results[results$qval <= .05, 'bpval'])[1]  # largest p-value below FDR
pthreshold = (lb+ub)/2
cat('-----> smallest p-value above FDR: ', ub, '\n')
cat('-----> largest p-value below FDR: ', lb, '\n')
cat('-----> min p-value threshold @ FDR 0.05: ', pthreshold, '\n', sep='')
results[, 'pval_nominal_threshold'] = signif(qbeta(pthreshold, results[, 'shape1'], results[, 'shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

sig_pheno = filter(results, qval<=.05)


# if(nrow(filter(sig_egenes, pval_nominal > pval_true_df)) >0) {
#     cat(nrow(filter(sig_egenes, pval_nominal > pval_true_df)), ' eGenes with df_corrected p-value < nominal p-value\n', sep='')
#     cat(nrow(filter(sig_egenes, pval_nominal > pval_nominal_threshold)), ' sig eGenes have nominal p > threshold. Modifying threshold so top eQTL will be sig\n', sep='')
#     sig_egenes = sig_egenes %>% mutate(pval_nominal_threshold=ifelse(pval_nominal_threshold >= pval_nominal, pval_nominal_threshold, pval_nominal))
# }

write.table(sig_pheno, paste0(resdir,'/chrAll_combined.FDR05'), quote=F, sep='\t', col.names=T, row.names=F)


# 5. Write counts
num_pheno = length(unique(sig_pheno$pid))
result_df = data.frame('num_pheno' = c(num_pheno))
write.table(result_df, paste0(resdir,'/chrAll_combined.FDR05_count'), quote=F, sep='\t', col.names=F, row.names=F)
