# Directory containing the fastQTL output
resdir=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_0HCPs

# Combine results for all chromosomes into a single file
cat $resdir/chr*/chr*_norm_perm1000.txt.gz.* > $resdir/chrAll_combined

# Calculate stats on sQTL results
Rscript calculateFQTLstats.R $resdir/

# Make files for downstream analyses 
python3 make_analysis_files-ourSNPs.py $resdir/

# Add SNP/phenotype gene annotations
../analysis/addQTLgenes/count_eGenes.lsf $resdir/ sqtlphenos_unique.bed qtls.txt pid
../analysis/addQTLgenes/count_eGenes.lsf $resdir/ sqtlSNPS_unique.bed qtls+pid_ensg.txt sid

# Add boolean column to indicate if SNP & phenotype are located in the same gene
../analysis/addQTLgenes/get_counts_add_samegeneTF.py $resdir/
