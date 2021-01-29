chrsizes=/sc/arion/work/belmoj01/lncrna_map/isoseq_map_analysis/analysis/track_hub/hg19.sizes
allphenosbed=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/Phenotype.bed
allsnpsbed=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY.snps.bed
windowsize=100000

./identify_SNPs_in_cis_window.sh /sc/arion/work/belmoj01/lncrna_map/isoseq_map_analysis/analysis/track_hub/hg19.sizes /sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/Phenotype.bed /sc/arion/projects/EPIASD/splicingQTL/intermediate_files/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY.snps.bed 100000

./identify_SNPs_in_cis_window.sh $chrsizes $allphenosbed $allsnpsbed $windowsize 
