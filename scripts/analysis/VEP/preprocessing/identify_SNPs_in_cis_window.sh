ml bedtools

chrsizes=/sc/arion/work/belmoj01/lncrna_map/isoseq_map_analysis/analysis/track_hub/hg19.sizes

#allphenosbed=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/Phenotype.bed
#allsnpsbed=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY.snps.bed
allphenosbed=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PEC/PEC_all_genes.bed
allsnpsbed=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PEC/SNP_Information_Table_with_Alleles_all+strand.bed

windowsize=1000000

phenossloppedbed=${allphenosbed%.*}_plusmin${windowsize}bp.bed
inwindowsnpsbed=${allsnpsbed%.*}_in${windowsize}bpwindow.bed

# Expand boundaries of phenotypes to include full cis windows tested by FastQTL
bedtools slop -i $allphenosbed -g $chrsizes -b $windowsize > $phenossloppedbed

# Intersect SNP list with slopped phenotypes & extract the IDs of all overlapping SNPs from the results. This is the list of SNPs fall within range to be tested by FastQTL.
bedtools intersect -a $allsnpsbed -b $phenossloppedbed -u  > $inwindowsnpsbed
