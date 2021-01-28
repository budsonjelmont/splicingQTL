#!/bin/bash
ml plink/1.90

invcf=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY.vcf
idout=${invcf%.*}

plink --vcf $invcf --double-id --write-snplist --make-bed --out $idout

# Make bed file of all snps 
awk 'BEGIN{OFS=FS="\t"};{print "chr"$1,$4-1,$4,$2}' $idout.bim > $idout.snps.bed

