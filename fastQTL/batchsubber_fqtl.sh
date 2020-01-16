#!/bin/bash

covPath='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/covar_wasp/'
covFile='leafcutter-input_covar.20genoPCs.idsync.10splicingPCs.deduped.minCovars.40HCPs_nogeno.txt'
#vcfPath='/sc/orga/projects/EPIASD/PsychENCODE_phase12_genotypes/Capstone4VCFreformat/'
vcfPath='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/geno/'
vcfFile='Capstone4.sel.idsync.2allele.maf5.mind1.geno1.hwe1e-5.deduped.vcf.gz'
bedPath='/sc/hydra/scratch/belmoj01/splicingQTL/'
bedFile='out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.05_perind.counts.idsync.deduped.gz'
outPath='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars_40HCPs/'

# List the contents of the directory passed as argument 1 and read all files there into an array 
chrs=($(seq 1 1 22))

# Loop through the array and ship off each command to bsub
for chr in "${chrs[@]}"
do
  outFolder="chr${chr}/"
  bedFile="out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.05_perind.counts.idsync.deduped.gz.qqnorm_chr${chr}.gz"
  Command="ml gsl; ml fastqtl;fastQTL --vcf ${vcfPath}${vcfFile} --bed ${bedPath}${bedFile} --cov ${covPath}${covFile} --out ${outPath}${outFolder}chr${chr}_norm_perm1000.txt.gz --permute 1000 --normal --commands 50 ${outPath}runFQTL_splitcommand_chr${chr}.sh"
  echo ${Command}
  mkdir -p ${outPath}${outFolder}
  bsub -J fQTL_${chr} -P acc_pintod02b -q express -n 50 -W 0:05 -R rusage[mem=8000] -o /sc/orga/projects/EPIASD/splicingQTL/lsf_output/fQTL_norm_perm1000_chr$chr ${Command}
done
