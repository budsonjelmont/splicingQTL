#!/bin/bash
# Iterate over n_hcps and call fastqtl to generate commands to run the QTL detection into 50 parallel jobs.

covpath=/sc/arion/projects/EPIASD/splicingQTL/output/covar_wasp/
covfilebase=leafcutter-input_covar_WASP.20genoPCs.idsync.deduped.minCovars+seqPC9+4genoPCs.

vcfpath=/sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/
vcffile=Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.vcf.gz

countspath=/sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp/

outpathbase=/sc/arion/projects/EPIASD/splicingQTL/output/permute/minCovars+seqPC9/4genoPCs/

permute=false
normal=true

# Parameters to iterate 
hcps=($(seq 0 5 100))
chrs=($(seq 1 1 22))

if [ $permute == 'true' ]
then
  permuteflag=" --permute 1000 10000 "
elif [ $permute == 'false' ]
then
  permuteflag=""
else
  echo "You didn't specify the permute (true/false) flag"
  exit 1
fi

if [ $normal == 'true' ]
then
  normalflag=" --normal "
elif [ $normal == 'false' ]
then
  normalflag=""
else
  echo "You didn't specify the normal (true/false) flag"
  exit 1
fi

# Loop through the array and ship off each command to bsub
for hcp in "${hcps[@]}"
do
  covfile=$covfilebase${hcp}HCPs_nogeno.txt
  outpath=$outpathbase${hcp}HCPs/
  for chr in "${chrs[@]}"
  do
    outFolder="chr${chr}/"
    countsfile="out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.deduped.gz.qqnorm_chr${chr}.gz"
    Command="ml gsl; ml fastqtl;fastQTL --vcf ${vcfpath}${vcffile} --bed ${countspath}${countsfile} --cov ${covpath}${covfile} --out ${outpath}${outFolder}chr${chr}_norm_perm1000.txt.gz${normalflag}${permuteflag}--window 1e5 --commands 50 ${outpath}runFQTL_splitcommand_chr${chr}.sh"
    echo ${Command}
    mkdir -p ${outpath}${outFolder}
    bsub -J fQTL_${chr} -P acc_pintod02b -q express -n 50 -W 0:05 -R rusage[mem=8000] -o /sc/arion/projects/EPIASD/splicingQTL/lsf_output/fQTL_norm_perm1000_chr$chr ${Command}
  done
done
