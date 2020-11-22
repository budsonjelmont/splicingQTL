#!/bin/bash
# Run this script after batchsubber_fastqtl.sh  

resdirbase=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp_nominal_normal/minCovars+seqPC79/20genoPCs/

# Parameters to iterate 
hcps=($(seq 0 5 100))
chrs=($(seq 1 1 22))

# Loop through the array and ship off each command to bsub
for hcps in "${hcps[@]}"
do
  resdir=$resdirbase/${hcps}HCPs/
  for chr in "${chrs[@]}"
  do
    Command="ml selfsched; ml fastqtl; mpirun selfsched < ${resdir}runFQTL_splitcommand_chr${chr}.sh"
    echo ${Command}
    bsub -J fQTL_pt2_chr${chr} -P acc_pintod02b -q premium -n 50 -W 0:50 -R rusage[mem=5000] -o /sc/arion/projects/EPIASD/splicingQTL/lsf_output/fQTL_pt2_norm_perm1000_chr$chr ${Command}
  done
done
