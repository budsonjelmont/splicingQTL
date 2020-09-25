#!/bin/bash

baseDir='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_55HCPs/'

# List the contents of the directory passed as argument 1 and read all files there into an array 
chrs=($(seq 1 1 22))

# Loop through the array and ship off each command to bsub
for chr in "${chrs[@]}"
do
  Command="ml selfsched; ml fastqtl; mpirun selfsched < ${baseDir}runFQTL_splitcommand_chr${chr}.sh"
  echo ${Command}
  bsub -J fQTL_pt2_chr${chr} -P acc_pintod02b -q premium -n 50 -W 4:00 -R rusage[mem=5000] -o /sc/arion/projects/EPIASD/splicingQTL/lsf_output/fQTL_pt2_norm_perm1000_chr$chr ${Command}
done
