# Results directory base
resdirBase=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/4genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_

# Parameters to iterate 
hcps=($(seq 0 5 100))

# Loop through the array and ship off each command to bsub
for nhcps in "${hcps[@]}"
do
  echo "$nhcps ..."
  resdir=$resdirBase${nhcps}HCPs
  command="ml python/3.7.3; ./post_fastQTL_wrapup.sh $resdir"
  echo $command
  submitjob 1 -P acc_pintod02b -q premium -c 1 -n 1 -m 10 -J post_fastQTL_wrapup_${nhcps}HCPs ${command} 
done
