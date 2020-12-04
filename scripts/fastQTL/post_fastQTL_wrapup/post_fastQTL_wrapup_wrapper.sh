# Results directory base
resdirbase=/sc/arion/projects/EPIASD/splicingQTL/output/sqtls/permute/minCovars+seqPC79/4genoPCs
#isnominal=TRUE # TRUE/FALSE--not used currently, calculateFQTLstats.R counts # columns to determine if results are from permutation or nominal pass
#fdrmethod=st # {'st','bh'} e.g storey
annotdir=/sc/arion/projects/EPIASD/splicingQTL/ext_data/other/gencode/hg19/v33lift37
annotbed=gencode.v33lift37.annotation.bed
annotidmap=gencode_v33_idmap.txt 

# Parameters to iterate 
hcps=($(seq 0 5 100))
declare -a fdr
fdr=( st bh )

# Loop through the array and ship off each command to bsub
for fdrmethod in "${fdr[@]}"
do
  for nhcps in "${hcps[@]}"
  do
    echo "$nhcps ..."
    resdir=$resdirbase/${nhcps}HCPs
    command="ml python/3.7.3; ./post_fastQTL_wrapup.sh $resdir $fdrmethod $annotdir $annotbed $annotidmap"
    echo $command
    submitjob 8 -P acc_pintod02b -q premium -c 1 -n 1 -m 60 -J post_fastQTL_wrapup_${nhcps}HCPs ${command} 
  done
done
