# Results directory base
resdirbase=/sc/arion/projects/EPIASD/splicingQTL/output/sqtls/minCovars+seqPC9/4genoPCs
#fdrmethod=st # {'st','bh'} e.g storey
annotdir=/sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp/
annotfile=out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.clustids.gcv33ids
genedir=/sc/arion/projects/EPIASD/splicingQTL/ext_data/other/gencode/hg19/v33lift37/
genebed=gencode.v33lift37.annotation.geneModel.bed

# Parameters to iterate 
hcps=($(seq 0 5 100))
declare -a fdr
fdr=( st bh ) # Not used

# Loop through the array and ship off each command to bsub
#for fdrmethod in "${fdr[@]}" # Currently does nothing
#do
  for nhcps in "${hcps[@]}"
  do
    echo "$nhcps ..."
    resdir=$resdirbase/${nhcps}HCPs
    command="ml python/3.7.3; ./post_fastQTL_wrapup.sh $resdir $annotdir $annotfile $genedir $genebed"
    echo $command
    submitjob 8 -P acc_pintod02b -q premium -c 1 -n 1 -m 60 -J post_fastQTL_wrapup_${nhcps}HCPs ${command}
  done
#done
