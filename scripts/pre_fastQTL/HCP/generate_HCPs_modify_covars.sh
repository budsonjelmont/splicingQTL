covarsfile=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/covar_wasp/leafcutter-input_covar_WASP.20genoPCs.idsync.deduped.minCovars+seqPC9+20genoPCs.txt
countsfile=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.deduped.gz.qqnorm_allCombined
hcpbasedir=/sc/arion/projects/EPIASD/splicingQTL/HCP/
hcpscript=do_hcp_outputHCPs_minCovars-noGenoPCs.R
appendcovarscript=add_HCPs_to_covars_minCovars-noGenoPCs.py

# Parameters to iterate 
hcps=($(seq 0 5 100))

# Run HCP to derive latent variables, then add HCs to covariate data frame (this needs to be done separately in 2 steps because R mangles the column names of the covariates)
for nhcps in "${hcps[@]}"
do
  echo "$nhcps ..."
  # Run GREGOR on newly created .conf file
  command="ml python/3.7.3; ml R; Rscript $hcpbasedir$hcpscript $countsfile $covarsfile $nhcps; python3 $hcpbasedir$appendcovarscript $covarsfile $nhcps"
  echo $command
  submitjob 1 -P acc_pintod02b -q premium -c 1 -n 1 -m 10 -J make_${nhcps}HCPs_covar ${command}
done
