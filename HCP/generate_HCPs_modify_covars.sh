covarsfile=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/covar_wasp/leafcutter-input_covar_WASP.20genoPCs.idsync.deduped.minCovars+seqPC9+20genoPCs.txt
countsfile=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.deduped.gz.qqnorm_allCombined
hcpbasedir=/sc/arion/projects/EPIASD/splicingQTL/HCP/
hcpscript=do_hcp_outputHCPs_minCovars-noGenoPCs.R
appendcovarscript=add_HCPs_to_covars_minCovars-noGenoPCs.py
nhcps=55

Rscript $hcpbasedir$hcpscript $countsfile $covarsfile $nhcps
# Add HCs to covariate data frame (this needs to be done separately because R mangles the column names of the covariates)
python3 $hcpbasedir$appendcovarscript $covarsfile $nhcps
