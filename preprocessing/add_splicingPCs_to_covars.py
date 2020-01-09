import pandas as pd

covarsfile = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/covar/leafcutter-input_covar.20genoPCs.idsync'
sqtlpcsfile = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/pheno/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.deduped.gz.PCs'
newcovarsfile = covarsfile + '.10splicingPCs.txt'

covars = pd.read_csv(covarsfile + '.txt', sep='\t')
sqtlpcs = pd.read_csv(sqtlpcsfile, sep='\t')

sqtlpcs['id'] = 'splicingPC' + sqtlpcs['id'].astype(str)

covars = pd.concat([covars, sqtlpcs], sort=False)

covars.to_csv(newcovarsfile,index=False,sep='\t')
