import pandas as pd

covarsfile = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/covar/leafcutter-input_covar.20genoPCs.idsync.10splicingPCs.txt'
dropfile = '/sc/orga/projects/EPIASD/splicingQTL/PCA/SamplesToExcludeForPCA.txt'
outfile = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/covar/leafcutter-input_covar.20genoPCs.idsync.10splicingPCs.deduped.txt'

covars = pd.read_csv(covarsfile,sep='\t')

dropme = pd.read_csv('/sc/orga/projects/EPIASD/splicingQTL/PCA/SamplesToExcludeForPCA.txt',sep=' ',header=None)

#dropme.set_index(0,inplace=True,append=False,drop=False)

covars.drop(columns=dropme[0],inplace=True)

covars.to_csv(outfile,sep='\t',index=False,header=True)

#covars.set_index('id',inplace=True, drop=False, append=False, verify_integrity=True)
#dropcols = ['splicingPC','SV'
