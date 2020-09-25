import sys
import pandas as pd

covarsfile=sys.argv[1] # e.g. '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/covar_wasp/leafcutter-input_covar.20genoPCs.idsync.10splicingPCs.deduped.minCovars+seqPC9'
nhcps=sys.argv[2] # e.g. 5

hcpfile = str(nhcps) + 'HCPs.txt'
newcovarsfile = covarsfile[:-4] + '.' + str(nhcps) + 'HCPs_nogeno.txt'

covars = pd.read_csv(covarsfile, sep='\t')
hcps = pd.read_csv(hcpfile, sep='\t')

# Fix HCP column names after R mangles them
hcps.columns = hcps.columns.str.replace('\.', '-')

covars = pd.concat([covars, hcps], sort=False)

covars.to_csv(newcovarsfile,index=False,sep='\t')                                           
