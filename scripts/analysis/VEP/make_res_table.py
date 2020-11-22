# Read the tab-delimited output files from VEP into a dataframe and generate a table of counts per region
import numpy as np
import pandas as pd

resdir='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_finished/'

resfiles = {'our_sSNPs':resdir+'our_sQTLSNPs_offlineRun',
  'Walker_sSNPs':resdir+'Walker_sQTLSNPs_offlineRun',
  'Raj_sSNPs':resdir+'Raj_sQTLSNPs_offlineRun',
  'PEC_eSNPs':resdir+'PEC_eQTLSNPs_offlineRun',
  'PEC_isoSNPs':resdir+'PEC_isoQTLSNPs_offlineRun',
  'our_nonQTLSNPs':resdir+'our_nonQTLSNPs_offlineRun'
}

# Column # containing the annotation to summarize
categorycol=6

resdf = {}
for name,file in resfiles.items():
  resdf[name]=pd.read_csv(file,sep='\t',comment='#',header=None)
  resdf[name]['dataset']=name

# Convert list of df -> df
combined=pd.concat(resdf,axis=0,ignore_index=True)

# Rename category column
combined.rename(columns={categorycol:'annotation'},inplace=True)

# Make summary table of event counts & add final row of totals
summary_table = combined.groupby(['annotation','dataset']).size().unstack()
summary_table.loc['total'] = summary_table.sum(skipna=True)

summary_table.to_csv(resdir+'QTLSNPs.VEPregionAnnotation.summary',sep='\t')
