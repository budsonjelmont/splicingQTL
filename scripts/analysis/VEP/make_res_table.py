# Read the tab-delimited output files from VEP into a dataframe and generate a table of counts per region
import numpy as np
import pandas as pd

resdir='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/'

resfiles = {'our_sSNPs':resdir+'OurStudy_sQTLSNPs',
  'our_nonQTLSNPs':resdir+'OurStudy_allnonsQTLSNPs_withinCisWindow',
  'Walker_sSNPs':resdir+'Walker_sQTLSNPs',
  'Raj_sSNPs':resdir+'Raj_sQTLSNPs',
  'PEC_eSNPs':resdir+'PEC_eQTLSNPs',
  'PEC_noneSNPs':resdir+'PEC_allnoneQTLSNPs_withinCisWindow',
  'PEC_isoSNPs':resdir+'PEC_isoQTLSNPs',
  'PEC_nonisoSNPs':resdir+'PEC_allnonisoQTLSNPs_withinCisWindow',
  'PEC_tSNPs':resdir+'PEC_tQTLSNPs',
  'PEC_nontSNPs':resdir+'PEC_allnontQTLSNPs_withinCisWindow'
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
