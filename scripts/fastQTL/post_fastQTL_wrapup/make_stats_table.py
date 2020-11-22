import pandas as pd

basepath='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/4genoPCs_nogenoInHCP/'
basename='deduped_mincovars+seqPC79_'

summarydf = pd.DataFrame(columns=['n_hcps','n_qtl_introns','n_introns','n_clusts','n_snps','n_sgenes','n_phenogenes'])

for n in [x for x in range(0,105,5)]:
  resdir=basepath + basename + str(n) + 'HCPs/'
  print(resdir)
  res=pd.read_csv(resdir + 'qtls+pid_ensg+sid_ensg.txt',sep='\t')
  summarydf = summarydf.append(pd.Series(
    {'n_hcps':n,
     'n_qtl_introns':res.drop_duplicates(['sid','pid']).shape[0],
     'n_introns':res['pid'].nunique(),
     'n_clusts':res['pid'].str.split(':',expand=True)[3].nunique(),
     'n_snps':res['sid'].nunique(),
     'n_sgenes':res['sid_ensg'].nunique(),
     'n_phenogenes':res['pid_ensg'].nunique()
    }),ignore_index=True)

summarydf.sort_values('n_hcps',ascending=False).to_csv(basepath + 'HCP_titration_summary_seqPC79.tsv',sep='\t',index=False)
