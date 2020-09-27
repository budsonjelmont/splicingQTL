library(data.table)
library(ggplot2)

basepath='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/'
basename='deduped_mincovars+seqPC9_'


summarydf = pd.DataFrame(columns=['n_hcps','n_qtls','n_introns','n_clusts','n_snps','n_sgenes','n_phenogenes'])

for n in [x for x in range(0,55,5)]:
  resdir=basename + str(n) + 'HCPs/'
  res=pd.read_csv(resdir + 'qtls+pid_ensg+pid_enst.txt')
  summarydf.append(pd.Series(
    {'n_hcps':n,
     'n_qtls',
     'n_introns','n_clusts','n_snps','n_sgenes','n_phenogenes'}))
