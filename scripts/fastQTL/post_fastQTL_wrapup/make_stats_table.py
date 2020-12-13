import argparse as ap
import pandas as pd
import os

parser = ap.ArgumentParser(description='Iterate over the results folders from the fastQTL HCP titration and create a table summarizing the yields from each number of hidden covariates.')
parser.add_argument('basepath', metavar='basepath', type=str, nargs=1, help='Path directory containing the results folders from the fastQTL HCP titration')
parser.add_argument('subdir', metavar='subdir', type=str, nargs=1, help='Directory inside each *HCPs folder containing the post-fastQTL wrapup results. A directory with this name will be created in the base path to store the output of this script')
args = parser.parse_args()

basepath=args.basepath[0]
subdir=args.subdir[0]

# Make results folder
outdir=basepath+'/'+subdir
os.makedirs(outdir, exist_ok=True)

summarydf = pd.DataFrame(columns=['n_hcps','n_qtl_introns','n_introns','n_clusts','n_snps','n_sgenes','n_phenogenes'])

for n in [x for x in range(0,105,5)]:
  try:
    resdir=basepath + '/' + str(n) + 'HCPs/' + subdir + '/'
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
  except FileNotFoundError as not_found:
    print('No qtls+pid_ensg+sid_ensg.txt file found in ' + basepath + '/' + str(n) + 'HCPs. Skipping...')

summarydf.sort_values('n_hcps',ascending=False).to_csv(outdir + '/HCP_titration_summary.tsv',sep='\t',index=False)
