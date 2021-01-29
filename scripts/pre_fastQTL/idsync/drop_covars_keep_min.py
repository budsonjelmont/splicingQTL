import pandas as pd

covarsfile = '/sc/arion/projects/EPIASD/splicingQTL/output/covars/leafcutter-input_covar_star_WASP.20genoPCs.idsync.deduped.txt'
ngenopcs=10 # n genotype PCs to retain in the output
set='minCovars+seqPC79' # Covariate set in the dictionary below to be retained in the output

covarsets={'minCovars':['sex','tissue','RIN','study'],'minCovars+seqPC9':['sex','tissue','RIN','study','seqPC9'],'minCovars+seqPC79':['sex','tissue','RIN','study','seqPC7','seqPC9'],'gt05var':['seqPC9','seqPC11','seqPC7','seqPC8','seqPC1','seqPC2','seqPC4','seqPC6','RIN','genotypePC1','sex','tissue','study']}
if ngenopcs==0:
  outfile = covarsfile[:-4]+'.'+set+'.txt'
else:
  outfile = covarsfile[:-4]+'.'+set+'+'+str(ngenopcs)+'genoPCs.txt'


covars = pd.read_csv(covarsfile,sep='\t')

covars.set_index('id', inplace=True, drop=False, append=False, verify_integrity=True)

tokeep = covarsets[set] 
tokeep = tokeep + ['genotypePC' + str(x) for x in range(1,(ngenopcs+1))]

covars = covars.loc[tokeep]

covars.to_csv(outfile,sep='\t',index=False,header=True)
