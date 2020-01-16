# Get all rows from the FastQTL output & use the Phenotype ID to select the rows from Phenotype.bed for which QTLs were mapped
import pandas as pd

resdir = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/'
resfile = 'chrAll_combined'
rescols = ['pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval', 'ppval', 'bpval','qval']

allphenodir = '/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/pheno/'
allphenofile = 'Phenotype.bed'

allres = pd.read_csv(resdir + resfile, sep='\s+', header=None, names=rescols)
allphenos = pd.read_csv(allphenodir + allphenofile, sep='\t', header=None)

allres.set_index('pid', drop=False, append=False, inplace=True, verify_integrity=False)
allphenos.set_index(3, drop=False, append=False, inplace=True, verify_integrity=True)

joined = allres.join(allphenos, how='inner')
joined['phenoId'] = joined.index

joined['sid'] = 'aaaaaaa'

joined.to_csv(resdir + 'allphenos.bed', columns=[0,1,2,'phenoId','sid',5], sep='\t', index=False, header=None)
