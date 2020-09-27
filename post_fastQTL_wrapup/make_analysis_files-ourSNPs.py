# Get all rows from the FastQTL output & use the Phenotype ID to select the rows from Phenotype.bed for which QTLs were mapped. Also write another file containing just the phenotypes in a qval-significant sQTL
import sys
import pandas as pd

# New WASP-treated data
resdir = sys.argv[1]
resfile = 'chrAll_combined+qval'
rescols = ['pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval', 'slope', 'ppval', 'bpval','qval']

allphenodir = '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/'
allphenofile = 'Phenotype.bed'

# Pre-WASP data
#resdir = '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/'
#resfile = 'chrAll_combined'
#rescols = ['pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval', 'slope', 'ppval', 'bpval','qval']

#allphenodir = '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno/'
#allphenofile = 'Phenotype.bed'

allres = pd.read_csv(resdir + resfile, sep='\t')
allphenos = pd.read_csv(allphenodir + allphenofile, sep='\t', header=None)

# Join PID to phenotype file
allres.set_index('pid', drop=False, append=False, inplace=True, verify_integrity=False)
allphenos.set_index(3, drop=False, append=False, inplace=True, verify_integrity=True)

joined = allres.join(allphenos, how='inner')
joined['phenoId'] = joined.index

joined['filler'] = 'aaaaaaa'

# Reset index before dropping rows
joined.reset_index(inplace=True,drop=True)

# Sort merged data frame on FDR so we get best hit when dropping duplicate SNPs/phenotypes
joined.sort_values('qval',ascending=True,inplace=True)

# Get rows where qval < 0.05
qval05 = joined['qval']<0.05

# Export phenotype bed files (all & qval<.05 only)
joined.astype(str).to_csv(resdir + 'allphenos.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)
joined.drop_duplicates(['pid'],inplace=False).to_csv(resdir + 'allphenos_unique.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)
joined.loc[qval05].astype(str).to_csv(resdir + 'sqtlphenos.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)
joined.loc[qval05].astype(str).drop_duplicates(['pid'],inplace=False).to_csv(resdir + 'sqtlphenos_unique.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)

# Get SNP coords
joined['sid_end'] = joined['sid'].str.split(':', expand=True)[1]
joined['sid_start'] = joined.loc[~pd.isna(joined['sid_end'])]['sid_end'].astype(int)-1
joined['sid_start'] = joined['sid_start'].astype(str).str.split('.', expand = True)[0]

sid_notna = ~pd.isna(joined['sid'])

# Export sQTL SNP bed files (all & qval<.05 only)
joined.loc[sid_notna].astype(str).to_csv(resdir + 'allSNPS.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
joined.loc[sid_notna].astype(str).drop_duplicates('sid').to_csv(resdir + 'allSNPS_unique.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
joined.loc[sid_notna & qval05].astype(str).to_csv(resdir + 'sqtlSNPS.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
joined.loc[sid_notna & qval05].astype(str).drop_duplicates('sid').to_csv(resdir + 'sqtlSNPS_unique.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
# Also export sQTL SNP bed file of just best hits for each pheno (THIS is the file to use for fastQTL)
joined.loc[sid_notna & qval05].astype(str).drop_duplicates('pid').to_csv(resdir + 'besthitSNPuniquephenos.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)

# Export standardized output file:
joined['pid_chr']=joined[0]
standardizedcoldict={'sid':'sid',0:'sid_chr','sid_start':'sid_start','sid_end':'sid_end','dist':'dist','pid':'pid','pid_chr':'pid_chr',1:'pid_start',2:'pid_end','phenoId':'phenoId','slope':'slope','pval':'ppval','qval':'qval'}
joined.rename(columns=standardizedcoldict,inplace=True)
joined.loc[sid_notna & qval05].to_csv(resdir + 'qtls.txt',columns=standardizedcoldict.values(),sep='\t',index=None)

# Export sQTL SNP text file of just best hits for each pheno (THIS is the file to use for GREGOR)
joined.loc[sid_notna & qval05].astype(str).drop_duplicates('pid').drop_duplicates('sid').to_csv(resdir + 'uniqueLeadSNPuniquephenos-chrpos.txt', sep='\t', columns=['sid'], header=None, index=False)
# Also export SNP IDs of SNPs that are NOT found in any significant QTLs
joined.loc[~joined.sid.isin(joined.loc[qval05,'sid']) & sid_notna].drop_duplicates('sid').to_csv(resdir + 'uniquenonsigSNPs-chrpos.txt', sep='\t', columns=['sid'], header=None, index=False)
# Make assoc file for GREGOR
joined['qval'] = joined['qval'].map(lambda x: '%.55f' % x) # Set the # of decimal places manually
joined.rename(columns={'sid':'SNP','qval':'P'}).loc[sid_notna & qval05].astype(str).drop_duplicates('pid').drop_duplicates('SNP').to_csv(resdir + 'uniqueLeadSNPuniquephenos-chrpos.assoc', sep='\t', columns=['SNP','P'], index=False)
