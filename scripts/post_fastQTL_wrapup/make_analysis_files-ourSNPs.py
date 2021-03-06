# Get all rows from the FastQTL output & use the Phenotype ID to select the rows from Phenotype.bed for which QTLs were mapped. Also write another file containing just the phenotypes in a qval-significant sQTL
import sys
import pandas as pd
import numpy as np
import scipy.stats as st

# New WASP-treated data
resdir = sys.argv[1] + '/'
resfile = 'chrAll_combined+qval'
rescols = ['pid', 'nvar', 'shape1', 'shape2', 'dummy', 'sid', 'dist', 'npval', 'slope', 'ppval', 'bpval','qval']

hg38bedfile = sys.argv[2] # e.g. '/sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/hg38_liftover/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.hg38.bed'
frqfile = sys.argv[3] # e.g. '/sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.frq'
allsnpsfile = sys.argv[4] # e.g. '/sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.frq'

try:
  hg38bed=pd.read_csv(hg38bedfile,sep='\t',header=None)
  liftover=True
except FileNotFoundError:
  liftover=False

try:
  frq=pd.read_csv(frqfile,sep='\t',header=None)
  addFrq=True
except FileNotFoundError:
  addFrq=False

allphenodir = '/sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp/'
allphenofile = 'Phenotype.bed_noChr'

allres = pd.read_csv(resdir + resfile, sep='\t')
allphenos = pd.read_csv(allphenodir + allphenofile, sep='\t', header=None)

# Join PID to phenotype file
allres.set_index('pid', drop=False, append=False, inplace=True, verify_integrity=False)
allphenos.set_index(3, drop=False, append=False, inplace=True, verify_integrity=True)

joined = allres.join(allphenos, how='inner')
joined['phenoId'] = joined.index

# Manual wrangling
joined['filler'] = 'aaaaaaa'
# Add sample size for coloc
joined['N'] = 1301

# Reset index before dropping rows
joined.reset_index(inplace=True,drop=True)

# Sort merged data frame on FDR so we get best hit when dropping duplicate SNPs/phenotypes
joined.sort_values('qval',ascending=True,inplace=True)

# If a .BED file of hg38 coordinates was passed, merge the bed file to the data frame to get the hg38 coordinates 
if liftover:
  joined.set_index('sid',drop=False,append=False,inplace=True,verify_integrity=False)
  hg38bed.set_index(3,drop=False,append=False,inplace=True,verify_integrity=True)
  joined['sid_chr_hg38'] = hg38bed[0]
  joined['sid_start_hg38'] = hg38bed[1]
  joined['sid_end_hg38'] = hg38bed[2]
  joined['sid_start_hg38'] = joined['sid_start_hg38'].astype('Int64')
  joined['sid_end_hg38'] = joined['sid_end_hg38'].astype('Int64')
  joined['sid_hg38'] = joined.apply(lambda x: np.nan if pd.isna(x['sid_chr_hg38']) else x['sid_chr_hg38'][3:] + ':' + str(x['sid_end_hg38']), axis=1)
  joined.reset_index(drop=True,inplace=True)

# Get hg19 SNP coords
joined['sid_end'] = joined['sid'].str.split(':', expand=True)[1]
joined['sid_start'] = joined.loc[~pd.isna(joined['sid_end'])]['sid_end'].astype(int)-1
joined['sid_start'] = joined['sid_start'].astype(str).str.split('.', expand = True)[0]

sid_notna = ~pd.isna(joined['sid'])

# Get rows where qval < 0.05
qval05 = joined['qval']<0.05

# Export phenotype bed files (all & qval<.05 only)
joined.astype(str).to_csv(resdir + 'allphenos.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)
joined.drop_duplicates(['pid'],inplace=False).to_csv(resdir + 'allphenos_unique.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)
joined.loc[qval05].astype(str).to_csv(resdir + 'sqtlphenos.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)
joined.loc[qval05].astype(str).drop_duplicates(['pid'],inplace=False).to_csv(resdir + 'sqtlphenos_unique.bed', columns=[0,1,2,'phenoId','filler',5], sep='\t', index=False, header=None)

# Export sQTL SNP bed files (all & qval<.05 only)
joined.loc[sid_notna].astype(str).to_csv(resdir + 'allSNPS.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
joined.loc[sid_notna].astype(str).drop_duplicates('sid').to_csv(resdir + 'allSNPS_unique.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
joined.loc[sid_notna & qval05].astype(str).to_csv(resdir + 'sqtlSNPS.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
joined.loc[sid_notna & qval05].astype(str).drop_duplicates('sid').to_csv(resdir + 'sqtlSNPS_unique.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
# Also export sQTL SNP bed file of just best hits for each pheno (THIS is the file to use for fastQTL)
joined.loc[sid_notna & qval05].astype(str).drop_duplicates('pid').drop_duplicates('sid').to_csv(resdir + 'besthitSNPuniquephenos.bed', sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)

# Export standardized output file:
joined['pid_chr']=joined[0]

# Column renamer
standardizedcoldict={'sid':'sid',0:'sid_chr','sid_start':'sid_start','sid_end':'sid_end','dist':'dist','pid':'pid','pid_chr':'pid_chr',1:'pid_start',2:'pid_end','phenoId':'phenoId','slope':'slope','bpval':'bpval','qval':'qval'} 
outfilename='qtls.txt'

joined.rename(columns=standardizedcoldict,inplace=True)
joined.loc[sid_notna & qval05].to_csv(resdir + outfilename,columns=standardizedcoldict.values(),sep='\t',index=None)

if liftover: 
  hg38outfilename='qtls+hg38.txt'
  hg38cols={'sid_hg38':'sid_hg38','sid_chr_hg38':'sid_chr_hg38','sid_start_hg38':'sid_start_hg38','sid_end_hg38':'sid_end_hg38'}
  standardizedcoldict.update(hg38cols)
  joined.rename(columns=standardizedcoldict,inplace=True)
  joined.loc[sid_notna & qval05].to_csv(resdir + hg38outfilename,columns=standardizedcoldict.values(),sep='\t',index=None)

# Make moloc input file
# Calculate SE from beta & pvalue. See https://github.com/stephenslab/gtexresults/commit/d1a14f3792844b121d0142c1d5c438ed91dfb258
def get_se(bhat, pval):
  if bhat < 0:
   z = st.norm.ppf(pval / 2)
  else:
   z = st.norm.ppf(1 - pval / 2)
  if z != 0:
    se = bhat/z
  else:
    se = np.nan
  return se

joined.loc[joined['bpval']==0,'bpval']=0.000000000000001
joined['se']=joined.apply(lambda x: get_se(x.slope,x.bpval),axis=1)

# Moloc
molocrenamerdict={'pid':'ProbeID','chr':'CHR','sid_chr_hg38':'CHR_hg38','sid_end':'POS','sid_end_hg38':'POS_hg38','slope':'BETA','bpval':'PVAL','se':'SE'}

joined.drop_duplicates('sid').rename(columns=molocrenamerdict).to_csv(resdir + 'allSNPs.moloc.dat', sep='\t', index=False, columns=molocrenamerdict.values())

# Coloc
colocrenamerdict={'sid':'snp','sid_hg38':'snp_hg38','chr':'CHR','sid_chr_hg38':'CHR_hg38','sid_end':'POS','sid_end_hg38':'POS_hg38','slope':'beta','bpval':'pvalues','se':'varbeta','N':'N'}

joined.drop_duplicates('sid').rename(columns=colocrenamerdict).to_csv(resdir + 'allSNPs_unique.coloc.dat', sep='\t', index=False, columns=colocrenamerdict.values())

# Export sQTL SNP text file of just best hits for each pheno (THIS is the file to use for GREGOR)
joined.loc[sid_notna & qval05].astype(str).drop_duplicates('pid').drop_duplicates('sid').to_csv(resdir + 'uniqueLeadSNPuniquephenos-chrpos.txt', sep='\t', columns=['sid'], header=None, index=False)
# Also export SNP IDs of SNPs that are NOT found in any significant QTLs
joined.loc[~joined.sid.isin(joined.loc[qval05,'sid']) & sid_notna].drop_duplicates('sid').to_csv(resdir + 'uniquenonsigSNPs-chrpos.txt', sep='\t', columns=['sid'], header=None, index=False)
# Make assoc file for GREGOR
joined['qval'] = joined['qval'].map(lambda x: '%.55f' % x) # Set the # of decimal places manually
joined.rename(columns={'sid':'SNP','qval':'P'}).loc[sid_notna & qval05].astype(str).drop_duplicates('pid').drop_duplicates('SNP').to_csv(resdir + 'uniqueLeadSNPuniquephenos-chrpos.assoc', sep='\t', columns=['SNP','P'], index=False)
