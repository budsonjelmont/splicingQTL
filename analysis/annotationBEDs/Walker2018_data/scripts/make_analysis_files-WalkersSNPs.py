# Make fastQTL input files from Walker 2018 data 
import pandas as pd

# New WASP-treated data
resdir = '/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/sQTLs/'
resfile = 'WalkerCell_significantSQTLS.txt'

res = pd.read_csv(resdir+resfile,sep='\t')

# Addtl cols that need to be added manually
res['sid_start']=res['pos']-1

res['sid_chr']=res['chr']
res['pid_chr']=res['chr']
res['score_dummy']='0'
res['strand_dummy']='+'

# Additional processing for existing columns
res['variant_id']=res['sid_chr'].astype(str)+':'+res['pos'].astype(str)

# Standardize columns
standardizedcoldict={
  'variant_id':'sid',
  'sid_chr':'sid_chr',
  'sid_start':'sid_start',
  'pos':'sid_end',
  'tss_distance':'dist',
  'intron_id':'pid',
  'pid_chr':'pid_chr',
  'intron_start':'pid_start',
  'intron_end':'pid_end',
  'slope':'slope',
  'pval':'pval_beta',
'qval':'qval'
}
res.rename(columns=standardizedcoldict,inplace=True)

# Reset index before dropping rows
# Sort merged data frame on FDR so we get best hit when dropping duplicate SNPs/phenotypes
res.sort_values('qval',ascending=True,inplace=True)

# Get rows where qval < 0.05
qval05 = res['qval']<0.05

# Export phenotype bed files (all & qval<.05 only)
phenobedcols=['pid_chr','pid_start','pid_end','pid','score_dummy','strand_dummy']
res.astype(str).to_csv(resdir + 'allphenos.bed', columns=phenobedcols, sep='\t', index=False, header=None)
res.drop_duplicates(['pid'],inplace=False).to_csv(resdir + 'allphenos_unique.bed', columns=phenobedcols, sep='\t', index=False, header=None)
res.loc[qval05].astype(str).to_csv(resdir + 'sqtlphenos.bed', columns=phenobedcols, sep='\t', index=False, header=None)
res.loc[qval05].astype(str).drop_duplicates(['pid'],inplace=False).to_csv(resdir + 'sqtlphenos_unique.bed', columns=phenobedcols, sep='\t', index=False, header=None)

# Get SNP coords
sid_notna = ~pd.isna(res['sid'])

# Export sQTL SNP bed files (all & qval<.05 only)
snpbedcols=['sid_chr','sid_start','sid_end','sid','score_dummy','strand_dummy']
res.loc[sid_notna].astype(str).to_csv(resdir + 'allSNPS.bed', sep='\t', columns=snpbedcols, header=None, index=False)
res.loc[sid_notna].astype(str).drop_duplicates('sid').to_csv(resdir + 'allSNPS_unique.bed', sep='\t', columns=snpbedcols, header=None, index=False)
res.loc[sid_notna & qval05].astype(str).to_csv(resdir + 'sqtlSNPS.bed', sep='\t', columns=snpbedcols, header=None, index=False)
res.loc[sid_notna & qval05].astype(str).drop_duplicates('sid').to_csv(resdir + 'sqtlSNPS_unique.bed', sep='\t', columns=snpbedcols, header=None, index=False)
# Export sQTL SNP bed file of just best hits for each pheno (THIS is the file to use for fastQTL)
res.loc[sid_notna & qval05].astype(str).drop_duplicates('pid').drop_duplicates('pid').to_csv(resdir + 'besthitSNPuniquephenos.bed', sep='\t', columns=snpbedcols, header=None, index=False)

# Also export sQTL SNP text file of just best hits for each pheno (THIS is the file to use for GREGOR)
res.loc[sid_notna & qval05].astype(str).drop_duplicates('pid').drop_duplicates('sid').to_csv(resdir + 'uniqueLeadSNPuniquephenos-chrpos.txt', sep='\t', columns=['sid'], header=None, index=False)
# Also export sQTL SNP text file of all SNPs 
res.loc[sid_notna & qval05].astype(str).drop_duplicates('sid').to_csv(resdir + 'uniqueSNPs-chrpos.txt', sep='\t', columns=['sid'], header=None, index=False)

# Make assoc file for GREGOR
res['qval'] = res['qval'].map(lambda x: '%.55f' % x) # Set the # of decimal places manually
res.loc[sid_notna & qval05].astype(str).rename(columns={'sid':'SNP','qval':'P'}).drop_duplicates('pid').drop_duplicates('SNP').to_csv(resdir + 'uniqueLeadSNPuniquephenos-chrpos.assoc', sep='\t', columns=['SNP','P'], index=False)

# Export standardized output file:
res.loc[sid_notna & qval05].to_csv(resdir + 'qtls.txt',columns=standardizedcoldict.values(),sep='\t',index=None)
