import pandas as pd
import sys

covarsfile = sys.argv[1] # '/sc/orga/projects/pintod02b/capstone/leafcutter/data-freeze4/leafcutter-input_covar.txt'
metadatfile = sys.argv[2] # '/sc/orga/projects/EPIASD/splicingQTL/meta_matchedIDs.csv'
pcafile = sys.argv[3] # '/sc/orga/projects/EPIASD/splicingQTL/PCA/plinkPCAresults/Capstone4.sel.idsync.2allele.maf5.mind1.geno1.hwe1e-5.deduped.highLDexcl.indep1500_150_.2.eigenvec'
relationfile = sys.argv[4] # '/sc/orga/projects/EPIASD/splicingQTL/PCA/Capstone4.sel.idsync.maf5.mind1.geno1.hwe1e-5.vcf.relatedness.genome'
dropfile = sys.argv[5] # '/sc/orga/projects/EPIASD/splicingQTL/PCA/SamplesToExcludeForPCA.txt'
outfile = sys.argv[6]

#countsfile = 'pheno_sampleIDs.txt'

covars = pd.read_csv(covarsfile,sep='\t')
covars.set_index('SampleID',inplace=True, drop=False, append=False, verify_integrity=True)

dropped = True

try:
  reln = pd.read_csv(relationfile, sep='\s+')
  drop = pd.read_csv(dropfile,header=None, sep='\s+')
except pd.errors.EmptyDataError:
  dropped = False

# Read mapping files
idmap = pd.read_csv(metadatfile)

# Process pca file w/ new sample names (& duplicate rows as needed)
pca = pd.read_csv(pcafile, sep='\s+')
#pca = pca.loc[pca['sample_id'].isin(idmap['individualID'])]

pca.set_index('FID', inplace=True, drop=False, append=False, verify_integrity=True)

#####################################################################
# If pca output uses the pre-idsync'ed IDs, fix them here
#idmap.set_index('individualID',inplace=True, drop=False, append=False, verify_integrity=False)

#for i in pca.index:
#  for j in idmap['PSIcountID'].loc[[i]]:
#    pca.loc[j] = pca.loc[i]

# Now drop the rows that don't match a PSI count ID (which should be all and only the original rows in the file
#pca = pca.loc[idmap['PSIcountID']]
#####################################################################

# Add rows to Peddy output for dropped samples
if dropped:
  reln = reln.loc[reln['PI_HAT'] > 0.2]
  relnlookup = reln[['FID1','FID2']].append(reln[['FID2','FID1']].rename(columns={'FID2':'FID1','FID1':'FID2'}))
  relnlookup.set_index('FID1',inplace=True,drop=False,append=False,verify_integrity=False)
  drop_rel = drop.merge(reln[['FID1','FID2']].rename(columns={'FID2':'relative'}),left_on=0,right_on='FID1',how='left')
  drop_rel = drop_rel.merge(reln[['FID1','FID2']].rename(columns={'FID1':'relative_2'}),left_on=0,right_on='FID2',how='left')
  drop_rel['relative'].fillna(drop_rel['relative_2'],inplace=True)
  drop_rel_PCA = drop_rel.merge(pca[['PC' + str(pc) for pc in range(1,21)]], left_on='relative',right_index=True)
  drop_rel_PCA = drop_rel_PCA.loc[~drop_rel_PCA[0].duplicated()]
  drop_rel_PCA.set_index(0,drop=False,inplace=True,verify_integrity=False)
  pca = pca.append(drop_rel_PCA.rename(columns={0:'FID'}))

# Add ancestry PCs from pca output
npcs = 20
covars[['genotypePC' + str(pc) for pc in range(1,(npcs+1))]] = pca[['PC' + str(pc) for pc in range(1,(npcs+1))]]

#covars['genotypePC1'] = pca['PC1']
#covars['genotypePC2'] = pca['PC2']
#covars['genotypePC3'] = pca['PC3']
#covars['genotypePC4'] = pca['PC4']

covars = covars.loc[idmap['PSIcountID']]

# Rename SampleID column
covars.rename(columns={'SampleID':'id'},inplace=True)

# Transpose covariate dataframe
covars_t = covars.transpose()

# Drop #Filename row and then insert id column
covars_t.drop(index='#Filename',inplace=True)
covars_t.insert(0,'id',covars_t.index)

# Write out all but 1st column
covars_t.to_csv(outfile,sep='\t',index=False,header=False)
