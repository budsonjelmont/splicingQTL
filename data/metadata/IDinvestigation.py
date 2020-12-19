# Script to determine which PEC samples have both genotype and phenotype data available and can be used for sQTL analysis. Reads in genotype and phenotype IDs from the VCF file and leafcutter output respectively. Merges all the 

import pandas as pd

geno = pd.read_csv('genotype_ids.txt',sep='\t',header=None) # All sample IDs in the Capstone4 VCF 
pheno = pd.read_csv('phenotype_ids.txt',sep='\t',header=None) # All sample IDs in the learfcutter counts file

# DP's merged metadata file
meta = pd.read_csv('metaData_RNAseqFreeze1+2.csv')
# PEC metadata file
PECmeta = pd.read_csv('PEC Capstone Collection_GenotypesMetadata_fromSynapse.tsv',sep='\t')
# DP's fixed PSIcounts file
fixedpsi = pd.read_csv('psicountsfile_sampleIDs-fixed.txt',sep='\t')

geno.set_index(0,drop=False,append=False,inplace=True,verify_integrity=True)
pheno.set_index(0,drop=False,append=False,inplace=True,verify_integrity=True)

genoids = geno[0]
phenoids = pheno[0]

genoids.index.name='genoids'
phenoids.index.name='phenoids'

fixedphenoids = fixedpsi['PSIcountID_fixed']
fixedphenoids = fixedphenoids.loc[~pd.isna(fixedphenoids)]

##################################################################################################################
# Rows where unnamed column A != specimenID
# meta_diffColASpecID = meta.loc[meta['Unnamed: 0']!=meta['specimenID']]
#
# # Fix study name
# meta['study'].replace({'LIBD_szControl':'LIBD__szControl'},inplace=True)
#
# # VCF ids that match the first metadata column
# genoids.loc[genoids.isin(meta['Unnamed: 0'])]
# # VCF ids that match the metadata 'specimenID' column
# genoids.loc[genoids.isin(meta['specimenID'])]
# # VCF ids that match the metadata 'individualID' column
# genoids.loc[genoids.isin(meta['individualID'])]
# # VCF ids that match the PEC files' 'genotypingID' column
# genoids.loc[genoids.isin(PECmeta['genotypingID'])]
#
#
#
# # make concat columns to match to phenotype file
# meta['study_specimen'] = meta['study'] + '_' + meta['specimenID']
# meta['study_colA'] = meta['study'] + '_' + meta['Unnamed: 0']
# meta['study_indiv'] = meta['study'] + '_' + meta['individualID']
#
# # Genotype - Metadata matching
# # Match colA to genotype file
# meta['colAinGeno'] = meta['Unnamed: 0'].isin(genoids)
# # Same for specimenID
# meta['specIDinGeno'] = meta['specimenID'].isin(genoids)
# # Same for individual ID
# meta['indivIDinGeno'] = meta['specimenID'].isin(genoids)
#
# # Phenotype - Metadata matching
# # Match colA to genotype file
# meta['study_colAinPheno'] = meta['study_colA'].isin(phenoids)
# # Same for specimenID
# meta['study_specIDinPheno'] = meta['study_specimen'].isin(phenoids)
# # Same for individual ID
# meta['study_indivIDinPheno'] = meta['study_indiv'].isin(phenoids)
#
# # Phenotype - Metadata matching w/ trimmed pheno IDs
# # Match colA to genotype file
# meta['colAinFixedPheno'] = meta['Unnamed: 0'].isin(fixedphenoids)
# # Same for specimenID
# meta['specIDinFixedPheno'] = meta['specimenID'].isin(fixedphenoids)
# # Same for individual ID
# meta['indivIDinFixedPheno'] = meta['individualID'].isin(fixedphenoids)
#
#
# # Total rows that have a match in both the Genotype and Phenotype files
# meta.query('(colAinGeno or specIDinGeno or indivIDinGeno) and (colAinPheno or specIDinPheno or indivIDinPheno)')
#
# # Total rows that have a match in both the Genotype and FIXED Phenotype files
# meta.query('(colAinGeno or specIDinGeno or indivIDinGeno) and (colAinFixedPheno or specIDinFixedPheno or indivIDinFixedPheno)')
#
# meta.to_csv('metaData_phase12_jmb_columnsadded.csv')
#
# ########################################################################################################
# # Do the same as above but for PECmeta
# # make concat columns to match to phenotype file
# PECmeta['study_geno'] = PECmeta['study'] + '_' + PECmeta['genotypingID']
# PECmeta['study_indiv'] = PECmeta['study'] + '_' + PECmeta['individualID']
#
# # Genotype - PECmetadata matching
# # Match colA to genotype file
# PECmeta['genoIDinGeno'] = PECmeta['genotypingID'].isin(genoids)
# # Same for individual ID
# PECmeta['indivIDinGeno'] = PECmeta['individualID'].isin(genoids)
#
# # Phenotype - PECmetadata matching
# # Match colA to genotype file
# PECmeta['study_genoIDinPheno'] = PECmeta['study_geno'].isin(phenoids)
# # Same for individual ID
# PECmeta['study_indivIDinPheno'] = PECmeta['study_indiv'].isin(phenoids)
#
# # Phenotype - PECmetadata matching w/ trimmed pheno IDs
# # Same for genoID
# PECmeta['genoIDinFixedPheno'] = PECmeta['genotypingID'].isin(fixedphenoids)
# # Same for individual ID
# PECmeta['indivIDinFixedPheno'] = PECmeta['individualID'].isin(fixedphenoids)
#
# # Total rows that have a match in both the Genotype and Phenotype files
# PECmeta.query('(genoIDinGeno or indivIDinGeno) and (study_genoIDinPheno or study_indivIDinPheno)')
#
# # Total rows that have a match in both the Genotype and FIXED Phenotype files
# PECmeta.query('(genoIDinGeno or indivIDinGeno) and (genoIDinFixedPheno or indivIDinFixedPheno)')
#
# PECmeta.to_csv('PECmetaData_phase12_jmb_columnsadded.csv')

##################################################################################################################
# VCF IDs = metadata file individualID
# PSI count IDs = metadata file column A

# meta.set_index('individualID', drop=False, append=False, inplace=True, verify_integrity=False)

fixedpsi = fixedpsi.loc[~pd.isna(fixedpsi['PSIcountID_fixed'])]
fixedpsi.set_index('PSIcountID_fixed',append=False,drop=False,inplace=True,verify_integrity=True)

# psimetamerge = fixedpsi.merge(meta[['Unnamed: 0','study','specimenID','individualID','RIN']],left_index=True,right_on='Unnamed: 0')
psimetamerge = fixedpsi.merge(meta,left_index=True,right_on='Unnamed: 0') # get ALL metadata columns
metapsimerge = meta.merge(fixedpsi,right_index=True,left_on='Unnamed: 0') # get ALL metadata columns
#psimetamerge = fixedpsi.merge(meta[['Unnamed: 0','study','specimenID','individualID']],left_index=True,right_on='specimenID') # different result b/c Unnamed: 0 != specimenID always
#psimetajoin = fixedpsi.join(meta[['Unnamed: 0','study','specimenID','individualID']], on='PSIcountID_fixed', how='inner')

# Join VCF IDs to the merged PSI counts-PEC metadata table
psimetamerge.set_index('individualID',append=False,drop=False,inplace=True,verify_integrity=False)
psimetamerge['VCFID'] = genoids

# Sanity check here: should return ONLY GTEX IDs
genoids.loc[~genoids.isin(meta['individualID'])]

# Rows in common between VCF and PSI counts files
genophenometamatched = psimetamerge[~pd.isna(psimetamerge['VCFID'])]

# VCF IDs that aren't found in the metadata-phenotype merge
genoids_notinmerged = genoids.loc[~genoids.isin(psimetamerge['individualID'])]

# Find dupes...
dupes = psimetamerge.loc[(psimetamerge['VCFID'].duplicated() & ~pd.isna(psimetamerge['VCFID']))].sort_values('VCFID')#.to_csv('PhenoSamplesMatchingToMultipleGenoIndivs.csv')
psimetamerge.loc[dupes['VCFID']]

# Check for rows that directly join between phenotypes and metadata
meta['study_specimenID'] = meta['study'].replace({'LIBD_szControl':'LIBD__szControl'}) + '_' + meta['specimenID']
# Subtract out 'study_' tag from study_sampleID columns 
studies = meta.study.unique()
studies = [x + '_' for x in studies]
meta['Unnamed: 0'].str.replace('|'.join(studies),'',regex=True)

meta['study_specimenID'] = meta['study'].replace({'LIBD_szControl':'LIBD__szControl'}) + '_' + meta['specimenID']

# Rows that don't match the fixed PSI ID file but DO match the metadata file on study_specimenID
phenoids.loc[(phenoids.isin(meta.study_specimenID)) & (~phenoids.isin(fixedpsi['PSIcountID']))]
# All rows that can be matched between the metadata file and the PSI counts with or without fixed PSI count file 
phenoids.loc[(phenoids.isin(meta.study_specimenID)) | (phenoids.isin(fixedpsi['PSIcountID']))]

### Find all phenotype IDs that can be joined to metadata
# First L join phenotype IDs to fixed PSI file
phenoids_ljoin_fixed = pd.DataFrame(phenoids).merge(fixedpsi, how='left', left_index=True, right_on='PSIcountID') 
phenoids_ljoin_fixed_ljoin_meta = phenoids_ljoin_fixed.merge(meta, how='left', left_on='PSIcountID_fixed', right_on='Unnamed: 0')

# Find rows that didn't join by way of the fixed PSI ID file but DO join on PSIcountID -> study_specimenID
phenoids_ljoin_fixed_ljoin_meta[phenoids_ljoin_fixed_ljoin_meta['Unnamed: 0'].isna()] = phenoids_ljoin_fixed.merge(meta, how='left', left_on='PSIcountID', right_on='study_specimenID')

# Now check how many of these match to VCF
phenoids_ljoin_fixed_ljoin_meta.set_index('individualID',append=False,drop=False,inplace=True,verify_integrity=False)
phenoids_ljoin_fixed_ljoin_meta['VCFID'] = genoids

# Get the unique samples from the 1310  
matched_all_3 = phenoids_ljoin_fixed_ljoin_meta[~phenoids_ljoin_fixed_ljoin_meta['VCFID'].isna() & ~phenoids_ljoin_fixed_ljoin_meta[0].isna()]

# Write metadata for all 1595 samples from 1310 unique subjects w/ RNAseq data
matched_all_3.to_csv('meta_matchedIDs_allSamplesWithGenotypes-total_1595.tsv',sep='\t',index=False, header=True)
# Write IDs of just 1310 unique subjects w/ their highest RIN value RNAseq sample
matched_all_3['VCFID2']=matched_all_3['VCFID']
matched_all_3['PSIcountID2']=matched_all_3['PSIcountID']

matched_all_3_best_sample_per_subject = matched_all_3.sort_values(['VCFID','RIN'],ascending=[True,False]).drop_duplicates('VCFID',keep='first')

matched_all_3_best_sample_per_subject.to_csv('meta_matchedIDs_highestRINSampleperSubject_total_1310.tsv', sep='\t', index=False, header=True)
# Also write out in plink --keep format
matched_all_3_best_sample_per_subject.to_csv('meta_matchedIDs_highestRINSampleperSubject_total_1310.idmap', columns=['VCFID','VCFID2','PSIcountID','PSIcountID2'], index=False, header=True,sep='\t')
#################################################################################
# List all 142 VCF IDs that DON'T have an RNAseq sample/aren't in the metadata
genoids.loc[~genoids.isin(matched_all_3.VCFID)]

# Find which subjects don't have RNAseq
genoids[~genoids.isin(matched_all_3.VCFID)]
# Find which samples in metadata don't have both genotypes and RNAseq 
meta.loc[~meta.individualID.isin(matched_all_3.individualID)].to_csv('meta_sample_in_metadata_without_either_genotypes_or_rnaseq-total_1275.tsv',sep='\t',index=False,header=True)
# Find which samples in metadata have genotypes, don't have RNAseq
meta.loc[~meta.individualID.isin(matched_all_3.individualID) & meta.individualID.isin(genoids)].to_csv('meta_sample_in_metadata_with_rnaseq_without_genotypes-total_121.tsv',sep='\t',index=False,header=True)

