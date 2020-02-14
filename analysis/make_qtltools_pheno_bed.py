# Make .BED file of ALL phenotypes (intron clusters) used as input to FastQTL
import pandas as pd
import glob

phenoDir='/sc/hydra/scratch/belmoj01/splicingQTL/'
phenoFile='out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.deduped.gz.qqnorm_allChr'
combojuncFile='/sc/hydra/projects/pintod02c/WASP_leafcutter/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01/junctions/all_juncs_combined.junc'

juncpath = '/sc/orga/projects/pintod02b/capstone/leafcutter/data-freeze4/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01/junctions'
all_juncs = glob.glob(juncpath + '/*.junc')

# Read and concatenate all .junc files
#li = []

#for filename in all_juncs:
#    df = pd.read_csv(filename, sep='\s+', header=None)
#    li.append(df)

#junc = pd.concat(li, axis=0, ignore_index=True)

# Read combined .junc file directly
junc = pd.read_csv(combojuncFile, sep='\s+', header=None)
#junc = pd.read_csv('/sc/orga/projects/pintod02b/capstone/leafcutter/data-freeze4/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01/junctions/all_juncs_combined_headn1000.junc', sep='\s+', header=None)

# Strip 'chr' from chromosome column
junc[0] = junc[0].str.replace('chr','').str.strip()

pheno = pd.read_csv(phenoDir + phenoFile, sep='\s+')

# Make first 4 BED columns (feature chr, start, end, ID)
bed = pheno[['start','end','ID']] 

#bed.insert(0, 'chr', pheno['#Chr'].str.split(':',expand=True)[1].str.strip())
bed.insert(0, 'chr', pheno['#Chr'])

# Remove any rows that are artifacts of file concatenation
bed = bed.loc[bed['chr']!='#Chr']

# Type conversion
bed['chr'] = bed['chr'].astype(str)
bed['start'] = bed['start'].astype(int)
junc[0] = junc[0].astype(str)

# 4th column is dummy column of top variants -- I use a single variant ID for all entries
bed[4] = 'whatever'

# The 'end' coordinate in the junction file seems to be 1 greater than the .junc file coordinate
bed['end-1'] = bed['end'].astype(int)-1
#bed['end-1'] = bed['end'].astype(int)

# Get strand orientation (column 5 in .junc file) by joining to FIRST MATCH in .junc file
leftjoined = bed.merge(junc, how='left', left_on=['chr','start','end-1'], right_on=[0,1,2])
bed = leftjoined.groupby(['chr','start','end'])[['ID','4_x',5]].first().reset_index()
# Report on how many rows don't match, then throw them out
print(str(bed.loc[~pd.isna(bed[5])].shape[0]) + ' rows out of ' + str(bed.shape[0]) + 'matched to a .junc file intron')
bed = bed.loc[~pd.isna(bed[5])]

bed.to_csv(phenoDir + 'Phenotype.bed', index=False, header=False, sep='\t')
