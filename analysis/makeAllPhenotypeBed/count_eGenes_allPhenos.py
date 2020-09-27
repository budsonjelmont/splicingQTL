import pandas as pd
import sys

overlapsfile = sys.argv[1] #'intron-gene_annotations.txt'
ensemblfile = sys.argv[2] #'/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/gencodeV19/gencode_ENSG_to_ENST.txt'
fqtloutfile = sys.argv[3] # If annotatefqtloutput=TRUE, path to fastqtl results file # '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/Phenotype.bed'

ol = pd.read_csv(overlapsfile, sep='\t', header=None)
ensembl = pd.read_csv(ensemblfile, sep='\t')

ol['Transcript stable ID'] = ol[9].str.split('\.',expand=True)[0]

ensembl.set_index('Transcript stable ID',inplace=True, verify_integrity=True)

# Column 9 of overlaps contains the ENST ID
#ol['gene'] = ol[15].str.extract('gene_id \"(ENSG[a-zA-Z0-9]*)\.[0-9]+\"', expand=True)
uniqueENST = ol['Transcript stable ID'].unique()
ensembl.loc[uniqueENST].dropna(inplace=True)

# Print unique genes
#print('# of unique ENSGs: ' + str(len(uniqueENST)))
print('# of unique ENSGs: ' + str(ensembl.loc[uniqueENST].dropna()['Gene stable ID'].nunique()))

# Read in fastQTL results
fqtlout = pd.read_csv(fqtloutfile,sep='\s+',header=None) 
# Use the phenotype id as the key to grab the (first) overlapping transcript ID...
fqtlout.set_index(3, inplace=True, drop=False, append=False, verify_integrity=False)
fqtlout['enst'] = ol.groupby(3)['Transcript stable ID'].first()
# ...And then use the transcript ID to get the gene ID
fqtlout.set_index('enst', inplace=True, drop=False, append=False, verify_integrity=False)
fqtlout['ensg'] = ensembl['Gene stable ID']

# Print unique genes
#print('# of unique ENSGs: ' + str(len(uniqueENST)))
print('# of unique ENSGs: ' + str(fqtlout['ensg'].dropna().nunique()))

# Write out the new file with the gene ID included
fqtlout.to_csv(fqtloutfile[:-4]+'+'+pidorsid+'_ensg.txt', sep='\t', index=False, header=True)
