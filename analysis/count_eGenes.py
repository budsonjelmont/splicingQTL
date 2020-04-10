import pandas as pd
import sys

overlapsfile = sys.argv[1] #'intron-gene_annotations.txt'
annotatefqtloutput = sys.argv[2] #T/F whether to append sGene ENSG to fastQTL results data frame
ensemblfile = '/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/biomart_ENSG_to_ENST.txt'

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
print('# of unique ENSGs: ' + str(len(ensembl.loc[uniqueENST].dropna()['Gene stable ID'].unique())))

# If you don't want to add the sGene column to the fastQTL results, exit here
if not annotatefqtloutput:
  sys.exit()

# Otherwise, read the fastQTL results file
fqtloutfile = sys.argv[3] # If annotatefqtloutput=TRUE, path to fastqtl results file # '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.csv'
  
# Read in fastQTL results
fqtlout = pd.read_csv(fqtloutfile,sep='\s+') #'/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.csv'
# Use the phenotype id as the key to grab the (first) overlapping transcript ID...
fqtlout.set_index('pid', inplace=True, drop=False, append=False, verify_integrity=False)
fqtlout['enst'] = ol.groupby(3)['Transcript stable ID'].first()
# ...And then use the transcript ID to get the gene ID
fqtlout.set_index('enst', inplace=True, drop=False, append=False, verify_integrity=False)
fqtlout['ensg'] = ensembl['Gene stable ID']


