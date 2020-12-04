import pandas as pd
import sys

overlapsfile = sys.argv[1] #'intron-gene_annotations.txt'
annotatefqtloutput = sys.argv[2] #T/F whether to append sGene ENSG to fastQTL results data frame
ensemblfile = sys.argv[3] #'/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/gencodeV19/gencode_ENSG_to_ENST.txt'
fqtloutfile = sys.argv[4] # If annotatefqtloutput=TRUE, path to fastqtl results file # '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.csv'
pidorsid = sys.argv[5] # {'pid','sid'} to determine if the overlap is being calculated for SNPs or phenotypes (intron)

if pidorsid not in ['pid','sid']:
  print('ERROR: Argument 5 must be one of \'pid\' or \'sid\'')
  exit()
else:
  tscrcol=pidorsid+'_enst'
  genecol=pidorsid+'_ensg'


ol = pd.read_csv(overlapsfile, sep='\t', header=None)
ensembl = pd.read_csv(ensemblfile, sep='\t')

ol['Transcript stable ID'] = ol[9].str.split('\.',expand=True)[0]

# Drop duplicate rows or rows where Transcript stable ID is NA before setting index
ensembl.dropna(axis='rows',subset=['Transcript stable ID'],inplace=True)
ensembl.drop_duplicates(subset='Transcript stable ID',inplace=True)

ensembl.set_index('Transcript stable ID',inplace=True, verify_integrity=True)

# Column 9 of overlaps contains the ENST ID
#ol['gene'] = ol[15].str.extract('gene_id \"(ENSG[a-zA-Z0-9]*)\.[0-9]+\"', expand=True)
uniqueENST = ol['Transcript stable ID'].unique()
ensembl.loc[uniqueENST].dropna(inplace=True)

# Read in fastQTL results
fqtlout = pd.read_csv(fqtloutfile,sep='\s+') #'/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.csv'
# Use the phenotype id as the key to grab the (first) overlapping transcript ID...
fqtlout.set_index(pidorsid, inplace=True, drop=False, append=False, verify_integrity=False)
fqtlout[tscrcol] = ol.groupby(3)['Transcript stable ID'].first()
# ...And then use the transcript ID to get the gene ID
fqtlout.set_index(tscrcol, inplace=True, drop=False, append=False, verify_integrity=False)
fqtlout[genecol] = ensembl['Gene stable ID']

# Print unique genes
#print('# of unique ENSGs: ' + str(len(uniqueENST)))
print('# of unique ENSGs: ' + str(len(fqtlout[genecol].dropna().unique())))

# Write out the new file with the gene ID included
fqtlout.to_csv(fqtloutfile[:-4]+'+'+pidorsid+'_ensg.txt', sep='\t', index=False, header=True)
