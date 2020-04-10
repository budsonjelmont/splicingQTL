# Reads in the GENCODE GTF, extracts transcript features, and creates a map of ENST->ENSG
# Emulates the format of this file: /sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/biomart_ENSG_to_ENST.txt, since many of the downstream analysis scripts here were written around that file 
import pandas as pd

gtffile = 'gencode.v19.annotation.gtf'

gtf = pd.read_csv(gtffile, sep='\t', skiprows=5, header=None)

# Extract transcripts
gtf = gtf.loc[gtf[2]=='transcript']
# Get transcript & gene IDs
# Headers: Gene stable ID	Gene stable ID version	Transcript stable ID	Transcript stable ID version
gtf['Gene stable ID version']=gtf[8].str.extract('gene_id \"(ENSG[a-zA-Z0-9]*\.[0-9]+)\"')
gtf['Transcript stable ID version']=gtf[8].str.extract('transcript_id \"(ENST[a-zA-Z0-9]*\.[0-9]+)\"')
gtf['Gene symbol']=gtf[8].str.extract('gene_name \"([a-zA-Z0-9\.\\-]*)\"')
gtf['Gene stable ID'] = gtf['Gene stable ID version'].str.split('\.',expand=True)[0]
gtf['Transcript stable ID'] = gtf['Transcript stable ID version'].str.split('\.',expand=True)[0]

gtf.to_csv('gencode_ENSG_to_ENST.txt', sep='\t', index=False, header=True, columns=['Gene stable ID','Gene stable ID version','Transcript stable ID','Transcript stable ID version','Gene symbol'])
