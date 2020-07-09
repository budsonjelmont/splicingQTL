# Make bed files of all & unique SNPs from the PEC isoQTL results files 
import pandas as pd
import numpy as np

#allsnpbedfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoqtls/humanSNPs/allSNPs.bed'
indir = '/sc/hydra/projects/pintod02c/datasets-external/PEC/'
infile = 'DER-10b_hg38_isoQTL.FPKM5.all.txt'
outdir = '/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoqtls/hg38_liftOver/'
#snpoutfile = 'DER-10b_hg38_isoQTL.FPKM5.all.SNPs.bed'
#uqsnpoutfile = 'DER-10b_hg38_isoQTL.FPKM5.all.uniqueSNPs.bed'

inbed = pd.read_csv(indir+infile,sep='\t')

renamer={'SNP_chr':'chr','SNP_start':'start','SNP_end':'end','SNP_id':'sid','transcript_id':'tscrid','strand':'strand'}
inbed.rename(columns=renamer,inplace=True)

# Modify transcript ID to drop version suffix
inbed['tscrid'] = inbed['tscrid'].str.split('.',expand=True)[0]

# Sort on FDR
inbed.sort_values('FDR',ascending=True,inplace=True)

# Write out BEDfile of all SNPs
inbed.to_csv(outdir+infile[:-4]+'.SNPs.bed',sep='\t',header=None,columns=[x for x in renamer.values()],index=False)

# Remove duplicate SNPs and write unique SNPs bedfile
inbed.drop_duplicates('sid',inplace=False).to_csv(outdir+infile[:-4]+'.uniqueSNPs.bed',sep='\t',header=None,columns=[x for x in renamer.values()],index=False)
# Remove duplicate Phenotypes and write unique SNPs bedfile
inbed.drop_duplicates('tscrid',inplace=False).to_csv(outdir+infile[:-4]+'.besthitSNPuniquephenos.bed',sep='\t',header=None,columns=[x for x in renamer.values()],index=False)
