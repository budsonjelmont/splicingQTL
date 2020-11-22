# After annotating QTL results file, report # unique genes found & add a boolean column stating if the SNP is found inside the same gene as the phenotype event
import sys
import pandas as pd

indir=sys.argv[1]
infile='qtls+pid_ensg+sid_ensg.txt'

res=pd.read_csv(indir+infile,sep='\t')

print("Unique genes (pid):" + str(res['pid_ensg'].nunique()))
print("Unique genes (sid):" + str(res['sid_ensg'].nunique()))

res['snp_in_sgene']=res['pid_ensg']==res['sid_ensg']

print("QTLs with both events in same gene:" + str(res.loc[res['snp_in_sgene']].shape[0]))

res.to_csv(indir+infile[:-4]+'+in_sgene.txt',index=False,header=True,sep='\t')
