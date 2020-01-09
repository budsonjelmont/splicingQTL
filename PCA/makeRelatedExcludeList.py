import pandas as pd
import sys

imissfile = sys.argv[1] #/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/plink_reports/Capstone4.sel.idsync.maf5.mind1.geno1.hwe1e-5.highLDexcl.indep15000_150_.2.missingStats.imiss
relfile = sys.argv[2] #/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/plink_reports/Capstone4.sel.idsync.maf5.mind1.geno1.hwe1e-5.highLDexcl.indep50_5_.5-relatedness-report.genome

imiss = pd.read_csv(imissfile,sep='\s+')
gen = pd.read_csv(relfile,sep='\s+')

imiss.set_index('FID',inplace=True,drop=False)

rel = gen.loc[gen.PI_HAT > 0.2] 

def buildExcludeList(rows):
  fid1 = rows['FID1']
  iid1 = rows['IID1']
  fid2 = rows['FID2']
  iid2 = rows['IID2']
  fid1_fmiss = imiss.loc[fid1,'F_MISS']
  fid2_fmiss = imiss.loc[fid2,'F_MISS']
  if fid1_fmiss == fid2_fmiss:
    return ' '.join((fid1,iid1)) 
  elif fid1_fmiss > fid2_fmiss:
    return ' '.join((fid1,iid1))
  elif fid1_fmiss < fid2_fmiss:
    return ' '.join((fid2,iid2))
  else:
    print('???')

exclude = rel.apply(buildExcludeList,axis=1) 

exclude.to_csv('SamplesToExcludeForPCA.txt',index=False,header=False)
