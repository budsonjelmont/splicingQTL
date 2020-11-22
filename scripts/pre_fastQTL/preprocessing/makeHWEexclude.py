import pandas as pd
import sys

plinkfile = sys.argv[1]
pthresh = float(sys.argv[2])

df = pd.read_csv(plinkfile, sep='\s+')

df.loc[df['P']<pthresh,['SNP']].to_csv('hwe_excludes.out',index=False,header=False)
