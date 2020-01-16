import pandas as pd
import sys

phenobedfile = sys.argv[1]
sqtlsfile = sys.argv[2]

outfile = sqtlsfile[:-3] + 'SNPs.bed'

pheno = pd.read_csv(phenobedfile, sep='\t', header=None)
sqtls = pd.read_csv(sqtlsfile, sep='\s+')

# Get SNP coords
sqtls['sid_end'] = sqtls['sid'].str.split(':', expand=True)[1].astype(int)
sqtls['sid_start'] = sqtls['sid_end']#-1

pheno.set_index(3, drop=False, append=False, inplace=True, verify_integrity=True)
sqtls.set_index('pid', drop=False, append=False, inplace=True, verify_integrity=False)

joined = sqtls.join(pheno)

joined.to_csv(outfile, sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
#pheno.loc[pheno[3].isin(sqtls['pid'])].to_csv(outfile, sep='\t', header=None, index=False)

# Now make the bed of uniqueSNPs only

outfile = sqtlsfile[:-3] + 'uniqueSNPs.bed' 

joined.drop_duplicates('sid', inplace=True)

joined.to_csv(outfile, sep='\t', columns=[0,'sid_start','sid_end','sid',3,5], header=None, index=False)
