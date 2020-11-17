# Combine the results of 2 VEP runs (e.g. QTL and non-QTL SNPs) into one data frame, then extract only SNPs falling with the fastQTL cis-window, leaving only those that can be sampled for bootstrapping.

import pandas as pd

# VEP files to subtract (file1-file2)
#file1='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/OurStudy_allSNPs'
#file2='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/OurStudy_sQTLSNPs'
file1='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/PEC_allQTLSNPs'
file2='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/PEC_tQTLSNPs'

# BED file containing IDs in the 4th column for ALL SNPs that can be used for resampling
#bedfile='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY.snps_in100000bpwindow.bed'
bedfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PEC/SNP_Information_Table_with_Alleles_all+strand_in100000bpwindow.bed'

#outfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/OurStudy_allnonsQTLSNPs_withinCisWindow'
outfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/PEC_allnontQTLSNPs_withinCisWindow'

vepres1=pd.read_csv(file1,sep='\t',comment='#',header=None)
vepres2=pd.read_csv(file2,sep='\t',comment='#',header=None)

# Combine VEP results
vepres=vepres1.loc[~vepres1[0].isin(vepres2[0])]

# Get IDs from BED file and select corresponding rows from the VEP results
bed = pd.read_csv(bedfile, sep='\t', header=None)

# Write out just the rows that match an ID in the bed file
vepres.loc[vepres[0].isin(bed[3])].to_csv(outfile, sep='\t', index=False, header=False)
