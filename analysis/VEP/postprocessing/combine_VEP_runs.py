# Combine the results of 2 VEP runs (e.g. QTL and non-QTL SNPs) into one data frame, then extract only SNPs falling with the fastQTL cis-window, leaving only those that can be sampled for bootstrapping.

import pandas as pd

# VEP files to combine
#file1='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111520/OurStudy_sQTLSNPs'
#file2='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111520/OurStudy_nonsQTLSNPs'
file1='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111520/OurStudy_allSNPs'

# BED file containing IDs in the 4th column for ALL SNPs that can be used for resampling
bedfile='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/geno_wasp/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY.snps_in100000bpwindow.bed'

outfile='/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111520/OurStudy_allSNPs_withinCisWindow'

vepres1=pd.read_csv(file1,sep='\t',comment='#',header=None)
vepres2=pd.read_csv(file2,sep='\t',comment='#',header=None)

# Combine VEP results
vepres=vepres1.append(vepres2)

# Get IDs from BED file and select corresponding rows from the VEP results
bed = pd.read_csv(bedfile, sep='\t', header=None)

# Write out just the rows that match an ID in the bed file
vepres.loc[vepres[0].isin(bed[3])].to_csv(outfile, sep='\t', index=False, header=False)
