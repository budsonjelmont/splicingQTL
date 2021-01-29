qtls*.txt file format:

sid:SNP id
sid_chr:SNP chr
sid_start:SNP start
sid_end:SNP end
dist:distance from SNP to phenotype
pid:Phenotype id
pid_chr:Phenotype chr
pid_start:Phenotype start
pid_end:Phenotype end
slope:Regression slope
npval:Nominal P value
bpval:Beta-permuted P value
qval:FDR-corrected P value

Addtl columns may include:
sid_enst/ensg: transcript/gene the SNP overlaps
pid_enst/ensg: transcript/gene the phenotype overlaps 
snp_in_sgene: (T/F) is the SNP located inside the same gene as the phenotype?

Other standardized QTL files (may or may not all be present depending on what data is available):
qtl_phenos.bed - 6 column bed of QTL phenotypes
qtl_phenos_unique.bed - same as above but retain only unique phenotypes (drop duplicates after sorting on qval ascending)
qtl_snps.bed - 6 column bed of QTL SNPs 
qtl_snps_unique.bed - same as above but retain only unique SNPs
qtl_leadsnps.bed - 6 column bed of lead QTL SNPs per phenotype 
qtl_leadsnps_unique.bed - Same as above but retain only unique SNPs 
qtl_snps_unique.txt - chr:pos of unique SNPs 
qtl_leadsnps_unique.txt - chr:pos of unique lead SNPs 
qtl_snps_unique.assoc - chr:pos of unique SNPs & test statistic of strongest assocation 
qtl_leadsnps_unique.assoc - chr:pos of unique lead SNPs & test statistic of strongest association 
nonqtl_snps_unique.txt - chr:pos of all SNPs that are NOT in any QTLs
