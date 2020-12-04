# Directory containing the fastQTL output
resdir=$1 # e.g./sc/arion/projects/EPIASD/splicingQTL/output/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_55HCPs
#isnominal=$2 # # {TRUE,FALSE}--not used currently, calculateFQTLstats.R counts # columns to determine if results are from permutation or nominal pass
fdrmethod=$2 # {'st','bh'}
annotdir=$3 # e.g. /sc/arion/projects/EPIASD/splicingQTL/ext_data/other/gencode/hg19/v33lift37
annotbed=$4 # e.g. gencode.v33lift37.annotation.bed
annotidmap=$5 # e.g. gencode_v33_idmap.txt 

# Subdir to place output files
subdir=${fdrmethod}_fdr/${annotbed%.*}_annotated

mkdir -p $resdir/$subdir

# Combine results for all chromosomes into a single file
cat $resdir/chr*/chr*_norm_perm1000.txt.gz.* > $resdir/$subdir/chrAll_combined

# Calculate stats on sQTL results
Rscript calculateFQTLstats.R $resdir/$subdir $fdrmethod 

# Make files for downstream analyses 
python3 make_analysis_files-ourSNPs.py $resdir/$subdir /sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/hg38_liftover/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.hg38.bed NA NA

# Add SNP/phenotype gene annotations
../addQTLgenes/count_eGenes.lsf $resdir/$subdir sqtlphenos_unique.bed qtls.txt pid $annotdir $annotbed $annotidmap
../addQTLgenes/count_eGenes.lsf $resdir/$subdir sqtlSNPS_unique.bed qtls+pid_ensg.txt sid $annotdir $annotbed $annotidmap

# Add boolean column to indicate if SNP & phenotype are located in the same gene
python3 ../addQTLgenes/get_counts_add_samegeneTF.py $resdir/$subdir
