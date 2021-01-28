# Directory containing the fastQTL output
resdir=$1 # e.g./sc/arion/projects/EPIASD/splicingQTL/output/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_55HCPs
#isnominal=$2 # # {TRUE,FALSE}--not used currently, calculateFQTLstats.R counts # columns to determine if results are from permutation or nominal pass
#fdrmethod=$2 # {'st','bh'}
annotdir=$2 # e.g. /sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp/ 
annotfile=$3 # e.g. out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.clustids.gcv33ids 
genedir=$4 # e.g. /sc/arion/projects/EPIASD/splicingQTL/ext_data/other/gencode/hg19/v33lift37/ 
genebed=$5 # e.g. gencode.v33lift37.annotation.geneModel.geneFeatOnly.bed 

# Subdir within each nHCPs/ to place output files
#subdir=${fdrmethod}_fdr/${annotbed%.*}_annotated
subdir=wrapup
# Dir to write the called QTLs (after nominal p-val threshold)
qtldir=qtls

modes=(permute nominal)

for mode in "${modes[@]}"
do
  
  # Combine results for all chromosomes into a single file
  mkdir -p $resdir/$mode/$subdir
  cat $resdir/$mode/chr*/chr*_norm_perm1000.txt.gz.* > $resdir/$mode/$subdir/chrAll_combined
  if [ $mode == 'permute' ]
  then 
   Rscript call_perm.R --resdir $resdir/$mode/$subdir/ 
  elif [ $mode == 'nominal' ]
  then
   Rscript call_nominal.R --resdir $resdir/$mode/$subdir/
  fi 

  # Calculate stats on sQTL results
  #Rscript calculateFQTLstats.R $resdir/$mode/$subdir $fdrmethod 

  # Make files for downstream analyses 
  #python3 make_analysis_files-ourSNPs.py $resdir/$subdir /sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/hg38_liftover/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.hg38.bed NA NA

  # Add SNP/phenotype gene annotations
  #addQTLgenes/count_eGenes.lsf $resdir/$subdir sqtlphenos_unique.bed qtls.txt pid $annotdir $annotbed $annotidmap
  #addQTLgenes/count_eGenes.lsf $resdir/$subdir sqtlSNPS_unique.bed qtls+pid_ensg.txt sid $annotdir $annotbed $annotidmap
  # Add

  # Add boolean column to indicate if SNP & phenotype are located in the same gene
  #python3 addQTLgenes/get_counts_add_samegeneTF.py $resdir/$subdir
done
# Now that we've called nominal & permuted results, use the nominal p-val threshold derived from the permuted results to call our QTLs
rm -f $resdir/$qtldir/qtl*
mkdir -p $resdir/$qtldir
python3 call_qtl_chunk.py --egenes_file $resdir/permute/$subdir/chrAll_combined.FDR05 --nominal_results $resdir/nominal/$subdir/chrAll_combined  --output_dir $resdir/$qtldir

# Make standardized output files
echo "python3 make_analysis_files-ourSNPs.py $resdir/$qtldir qtl_chunk_concat.txt /sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/hg38_liftover/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.hg38.bed NA NA"
python3 make_analysis_files-ourSNPs.py $resdir/$qtldir qtl_chunk_concat.txt /sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/hg38_liftover/Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.hg38.bed NA NA

# Annotate sGenes
echo "python3 addQTLgenes_leafcutter/annotate_sGenes.py $resdir/$qtldir/qtls.txt $annotdir/$annotfile" 
python3 addQTLgenes_leafcutter/annotate_sGenes.py $resdir/$qtldir/qtls.txt $annotdir/$annotfile 
# Annotate genes overlapping SNPs
echo "addQTLgenes/count_eGenes.lsf $resdir/$qtldir sqtlSNPS_unique.bed qtls+pid_ensg.txt sid $genedir $genebed"
addQTLgenes/count_eGenes.lsf $resdir/$qtldir sqtlSNPS_unique.bed qtls+pid_ensg.txt sid $genedir $genebed
# Annotate sGenes (bedtools)
echo "addQTLgenes/count_eGenes.lsf $resdir/$qtldir sqtlSNPS_unique.bed qtls+pid_ensg+sid_ensg.txt pid $genedir $genebed"
addQTLgenes/count_eGenes.lsf $resdir/$qtldir allphenos_unique.bed qtls+pid_ensg+sid_ensg_jmb.txt pid $genedir $genebed
