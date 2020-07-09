# After liftover to hg38, remove sQTL SNPs with phenotypes that did not lift over since they will throw errors in QTLtools' fenrich

#indir='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/hg38_liftOver/'
#phenofile=allphenos_unique_hg38.bed
#snpfiles=(sqtlSNPS_unique_hg38.bed sqtlSNPS_hg38.bed)

indir='/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoqtls/hg38_liftOver/'
phenofile=isoqtl_pheno_tscript_hg38.bed
snpfiles=(DER-10b_hg38_isoQTL.FPKM5.all.uniqueSNPs.bed DER-10b_hg38_isoQTL.FPKM5.all.besthitSNPuniquephenos.bed)

for snps in ${snpfiles[@]}
do 
  awk 'BEGIN{OFS=FS="\t"}; NR==FNR{F1[$4];next}{if($5 in F1)print $0;}' ${indir}${phenofile} ${indir}${snps} > ${indir}$(basename ${snps} .bed)_pruned.bed
done

