module load ucsc-utils/2015-04-07

chainfile=/hpc/packages/minerva-centos7/liftover/09-Jul-2019/hg19ToHg38.over.chain

# Liftover our sQTLs
#indir=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/
#infiles=(allphenos.bed allphenos_unique.bed allSNPS.bed allSNPS_unique.bed sqtlSNPS_unique.bed sqtlSNPS.bed sqtlphenos_unique.bed sqtlphenos.bed)
#outdir=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/hg38_liftOver/
# Liftover isoQTL transcript bed
indir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoqtls/
infiles=(isoqtl_pheno_tscript.bed)
outdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoqtls/hg38_liftOver/

mkdir -p ${outdir}

for file in ${infiles[@]}
do
  sed 's/^/chr/g' ${indir}${file} > ${outdir}${file}_hg19withChr
  liftOver -bedPlus=6 ${outdir}${file}_hg19withChr ${chainfile} ${outdir}${file%.bed}_hg38.bed ${outdir}${file%.bed}_unmapped.bed 
done
