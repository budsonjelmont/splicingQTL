module load ucsc-utils/2015-04-07

chainfile=/hpc/packages/minerva-centos7/liftover/09-Jul-2019/hg19ToHg38.over.chain

# Liftover our sQTLs
indir=/sc/arion/projects/EPIASD/splicingQTL/output/vcf/
infiles=(Capstone4.sel.hasPhenosOnly.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.snps.bed)
outdir=/sc/arion/projects/EPIASD/splicingQTL/output/hg38_liftOver/

mkdir -p ${outdir}

for file in ${infiles[@]}
do
  ./lifthg19to38.sh  $indir $file $outdir false
done
